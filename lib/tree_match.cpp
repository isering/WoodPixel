/*
WoodPixel - Supplementary code for Computational Parquetry:
            Fabricated Style Transfer with Wood Pixels
            ACM Transactions on Graphics 39(2), 2020

Copyright (C) 2020 Julian Iseringhausen <opensource@iseringhausen.graphics>
Copyright (C) 2020 Matthias Hullin, University of Bonn <hullin@cs.uni-bonn.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tree_match.hpp"

#include <deque>
#include <functional>
#include <numeric>
#include <set>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/tokenizer.hpp>
#include <omp.h>
#include <opencv2/saliency.hpp>
#include <opencv2/saliency/saliencySpecializedClasses.hpp>

#include "adaptive_patch.hpp"
#include "affine_transformation.hpp"
#include "bezier_ransac.hpp"
#include "convex_hull.hpp"
#include "generate_patches.hpp"
#include "histogram.hpp"
#include "line_segment_detect.hpp"
#include "opencv_extra.hpp"
#include "texture.hpp"
#include "timer.hpp"
#include "tree_match.hpp"

namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

TreeMatch::TreeMatch(int min_patch_size, int patch_levels, double patch_quality_factor, int filter_resolution, double frequency_octaves, int num_filter_directions) :
m_patch_quality_factor(patch_quality_factor),
m_subpatch_size(min_patch_size/4, min_patch_size/4),
m_filter_bank(filter_resolution, frequency_octaves, num_filter_directions)
{
  cv::Point boundary_size(min_patch_size/4, min_patch_size/4);
  cv::Point current_patch_size(min_patch_size, min_patch_size);
  for (int i = 0; i < patch_levels; ++i)
  {
    m_patch_sizes.emplace_back(current_patch_size);
    current_patch_size = 2 * current_patch_size - boundary_size;
  }
}

void TreeMatch::add_target(const boost::filesystem::path& path, double dpi, double scale)
{
  /*
  * Load target texture
  */
  m_targets.emplace_back(Texture(path, dpi, scale));
  m_target_images.push_back(cv::imread(path.string(), cv::IMREAD_COLOR));
}

cv::Mat TreeMatch::compute_priority_map(const cv::Mat& texture_in)
{
  cv::Mat texture;
  texture_in.convertTo(texture, CV_32FC3, 1.0 / 65535.0);
  cv::cvtColor(texture, texture, cv::COLOR_BGR2GRAY);

  cv::GaussianBlur(texture, texture, cv::Size(3, 3), 0.0);
  cv::Laplacian(texture, texture, CV_32FC1);
  cv::normalize(cv::abs(texture), texture, 1.0, 0.0, cv::NORM_MINMAX);
  return texture;
}

void TreeMatch::add_texture(const boost::filesystem::path& path, double dpi, double scale, int num_rotations, const TextureMarker& marker, const std::string id)
{
  m_textures.emplace_back();
  m_textures.back().emplace_back(path, dpi, scale, marker, id);
  for (int i = 1; i < num_rotations; ++i)
  {
    m_textures.back().push_back(m_textures.back()[0].rotate(-i * 2.0f * boost::math::float_constants::pi / num_rotations));
  }
}

void TreeMatch::add_texture(const boost::filesystem::path& path, const boost::filesystem::path& mask, double dpi, double scale, int num_rotations, const TextureMarker& marker, const std::string id)
{
  if (mask.empty() || !fs::exists(mask) || !fs::is_regular_file(mask))
  {
    add_texture(path, dpi, scale, num_rotations, marker, id);
  }
  else
  {
    m_textures.emplace_back();
    m_textures.back().emplace_back(path, mask, dpi, scale, marker, id);
    for (int i = 1; i < num_rotations; ++i)
    {
      m_textures.back().push_back(m_textures.back()[0].rotate(-i * 2.0f * boost::math::float_constants::pi / num_rotations));
    }
  }
}

void TreeMatch::compute_responses(double weight_intensity, double weight_sobel, double weight_gabor, double histogram_matching_factor)
{
  FeatureEvaluator evaluator(weight_intensity, weight_sobel, weight_gabor, m_filter_bank);
  cv::Mat kernel = cv::Mat::ones(evaluator.max_filter_size(), CV_8UC1);

  if (m_textures.empty())
  {
    std::cerr << "WARNING: No textures supplied to TreeMatch::compute_responses." << std::endl;
    return;
  }

  std::cout << "Compute target responses..." << std::endl;
  for (Texture& target : m_targets)
  {
    if (histogram_matching_factor != 0.0)
    {
      std::vector<Texture> textures_histogram;
      for (const std::vector<Texture>& textures : m_textures)
      {
        textures_histogram.push_back(textures.front());
      }
      target.response = evaluator.evaluate_with_histogram_matching(target.texture, textures_histogram, target.mask_rotation, histogram_matching_factor);
    }
    else
    {
      target.response = evaluator.evaluate(target.texture, target.mask());
    }

    cv::erode(target.mask_rotation, target.mask_rotation, kernel, cv::Point(-1, -1), 1, cv::BORDER_CONSTANT, cv::Scalar(0));
  }
  std::cout << "done" << std::endl;

  int num_textures = 0;
  for (std::vector<Texture>& textures_rot : m_textures)
  {
    num_textures += static_cast<int>(textures_rot.size());
  }

  int i = 0;
  std::cout << "Compute " << num_textures << " texture features..." << std::endl;
  for (std::vector<Texture>& textures_rot : m_textures)
  {
    for (Texture& t : textures_rot)
    {
      i++;
      if ((i % 10) == 0)
      {
        std::cout << "(" << i << ")" << std::flush;
      }
      std::cout << "." << std::flush;

      t.response = evaluator.evaluate(t.texture, t.mask());
      cv::erode(t.mask_rotation, t.mask_rotation, kernel, cv::Point(-1, -1), 1, cv::BORDER_CONSTANT, cv::Scalar(0));
    }
  }
  std::cout << "done" << std::endl;
}

static bool is_valid_point(cv::Point p, cv::Mat mask)
{
  return p.x >= 0 && p.x < mask.cols &&
    p.y >= 0 && p.y < mask.rows &&
    mask.at<unsigned char>(p);
}

static cv::Mat rotated_rect_mask(cv::Size size, float angle)
{
  cv::RotatedRect rotated_rect(cv::Point2f(0.0f, 0.0f), size, angle);
  rotated_rect.center -= rotated_rect.boundingRect2f().tl();
  cv::Mat mask = cv::Mat::zeros(rotated_rect.boundingRect().size(), CV_8UC1);
  std::vector<cv::Point2f> points2f(4);
  rotated_rect.points(points2f.data());
  std::vector<cv::Point> points(points2f.begin(), points2f.end());
  cv::fillConvexPoly(mask, points, cv::Scalar(1));
  cv::imshow("Mask", 255*mask);
  cv::waitKey(0);
  return mask;
}

cv::Mat TreeMatch::fit_single_patch(const boost::filesystem::path& filename)
{
  if (m_textures.empty())
  {
    throw(std::invalid_argument("TreeMatch::match_patch_impl called but no textures supplied."));
  }

  FeatureEvaluator evaluator(0.5, 0.5, 0.0, m_filter_bank);
  Texture patch(filename, 300.0, 1.0);
  patch.response = evaluator.evaluate(patch.texture, patch.mask());
  
  struct MatchPatchResult
  {
    MatchPatchResult(int texture_index, int texture_rot) :
      texture_index(texture_index),
      texture_rot(texture_rot),
      cost(std::numeric_limits<double>::max()),
      texture_pos(-1, -1)
    {
    }

    int texture_index;
    int texture_rot;
    double cost;
    cv::Point texture_pos;
    cv::Mat error;
  };

  std::vector<MatchPatchResult> results;
  for (size_t i = 0; i < m_textures.size(); ++i)
  {
    for (size_t j = 0; j < m_textures[i].size(); ++j)
    {
      results.emplace_back(static_cast<int>(i), static_cast<int>(j));
    }
  }

#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(results.size()); ++i)
  {
    cv::Mat texture_mask;
    cv::erode(m_textures[results[i].texture_index][results[i].texture_rot].mask(), texture_mask, patch.mask(), cv::Point(0, 0), 1, cv::BORDER_CONSTANT, 0);
    texture_mask = texture_mask(cv::Rect(0, 0, texture_mask.cols - patch.response.cols() + 1, texture_mask.rows - patch.response.rows() + 1));

    if (cv::countNonZero(texture_mask) > 0)
    {
      cv::Mat match = m_textures[results[i].texture_index][results[i].texture_rot].template_match(patch);
      cv::minMaxLoc(match, &results[i].cost, 0, &results[i].texture_pos, 0, texture_mask);
    }
  }

  MatchPatchResult result_min = *std::min_element(results.begin(), results.end(), [](const MatchPatchResult& lhs, const MatchPatchResult& rhs){return lhs.cost < rhs.cost; });

  if (result_min.cost == std::numeric_limits<double>::max())
  {
    throw(std::runtime_error("No texture samples available."));
  }

  cv::Mat result = m_textures[result_min.texture_index][result_min.texture_rot].texture(cv::Rect(result_min.texture_pos, patch.texture.size()));
  result.convertTo(result, CV_8UC3, 1.0/255.0);

  return result;
}

bool TreeMatch::find_next_patch()
{
  if (m_reconstruction_regions.empty())
  {
    return false;
  }
  
  const PatchRegion region = m_reconstruction_regions.front();
  m_reconstruction_regions.pop_front();
 
  std::vector<Patch> patches = match_patch(region);

  for (const Patch& patch : patches)
  {
    if (patch.cost() == std::numeric_limits<double>::max())
    {
      return false;
    }
  }

  return true;
}

bool TreeMatch::find_next_patch_adaptive()
{
  if (m_patch_sizes.size() == 1)
  {
    return find_next_patch();
  }

  if (m_reconstruction_regions.empty())
  {
    return false;
  }

  const PatchRegion region = m_reconstruction_regions.front();
  m_reconstruction_regions.pop_front();

  /*
   * Hierarchical patch matching.
   */
  std::deque<AdaptivePatch> queue(1, AdaptivePatch(match_patch(region), static_cast<int>(m_patch_sizes.size()-1)));
  mask_patch_resources(queue.front());

  while (!queue.empty())
  {
    AdaptivePatch patch_coarse = queue.front();
    queue.pop_front();

    if (patch_coarse.cost() == std::numeric_limits<double>::max())
    {
      return false;
    }

    if (patch_coarse.level == 0)
    {
      std::cout << "add_patch_level: 0" << std::endl;
      for (const Patch& p : patch_coarse.patches)
      {
        add_patch(p);
      }
    }
    else
    {
      /*
       * Release temporary hold of patch area.
       */
      unmask_patch_resources(patch_coarse);

      /*
       * Get cost for finer solution.
       */
      double cost_fine = 0.0;
      int num_pixels_fine = 0;

      const int level_fine = patch_coarse.level - 1;
      std::vector<AdaptivePatch> patches_fine;
      cv::Size patch_size_fine = m_patch_sizes[level_fine];
      for (int y = 0; y < 2; ++y)
      {
        for (int x = 0; x < 2; ++x)
        {
          int x_fine = patch_coarse.patches[0].anchor_target().x + x * (patch_size_fine.width - m_subpatch_size.width);
          int y_fine = patch_coarse.patches[0].anchor_target().y + y * (patch_size_fine.height - m_subpatch_size.height);
          cv::Rect region_fine(cv::Point(x_fine, y_fine), patch_size_fine);
          
          AdaptivePatch patch_fine(match_patch(PatchRegion(patch_coarse.target_index(), patch_coarse.coordinate(), region_fine)), level_fine);

          if (patch_fine.cost() == std::numeric_limits<double>::max())
          {
            return false;
          }

          patches_fine.push_back(patch_fine);
          mask_patch_resources(patch_fine);
          cost_fine += patch_fine.cost();
          num_pixels_fine += patch_fine.num_pixels();
        }
      }

      /*
      * Decide which solution to accept.
      */
      const double cost_factor = (patch_coarse.cost() / patch_coarse.num_pixels()) / (cost_fine / num_pixels_fine);
      std::cout << "cost_factor: " << cost_factor << std::endl;

      if (cost_factor < m_patch_quality_factor)
      {
        std::cout << "add_patch_level: " << patch_coarse.level << std::endl;
        for (const AdaptivePatch& p : patches_fine)
        {
          unmask_patch_resources(p);
        }
        for (const Patch& p : patch_coarse.patches)
        {
          add_patch(p);
        }
      }
      else
      {
        queue.insert(queue.end(), patches_fine.begin(), patches_fine.end());
      }
    }
  }

  return true;
}

Patch TreeMatch::match_patch_impl(const PatchRegion& region, cv::Mat mask)
{
  if (m_textures.empty())
  {
    throw(std::invalid_argument("TreeMatch::match_patch_impl called but no textures supplied."));
  }

  struct MatchPatchResult
  {
    MatchPatchResult(int texture_index, int texture_rot) :
      texture_index(texture_index),
      texture_rot(texture_rot),
      cost(std::numeric_limits<double>::max()),
      texture_pos(-1, -1)
    {}

    int texture_index;
    int texture_rot;
    double cost;
    cv::Point texture_pos;
    cv::Mat error;
  };

  const bool is_rectangular = (mask.empty() || cv::countNonZero(mask) == mask.rows*mask.cols);
 

  Texture kernel = m_targets[region.target_index()](region.bounding_box());

  std::vector<MatchPatchResult> results;
  for (size_t i = 0; i < m_textures.size(); ++i)
  {
    for (size_t j = 0; j < m_textures[i].size(); ++j)
    {
      results.emplace_back(static_cast<int>(i), static_cast<int>(j));
    }
  }
   
  #pragma omp parallel for
  for (int i = 0; i < static_cast<int>(results.size()); ++i)
  {
    cv::Mat texture_mask;
    cv::erode(m_textures[results[i].texture_index][results[i].texture_rot].mask(), texture_mask, region.mask(), cv::Point(0, 0), 1, cv::BORDER_CONSTANT, 0);
    texture_mask = texture_mask(cv::Rect(0, 0, texture_mask.cols - kernel.response.cols() + 1, texture_mask.rows - kernel.response.rows() + 1));

    if (cv::countNonZero(texture_mask) > 0)
    {
      cv::Mat match;
      if (is_rectangular)
      {
        match = m_textures[results[i].texture_index][results[i].texture_rot].template_match(kernel);
      }
      else
      {
        match = m_textures[results[i].texture_index][results[i].texture_rot].template_match(kernel, region.mask());
      }
      cv::minMaxLoc(match, &results[i].cost, 0, &results[i].texture_pos, 0, texture_mask);
    }
  }

  MatchPatchResult result_min = *std::min_element(results.begin(), results.end(), [](const MatchPatchResult& lhs, const MatchPatchResult& rhs){return lhs.cost < rhs.cost; });


  if (result_min.cost == std::numeric_limits<double>::max())
  {
    std::cout << "Finished. No more texture samples available." << std::endl;
  }

  const cv::Mat error_mat = m_targets[region.target_index()].response(region.bounding_box()).dist_sqr_mat(m_textures[result_min.texture_index][result_min.texture_rot].response(cv::Rect(result_min.texture_pos, region.bounding_box().size())));

  Patch patch(region, result_min.texture_pos, m_textures[result_min.texture_index][result_min.texture_rot].transformation_matrix, result_min.texture_index, result_min.texture_rot, error_mat, result_min.cost);
  add_patch(patch);

  return patch;
}


std::vector<Patch> TreeMatch::match_patch(const PatchRegion& region)
{
  if (region.has_sub_regions())
  {
    std::vector<Patch> sub_patches;
    for (const PatchRegion& sub_region : region.sub_regions())
    {
      if (!sub_region.valid())
      {
        throw(std::invalid_argument((boost::format("TreeMatch::match_patch encountered invalid subpatch as input: %d %d %d") %
          sub_region.target_index() % sub_region.coordinate().x % sub_region.coordinate().y).str()));
      }

      sub_patches.push_back(match_patch_impl(sub_region, sub_region.mask()));

      if (!sub_patches.back().region_target.valid())
      {
        throw(std::runtime_error((boost::format("TreeMatch::match_patch generated invalid subpatch: %d %d %d") %
          sub_patches.back().region_target.target_index() %
          sub_patches.back().region_target.coordinate().x %
          sub_patches.back().region_target.coordinate().y).str()));
      }
    }
    return sub_patches;
  }
  else
  {
    Patch p = match_patch_impl(region, region.mask());

    /*
    if (!p.region_target.valid())
    {
      throw(std::runtime_error((boost::format("TreeMatch::match_patch generated invalid patch: %d %d %d") %
        p.region_target.target_index() %
        p.region_target.coordinate().x %
        p.region_target.coordinate().y).str()));
    }
    */

    return std::vector<Patch>(1, p);
  }
}

struct point_less
{
  bool operator()(const cv::Point& lhs, const cv::Point& rhs)
  {
    return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
  }
};

static std::vector<cv::Point> get_queue(cv::Mat mask)
{
  std::vector<cv::Point> queue;
  cv::Mat mask_new;
  cv::dilate(mask, mask_new, cv::Mat::ones(3, 3, CV_8UC1));
  mask_new -= mask;
  cv::findNonZero(mask_new, queue);
  return std::move(queue);
}

static std::vector<cv::Point> get_convex_hull(cv::Mat mask)
{
  std::vector<cv::Point> convex_hull, convex_hull_points;
  cv::findNonZero(mask, convex_hull_points);
  cv::convexHull(convex_hull_points, convex_hull);
  return std::move(convex_hull);
}

static std::vector<cv::Point> update_convex_hull(std::vector<cv::Point> convex_hull, cv::Point point_new)
{
  convex_hull.push_back(point_new);
  std::vector<cv::Point> convex_hull_new;
  cv::convexHull(convex_hull, convex_hull_new);
  return std::move(convex_hull_new);
}

bool is_valid_mask(const cv::Mat& mask, cv::Point p_patch, const cv::Mat& mask_target, cv::Point p_target, const cv::Mat& mask_source, cv::Point p_source)
{
  int y_target = p_target.y - p_patch.y;
  int y_source = p_source.y - p_patch.y;

  for (int y = 0; y < mask.rows; ++y, ++y_target, ++y_source)
  {
    const unsigned char* ptr_patch = mask.ptr(y);
    const unsigned char* ptr_target = mask_target.ptr(y_target);
    const unsigned char* ptr_source = mask_source.ptr(y_source);

    int x_target = p_target.x - p_patch.x;
    int x_source = p_source.x - p_patch.x;

    for (int x = 0; x < mask.cols; ++x, ++x_target, ++x_source)
    {
      if (ptr_patch[x] && !(ptr_target[x_target] && ptr_source[x_source]))
      {
        return false;
      }
    }
  }

  return true;
}

float compute_cost(const cv::Mat& mask, cv::Mat& cost_mat, cv::Point p_patch, const FeatureVector& response_target, cv::Point p_target, const FeatureVector& response_source, cv::Point p_source)
{
  float cost = 0.0f;
  int num_pixel = 0;

  p_target -= p_patch;
  p_source -= p_patch;

  for (int y = 0; y < mask.rows; ++y, ++p_target.y, ++p_source.y)
  {
    const unsigned char* ptr_patch = mask.ptr(y);
    float* ptr_cost = reinterpret_cast<float*>(cost_mat.ptr(y));

    cv::Point p_target_temp = p_target;
    cv::Point p_source_temp = p_source;

    for (int x = 0; x < mask.cols; ++x, ++p_target_temp.x, ++p_source_temp.x)
    {
      if (ptr_patch[x])
      {
        if (ptr_cost[x] < 0.0f)
        {
          ptr_cost[x] = response_target.dist(p_target_temp, response_source, p_source_temp);
        }
        cost += ptr_cost[x];
        ++num_pixel;
      }
    }
  }

  return cost / num_pixel;
}

static std::pair<std::vector<cv::Point>, std::vector<float>> mask_to_vec(const cv::Mat& mask, const cv::Mat& cost_mat, cv::Point p_patch)
{
  std::vector<cv::Point> patch;
  std::vector<float> cost_vec;

  cv::Point p_target = -p_patch;
  for (int y = 0; y < mask.rows; ++y, ++p_target.y)
  {
    const unsigned char* ptr = mask.ptr(y);
    const float* ptr_cost = reinterpret_cast<const float*>(cost_mat.ptr(y));
    cv::Point p_target_temp = p_target;
    for (int x = 0; x < mask.cols; ++x, ++p_target_temp.x)
    {
      if (ptr[x])
      {
        patch.push_back(p_target_temp);
        cost_vec.push_back(ptr_cost[x]);
      }
    }
  }
  return std::make_pair(patch, cost_vec);
}

cv::Mat TreeMatch::draw(int target_index, bool draw_target) const
{
  cv::Mat image;
  if (draw_target)
  {
    image = m_targets[target_index].texture.clone();
    image.convertTo(image, CV_8UC3, 1.0/255.0);
  }
  else
  {
    image = cv::Mat::zeros(m_targets[target_index].texture.size(), CV_8UC3);
  }

  for (const Patch& p : m_patches)
  {
    if (p.region_target.target_index() == target_index)
    {
      p.draw(image, 1.0, m_textures[p.source_index][0], 1.0);
    }
  }

  return image;
}

cv::Mat TreeMatch::draw_patch(const Patch& patch) const
{
  cv::Mat image = m_textures[patch.source_index][patch.source_rot].texture(cv::Rect(patch.anchor_source, patch.size()));
  image.convertTo(image, CV_8UC3, 1.0/255.0);
  return image;
}

static cv::Mat warp_nn(cv::Mat mat, cv::Mat transformation, cv::Size size)
{
  cv::Mat result;
  cv::warpAffine(mat, result, transformation, size, cv::INTER_NEAREST);
  return result;
}

static cv::Mat combine_with_mask(cv::Mat mat_1, cv::Mat mask, int type, cv::Mat mat_2)
{
  int num_channels = mat_1.channels();
  mask.convertTo(mask, type);
  std::vector<cv::Mat> mask_vec(num_channels, mask);
  cv::merge(mask_vec, mask);
  return mat_1.mul(mask) + mat_2.mul(cv::Scalar::all(1) - mask);
}

void TreeMatch::mask_patch_resources(const Patch& patch)
{
  mask_patch_resources(patch, patch.mask());
}

void TreeMatch::mask_patch_resources(const AdaptivePatch& adaptive_patch)
{
  for (const Patch& p : adaptive_patch.patches)
  {
    mask_patch_resources(p);
  }
}

void TreeMatch::mask_patch_resources(const Patch& patch, const cv::Mat& mask)
{
  const Texture& texture_source = m_textures[patch.source_index][patch.source_rot];
  cv::Mat mask_texture_new = cv::Mat::ones(texture_source.mask_done.size(), CV_8UC1);

  if (mask.empty())
  {
    mask_texture_new(cv::Rect(patch.anchor_source, patch.size())) = 0;
  }
  else
  {
#pragma omp parallel for
    for (int y = 0; y < mask.rows; ++y)
    {
      unsigned char* ptr_mask_new = reinterpret_cast<unsigned char*>(mask_texture_new.ptr(patch.anchor_source.y+y));
      const unsigned char* ptr_mask = reinterpret_cast<const unsigned char*>(mask.ptr(y));
      for (int x = 0; x < mask.cols; ++x)
      {
        if (ptr_mask[x])
        {
          ptr_mask_new[patch.anchor_source.x+x] = 0;
        }
      }
    }
  }

  for (size_t i = 0; i < m_textures[patch.source_index].size(); ++i)
  {
    cv::Mat transform_mask = AffineTransformation::concat(m_textures[patch.source_index][i].transformation_matrix, texture_source.transformation_matrix_inv);
    cv::Mat mask_rotated;
    cv::warpAffine(mask_texture_new, mask_rotated, transform_mask, m_textures[patch.source_index][i].mask_done.size(), cv::INTER_NEAREST, cv::BORDER_CONSTANT, cv::Scalar(1));
    cv::erode(mask_rotated, mask_rotated, cv::Mat::ones(3, 3, CV_8UC1));
    m_textures[patch.source_index][i].mask_done = cv::min(m_textures[patch.source_index][i].mask_done, mask_rotated);
  }
}

void TreeMatch::unmask_patch_resources(const Patch& patch)
{
  unmask_patch_resources(patch, patch.region_target.mask());
}

void TreeMatch::unmask_patch_resources(const AdaptivePatch& adaptive_patch)
{
  for (const Patch& p : adaptive_patch.patches)
  {
    unmask_patch_resources(p);
  }
}

void TreeMatch::unmask_patch_resources(const Patch& patch, const cv::Mat& mask)
{
  const Texture& texture_source = m_textures[patch.source_index][patch.source_rot];
  cv::Mat mask_texture_new = cv::Mat::ones(texture_source.mask_done.size(), CV_8UC1);

  if (mask.empty())
  {
    mask_texture_new(cv::Rect(patch.anchor_source, patch.size())) = 0;
  }
  else
  {
#pragma omp parallel for
    for (int y = 0; y < mask.rows; ++y)
    {
      unsigned char* ptr_mask_new = reinterpret_cast<unsigned char*>(mask_texture_new.ptr(patch.anchor_source.y+y));
      const unsigned char* ptr_mask = reinterpret_cast<const unsigned char*>(mask.ptr(y));
      for (int x = 0; x < mask.cols; ++x)
      {
        if (ptr_mask[x])
        {
          ptr_mask_new[patch.anchor_source.x+x] = 0;
        }
      }
    }
  }

  for (size_t i = 0; i < m_textures[patch.source_index].size(); ++i)
  {
    cv::Mat transform_mask = AffineTransformation::concat(m_textures[patch.source_index][i].transformation_matrix, texture_source.transformation_matrix_inv);
    cv::Mat mask_rotated;
    cv::warpAffine(mask_texture_new, mask_rotated, transform_mask, m_textures[patch.source_index][i].mask_done.size(), cv::INTER_NEAREST, cv::BORDER_CONSTANT, cv::Scalar(1));
    cv::erode(mask_rotated, mask_rotated, cv::Mat::ones(3, 3, CV_8UC1));
    mask_rotated = 1 - mask_rotated;
    m_textures[patch.source_index][i].mask_done = cv::max(m_textures[patch.source_index][i].mask_done, mask_rotated);
  }
}

void TreeMatch::add_patch(const Patch& patch)
{
  mask_patch_resources(patch);
  m_patches.push_back(patch);
}

cv::Mat TreeMatch::draw_matched_target(int target_index, double histogram_matching_factor) const
{
  FeatureEvaluator evaluator(1.0, 0.0, 0.0, m_filter_bank);
  FeatureVector response;

  if (histogram_matching_factor != 0.0)
  {
    std::vector<Texture> textures_histogram;
    for (const std::vector<Texture>& textures : m_textures)
    {
      textures_histogram.push_back(textures.front());
    }
    response = evaluator.evaluate_with_histogram_matching(m_targets[target_index].texture, textures_histogram, m_targets[target_index].mask_rotation, histogram_matching_factor);
  }
  else
  {
    response = evaluator.evaluate(m_targets[target_index].texture, m_targets[target_index].mask());
  }

  return response[0];
}

cv::Mat TreeMatch::draw_masked_target(int target_index) const
{
  cv::Mat image = m_targets[target_index].texture.clone();
  image.convertTo(image, CV_8UC3, 1.0/255.0);

  cv::Mat mask = m_targets[target_index].mask();
  cv::merge(std::vector<cv::Mat>(3, mask), mask);

  return cv::min(image, mask);
}

std::vector<cv::Mat> TreeMatch::draw_masked_textures() const
{
  std::vector<cv::Mat> textures;
  for (const std::vector<Texture>& textures_rot : m_textures)
  {
    cv::Mat texture;
    textures_rot[0].texture.convertTo(texture, CV_8UC3, 1.0/255.0);

    cv::Mat mask = textures_rot[0].mask();
    mask.setTo(255, mask > 0);
    cv::merge(std::vector<cv::Mat>(3, mask), mask);

    textures.push_back(cv::min(texture, mask));
  }
  return textures;
}

std::vector<cv::Mat> TreeMatch::draw_masked_textures_patch(const Patch& patch, cv::Scalar color, double alpha) const
{
  return draw_masked_textures_patch(std::vector<Patch>(1, patch), color, alpha);
}

std::vector<cv::Mat> TreeMatch::draw_masked_textures_patch(const std::vector<Patch>& patches, cv::Scalar color, double alpha) const
{
  return draw_masked_textures_patch(patches, std::vector<cv::Scalar>(patches.size(), color), std::vector<double>(patches.size(), alpha));
}

std::vector<cv::Mat> TreeMatch::draw_masked_textures_patch(const std::vector<Patch>& patches, const std::vector<cv::Scalar>& color, const std::vector<double>& alpha) const
{
  // Draw masked textures.
  std::vector<cv::Mat> textures;
  for (const std::vector<Texture>& textures_rot : m_textures)
  {
    cv::Mat texture;
    textures_rot[0].texture.convertTo(texture, CV_32FC3, 1.0/65535.0);
    textures.push_back(texture);
  }

  for (size_t i = 0; i < patches.size(); ++i)
  {
    // Generate mask from patch.
    const Texture& texture_patch = m_textures[patches[i].source_index][patches[i].source_rot];
    cv::Mat mask_patch = cv::Mat::zeros(texture_patch.mask_done.size(), CV_8UC1);
    mask_patch(cv::Rect(patches[i].anchor_source, patches[i].size())) = 255;

    // Compute transformation to unrotated texture.
    cv::Mat transform = AffineTransformation::concat(
      m_textures[patches[i].source_index][0].transformation_matrix,
      texture_patch.transformation_matrix_inv);

    // Warp mask onto unrotated texture.
    cv::Mat mask_rotated;
    cv::warpAffine(mask_patch, mask_rotated, transform, textures[patches[i].source_index].size(), cv::INTER_NEAREST, cv::BORDER_CONSTANT, cv::Scalar(0));
    cv::dilate(mask_rotated, mask_rotated, cv::Mat::ones(3, 3, CV_8UC1));

    cvx::multiplywithMask(textures[patches[i].source_index], cv::Scalar::all(1.0-alpha[i]), textures[patches[i].source_index], mask_rotated);
    cvx::addWithMask(textures[patches[i].source_index], alpha[i]*color[i], textures[patches[i].source_index], mask_rotated);
  }

  for (size_t i = 0; i < textures.size(); ++i)
  {
    textures[i].convertTo(textures[i], CV_8UC3, 255.0);
  }

  return textures;
}

static void mask_patch(cv::Mat mask, const Patch& p, cv::Mat transform)
{
  const float x1 = static_cast<float>(p.anchor_source.x);
  const float y1 = static_cast<float>(p.anchor_source.y);

  const float x2 = static_cast<float>(p.anchor_source.x + p.size().width);
  const float y2 = static_cast<float>(p.anchor_source.y + p.size().height);

  const std::vector<cv::Point2f> points = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
  const std::vector<cv::Point> points_transformed = AffineTransformation::transform<float, int>(transform, points);

  cv::fillConvexPoly(mask, points_transformed, 255);
}

std::vector<cv::Mat> TreeMatch::draw_masked_textures_patch_last(const std::vector<Patch>& patches, cv::Scalar color_1, double alpha_1, cv::Scalar color_2, double alpha_2, double scale) const
{
  std::vector<cv::Mat> textures, masks;
  for (const std::vector<Texture>& textures_rot : m_textures)
  {
    cv::Mat texture;
    textures_rot[0].texture.convertTo(texture, CV_32FC3, 1.0/65535.0);
    cv::resize(texture, texture, cv::Size(), scale, scale);
    textures.push_back(texture);

    masks.push_back(cv::Mat::zeros(texture.size(), CV_8UC1));
  }

  for (size_t i = 0; i < patches.size()-1; ++i)
  {
    const Patch& p = patches[i];
    const Texture& texture_patch = m_textures[p.source_index][p.source_rot];

    // Compute transformation to unrotated texture.
    cv::Mat transform = AffineTransformation::concat(
      m_textures[p.source_index][0].transformation_matrix,
      texture_patch.transformation_matrix_inv);
    transform = AffineTransformation::concat(AffineTransformation::T_scale(scale, scale), transform);

    mask_patch(masks[p.source_index], p, transform);
  }

  for (size_t i = 0; i < textures.size(); ++i)
  {
    cvx::multiplywithMask(textures[i], cv::Scalar::all(1.0-alpha_1), textures[i], masks[i]);
    cvx::addWithMask(textures[i], alpha_1*color_1, textures[i], masks[i]);

    cv::Mat mask_border;
    cv::erode(masks[i], mask_border, cv::Mat::ones(5, 5, CV_8UC1));
    mask_border = masks[i] - mask_border;
    textures[i].setTo(color_1, mask_border);
  }

  // Clear mask
  for (cv::Mat& mask : masks)
  {
    mask.setTo(0);
  }

  // Generate mask from patch.
  const Patch& p = patches.back();
  const Texture& texture_patch = m_textures[p.source_index][p.source_rot];

  // Compute transformation to unrotated texture.
  cv::Mat transform = AffineTransformation::concat(
    m_textures[p.source_index][0].transformation_matrix,
    texture_patch.transformation_matrix_inv);
  transform = AffineTransformation::concat(AffineTransformation::T_scale(scale, scale), transform);

  // Warp mask onto unrotated texture.
  mask_patch(masks[p.source_index], p, transform);

  for (size_t i = 0; i < textures.size(); ++i)
  {
    cvx::multiplywithMask(textures[i], cv::Scalar::all(1.0-alpha_2), textures[i], masks[i]);
    cvx::addWithMask(textures[i], alpha_2*color_2, textures[i], masks[i]);

    cv::Mat mask_border;
    cv::erode(masks[i], mask_border, cv::Mat::ones(5, 5, CV_8UC1));
    mask_border = masks[i] - mask_border;
    textures[i].setTo(color_2, mask_border);

    textures[i].convertTo(textures[i], CV_8UC3, 255.0);
  }

  return textures;
}

void TreeMatch::generate_patches_square(int target_index)
{
  cv::Size filter_kernel_size = FeatureEvaluator(1.0f, 1.0f, 1.0f, m_filter_bank).max_filter_size();
  std::vector<PatchRegion> patches = ::generate_patches_square(target_index, m_targets[target_index].texture.size(), m_patch_sizes.back(), filter_kernel_size);
  m_reconstruction_regions.insert(m_reconstruction_regions.end(), patches.begin(), patches.end());
}

void TreeMatch::generate_patches(int target_index, const Grid& morphed_grid, cv::Mat edge_image)
{
  std::vector<PatchRegion> patches;

  if (morphed_grid.empty() || edge_image.empty())
  {
    cv::Size filter_kernel_size = FeatureEvaluator(1.0f, 1.0f, 1.0f, m_filter_bank).max_filter_size();
    patches = ::generate_patches_square(target_index, m_targets[target_index].texture.size(), m_patch_sizes.back(), filter_kernel_size);
  }
  else
  {
    patches = ::generate_patches(target_index, morphed_grid, edge_image);
  }

  m_reconstruction_regions.insert(m_reconstruction_regions.end(), patches.begin(), patches.end());
}

void TreeMatch::add_patches(int target_index, const std::vector<PatchRegion>& patch_regions, double scale)
{
  m_reconstruction_regions.insert(m_reconstruction_regions.begin(), patch_regions.begin(), patch_regions.end());
  for (std::deque<PatchRegion>::iterator iter = m_reconstruction_regions.begin(); iter != m_reconstruction_regions.begin()+patch_regions.size(); ++iter)
  {
    iter->set_target_index(target_index);
    iter->scale(scale);
  }
}

void TreeMatch::downsample(int factor)
{
  if (factor == 1)
  {
    return;
  }

  for (Texture& target : m_targets)
  {
    target.downsample_nn(factor);
  }

  for (std::vector<Texture>& texture_vec : m_textures)
  {
    for (Texture& t : texture_vec)
    {
      t.downsample_nn(factor);
      t.transformation_matrix = Texture::compute_transformation_matrix(texture_vec[0].texture, t.angle_rad);
      cv::invertAffineTransform(t.transformation_matrix, t.transformation_matrix_inv);
    }
  }

  std::for_each(
    m_reconstruction_regions.begin(), m_reconstruction_regions.end(),
    std::bind(&PatchRegion::scale, std::placeholders::_1, 1.0/factor));
}

static int get_index(std::vector<std::string>& filenames, std::string filename)
{
  auto iter = std::find(filenames.begin(), filenames.end(), filename);
  if (iter != filenames.end())
  {
    return static_cast<int>(std::distance(filenames.begin(), iter));
  }
  else
  {
    filenames.push_back(filename);
    return static_cast<int>(filenames.size()) - 1;
  }
}

void TreeMatch::save(int target_index, boost::filesystem::path path) const
{
  path /= (boost::format("%04d") % target_index).str();

  if (!fs::exists(path))
  {
    fs::create_directories(path);
  }

  cv::Mat image = draw(target_index, true);
  cv::imwrite((path / "render.png").string(), image);
  m_targets[target_index].response.save(path / "responses" / "target");

  for (size_t i = 0; i < m_textures.size(); ++i)
  {
    m_textures[i][0].response.save(path / "responses" / (boost::format("texture_%04d") % i).str());
  }

  pt::ptree root;

  pt::ptree tree_patches;
  for (const Patch& patch : m_patches)
  {
    if (patch.region_target.target_index() == target_index)
    {
      tree_patches.push_back(std::make_pair("", patch.save(path, "patches")));
    }
  }
  root.add_child("patches", tree_patches);

  pt::ptree tree_textures;
  for (size_t i = 0; i < m_textures.size(); ++i)
  {
    tree_textures.push_back(std::make_pair("", m_textures[i][0].save(path, fs::path("source") / (boost::format("%04d") % i).str())));
  }
  root.add_child("textures_source", tree_textures);
  root.add_child("texture_target", m_targets[target_index].save(path, "target"));

  root.add<int>("subpatch_size", m_subpatch_size.width);

  pt::write_json((path / "result.json").string(), root);
}

void TreeMatch::find_markers(double marker_size_mm, int num_marker)
{
  for (std::vector<Texture>& texture : m_textures)
  {
    for (Texture& t : texture)
    {
      t.find_markers(marker_size_mm, num_marker);
    }
  }
}

cv::Size TreeMatch::max_filter_size(double weight_intensity, double weight_sobel, double weight_gabor) const
{
  cv::Size filter_size(1, 1);
  if (weight_sobel > 0.0)
  {
    filter_size = cv::Size(3, 3);
  }
  if (weight_gabor > 0.0)
  {
    cv::Size gabor_size = m_filter_bank.max_filter_size();
    filter_size.width = std::max(filter_size.width, gabor_size.width);
    filter_size.height = std::max(filter_size.height, gabor_size.height);
  }
  return filter_size;
}

void TreeMatch::sort_patches_by_saliency()
{
  struct SaliencySortData
  {
    PatchRegion region;
    double score;
  };

  std::vector<SaliencySortData> sort_data;
  cv::Ptr<cv::saliency::StaticSaliencyFineGrained> saliency_alg = cv::saliency::StaticSaliencyFineGrained::create();
  double max_saliency;

  std::vector<cv::Mat> saliency_maps;
  for (Texture& target : m_targets)
  {
    cv::Mat image, saliency_map;
    target.texture.convertTo(image, CV_8UC3, 1.0/255.0);
    saliency_alg->computeSaliency(image, saliency_map);
    saliency_maps.push_back(saliency_map);

    //cv::imshow("Saliency map", saliency_map);
    //cv::waitKey(0);
  }

  for (PatchRegion& region : m_reconstruction_regions)
  {
    cv::Mat local_saliency = saliency_maps[region.target_index()](region.bounding_box()).clone();
    local_saliency.setTo(0.0f, region.mask() == 0);
    cv::minMaxLoc(local_saliency, 0, &max_saliency);
    sort_data.push_back({region, max_saliency});
  }

  std::sort(sort_data.begin(), sort_data.end(), [](const SaliencySortData& lhs, const SaliencySortData& rhs){return lhs.score > rhs.score; });

  for (size_t i = 0; i < m_reconstruction_regions.size(); ++i)
  {
    m_reconstruction_regions[i] = sort_data[i].region;
  }
}

struct SortCenter
{
  SortCenter(const std::vector<cv::Point2d>& p_center) :
    p_center(p_center)
  {}

  bool operator()(const std::reference_wrapper<PatchRegion>& lhs, const std::reference_wrapper<PatchRegion>& rhs)
  {
    const double dist_1 = cv::norm(0.5 * cv::Point2d(lhs.get().bounding_box().tl() + lhs.get().bounding_box().br()) - p_center[lhs.get().target_index()]);
    const double dist_2 = cv::norm(0.5 * cv::Point2d(rhs.get().bounding_box().tl() + rhs.get().bounding_box().br()) - p_center[rhs.get().target_index()]);
    return dist_1 < dist_2;
  }

  std::vector<cv::Point2d> p_center;
};

void TreeMatch::sort_patches_by_center_distance()
{
  std::vector<cv::Point2d> p_center;
  for (const Texture& target : m_targets)
  {
    p_center.emplace_back(0.5 * target.texture.cols, 0.5 * target.texture.rows);
  }

  std::sort(m_reconstruction_regions.begin(), m_reconstruction_regions.end(), SortCenter(p_center));
}

std::pair<cv::Mat, cv::Mat> TreeMatch::draw_saliency(int target_index) const
{
  const double transparency = 0.6;
  const int num_regions = static_cast<int>(m_reconstruction_regions.size());

  cv::Mat color_map(1, 256, CV_8UC1);
  unsigned char* ptr_color_map_raw = color_map.ptr(0);
  std::iota(ptr_color_map_raw, ptr_color_map_raw+256, 0);
  std::reverse(ptr_color_map_raw, ptr_color_map_raw+256);
  cv::applyColorMap(color_map, color_map, cv::COLORMAP_JET);
  const cv::Vec3b* ptr_color_map = reinterpret_cast<const cv::Vec3b*>(color_map.ptr(0));

  cv::Mat image = m_target_images[target_index].clone();
  for (size_t i = 0; i < m_reconstruction_regions.size(); ++i)
  {
    const PatchRegion region = m_reconstruction_regions[i].scaled(1.0/m_targets[target_index].scale);

    const int color_map_index = static_cast<int>((static_cast<double>(i) / m_reconstruction_regions.size()) * 256);
    const cv::Vec3b overlay_color = ptr_color_map[color_map_index];
    cv::Mat mask = region.mask();

    for (int y = 0; y < region.size().height; ++y)
    {
      cv::Vec3b* image_ptr = reinterpret_cast<cv::Vec3b*>(image.ptr(region.bbox_tl().y + y));
      const unsigned char* mask_ptr = reinterpret_cast<const unsigned char*>(mask.ptr(y));

      for (int x = 0; x < region.size().width; ++x)
      {
        if (mask_ptr[x] > 0)
        {
          image_ptr[region.bbox_tl().x+x] = transparency * overlay_color + (1.0 - transparency) * image_ptr[region.bbox_tl().x+x];
        }
      }
    }    
  }

  return std::make_pair(image, color_map);
}

TreeMatch TreeMatch::load(const boost::filesystem::path& path, bool load_textures)
{
  pt::ptree root;
  pt::read_json(path.string(), root);

  fs::path path_json = path.parent_path();

  typedef struct
  {
    fs::path path;
    double scale;
    double dpi;
    std::vector<PatchRegion> reconstruction_regions;
  } target_json_t;

  typedef struct
  {
    std::string id;
    fs::path path_texture;
    fs::path path_mask;
    double scale;
    double dpi;
    TextureMarker markers;
  } texture_json_t;

  std::vector<target_json_t> targets_json;
  std::vector<texture_json_t> sources_json;

  int num_source_rotations;
  int min_patch_size;
  int patch_levels;
  double patch_quality_factor;

  double weight_intensity;
  double weight_sobel;
  double weight_gabor;

  double histogram_matching;
  int filter_resolution;
  int num_filter_directions;
  double filter_bandwidth_octaves;
  bool sort_patches_saliency;

  try
  {
    for (const auto& tree_source : root.get_child("source_textures"))
    {
      texture_json_t t;

      if (tree_source.second.count("json"))
      {
        pt::ptree tree_texture_json;
        fs::path path_texture_json = tree_source.second.get<fs::path>("json");

        if (!fs::exists(path_texture_json))
        {
          path_texture_json = path_json / path_texture_json;
        }

        pt::read_json(path_texture_json.string(), tree_texture_json);

        if (tree_texture_json.count("id") > 0)
        {
          t.id = tree_texture_json.get<std::string>("id");
        }
        else
        {
          t.id = "";
        }
        
        t.path_texture = tree_texture_json.get<fs::path>("texture");
        t.path_mask = tree_texture_json.get<fs::path>("mask");
        t.markers.load(path_texture_json.parent_path(), tree_texture_json);
        t.dpi = tree_texture_json.get<double>("dpi");
        t.scale = tree_source.second.get<double>("scale");
      
        if (!fs::exists(t.path_texture))
        {
          t.path_texture = path_texture_json.parent_path() / t.path_texture;
        }
        
        if (!fs::exists(t.path_mask))
        {
          t.path_mask = path_texture_json.parent_path() / t.path_mask;
        }
      }
      else
      {
        t.path_texture = tree_source.second.get<fs::path>("texture");
        
        if (tree_source.second.count("mask"))
        {
          t.path_mask = tree_source.second.get<fs::path>("mask");  
        }

        if (tree_source.second.count("markers_pixel"))
        {
          t.markers.load(path_json, tree_source.second);
        }

        if (tree_source.second.count("id"))
        {
          t.id = tree_source.second.get<std::string>("id");  
        }

        t.dpi = tree_source.second.get<double>("dpi");
        t.scale = tree_source.second.get<double>("scale");
      
        if (!fs::exists(t.path_texture))
        {
          t.path_texture = path_json / t.path_texture;
        }

        if (!fs::exists(t.path_mask))
        {
          t.path_mask = path_json / t.path_mask;
        }
      }

      sources_json.push_back(t);
    }

    num_source_rotations = root.get<int>("num_source_rotations");

    for (const auto& tree_target : root.get_child("targets"))
    {
      target_json_t t;
      t.path = tree_target.second.get<fs::path>("filename");
      if (!fs::exists(t.path))
      {
        t.path = path_json / t.path;
      }

      t.scale = tree_target.second.get<double>("scale");
      t.dpi = tree_target.second.get<double>("dpi");

      if (tree_target.second.count("morphed_grid"))
      {
        fs::path path_patches = tree_target.second.get<fs::path>("morphed_grid");
        if (!fs::exists(path_patches))
        {
          path_patches = path_json / path_patches;
        }

        if (!fs::exists(path_patches) || !fs::is_regular_file(path_patches))
        {
          std::cerr <<  "Could not find patches: " << path_patches << std::endl;
          std::exit(EXIT_FAILURE);
        }

        pt::ptree tree_patches;
        pt::read_json(path_patches.string(), tree_patches);
        for (const auto& tree_val : tree_patches.get_child("patches"))
        {
          PatchRegion p;
          p.load(path_patches.parent_path(), tree_val.second.get_child(""));
          t.reconstruction_regions.push_back(p);
        }
      }

      targets_json.push_back(t);
    }

    min_patch_size = root.get<int>("min_patch_size");
    patch_levels = root.get<int>("patch_levels");
    patch_quality_factor = root.get<double>("patch_quality_factor");

    weight_intensity = root.get<double>("weight_intensity");
    weight_sobel = root.get<double>("weight_sobel");
    weight_gabor = root.get<double>("weight_gabor");

    histogram_matching = root.get<double>("histogram_matching");
    filter_resolution = root.get<int>("filter_resolution");
    num_filter_directions = root.get<int>("num_filter_directions");
    filter_bandwidth_octaves = root.get<double>("filter_bandwidth_octaves");
    sort_patches_saliency = root.get<bool>("sort_patches_saliency");

    if (min_patch_size % 8 != 0)
    {
      std::cerr << "min_patch_size should be a multiple of 8." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    std::cout << "Settings:" << std::endl;
    for (const texture_json_t& t : sources_json)
    {
      std::cout << "source panel ID: " << t.id << std::endl
        << "  source_filename: " << t.path_texture << std::endl
        << "  source_filename_mask: " << t.path_mask << std::endl
        << "  source_scale: " << t.scale << std::endl
        << "  source_dpi: " << t.dpi << std::endl
        << "  num_markers: " << t.markers.markers_pix.size() << std::endl;
    }

    for (const target_json_t& t : targets_json)
    {
      std::cout << "target_filename: " << t.path << std::endl
        << "target_scale: " << t.scale << std::endl
        << "target_dpi: " << t.dpi << std::endl
        << "target_morphed_grid: " << (t.reconstruction_regions.empty() ? "no" : "yes") << std::endl;
    }

    std::cout << "min_patch_size: " << min_patch_size << std::endl
      << "patch_levels: " << patch_levels << std::endl
      << "patch_quality_factor: " << patch_quality_factor << std::endl
      << "weight_intensity: "  << weight_intensity << std::endl
      << "weight_sobel: " << weight_sobel << std::endl
      << "weight_gabor: " << weight_gabor << std::endl
      << "histogram_matching: " << histogram_matching << std::endl
      << "filter_resolution: " << filter_resolution << std::endl
      << "num_filter_directions: " << num_filter_directions << std::endl
      << "filter_bandwidth_octaves: " << filter_bandwidth_octaves << std::endl;

    if (root.count("downsample"))
    {
      std::cerr << "******************************************************" << std::endl
        << "WARNING: downsample is deprecated and will be ignored!" << std::endl
        << "******************************************************" << std::endl;
    }

  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TreeMatch matcher(min_patch_size, patch_levels, patch_quality_factor, filter_resolution, filter_bandwidth_octaves, num_filter_directions);

  for (const target_json_t& t : targets_json)
  {
    matcher.add_target(t.path, t.dpi, t.scale);
    if (t.reconstruction_regions.empty())
    {
      matcher.generate_patches_square(matcher.num_targets()-1);
    }
    else
    {
      matcher.add_patches(matcher.num_targets()-1, t.reconstruction_regions, t.scale);
    }
  }

  if (load_textures)
  {
    for (const texture_json_t& t : sources_json)
    {
      matcher.add_texture(t.path_texture, t.path_mask, t.dpi, t.scale, num_source_rotations, t.markers, t.id);
    }
    matcher.compute_responses(weight_intensity, weight_sobel, weight_gabor, histogram_matching);
  }

  if (sort_patches_saliency)
  {
    matcher.sort_patches_by_saliency();
  }
  else
  {
    matcher.sort_patches_by_center_distance();
  }

  return matcher;
}
