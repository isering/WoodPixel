/*
WoodPixel - Supplementary code for Computational Parquetry:
            Fabricated Style Transfer with Wood Pixels
            ACM Transactions on Graphics 39(2), 2020

Copyright (C) 2020  Julian Iseringhausen, University of Bonn, <iseringhausen@cs.uni-bonn.de>
Copyright (C) 2020  Matthias Hullin, University of Bonn, <hullin@cs.uni-bonn.de>

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

#include <random>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <opencv2/opencv.hpp>

#include "bezier_curve.hpp"
#include "cut_saver.hpp"
#include "mat.hpp"
#include "merge_patch.hpp"
#include "patch.hpp"
#include "texture.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

const int curve_degree = 3;

const double table_size_x_mm = 420.0;
const double table_size_y_mm = 297.0;


static cv::Vec3b hsv_to_bgr(int hue, int saturation, int value)
{
  cv::Mat hsv_mat(1, 1, CV_8UC3, cv::Vec3b(hue, saturation, value));
  cv::Mat rgb_mat;
  cv::cvtColor(hsv_mat, rgb_mat, cv::COLOR_HSV2BGR);
  return cv::Vec3b(rgb_mat.data[0], rgb_mat.data[1], rgb_mat.data[2]);
}

static cv::Mat draw_cut_matrix(const std::vector<MergePatch>& patches, cv::Size merged_size, int subpatch_size, double factor = 1.0)
{
  cv::Mat image = cv::Mat::zeros(merged_size, CV_8UC3);

  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, 180);

  for (const MergePatch& patch : patches)
  {
    cv::Vec3b color = hsv_to_bgr(distribution(generator), 128, 240);
    cv::Mat image_patch = image(cv::Rect(patch.region_subpatch.x * subpatch_size, patch.region_subpatch.y * subpatch_size, patch.region_subpatch.width * subpatch_size, patch.region_subpatch.height * subpatch_size));
    for (int y = 0; y < patch.active_pixel.rows; ++y)
    {
      const uint8_t* ptr_active = reinterpret_cast<const uint8_t*>(patch.active_pixel.ptr(y));
      cv::Vec3b* ptr_image = reinterpret_cast<cv::Vec3b*>(image_patch.ptr(y));
      for (int x = 0; x < patch.active_pixel.cols; ++x)
      {
        if (ptr_active[x])
        {
          ptr_image[x] = color;
        }
      }
    }
  }

  cv::resize(image, image, cv::Size(), factor, factor, cv::INTER_NEAREST);

  /*
  for (const MergePatch& patch : patches)
  {
    std::string str = (boost::format("(%02d,%02d)") % patch.pos_patch.x % patch.pos_patch.y).str();
    cv::putText(image, str, cv::Point(factor * patch.anchor_target.x, factor * patch.anchor_target.y), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0));
  }
  */

  return image;
}

static std::vector<double> smooth_cut_coordinates(const std::vector<int>& coordinates_int)
{
  std::vector<double> coordinates;
  coordinates.reserve(coordinates_int.size());

  if (coordinates_int.empty())
  {
    return coordinates;
  }

  coordinates.push_back(static_cast<double>(coordinates_int[0]));
  for (size_t i = 1; i < coordinates_int.size(); ++i)
  {
    if (coordinates_int[i - 1] == coordinates_int[i])
    {
      coordinates.push_back(static_cast<double>(coordinates_int[i]));
    }
    else if (coordinates_int[i - 1] > coordinates_int[i])
    {
      coordinates.push_back(coordinates_int[i] + 0.5);
    }
    else
    {
      coordinates.push_back(coordinates_int[i] - 0.5);
    }
  }

  return coordinates;
}

static int get_position_code_x(int x, const MergePatch& patch)
{
  if (x == patch.region_subpatch.x)
  {
    return 1;
  }
  else if (x == patch.region_subpatch.x + patch.region_subpatch.width - 1)
  {
    return 2;
  }
  return 4;
}

static int get_position_code_y(int y, const MergePatch& patch)
{
  if (y == patch.region_subpatch.y)
  {
    return 1;
  }
  else if (y == patch.region_subpatch.y + patch.region_subpatch.height - 1)
  {
    return 2;
  }
  return 4;
}

static bool is_y_cut(int x, const std::vector<MergePatch>& patches, const std::vector<int> index_vec)
{
  std::vector<int> position_codes_x(index_vec.size());
  std::transform(index_vec.begin(), index_vec.end(), position_codes_x.begin(), [x, &patches](int i) {return get_position_code_x(x, patches[i]); });
  return std::accumulate(position_codes_x.begin(), position_codes_x.end(), 0) == 7;
}

static bool is_x_cut(int y, const std::vector<MergePatch>& patches, const std::vector<int> index_vec)
{
  std::vector<int> position_codes_y(index_vec.size());
  std::transform(index_vec.begin(), index_vec.end(), position_codes_y.begin(), [y, &patches](int i) {return get_position_code_y(y, patches[i]); });
  return std::accumulate(position_codes_y.begin(), position_codes_y.end(), 0) == 7;
}

static int get_start(cv::Mat active_pixel)
{
  int index = 0;
  for (auto iter = active_pixel.begin<uint8_t>(); iter != active_pixel.end<uint8_t>(); ++iter, ++index)
  {
    if (*iter == 0)
    {
      break;
    }
  }
  return index;
}

static bool is_valid_y_cut(const std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat, cv::Mat is_finished_mat, int y, int x)
{
  const std::vector<int>& index_vec = subpatch_index_mat(y, x);
  const int num_patches = static_cast<int>(index_vec.size());
  return (
    x >= 0 &&
    x < subpatch_index_mat.width() &&
    !is_finished_mat.at<uint8_t>(y, x) &&
    (num_patches == 2 || (num_patches == 3 && is_y_cut(x, patches, index_vec))));
}

static bool is_valid_x_cut(const std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat, cv::Mat is_finished_mat, int y, int x)
{
  const std::vector<int>& index_vec = subpatch_index_mat(y, x);
  const int num_patches = static_cast<int>(index_vec.size());
  return (
    y >= 0 &&
    y < subpatch_index_mat.height() &&
    !is_finished_mat.at<uint8_t>(y, x) &&
    (num_patches == 2 || (num_patches == 3 && is_x_cut(y, patches, index_vec))));
}

typedef struct
{
  std::vector<MergePatch> merge_patches;
  std::vector<MergePatch> merge_patches_trimmed;
  std::vector<Texture> textures_source;
  Texture texture_target;
  cv::Size merged_size;
  int subpatch_size;
  int x_min, y_min;
  boost::filesystem::path base_path;
} MergeResult;

MergeResult compute_patches(const fs::path& path_in, int target_id, double w_reg, double w_slope)
{
  pt::ptree root;
  pt::read_json(path_in.string(), root);

  const fs::path base_path = path_in.parent_path();

  const int subpatch_size = root.get<int>("subpatch_size");

  std::vector<Patch> patches = Serializable::deserialize_vec<Patch>(root, "patches", base_path);
  std::vector<MergePatch> merge_patches = MergePatch::from_patches(patches, subpatch_size, target_id);

  std::vector<Texture> textures_source = Serializable::deserialize_vec<Texture>(root, "textures_source", base_path);

  Texture texture_target;
  texture_target.load(base_path, root.get_child("texture_target"));

  mat<std::vector<int>> subpatch_index_mat = MergePatch::get_subpatch_index_mat(merge_patches, subpatch_size);
  cv::Mat is_finished_mat = cv::Mat::zeros(subpatch_index_mat.height(), subpatch_index_mat.width(), CV_8UC1);

  cv::Size merged_size(subpatch_index_mat.width() * subpatch_size, subpatch_index_mat.height() * subpatch_size);

  /*
  * Get top left corner in pixels.
  */
  int x_min = std::numeric_limits<int>::max();
  int y_min = std::numeric_limits<int>::max();
  for (const Patch& p : patches)
  {
    x_min = std::min(x_min, p.anchor_target().x);
    y_min = std::min(y_min, p.anchor_target().y);
  }

  /*
  * Merge parts with three neighboring patches first.
  */
  for (int y = 0; y < subpatch_index_mat.height(); ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width(); ++x)
    {
      if (subpatch_index_mat(y, x).size() == 3)
      {
        /*
        * Decide whether this is a vertical or horizontal cut.
        */
        const std::vector<int>& index_vec = subpatch_index_mat(y, x);

        std::vector<int> position_codes_x(index_vec.size());
        std::transform(index_vec.begin(), index_vec.end(), position_codes_x.begin(), [x, &merge_patches](int i) {return get_position_code_x(x, merge_patches[i]); });
        const int code_sum_x = std::accumulate(position_codes_x.begin(), position_codes_x.end(), 0);

        std::vector<int> position_codes_y(index_vec.size());
        std::transform(index_vec.begin(), index_vec.end(), position_codes_y.begin(), [y, &merge_patches](int i) {return get_position_code_y(y, merge_patches[i]); });
        const int code_sum_y = std::accumulate(position_codes_y.begin(), position_codes_y.end(), 0);

        if ((code_sum_x == 7 && code_sum_y == 7) || (code_sum_x != 7 && code_sum_y != 7))
        {
          throw(std::runtime_error("3 neighbor merge failed."));
        }
        if (code_sum_x == 7)
        {
          // This is a cut in x-direction.
          const int index_left = position_codes_x[0] == 2 ? index_vec[0] : position_codes_x[1] == 2 ? index_vec[1] : index_vec[2];
          const int index_right = position_codes_x[0] == 1 ? index_vec[0] : position_codes_x[1] == 1 ? index_vec[1] : index_vec[2];
          MergePatch::merge_patches_x(merge_patches[index_left], merge_patches[index_right], y, x);
        }
        else
        {
          // This is a cut in y-direction.
          const int index_top = position_codes_y[0] == 2 ? index_vec[0] : position_codes_y[1] == 2 ? index_vec[1] : index_vec[2];
          const int index_bottom = position_codes_y[0] == 1 ? index_vec[0] : position_codes_y[1] == 1 ? index_vec[1] : index_vec[2];
          MergePatch::merge_patches_y(merge_patches[index_top], merge_patches[index_bottom], y, x);
        }
      }
    }
  }

  /*
  * Merge parts with two neighboring patches.
  */
  for (int y = 0; y < subpatch_index_mat.height() - 1; ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width() - 1; ++x)
    {
      if (!is_finished_mat.at<uint8_t>(y, x) && subpatch_index_mat(y, x).size() == 2)
      {
        /*
        * Check whether this is the beginning of a horizontal or vertical merge.
        */
        if (subpatch_index_mat(y, x + 1).size() == 2)
        {
          int x_start = x;
          int x_end = x;
          while (is_valid_y_cut(merge_patches, subpatch_index_mat, is_finished_mat, y, x_end))
          {
            is_finished_mat.at<uint8_t>(y, x_end) = 1;
            ++x_end;
          }

          int y_left = -1;
          if (x_start > 0 && subpatch_index_mat(y, x_start - 1).size() == 3)
          {
            if (is_x_cut(y, merge_patches, subpatch_index_mat(y, x_start - 1)))
            {
              for (int i : subpatch_index_mat(y, x_start - 1))
              {
                if (merge_patches[i].region_subpatch.y + merge_patches[i].region_subpatch.height - 1 == y)
                {
                  cv::Mat active = merge_patches[i].get_active_pixel(y, x_start - 1);
                  active = active.col(active.cols - 1);
                  y_left = get_start(active);
                }
              }
            }
          }

          int y_right = -1;
          if (x_end < subpatch_index_mat.width() && subpatch_index_mat(y, x_end).size() == 3)
          {
            if (is_x_cut(y, merge_patches, subpatch_index_mat(y, x_end)))
            {
              for (int i : subpatch_index_mat(y, x_end))
              {
                if (merge_patches[i].region_subpatch.y + merge_patches[i].region_subpatch.height - 1 == y)
                {
                  cv::Mat active = merge_patches[i].get_active_pixel(y, x_end);
                  active = active.col(0);
                  y_right = get_start(active);
                }
              }
            }
          }

          MergePatch::merge_patches_y(merge_patches, subpatch_index_mat, y, x_start, x_end, subpatch_size, y_left, y_right);
        }

        if (subpatch_index_mat(y + 1, x).size() == 2)
        {
          int y_start = y;
          int y_end = y;
          while (is_valid_x_cut(merge_patches, subpatch_index_mat, is_finished_mat, y_end, x))
          {
            is_finished_mat.at<uint8_t>(y_end, x) = 1;
            ++y_end;
          }

          int x_top = -1;
          if (y_start > 0 && subpatch_index_mat(y_start - 1, x).size() == 3)
          {
            if (is_y_cut(x, merge_patches, subpatch_index_mat(y_start - 1, x)))
            {
              for (int i : subpatch_index_mat(y_start - 1, x))
              {
                if (merge_patches[i].region_subpatch.x + merge_patches[i].region_subpatch.width - 1 == x)
                {
                  cv::Mat active = merge_patches[i].get_active_pixel(y_start - 1, x);
                  active = active.row(active.rows - 1);
                  x_top = get_start(active);
                }
              }
            }
          }

          int x_bottom = -1;
          if (y_end < subpatch_index_mat.height() && subpatch_index_mat(y_end, x).size() == 3)
          {
            if (is_y_cut(x, merge_patches, subpatch_index_mat(y_end, x)))
            {
              for (int i : subpatch_index_mat(y_end, x))
              {
                if (merge_patches[i].region_subpatch.x + merge_patches[i].region_subpatch.width - 1 == x)
                {
                  cv::Mat active = merge_patches[i].get_active_pixel(y_end, x);
                  active = active.row(0);
                  x_bottom = get_start(active);
                }
              }
            }
          }

          MergePatch::merge_patches_x(merge_patches, subpatch_index_mat, x, y_start, y_end, subpatch_size, x_top, x_bottom);
        }
      }
    }
  }

  /*
  * Fix isolated pixels in triple corners.
  */
  for (int y = 0; y < subpatch_index_mat.height(); ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width(); ++x)
    {
      const std::vector<int>& index_vec = subpatch_index_mat(y, x);
      std::vector<cv::Mat> active_mat_vec(index_vec.size());
      std::transform(index_vec.begin(), index_vec.end(), active_mat_vec.begin(), [&merge_patches, x, y](int i) {return merge_patches[i].get_active_pixel(y, x); });
      for (int i = 0; i < static_cast<int>(active_mat_vec.size()); ++i)
      {
        for (int y_pix = 0; y_pix < subpatch_size; ++y_pix)
        {
          for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
          {
            if (active_mat_vec[i].at<uint8_t>(y_pix, x_pix) == 1)
            {
              const bool is_connected = (
                (x_pix == 0 && x > merge_patches[i].region_subpatch.x) ||
                (x_pix == subpatch_size - 1 && x < merge_patches[i].region_subpatch.x + merge_patches[i].region_subpatch.width - 1) ||
                (y_pix == 0 && y > merge_patches[i].region_subpatch.y) ||
                (y_pix == subpatch_size - 1 && y < merge_patches[i].region_subpatch.y + merge_patches[i].region_subpatch.height - 1) ||
                (x_pix > 0 && active_mat_vec[i].at<uint8_t>(y_pix, x_pix - 1) == 1) ||
                (x_pix < subpatch_size - 1 && active_mat_vec[i].at<uint8_t>(y_pix, x_pix + 1) == 1) ||
                (y_pix > 0 && active_mat_vec[i].at<uint8_t>(y_pix - 1, x_pix) == 1) ||
                (y_pix < subpatch_size - 1 && active_mat_vec[i].at<uint8_t>(y_pix + 1, x_pix) == 1));
              if (!is_connected)
              {
                active_mat_vec[i].at<uint8_t>(y_pix, x_pix) = 0;
                active_mat_vec[(i + 1) % active_mat_vec.size()].at<uint8_t>(y_pix, x_pix) = 1;
              }
            }
          }
        }
      }
    }
  }

  /*
  * Merge cross sections (4 neighboring patches).
  */
  is_finished_mat = cv::Mat::zeros(subpatch_index_mat.height(), subpatch_index_mat.width(), CV_8UC1);
  for (int y = 0; y < subpatch_index_mat.height(); ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width(); ++x)
    {
      if (!is_finished_mat.at<uint8_t>(y, x) && subpatch_index_mat(y, x).size() == 4)
      {
        MergePatch::merge_patches_cross(merge_patches, subpatch_index_mat(y, x), subpatch_size);
      }
    }
  }

  /*
  * Remove temporary patches
  */
  merge_patches.erase(std::remove_if(merge_patches.begin(), merge_patches.end(), [](const MergePatch& p) {return p.pos_patch.x < 0; }), merge_patches.end());
  for (int y = 0; y < subpatch_index_mat.height(); ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width(); ++x)
    {
      std::vector<int>& indices = subpatch_index_mat(y, x);
      indices.erase(std::remove_if(indices.begin(), indices.end(), [&merge_patches](int i) {return i >= merge_patches.size(); }), indices.end());
    }
  }

  /*
  * Fit horizontal Bezier curves.
  */
  is_finished_mat.setTo(0);
  for (int y = 0; y < subpatch_index_mat.height(); ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width() - 1; ++x)
    {
      if (!is_finished_mat.at<uint8_t>(y, x) && subpatch_index_mat(y, x).size() > 1 && subpatch_index_mat(y, x + 1).size() > 1)
      {
        int x_cut = x;
        std::vector<int> cut_coordinates_int;
        while (x_cut < subpatch_index_mat.width() && subpatch_index_mat(y, x_cut).size() > 1)
        {
          // Find cut y coordinates.
          std::vector<int> cut_coordinates_local = MergePatch::get_cut_coordinates_horiz(merge_patches, subpatch_index_mat, y, x_cut, subpatch_size);
          cut_coordinates_int.insert(cut_coordinates_int.end(), cut_coordinates_local.begin(), cut_coordinates_local.end());
          is_finished_mat.at<uint8_t>(y, x_cut) = 1;
          ++x_cut;
        }

        // Smooth cut coordinates and fit Bezier curves.
        std::vector<double> cut_coordinates = smooth_cut_coordinates(cut_coordinates_int);
        cut_coordinates.push_back(cut_coordinates.back());
        std::vector<BezierCurve> curves = BezierCurve::fit(cut_coordinates, curve_degree, subpatch_size, w_reg, w_slope);

        // Sort curves back into patches.
        for (int i = 0; i < static_cast<int>(curves.size()); ++i)
        {
          for (int j : subpatch_index_mat(y, x + i))
          {
            const int position_code = get_position_code_y(y, merge_patches[j]);
            if (subpatch_index_mat(y, x + i).size() != 3 || position_code == 1 || position_code == 2)
            {
              merge_patches[j].add_curve_horiz(y, x + i, curves[i] + cv::Point2d(x*subpatch_size+x_min, y*subpatch_size + y_min));
            }
          }
        }
      }
    }
  }

  /*
  * Fit vertical Bezier curves.
  */
  is_finished_mat.setTo(0);
  for (int y = 0; y < subpatch_index_mat.height() - 1; ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width(); ++x)
    {
      if (!is_finished_mat.at<uint8_t>(y, x) && subpatch_index_mat(y, x).size() > 1 && subpatch_index_mat(y + 1, x).size() > 1)
      {
        int y_cut = y;
        std::vector<int> cut_coordinates_int;
        while (y_cut < subpatch_index_mat.height() && subpatch_index_mat(y_cut, x).size() > 1)
        {
          // Find cut x coordinates.
          std::vector<int> cut_coordinates_local = MergePatch::get_cut_coordinates_vert(merge_patches, subpatch_index_mat, y_cut, x, subpatch_size);
          cut_coordinates_int.insert(cut_coordinates_int.end(), cut_coordinates_local.begin(), cut_coordinates_local.end());
          is_finished_mat.at<uint8_t>(y_cut, x) = 1;
          ++y_cut;
        }

        // Smooth cut coordinates and fit Bezier curves.
        std::vector<double> cut_coordinates = smooth_cut_coordinates(cut_coordinates_int);
        cut_coordinates.push_back(cut_coordinates.back());
        std::vector<BezierCurve> curves = BezierCurve::fit(cut_coordinates, curve_degree, subpatch_size, w_reg, w_slope);
        std::for_each(curves.begin(), curves.end(), std::bind(&BezierCurve::swap_x_y, std::placeholders::_1));

        // Sort curves back into patches.
        for (int i = 0; i < static_cast<int>(curves.size()); ++i)
        {
          for (int j : subpatch_index_mat(y + i, x))
          {
            const int position_code = get_position_code_x(x, merge_patches[j]);
            if (subpatch_index_mat(y + i, x).size() != 3 || position_code == 1 || position_code == 2)
            {
              merge_patches[j].add_curve_vert(y + i, x, curves[i] + cv::Point2d(x*subpatch_size+x_min, y*subpatch_size+y_min));
            }
          }
        }
      }
    }
  }

 /*
  * Add image edge curves.
  */
  for (int y = 0; y < subpatch_index_mat.height(); ++y)
  {
    for (int x = 0; x < subpatch_index_mat.width(); ++x)
    {
      if (subpatch_index_mat(y, x).empty())
      {
        continue;
      }
      if (x == 0 || subpatch_index_mat(y, x - 1).empty())
      {
        for (int i : subpatch_index_mat(y, x))
        {
          cv::Point2f p_start(static_cast<float>(x_min + (x-0.1f)*subpatch_size), static_cast<float>(y_min + y*subpatch_size));
          cv::Point2f p_end(static_cast<float>(x_min + (x-0.1f)*subpatch_size), static_cast<float>(y_min + (y + 1)*subpatch_size));
          merge_patches[i].add_curve_vert(y, x, BezierCurve(p_start, p_end, 3));
        }
      }
      if (x == subpatch_index_mat.width() - 1 || subpatch_index_mat(y, x + 1).empty())
      {
        for (int i : subpatch_index_mat(y, x))
        {
          cv::Point2f p_start(static_cast<float>(x_min + (x + 1.1f)*subpatch_size), static_cast<float>(y_min + y*subpatch_size));
          cv::Point2f p_end(static_cast<float>(x_min + (x + 1.1f)*subpatch_size), static_cast<float>(y_min + (y + 1)*subpatch_size));
          merge_patches[i].add_curve_vert(y, x, BezierCurve(p_start, p_end, 3));
        }
      }
      if (y == 0 || subpatch_index_mat(y - 1, x).empty())
      {
        for (int i : subpatch_index_mat(y, x))
        {
          cv::Point2f p_start(static_cast<float>(x_min + x*subpatch_size), static_cast<float>(y_min + (y-0.1f)*subpatch_size));
          cv::Point2f p_end(static_cast<float>(x_min + (x + 1)*subpatch_size), static_cast<float>(y_min + (y-0.1f)*subpatch_size));
          merge_patches[i].add_curve_horiz(y, x, BezierCurve(p_start, p_end, 3));
        }
      }
      if (y == subpatch_index_mat.height() - 1 || subpatch_index_mat(y + 1, x).empty())
      {
        for (int i : subpatch_index_mat(y, x))
        {
          cv::Point2f p_start(static_cast<float>(x_min + x*subpatch_size), static_cast<float>(y_min + (y + 1.1f)*subpatch_size));
          cv::Point2f p_end(static_cast<float>(x_min + (x + 1)*subpatch_size), static_cast<float>(y_min + (y + 1.1f)*subpatch_size));
          merge_patches[i].add_curve_horiz(y, x, BezierCurve(p_start, p_end, 3));
        }
      }
    }
  }

  /*
  * Fix curve corners.
  */
  //MergePatch::fix_curve_corners(merge_patches, subpatch_index_mat);

  /*
  * Trim Bezier curves.
  */
  std::vector<MergePatch> merge_patches_trimmed = merge_patches;
  for (MergePatch& p : merge_patches_trimmed)
  {
    p.trim_bezier_curves();
  }

  return {merge_patches, merge_patches_trimmed, textures_source, texture_target, merged_size, subpatch_size, x_min, y_min, base_path};
}

int main(int argc, char* argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show this help message")
    ("in,i", po::value<std::vector<fs::path>>(), "Input JSON file")
    ("out,o", po::value<fs::path>(), "Output directory")
    ("scale,s", po::value<double>(), "Output scale")
    ("target_id,t", po::value<int>(), "Target image ID")
    ("w_reg", po::value<double>(), "Dynamic programming regularization penalty")
    ("w_slope", po::value<double>(), "Dynamic programming slope penalty")
    ("render_only", "Don't output cut patterns");

  std::vector<fs::path> paths_in;
  fs::path path_out;
  double scale_output = 4.0;
  int target_id = 0;
  double w_reg = 10.0;
  double w_slope = 1.0;
  bool render_only = false;

  /*
  std::vector<cv::Point2d> control_points_1({
    cv::Point2d(208.00001340046063, 424.00000000000000),
    cv::Point2d(208.00138918108127, 439.99999994085084),
    cv::Point2d(208.00276496170190, 455.99999988170174),
    cv::Point2d(208.00414074232253, 471.99999982255258)});

  std::vector<cv::Point2d> control_points_2({
    cv::Point2d(202.00000000000000, 424.59997558593750),
    cv::Point2d(204.00000000000000, 424.59997558593750),
    cv::Point2d(206.00000000000000, 424.59997558593750),
    cv::Point2d(208.00000000000000, 424.59997558593750)});

  std::vector<cv::Point2d> control_points_3({
    cv::Point2d(208.00000000000000, 424.59997558593750),
    cv::Point2d(210.00000000000000, 424.59997558593750),
    cv::Point2d(212.00000000000000, 424.59997558593750),
    cv::Point2d(214.00000000000000, 424.59997558593750)});

  

  BezierCurve left(control_points_1);
  BezierCurve bot_1(control_points_2);
  BezierCurve bot_2(control_points_3);

  double t1, t2;
  bool result = intersect_cubic_test(left, bot_2, t1, t2);
  std::cout << result << std::endl
    <<t1 << " " << t2 << std::endl;

  cv::Rect bbox_1 = left.bounding_box();
  cv::Rect bbox_2 = bot_2.bounding_box();
  cv::Rect bbox = bbox_1 | bbox_2;
  cv::Mat im = cv::Mat::zeros(bbox.size(), CV_8UC3);
  left.draw(im, cv::Vec3b(255, 255, 0), bbox.tl());
  bot_2.draw(im, cv::Vec3b(255, 0, 255), bbox.tl());
  cv::resize(im, im, cv::Size(), 4.0, 4.0, cv::INTER_NEAREST);
  cv::imshow("Im", im);
  cv::waitKey();




  std::cout << result << std::endl;
  return 0;
  */

  /*
  std::vector<cv::Point2d> control_points({
    cv::Point2d(100.0, 100.0),
    cv::Point2d(233.0, 25.0),
    cv::Point2d(250.0, 150.0),
    cv::Point2d(200.0, 200.0)});
  BezierCurve c(control_points);
  cv::Rect2d bbox = c.bounding_box_cubic();
  bbox.width += 1.0;
  bbox.height += 1.0;

  cv::Mat image = cv::Mat::zeros(300, 300, CV_8UC3);
  c.draw(image, cv::Vec3b(255, 255, 255));
  cv::rectangle(image, bbox, cv::Vec3b(255, 255, 0));
  for (int i = 0; i < 4; ++i)
  {
    cv::circle(image, c.control_point(i), 3, cv::Scalar(0, 255, 255), -1);
  }

  cv::imshow("Image", image);
  cv::waitKey();
  return 0;
  */

  try
  {
    po::variables_map vm;    
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
      std::cout << desc << std::endl;
      return 0;
    }

    if (vm.count("in"))
    {
      paths_in = vm["in"].as<std::vector<fs::path>>();
      for (const fs::path& p : paths_in)
      {
        if (!fs::exists(p) || !fs::is_regular_file(p))
        {
          std::cerr << "Input JSON file does not exist or is no regular file: " << p << std::endl;
          return -1;
        }
      }
    }
    else
    {
      std::cerr << "No input JSON file specified." << std::endl
        << desc << std::endl;
      return -1;
    }
    
    if (vm.count("out"))
    {
      path_out = vm["out"].as<fs::path>();
      if (!fs::exists(path_out))
      {
        fs::create_directories(path_out);
      }
      if (!fs::exists(path_out) || !fs::is_directory(path_out))
      {
        std::cerr << "Unable to create output directory." << std::endl
          << desc << std::endl;
        return -1;
      }
    }
    
    if (vm.count("scale"))
    {
      scale_output = vm["scale"].as<double>();
      if (scale_output <= 0.0)
      {
        std::cout << "Output scale must be greater zero." << std::endl
          << desc << std::endl;
        return -1;
      }
    }
    
    if (vm.count("target_id"))
    {
      target_id = vm["target_id"].as<int>();
    }
    
    if (vm.count("w_reg"))
    {
      w_reg = vm["w_reg"].as<double>();
    }
    
    if (vm.count("w_slope"))
    {
      w_slope = vm["w_slope"].as<double>();
    }

    if (vm.count("render_only"))
    {
      render_only = true;
    }
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl
      << desc << std::endl;
    return -1;
  }

  try
  {
    /*
     * Merge results.
     */
    std::vector<MergeResult> merge_results;
    for (size_t i = 0; i < paths_in.size(); ++i)
    {
      merge_results.push_back(compute_patches(paths_in[i], static_cast<int>(i), w_reg, w_slope));
    }

    std::vector<cv::Mat> images;
    for (const Texture& t : merge_results.front().textures_source)
    {
      cv::Mat image = t.texture.clone();
      image.convertTo(image, CV_8UC3, 1.0 / 255.0);
      images.push_back(image);
    }

    for (const MergeResult& r : merge_results)
    {
      MergePatch::draw_bezier_curves_source(images, r.merge_patches_trimmed, 1.0);
    }

    for (size_t i = 0; i < images.size(); ++i)
    {
      const fs::path path_source = path_out / (boost::format("source_image_bezier_%04d.png") % i).str();
      cv::imwrite(path_source.string(), images[i]);
    }

    /*
     * Write results to disk.
     */
    std::cout << "Write results to disk" << std::endl;
    /*
    std::vector<SVGSaver> saver_sources_svg;
    for (const MergeResult& r : merge_results)
    {
      CutSaver::add_sources_bezier_svg(saver_sources_svg, path_out, r.textures_source, r.merge_patches);
    }
    CutSaver::save<SVGSaver>(saver_sources_svg, path_out, "source_bezier");
    CutSaver::save_calibrated_front<SVGSaver>(saver_sources_svg, path_out, "source_bezier_calibrated_front", table_size_x_mm, table_size_y_mm);
    CutSaver::save_calibrated_back<SVGSaver>(saver_sources_svg, path_out, "source_bezier_calibrated_back", table_size_x_mm, table_size_y_mm);
    */

    if (!render_only)
    {
      //FIXME
      /*
      std::vector<SVGSaver> saver_sources_rect_svg;
      for (const MergeResult& r : merge_results)
      {
        CutSaver::add_sources_rect_svg(saver_sources_rect_svg, path_out, r.textures_source, r.merge_patches_trimmed);
      }
      CutSaver::save<SVGSaver>(saver_sources_rect_svg, path_out, "source_rect");

      std::vector<EPSSaver> saver_sources_eps;
      for (const MergeResult& r : merge_results)
      {
        CutSaver::add_sources_bezier_eps(saver_sources_eps, path_out, r.textures_source, r.merge_patches_trimmed);
      }
      CutSaver::save<EPSSaver>(saver_sources_eps, path_out, "source_bezier");

      std::vector<EPSSaver> saver_sources_rect_eps;
      for (const MergeResult& r : merge_results)
      {
        CutSaver::add_sources_rect_eps(saver_sources_rect_eps, path_out, r.textures_source, r.merge_patches_trimmed);
      }
      CutSaver::save<EPSSaver>(saver_sources_rect_eps, path_out, "source_rect");

      std::vector<SVGSaver> saver_target_svg;
      for (const MergeResult& r : merge_results)
      {
        CutSaver::add_target_bezier_svg(saver_target_svg, path_out, r.texture_target, r.merge_patches_trimmed);
      }
      CutSaver::save<SVGSaver>(saver_target_svg, path_out, "target_bezier");

      std::vector<SVGSaver> saver_target_rect_svg;
      for (const MergeResult& r : merge_results)
      {
        CutSaver::add_target_rect_svg(saver_target_rect_svg, path_out, r.texture_target, r.merge_patches_trimmed);
      }
      CutSaver::save<SVGSaver>(saver_target_rect_svg, path_out, "target_rect");

      std::vector<EPSSaver> saver_target_eps;
      for (const MergeResult& r : merge_results)
      {
        CutSaver::add_target_bezier_eps(saver_target_eps, path_out, r.texture_target, r.merge_patches_trimmed);
      }
      CutSaver::save<EPSSaver>(saver_target_eps, path_out, "target_bezier");

      std::vector<EPSSaver> saver_target_rect_eps;
      for (const MergeResult& r : merge_results)
      {
        CutSaver::add_target_rect_eps(saver_target_rect_eps, path_out, r.texture_target, r.merge_patches_trimmed);
      }
      CutSaver::save<EPSSaver>(saver_target_rect_eps, path_out, "target_rect");
      */
    }

    for (size_t i = 0; i < merge_results.size(); ++i)
    {
      const MergeResult& r = merge_results[i];

      cv::Mat cut_image = draw_cut_matrix(r.merge_patches_trimmed, r.merged_size, r.subpatch_size, scale_output);
      cv::imwrite((path_out / (boost::format("cut_image_%04d.png") % i).str()).string(), cut_image);

      cv::Mat curve_image = cut_image.clone();
      MergePatch::draw_bezier_curves(curve_image, r.merge_patches_trimmed, scale_output, cv::Point2d(static_cast<double>(-r.x_min), static_cast<double>(-r.y_min)));
      cv::imwrite((path_out / (boost::format("curve_image_%04d.png") % i).str()).string(), curve_image);

      cv::Mat image = MergePatch::draw_fullres(r.merge_patches_trimmed, r.base_path, r.textures_source, scale_output, false);
      cv::imwrite((path_out / (boost::format("image_%04d.jpg") % i).str()).string(), image);

      cv::Mat image_boundary = MergePatch::draw_fullres(r.merge_patches_trimmed, r.base_path, r.textures_source, scale_output, true);
      cv::imwrite((path_out / (boost::format("image_boundary_%04d.jpg") % i).str()).string(), image_boundary);

      cv::Mat image_rect = MergePatch::draw_rect_fullres(r.merge_patches_trimmed, r.base_path, r.textures_source, scale_output, r.subpatch_size, false);
      cv::imwrite((path_out / (boost::format("image_rect_%04d.jpg") % i).str()).string(), image_rect);

      cv::Mat image_rect_boundary = MergePatch::draw_rect_fullres(r.merge_patches_trimmed, r.base_path, r.textures_source, scale_output, r.subpatch_size, true);
      cv::imwrite((path_out / (boost::format("image_rect_boundary_%04d.jpg") % i).str()).string(), image_rect_boundary);

      cv::Mat image_half_rect_boundary = MergePatch::draw_half_rect_fullres(r.merge_patches, r.base_path, r.textures_source, scale_output, r.subpatch_size, true);
      cv::imwrite((path_out / (boost::format("image_half_rect_boundary_%04d.jpg") % i).str()).string(), image_half_rect_boundary);

      cv::Mat image_baseline;
      cv::Rect region_baseline(r.subpatch_size, r.subpatch_size, image_rect.cols - 2 * r.subpatch_size, image_rect.rows - 2 * r.subpatch_size);
      cv::Mat image_center = image_rect(region_baseline);
      cv::resize(image_center, image_baseline, image_center.size() / static_cast<int>((3 * r.subpatch_size * scale_output)), 0.0, 0.0, cv::INTER_LINEAR);
      cv::resize(image_baseline, image_baseline, image_center.size(), 0.0, 0.0, cv::INTER_NEAREST);
      cv::imwrite((path_out / (boost::format("image_baseline_%04d.jpg") % i).str()).string(), image_baseline);
    }
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  return 0;
}
