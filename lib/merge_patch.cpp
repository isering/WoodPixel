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

#include "affine_transformation.hpp"
#include "eps_saver.hpp"
#include "merge_patch.hpp"
#include "svg_saver.hpp"

const std::string MergePatch::output_characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+-*#$%!<=>?~^";

std::vector<MergePatch> MergePatch::from_patches(std::vector<Patch>& patches_in, int subpatch_size, int target_id)
{
  /*
   * Get offset.
   */
  int x_min = std::numeric_limits<int>::max();
  int y_min = std::numeric_limits<int>::max();
  for (const Patch& p : patches_in)
  {
    x_min = std::min(x_min, p.anchor_target().x);
    y_min = std::min(y_min, p.anchor_target().y);
  }

  /*
   * Find patch position in global grid.
   */
  std::set<int> x_coords, y_coords;
  for (const Patch& patch : patches_in)
  {
    x_coords.insert(patch.anchor_target().x);
    y_coords.insert(patch.anchor_target().y);
  }

  /*
   * Generate merge patches.
   */
  std::vector<MergePatch> patches;
  for (const Patch& p : patches_in)
  {
    cv::Rect region_subpatch;
    region_subpatch.x = (p.anchor_target().x - x_min) / subpatch_size;
    region_subpatch.y = (p.anchor_target().y - y_min) / subpatch_size;
    region_subpatch.width = p.size().width / subpatch_size;
    region_subpatch.height = p.size().height / subpatch_size;

    int pos_patch_x = static_cast<int>(std::distance(x_coords.begin(), std::find(x_coords.begin(), x_coords.end(), p.anchor_target().x)));
    int pos_patch_y = static_cast<int>(std::distance(y_coords.begin(), std::find(y_coords.begin(), y_coords.end(), p.anchor_target().y)));

    patches.emplace_back(p, region_subpatch, cv::Point(pos_patch_x, pos_patch_y), subpatch_size, target_id);
  }

  return patches;
}

mat<std::vector<int>> MergePatch::get_subpatch_index_mat(std::vector<MergePatch>& patches, int subpatch_size)
{
  /*
   * Get size of index matrix and create map from coordinates to index.
   */
  int rows = 0;
  int cols = 0;
  int anchor_x_min = std::numeric_limits<int>::max();
  int anchor_y_min = std::numeric_limits<int>::max();
  for (const MergePatch& p : patches)
  {
    rows = std::max(rows, p.region_subpatch.y + p.region_subpatch.height);
    cols = std::max(cols, p.region_subpatch.x + p.region_subpatch.width);
    anchor_x_min = std::min(anchor_x_min, p.anchor_target.x);
    anchor_y_min = std::min(anchor_y_min, p.anchor_target.y);
  }

  cv::Mat subpatch_mat_filled = cv::Mat::zeros(rows, cols, CV_8UC1);

  mat<std::vector<int>> index_mat(rows, cols);
  for (size_t i = 0; i < patches.size(); ++i)
  {
    for (int y = 0; y < patches[i].region_subpatch.height; ++y)
    {
      for (int x = 0; x < patches[i].region_subpatch.width; ++x)
      {
        index_mat(patches[i].region_subpatch.y + y, patches[i].region_subpatch.x + x).push_back(static_cast<int>(i));
        subpatch_mat_filled.at<unsigned char>(patches[i].region_subpatch.y + y, patches[i].region_subpatch.x + x) = 255;
      }
    }
  }

  for (int y = 0; y < rows-1; y += 3)
  {
    for (int x = 0; x < cols-1; x += 3)
    {
      if (subpatch_mat_filled.at<unsigned char>(y+1, x+1) == 0)
      {
        const int patch_index = static_cast<int>(patches.size());
        patches.emplace_back(cv::Rect(x, y, 4, 4), cv::Point(anchor_x_min+x*subpatch_size, anchor_y_min+y*subpatch_size), subpatch_size);

        for (int y_delta = 0; y_delta < 4; ++y_delta)
        {
          for (int x_delta = 0; x_delta < 4; ++x_delta)
          {
            index_mat(y+y_delta, x+x_delta).push_back(patch_index);
          }
        }
      }
    }
  }

  return index_mat;
}

enum directions
{
  left = 0,
  top = 0,
  right = 1,
  bottom = 1,
  ignore = 254,
  uninitialized = 255
};

void MergePatch::merge_patches_x(MergePatch& patch_left, cv::Rect region_left, MergePatch& patch_right, cv::Rect region_right, int start_top, int start_bottom)
{
  if (region_left.size() != region_right.size())
  {
    throw(std::invalid_argument("Region sizes differ."));
  }

  cv::Mat active = merge_patches_x(patch_left.error(region_left), patch_right.error(region_right), start_top, start_bottom);
  cv::Mat active_left = patch_left.active_pixel(region_left);
  cv::Mat active_right = patch_right.active_pixel(region_right);

  for (int y = 0; y < active.rows; ++y)
  {
    const uint8_t* ptr = reinterpret_cast<const uint8_t*>(active.ptr(y));
    uint8_t* ptr_left = reinterpret_cast<uint8_t*>(active_left.ptr(y));
    uint8_t* ptr_right = reinterpret_cast<uint8_t*>(active_right.ptr(y));
    for (int x = 0; x < active.cols; ++x)
    {
      ptr_left[x] = ptr[x] == left ? 1 : 0;
      ptr_right[x] = 1 - ptr_left[x];
    }
  }
}

cv::Mat MergePatch::merge_patches_x(cv::Mat error_left, cv::Mat error_right, int start_top, int start_bottom)
{
  if (error_left.size() != error_right.size())
  {
    throw(std::invalid_argument("Region sizes differ."));
  }

  const int h = error_left.rows;
  const int w = error_left.cols;

  cv::Mat error = cv::Mat::zeros(h, w + 1, CV_32FC1);

  for (int y = 0; y < h; ++y)
  {
    const float* ptr_error_top = reinterpret_cast<const float*>(error_left.ptr(y));
    const float* ptr_error_bottom = reinterpret_cast<const float*>(error_right.ptr(y));
    float* ptr_error = reinterpret_cast<float*>(error.ptr(y));

    for (int x = 0; x < w + 1; ++x)
    {
      for (int x_left = 0; x_left < x; ++x_left)
      {
        ptr_error[x] += ptr_error_top[x_left];
      }
      for (int x_right = x; x_right < w; ++x_right)
      {
        ptr_error[x] += ptr_error_bottom[x_right];
      }
    }
  }

  if (start_top >= 0)
  {
    float* ptr_error = reinterpret_cast<float*>(error.ptr(0));
    for (int x = 0; x < start_top - 1; ++x)
    {
      ptr_error[x] = std::numeric_limits<float>::infinity();
    }
    for (int x = start_top + 2; x < w + 1; ++x)
    {
      ptr_error[x] = std::numeric_limits<float>::infinity();
    }
    //ptr_error[0] = std::numeric_limits<float>::infinity();
    //ptr_error[w] = std::numeric_limits<float>::infinity();
  }

  if (start_bottom >= 0)
  {
    float* ptr_error = reinterpret_cast<float*>(error.ptr(h - 1));
    for (int x = 0; x < start_bottom - 1; ++x)
    {
      ptr_error[x] = std::numeric_limits<float>::infinity();
    }
    for (int x = start_bottom + 2; x < w + 1; ++x)
    {
      ptr_error[x] = std::numeric_limits<float>::infinity();
    }
    //ptr_error[0] = std::numeric_limits<float>::infinity();
    //ptr_error[w] = std::numeric_limits<float>::infinity();
  }

  // Accumulate errors / forward pass
  for (int y = 1; y < h; ++y)
  {
    const float* ptr_last = reinterpret_cast<const float*>(error.ptr(y - 1));
    float* ptr_error = reinterpret_cast<float*>(error.ptr(y));

    ptr_error[0] += std::min(ptr_last[0], ptr_last[1]);
    ptr_error[w] += std::min(ptr_last[w - 1], ptr_last[w]);

    for (int x = 1; x < w; ++x)
    {
      ptr_error[x] += std::min({ ptr_last[x - 1], ptr_last[x], ptr_last[x + 1] });
    }
  }

  // Backward pass
  std::vector<int> x_vals;
  {
    const float* ptr_error = reinterpret_cast<const float*>(error.ptr(h - 1));
    x_vals.push_back(static_cast<int>(std::distance(ptr_error, std::min_element(ptr_error, ptr_error + w + 1))));
  }

  for (int y = h - 2; y >= 0; --y)
  {
    const float* ptr_error = reinterpret_cast<const float*>(error.ptr(y));
    int x_min = x_vals.back() == 0 ? 0 : x_vals.back() - 1;
    int x_max = x_vals.back() == w ? w + 1 : x_vals.back() + 2;
    x_vals.push_back(static_cast<int>(std::distance(ptr_error, std::min_element(ptr_error + x_min, ptr_error + x_max))));
  }

  std::reverse(x_vals.begin(), x_vals.end());

  cv::Mat active(h, w, CV_8UC1);
  for (int y = 0; y < h; ++y)
  {
    uint8_t* ptr = reinterpret_cast<uint8_t*>(active.ptr(y));
    for (int x = 0; x < x_vals[y]; ++x)
    {
      ptr[x] = left;
    }
    for (int x = x_vals[y]; x < w; ++x)
    {
      ptr[x] = right;
    }
  }

  return active;
}

void MergePatch::merge_patches_y(MergePatch& patch_top, cv::Rect region_top, MergePatch& patch_bottom, cv::Rect region_bottom, int start_left, int start_right)
{
  if (region_top.size() != region_bottom.size())
  {
    throw(std::invalid_argument("Region sizes differ."));
  }

  cv::Mat active = merge_patches_y(patch_top.error(region_top), patch_bottom.error(region_bottom), start_left, start_right);
  cv::Mat active_top = patch_top.active_pixel(region_top);
  cv::Mat active_bottom = patch_bottom.active_pixel(region_bottom);

  for (int y = 0; y < active.rows; ++y)
  {
    const uint8_t* ptr = reinterpret_cast<const uint8_t*>(active.ptr(y));
    uint8_t* ptr_top = reinterpret_cast<uint8_t*>(active_top.ptr(y));
    uint8_t* ptr_bottom = reinterpret_cast<uint8_t*>(active_bottom.ptr(y));

    for (int x = 0; x < active.cols; ++x)
    {
      ptr_top[x] = ptr[x] == top ? 1 : 0;
      ptr_bottom[x] = 1 - ptr_top[x];
    }
  }
}

cv::Mat MergePatch::merge_patches_y(cv::Mat error_top, cv::Mat error_bottom, int start_left, int start_right)
{
  if (error_top.size() != error_bottom.size())
  {
    throw(std::invalid_argument("Region sizes differ."));
  }

  const int h = error_top.rows;
  const int w = error_bottom.cols;

  cv::Mat error = cv::Mat::zeros(h+1, w, CV_32FC1);

  for (int x = 0; x < w; ++x)
  {
    for (int y = 0; y < h + 1; ++y)
    {
      float& val = error.at<float>(y, x);
      for (int y_left = 0; y_left < y; ++y_left)
      {
        val += error_top.at<float>(y_left, x);
      }
      for (int y_right = y; y_right < h; ++y_right)
      {
        val += error_bottom.at<float>(y_right, x);
      }
    }
  }

  if (start_left >= 0)
  {
    for (int y = 0; y < start_left - 1; ++y)
    {
      error.at<float>(y, 0) = std::numeric_limits<float>::infinity();
    }
    for (int y = start_left + 2; y < h + 1; ++y)
    {
      error.at<float>(y, 0) = std::numeric_limits<float>::infinity();
    }
  }

  if (start_right >= 0)
  {
    for (int y = 0; y < start_right - 1; ++y)
    {
      error.at<float>(y, w - 1) = std::numeric_limits<float>::infinity();
    }
    for (int y = start_right + 2; y < h + 1; ++y)
    {
      error.at<float>(y, w - 1) = std::numeric_limits<float>::infinity();
    }
  }

  // Accumulate errors / forward pass
  for (int x = 1; x < w; ++x)
  {
    error.at<float>(0, x) += std::min(error.at<float>(0, x - 1), error.at<float>(1, x - 1));
    error.at<float>(h, x) += std::min(error.at<float>(h - 1, x - 1), error.at<float>(h, x - 1));
    for (int y = 1; y < h; ++y)
    {
      error.at<float>(y, x) += std::min({ error.at<float>(y - 1, x - 1), error.at<float>(y, x - 1), error.at<float>(y + 1, x - 1) });
    }
  }

  // Backward pass
  std::vector<int> y_vals;
  {
    float min_val = std::numeric_limits<float>::infinity();
    int min_index = -1;
    for (int y = 0; y < h + 1; ++y)
    {
      if (error.at<float>(y, w - 1) <= min_val)
      {
        min_val = error.at<float>(y, w - 1);
        min_index = y;
      }
    }
    y_vals.push_back(min_index);
  }

  for (int x = w - 2; x >= 0; --x)
  {
    int val_old = y_vals.back();

    int y_min = y_vals.back() == 0 ? 0 : y_vals.back() - 1;
    int y_max = y_vals.back() == h ? h + 1 : y_vals.back() + 2;

    float min_val = std::numeric_limits<float>::infinity();
    int min_index = -1;
    for (int y = y_min; y < y_max; ++y)
    {
      if (error.at<float>(y, x) <= min_val)
      {
        min_val = error.at<float>(y, x);
        min_index = y;
      }
    }

    y_vals.push_back(min_index);
  }

  std::reverse(y_vals.begin(), y_vals.end());

  cv::Mat active(h, w, CV_8UC1);
  for (int y = 0; y < h; ++y)
  {
    uint8_t* ptr = reinterpret_cast<uint8_t*>(active.ptr(y));
    
    for (int x = 0; x < w; ++x)
    {
      ptr[x] = y < y_vals[x] ? top : bottom;
    }
  }

  return active;
}

int get_start(const cv::Mat& active_pixel)
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

void MergePatch::merge_patches_cross(MergePatch& patch_tl, cv::Rect region_tl, MergePatch& patch_bl, cv::Rect region_bl, MergePatch& patch_tr, cv::Rect region_tr, MergePatch& patch_br, cv::Rect region_br)
{
  if (region_tl.size() != region_bl.size() || region_bl.size() != region_tr.size() || region_tr.size() != region_br.size())
  {
    throw(std::invalid_argument("Bounding box sizes differ."));
  }

  if (region_tl.width != region_tl.height)
  {
    throw(std::invalid_argument("Only quadratic bounding boxes supported."));
  }

  int num_valid = (patch_tl.pos_patch.x >= 0 ? 1 : 0) +
                  (patch_bl.pos_patch.x >= 0 ? 1 : 0) +
                  (patch_tr.pos_patch.x >= 0 ? 1 : 0) +
                  (patch_br.pos_patch.x >= 0 ? 1 : 0);

  if (num_valid == 1)
  {
    if (patch_tl.pos_patch.x < 0)
    {
      patch_tl.active_pixel(region_tl) = 0;
    }
    if (patch_bl.pos_patch.x < 0)
    {
      patch_bl.active_pixel(region_bl) = 0;
    }
    if (patch_tr.pos_patch.x < 0)
    {
      patch_tr.active_pixel(region_tr) = 0;
    }
    if (patch_br.pos_patch.x < 0)
    {
      patch_br.active_pixel(region_br) = 0;
    }
    return;
  }

  const int h = region_tl.height;
  const int w = region_tl.width;

  cv::Mat error[2][2];
  error[0][0] = patch_tl.error(region_tl);
  error[0][1] = patch_tr.error(region_tr);
  error[1][0] = patch_bl.error(region_bl);
  error[1][1] = patch_br.error(region_br);

  cv::Mat active[2][2];
  active[0][0] = patch_tl.active_pixel(region_tl);
  active[0][1] = patch_tr.active_pixel(region_tr);
  active[1][0] = patch_bl.active_pixel(region_bl);
  active[1][1] = patch_br.active_pixel(region_br);

  int start_left = get_start(patch_tl.active_pixel(cv::Rect(region_tl.x-1, region_tl.y, 1, region_tl.height)));
  int start_top = get_start(patch_tl.active_pixel(cv::Rect(region_tl.x, region_tl.y - 1, region_tl.width, 1)));
  int start_right = get_start(patch_tr.active_pixel(cv::Rect(region_tr.x + region_tr.width, region_tr.y, 1, region_tr.height)));
  int start_bottom = get_start(patch_bl.active_pixel(cv::Rect(region_bl.x, region_bl.y + region_bl.height, region_bl.width, 1)));

  cv::Mat active_mat(region_tl.size(), CV_8UC2);
  for (int y = 0; y < active_mat.rows; ++y)
  {
    cv::Vec2b* ptr = reinterpret_cast<cv::Vec2b*>(active_mat.ptr(y));
    for (int x = 0; x < active_mat.cols; ++x)
    {
      ptr[x][0] = y < active_mat.rows / 2 ? top : bottom;
      ptr[x][1] = x < active_mat.cols / 2 ? left : right;
    }
  }

  cv::Mat error_1(h, w, CV_32FC1);
  cv::Mat error_2(h, w, CV_32FC1);

  for (int i = 0; i < 3; ++i)
  {
    /* 
     * Merge vertically neighboring patches.
     */
    for (int y = 0; y < h; ++y)
    {
      const cv::Vec2b* ptr_active = reinterpret_cast<const cv::Vec2b*>(active_mat.ptr(y));
      const float* ptr_error_tl = reinterpret_cast<const float*>(error[0][0].ptr(y));
      const float* ptr_error_tr = reinterpret_cast<const float*>(error[0][1].ptr(y));
      const float* ptr_error_bl = reinterpret_cast<const float*>(error[1][0].ptr(y));
      const float* ptr_error_br = reinterpret_cast<const float*>(error[1][1].ptr(y));
      float* ptr_error_1 = reinterpret_cast<float*>(error_1.ptr(y));
      float* ptr_error_2 = reinterpret_cast<float*>(error_2.ptr(y));
      for (int x = 0; x < w; ++x)
      {
        ptr_error_1[x] = ptr_active[x][1] == left ? ptr_error_tl[x] : ptr_error_tr[x];
        ptr_error_2[x] = ptr_active[x][1] == left ? ptr_error_bl[x] : ptr_error_br[x];
      }
    }

    cv::Mat active_y = merge_patches_y(error_1, error_2, start_left, start_right);

    for (int y = 0; y < h; ++y)
    {
      cv::Vec2b* ptr_active = reinterpret_cast<cv::Vec2b*>(active_mat.ptr(y));
      const uint8_t* ptr_active_y = reinterpret_cast<const uint8_t*>(active_y.ptr(y));
      for (int x = 0; x < w; ++x)
      {
        ptr_active[x][0] = ptr_active_y[x];
      }
    }

    /*
     * Merge horizontally neighboring patches.
     */
    for (int y = 0; y < h; ++y)
    {
      const cv::Vec2b* ptr_active = reinterpret_cast<const cv::Vec2b*>(active_mat.ptr(y));
      const float* ptr_error_tl = reinterpret_cast<const float*>(error[0][0].ptr(y));
      const float* ptr_error_tr = reinterpret_cast<const float*>(error[0][1].ptr(y));
      const float* ptr_error_bl = reinterpret_cast<const float*>(error[1][0].ptr(y));
      const float* ptr_error_br = reinterpret_cast<const float*>(error[1][1].ptr(y));
      float* ptr_error_1 = reinterpret_cast<float*>(error_1.ptr(y));
      float* ptr_error_2 = reinterpret_cast<float*>(error_2.ptr(y));
      for (int x = 0; x < w; ++x)
      {
        ptr_error_1[x] = ptr_active[x][0] == top ? ptr_error_tl[x] : ptr_error_bl[x];
        ptr_error_2[x] = ptr_active[x][0] == top ? ptr_error_tr[x] : ptr_error_br[x];
      }
    }

    cv::Mat active_x = merge_patches_x(error_1, error_2, start_top, start_bottom);

    for (int y = 0; y < h; ++y)
    {
      cv::Vec2b* ptr_active = reinterpret_cast<cv::Vec2b*>(active_mat.ptr(y));
      const uint8_t* ptr_active_x = reinterpret_cast<const uint8_t*>(active_x.ptr(y));
      for (int x = 0; x < w; ++x)
      {
        ptr_active[x][1] = ptr_active_x[x];
      }
    }
  }

  active[0][0].setTo(0);
  active[0][1].setTo(0);
  active[1][0].setTo(0);
  active[1][1].setTo(0);
  
  for (int y = 0; y < h; ++y)
  {
    const cv::Vec2b* ptr_active = reinterpret_cast<const cv::Vec2b*>(active_mat.ptr(y));
    uint8_t* ptr_active_tl = reinterpret_cast<uint8_t*>(active[0][0].ptr(y));
    uint8_t* ptr_active_tr = reinterpret_cast<uint8_t*>(active[0][1].ptr(y));
    uint8_t* ptr_active_bl = reinterpret_cast<uint8_t*>(active[1][0].ptr(y));
    uint8_t* ptr_active_br = reinterpret_cast<uint8_t*>(active[1][1].ptr(y));
    for (int x = 0; x < w; ++x)
    {
      if (ptr_active[x][0] == top)
      {
        if (ptr_active[x][1] == left)
        {
          ptr_active_tl[x] = 1;
        }
        else
        {
          ptr_active_tr[x] = 1;
        }
      }
      else
      {
        if (ptr_active[x][1] == left)
        {
          ptr_active_bl[x] = 1;
        }
        else
        {
          ptr_active_br[x] = 1;
        }
      }
    }
  }
}

static cv::Mat get_patch_boundary(const MergePatch& patch, int overscan, double scale)
{
  const double delta = 1.0 / scale;

  std::vector<double> x_left, x_right;
  std::vector<double> y_top, y_bottom;

  cv::Point p, p_last;
  const cv::Point p_invalid(-1, -1);

  cv::Mat mask = cv::Mat::zeros(static_cast<int>(patch.size.height * scale) + 2 * overscan, static_cast<int>(patch.size.width * scale) + 2 * overscan, CV_8UC1);

  cv::Point top_first = p_invalid;
  cv::Point top_last = p_invalid;
  p_last = p_invalid;
  for (int x = 0; x < mask.cols; ++x)
  {
    double pos_x = patch.anchor_target.x + (x - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_top())
    {
      const double t = b.find_x(pos_x);
      p = p_invalid;
      if (t >= 0.0 && t <= 1.0)
      {
        int y = static_cast<int>((b.eval_y(t) - patch.anchor_target.y) * scale + overscan);
        y = std::max(y, 0);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (top_first == p_invalid)
      {
        top_first = p;
      }
      top_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::Point bottom_first = p_invalid;
  cv::Point bottom_last = p_invalid;
  p_last = p_invalid;
  for (int x = 0; x < mask.cols; ++x)
  {
    double pos_x = patch.anchor_target.x + (x - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_bottom())
    {
      const double t = b.find_x(pos_x);
      p = p_invalid;
      if (t >= 0.0 && t <= 1.0)
      {
        int y = static_cast<int>((b.eval_y(t) - patch.anchor_target.y) * scale + overscan);
        y = std::max(y, 0);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (bottom_first == p_invalid)
      {
        bottom_first = p;
      }
      bottom_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::Point left_first = p_invalid;
  cv::Point left_last = p_invalid;
  p_last = p_invalid;
  for (int y = 0; y < mask.rows; ++y)
  {
    double pos_y = patch.anchor_target.y + (y - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_left())
    {
      const double t = b.find_y(pos_y);
      if (t >= 0.0 && t <= 1.0)
      {
        int x = static_cast<int>((b.eval_x(t) - patch.anchor_target.x) * scale + overscan);
        x = std::max(x, 0);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (left_first == p_invalid)
      {
        left_first = p;
      }
      left_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::Point right_first = p_invalid;
  cv::Point right_last = p_invalid;
  p_last = p_invalid;
  for (int y = 0; y < mask.rows; ++y)
  {
    double pos_y = patch.anchor_target.y + (y - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_right())
    {
      const double t = b.find_y(pos_y);
      if (t >= 0.0 && t <= 1.0)
      {
        int x = static_cast<int>((b.eval_x(t) - patch.anchor_target.x) * scale + overscan);
        x = std::min(x, mask.cols-1);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (right_first == p_invalid)
      {
        right_first = p;
      }
      right_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::line(mask, top_last, right_first, cv::Scalar(255));
  cv::line(mask, right_last, bottom_last, cv::Scalar(255));
  cv::line(mask, bottom_first, left_last, cv::Scalar(255));
  cv::line(mask, left_first, top_first, cv::Scalar(255));

  return mask;
}

static cv::Mat get_source_patch_boundary(const MergePatch& patch, int overscan, double scale)
{
  const double delta = 1.0 / scale;

  std::vector<double> x_left, x_right;
  std::vector<double> y_top, y_bottom;

  cv::Point p, p_last;
  const cv::Point p_invalid(-1, -1);

  cv::Mat mask = cv::Mat::zeros(static_cast<int>(patch.size.height * scale) + 2 * overscan, static_cast<int>(patch.size.width * scale) + 2 * overscan, CV_8UC1);

  cv::Point top_first = p_invalid;
  cv::Point top_last = p_invalid;
  p_last = p_invalid;
  for (int x = 0; x < mask.cols; ++x)
  {
    double pos_x = patch.anchor_target.x + (x - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_top())
    {
      const double t = b.find_x(pos_x);
      p = p_invalid;
      if (t >= 0.0 && t <= 1.0)
      {
        int y = static_cast<int>((b.eval_y(t) - patch.anchor_target.y) * scale + overscan);
        y = std::max(y, 0);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (top_first == p_invalid)
      {
        top_first = p;
      }
      top_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::Point bottom_first = p_invalid;
  cv::Point bottom_last = p_invalid;
  p_last = p_invalid;
  for (int x = 0; x < mask.cols; ++x)
  {
    double pos_x = patch.anchor_target.x + (x - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_bottom())
    {
      const double t = b.find_x(pos_x);
      p = p_invalid;
      if (t >= 0.0 && t <= 1.0)
      {
        int y = static_cast<int>((b.eval_y(t) - patch.anchor_target.y) * scale + overscan);
        y = std::max(y, 0);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (bottom_first == p_invalid)
      {
        bottom_first = p;
      }
      bottom_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::Point left_first = p_invalid;
  cv::Point left_last = p_invalid;
  p_last = p_invalid;
  for (int y = 0; y < mask.rows; ++y)
  {
    double pos_y = patch.anchor_target.y + (y - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_left())
    {
      const double t = b.find_y(pos_y);
      if (t >= 0.0 && t <= 1.0)
      {
        int x = static_cast<int>((b.eval_x(t) - patch.anchor_target.x) * scale + overscan);
        x = std::max(x, 0);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (left_first == p_invalid)
      {
        left_first = p;
      }
      left_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::Point right_first = p_invalid;
  cv::Point right_last = p_invalid;
  p_last = p_invalid;
  for (int y = 0; y < mask.rows; ++y)
  {
    double pos_y = patch.anchor_target.y + (y - overscan) * delta;

    p = p_invalid;
    for (const BezierCurve& b : patch.target_curves_right())
    {
      const double t = b.find_y(pos_y);
      if (t >= 0.0 && t <= 1.0)
      {
        int x = static_cast<int>((b.eval_x(t) - patch.anchor_target.x) * scale + overscan);
        x = std::min(x, mask.cols-1);
        p.x = x;
        p.y = y;
        break;
      }
    }

    if (p != p_invalid)
    {
      if (right_first == p_invalid)
      {
        right_first = p;
      }
      right_last = p;
      if (p_last != p_invalid)
      {
        cv::line(mask, p_last, p, cv::Scalar(255));
      }
      p_last = p;
    }
  }

  cv::line(mask, top_last, right_first, cv::Scalar(255));
  cv::line(mask, right_last, bottom_last, cv::Scalar(255));
  cv::line(mask, bottom_first, left_last, cv::Scalar(255));
  cv::line(mask, left_first, top_first, cv::Scalar(255));

  return mask;
}

static cv::Mat get_patch_mask(const MergePatch& patch, int overscan, double scale)
{
  /*
   * Get patch boundary.
   */
  cv::Mat mask = get_patch_boundary(patch, overscan, scale);

  /* 
   * Fill holes
   */
  std::deque<cv::Point> stack;

  for (int x = 0; x < mask.cols; ++x)
  {
    stack.emplace_back(x, 0);
    stack.emplace_back(x, mask.rows-1);
  }

  for (int y = 0; y < mask.rows; ++y)
  {
    stack.emplace_back(0, y);
    stack.emplace_back(mask.cols-1, y);
  }

  while (!stack.empty())
  {
    const cv::Point p = stack.front();
    stack.pop_front();

    if (p.x >= 0 && p.x < mask.cols &&
      p.y >= 0 && p.y < mask.rows &&
      !mask.at<unsigned char>(p))
    {
      mask.at<unsigned char>(p) = 255;
      stack.emplace_back(p.x, p.y + 1);
      stack.emplace_back(p.x, p.y - 1);
      stack.emplace_back(p.x - 1, p.y);
      stack.emplace_back(p.x + 1, p.y);
    }
  }

  mask = 255 - mask;
  cv::Mat dilate_mask = (cv::Mat_<unsigned char>(3, 3) << 0, 1, 0, 1, 1, 1, 0, 1, 0);
  cv::dilate(mask, mask, cv::Mat::ones(3, 3, CV_8UC1));

  return mask;
}

static void get_patch_boundaries(const MergePatch& patch, float dx, float dy, std::vector<float>& x_vals, std::vector<float>& y_vals, std::vector<float>& x_start, std::vector<float>& x_end, std::vector<float>& y_start, std::vector<float>& y_end)
{
  const std::vector<BezierCurve> target_curves_top = patch.target_curves_top();
  const std::vector<BezierCurve> target_curves_bottom = patch.target_curves_bottom();
  const std::vector<BezierCurve> target_curves_left = patch.target_curves_left();
  const std::vector<BezierCurve> target_curves_right = patch.target_curves_right();

  float x, y;

  x_vals.clear();
  y_vals.clear();
  x_start.clear();
  x_end.clear();
  y_start.clear();
  y_end.clear();

  x_vals.push_back(static_cast<float>(patch.anchor_target.x));
  y_start.push_back(static_cast<float>(std::min(target_curves_top[0].eval_y(0.0), target_curves_top[0].eval_y(1.0))));
  y_end.push_back(static_cast<float>(std::max(target_curves_bottom[0].eval_y(0.0), target_curves_bottom[0].eval_y(1.0))));
  for (x = patch.anchor_target.x + dx; x < patch.anchor_target.x + patch.size.width; x += dx)
  {
    x_vals.push_back(x);
    bool success = false;
    for (const BezierCurve& b : target_curves_top)
    {
      double t = b.find_x(x);
      if (t >= 0.0 && t <= 1.0)
      {
        y_start.push_back(static_cast<float>(b.eval_y(t)));
        success = true;
        break;
      }
    }
    if (!success)
    {
      y_start.push_back(static_cast<float>(patch.anchor_target.y));
    }

    success = false;
    for (const BezierCurve& b : target_curves_bottom)
    {
      double t = b.find_x(x);
      if (t >= 0.0 && t <= 1.0)
      {
        y_end.push_back(static_cast<float>(b.eval_y(t)));
        success = true;
        break;
      }
    }
    if (!success)
    {
      y_end.push_back(static_cast<float>(patch.anchor_target.y + patch.size.height));
    }
  }

  y_vals.push_back(static_cast<float>(patch.anchor_target.y));
  x_start.push_back(static_cast<float>(std::min(target_curves_left[0].eval_x(0.0), target_curves_left[0].eval_x(1.0))));
  x_end.push_back(static_cast<float>(std::max(target_curves_right[0].eval_x(0.0), target_curves_right[0].eval_x(1.0))));
  for (y = patch.anchor_target.y + dy; y < patch.anchor_target.y + patch.size.height; y += dy)
  {
    y_vals.push_back(y);
    bool success = false;
    for (const BezierCurve& b : target_curves_left)
    {
      double t = b.find_y(y);
      if (t >= 0.0 && t <= 1.0)
      {
        x_start.push_back(static_cast<float>(b.eval_x(t)));
        success = true;
        break;
      }
    }
    if (!success)
    {
      x_start.push_back(static_cast<float>(patch.anchor_target.x));
    }

    success = false;
    for (const BezierCurve& b : target_curves_right)
    {
      double t = b.find_y(y);
      if (t >= 0.0 && t <= 1.0)
      {
        x_end.push_back(static_cast<float>(b.eval_x(t)));
        success = true;
        break;
      }
    }
    if (!success)
    {
      x_end.push_back(static_cast<float>(patch.anchor_target.x + patch.size.width));
    }
  }
}

cv::Mat MergePatch::draw(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, double scale, bool draw_boundaries)
{
  if (patches.empty())
  {
    return cv::Mat();
  }

  cv::Rect bbox(patches[0].anchor_target, patches[0].size);
  for (const MergePatch& patch : patches)
  {
    bbox = cv::boundingRect(std::vector<cv::Point>({bbox.tl(), bbox.br() - cv::Point(1, 1), patch.anchor_target, patch.anchor_target + cv::Point(patch.size)}));
  }

  cv::Mat image(static_cast<int>((bbox.height-1) * scale), static_cast<int>((bbox.width-1) * scale), CV_8UC3, cv::Scalar(255, 255, 255));
  cv::Mat boundary_mask;
  if (draw_boundaries)
  {
    boundary_mask = cv::Mat::zeros(image.size(), CV_32FC1);
  }

  for (const MergePatch& patch : patches)
  {
    const int overscan = static_cast<int>(2 * scale);
    cv::Mat mask = get_patch_mask(patch, overscan, scale);

    cv::Rect rect_roi(static_cast<int>(scale * (patch.anchor_target.x - bbox.x) - overscan), static_cast<int>(scale * (patch.anchor_target.y - bbox.y) - overscan), mask.cols, mask.rows);
    int x_beg = 0;
    int x_end = mask.cols;

    int y_beg = 0;
    int y_end = mask.rows;

    if (rect_roi.x < 0)
    {
      x_beg = -rect_roi.x;
      rect_roi.width += rect_roi.x;
      rect_roi.x = 0;
    }

    if (rect_roi.y < 0)
    {
      y_beg = -rect_roi.y;
      rect_roi.height += rect_roi.y;
      rect_roi.y = 0;
    }

    if (rect_roi.x + rect_roi.width > image.cols)
    {
      x_end -= rect_roi.x + rect_roi.width - image.cols;
      rect_roi.width -= rect_roi.x + rect_roi.width - image.cols;
    }

    if (rect_roi.y + rect_roi.height > image.rows)
    {
      y_end -= rect_roi.y + rect_roi.height - image.rows;
      rect_roi.height -= rect_roi.y + rect_roi.height - image.rows;
    }

    cv::Mat image_roi = image(rect_roi);
    for (int y = y_beg; y < y_end; ++y)
    {
      const unsigned char* mask_ptr = reinterpret_cast<const unsigned char*>(mask.ptr(y));
      cv::Vec3b* image_ptr = reinterpret_cast<cv::Vec3b*>(image_roi.ptr(y-y_beg));
      for (int x = x_beg; x < x_end; ++x)
      {
        if (mask_ptr[x])
        {
          cv::Point2f p_source = AffineTransformation::transform(patch.transformation_source_inv,
            cv::Point2f(static_cast<float>(patch.anchor_source.x + (x - overscan) / scale), static_cast<float>(patch.anchor_source.y + (y - overscan) / scale)));
          image_ptr[x-x_beg] = textures[patch.source_index].interpolate_texture(p_source);
        }
      }
    }

    if (draw_boundaries)
    {
      cv::Mat patch_boundary = get_patch_boundary(patch, overscan, scale);
      cv::Mat boundary_mask_roi = boundary_mask(rect_roi);
      for (int y = y_beg; y < y_end; ++y)
      {
        const unsigned char* mask_ptr = reinterpret_cast<const unsigned char*>(patch_boundary.ptr(y));
        float* boundary_mask_ptr = reinterpret_cast<float*>(boundary_mask_roi.ptr(y-y_beg));
        for (int x = x_beg; x < x_end; ++x)
        {
          if (mask_ptr[x])
          {
            boundary_mask_ptr[x-x_beg] = 1.0f;
          }
        }
      }
    }
  }

  if (draw_boundaries)
  {
    const cv::Vec3b dark_brown(14, 29, 43);
    cv::GaussianBlur(boundary_mask, boundary_mask, cv::Size(5, 5), 0.0);
    
    for (int y = 0; y < image.rows; ++y)
    {
      cv::Vec3b* ptr_image = reinterpret_cast<cv::Vec3b*>(image.ptr(y));
      const float* ptr_mask = reinterpret_cast<const float*>(boundary_mask.ptr(y));
      for (int x = 0; x < image.cols; ++x)
      {
        ptr_image[x] = ptr_mask[x] * dark_brown + (1.0f - ptr_mask[x]) * ptr_image[x];
      }
    }
  }

  return image;
}

cv::Mat MergePatch::draw_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale, bool draw_boundaries)
{
  std::vector<Texture> textures_fullres;
  for (const Texture& t : textures)
  {
    textures_fullres.emplace_back(base_path / t.filename, t.dpi, 1.0);
  }
  return draw_fullres(patches, textures, textures_fullres, scale, draw_boundaries);
}

cv::Mat MergePatch::draw_fullres(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale, bool draw_boundaries)
{
  if (patches.empty())
  {
    return cv::Mat();
  }

  cv::Rect bbox(patches[0].anchor_target, patches[0].size);
  for (const MergePatch& patch : patches)
  {
    bbox = cv::boundingRect(std::vector<cv::Point>({bbox.tl(), bbox.br() - cv::Point(1, 1), patch.anchor_target, patch.anchor_target + cv::Point(patch.size)}));
  }

  const float dx = static_cast<float>(1.0f / scale);
  const float dy = static_cast<float>(1.0f / scale);

  cv::Mat image(static_cast<int>((bbox.height-1) * scale), static_cast<int>((bbox.width-1) * scale), CV_8UC3, cv::Scalar(255, 255, 255));
  cv::Mat boundary_mask;
  if (draw_boundaries)
  {
    boundary_mask = cv::Mat::zeros(image.size(), CV_32FC1);
  }

  for (const MergePatch& patch : patches)
  {
    const double texture_scale = 1.0 / textures[patch.source_index].scale;
    const cv::Mat T_source = AffineTransformation::concat(AffineTransformation::T_scale(texture_scale, texture_scale), patch.transformation_source_inv);

    const int overscan = static_cast<int>(2 * scale);
    cv::Mat mask = get_patch_mask(patch, overscan, scale);

    cv::Rect rect_roi(static_cast<int>(scale * (patch.anchor_target.x - bbox.x) - overscan), static_cast<int>(scale * (patch.anchor_target.y - bbox.y) - overscan), mask.cols, mask.rows);
    int x_beg = 0;
    int x_end = mask.cols;

    int y_beg = 0;
    int y_end = mask.rows;

    if (rect_roi.x < 0)
    {
      x_beg = -rect_roi.x;
      rect_roi.width += rect_roi.x;
      rect_roi.x = 0;
    }

    if (rect_roi.y < 0)
    {
      y_beg = -rect_roi.y;
      rect_roi.height += rect_roi.y;
      rect_roi.y = 0;
    }

    if (rect_roi.x + rect_roi.width > image.cols)
    {
      x_end -= rect_roi.x + rect_roi.width - image.cols;
      rect_roi.width -= rect_roi.x + rect_roi.width - image.cols;
    }

    if (rect_roi.y + rect_roi.height > image.rows)
    {
      y_end -= rect_roi.y + rect_roi.height - image.rows;
      rect_roi.height -= rect_roi.y + rect_roi.height - image.rows;
    }

    cv::Mat image_roi = image(rect_roi);

    for (int y = y_beg; y < y_end; ++y)
    {
      const unsigned char* mask_ptr = reinterpret_cast<const unsigned char*>(mask.ptr(y));
      cv::Vec3b* image_ptr = reinterpret_cast<cv::Vec3b*>(image_roi.ptr(y-y_beg));
      for (int x = x_beg; x < x_end; ++x)
      {
        if (mask_ptr[x])
        {
          cv::Point2f p_source = AffineTransformation::transform(T_source,
            cv::Point2f(static_cast<float>(patch.anchor_source.x + (x - overscan) / scale), static_cast<float>(patch.anchor_source.y + (y - overscan) / scale)));
          image_ptr[x-x_beg] = textures_fullres[patch.source_index].interpolate_texture(p_source);
        }
      }
    }

    if (draw_boundaries)
    {
      cv::Mat patch_boundary = get_patch_boundary(patch, overscan, scale);
      cv::Mat boundary_mask_roi = boundary_mask(rect_roi);
      for (int y = y_beg; y < y_end; ++y)
      {
        const unsigned char* mask_ptr = reinterpret_cast<const unsigned char*>(patch_boundary.ptr(y));
        float* boundary_mask_ptr = reinterpret_cast<float*>(boundary_mask_roi.ptr(y-y_beg));
        for (int x = x_beg; x < x_end; ++x)
        {
          if (mask_ptr[x])
          {
            boundary_mask_ptr[x-x_beg] = 1.0f;
          }
        }
      }
    }
  }

  if (draw_boundaries)
  {
    const cv::Vec3b dark_brown(14, 29, 43);
    cv::GaussianBlur(boundary_mask, boundary_mask, cv::Size(5, 5), 0.0);

    for (int y = 0; y < image.rows; ++y)
    {
      cv::Vec3b* ptr_image = reinterpret_cast<cv::Vec3b*>(image.ptr(y));
      const float* ptr_mask = reinterpret_cast<const float*>(boundary_mask.ptr(y));
      for (int x = 0; x < image.cols; ++x)
      {
        ptr_image[x] = ptr_mask[x] * dark_brown + (1.0f - ptr_mask[x]) * ptr_image[x];
      }
    }
  }

  return image;
}

std::vector<cv::Mat> MergePatch::draw_source_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale)
{
  std::vector<Texture> textures_fullres;
  for (const Texture& t : textures)
  {
    textures_fullres.emplace_back(base_path / t.filename, t.dpi, 1.0);
  }
  return draw_source_fullres(patches, textures, textures_fullres, scale);
}

std::vector<cv::Mat> MergePatch::draw_source_fullres(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale)
{
  std::vector<cv::Mat> images;
  for (const Texture& t : textures)
  {
    cv::Mat image = t.texture.clone();
    image.convertTo(image, CV_8UC3, 1.0 / 255.0);
    images.push_back(image);
  }

  draw_bezier_curves_source(images, patches, scale);

  for (const cv::Mat& image : images)
  {
    cv::imshow("Image", image);
    cv::waitKey(0);
  }
  std::exit(0);
}

void MergePatch::merge_patches_x(std::vector<MergePatch>& patches, const mat<std::vector<int>>& patch_indices, int x, int y_start, int y_end, int subpatch_size, int start_top, int start_bottom)
{
  /*
  * Build error matrix.
  */
  cv::Mat error_left((y_end - y_start) * subpatch_size, subpatch_size, CV_32FC1, std::numeric_limits<float>::infinity());
  cv::Mat error_right((y_end - y_start) * subpatch_size, subpatch_size, CV_32FC1, std::numeric_limits<float>::infinity());

  std::vector<std::vector<int>> indices_left_vec, indices_right_vec;
  std::vector<std::vector<cv::Mat>> active_left_vec, active_right_vec;
  std::vector<std::vector<cv::Mat>> error_left_vec, error_right_vec;

  for (int y = y_start; y < y_end; ++y)
  {
    const std::vector<int>& indices_coord = patch_indices(y, x);

    indices_left_vec.emplace_back();
    indices_right_vec.emplace_back();
    active_left_vec.emplace_back();
    active_right_vec.emplace_back();
    error_left_vec.emplace_back();
    error_right_vec.emplace_back();

    for (int i : indices_coord)
    {
      if (patches[i].region_subpatch.x < x)
      {
        indices_left_vec.back().push_back(i);
        active_left_vec.back().push_back(patches[i].get_active_pixel(y, x));
        error_left_vec.back().push_back(patches[i].get_error_mat(y, x));
      }
      else
      {
        indices_right_vec.back().push_back(i);
        active_right_vec.back().push_back(patches[i].get_active_pixel(y, x));
        error_right_vec.back().push_back(patches[i].get_error_mat(y, x));
      }
    }

    cv::Mat error_left_subpatch = error_left(cv::Rect(0, (y - y_start) * subpatch_size, subpatch_size, subpatch_size));
    cv::Mat error_right_subpatch = error_right(cv::Rect(0, (y - y_start) * subpatch_size, subpatch_size, subpatch_size));

    for (int y_pix = 0; y_pix < subpatch_size; ++y_pix)
    {
      float* ptr_error_left = reinterpret_cast<float*>(error_left_subpatch.ptr(y_pix));
      for (size_t i = 0; i < indices_left_vec.back().size(); ++i)
      {
        const uint8_t* ptr_active = reinterpret_cast<const uint8_t*>(active_left_vec.back()[i].ptr(y_pix));
        const float* ptr_error = reinterpret_cast<const float*>(error_left_vec.back()[i].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix])
          {
            ptr_error_left[x_pix] = ptr_error[x_pix];
          }
        }
      }

      float* ptr_error_right = reinterpret_cast<float*>(error_right_subpatch.ptr(y_pix));
      for (size_t i = 0; i < indices_right_vec.back().size(); ++i)
      {
        const uint8_t* ptr_active = reinterpret_cast<const uint8_t*>(active_right_vec.back()[i].ptr(y_pix));
        const float* ptr_error = reinterpret_cast<const float*>(error_right_vec.back()[i].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix])
          {
            ptr_error_right[x_pix] = ptr_error[x_pix];
          }
        }
      }
    }
  }

  cv::Mat active = merge_patches_x(error_left, error_right, start_top, start_bottom);
  for (int i = 0; i < static_cast<int>(indices_left_vec.size()); ++i)
  {
    cv::Mat active_subpatch = active(cv::Rect(0, i * subpatch_size, subpatch_size, subpatch_size));
    for (int y_pix = 0; y_pix < subpatch_size; ++y_pix)
    {
      const uint8_t* ptr_active = reinterpret_cast<const uint8_t*>(active_subpatch.ptr(y_pix));
      for (size_t j = 0; j < indices_left_vec[i].size(); ++j)
      {
        uint8_t* ptr_active_left = reinterpret_cast<uint8_t*>(active_left_vec[i][j].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix] == right)
          {
            ptr_active_left[x_pix] = 0;
          }
        }
      }
      for (size_t j = 0; j < indices_right_vec[i].size(); ++j)
      {
        uint8_t* ptr_active_right = reinterpret_cast<uint8_t*>(active_right_vec[i][j].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix] == left)
          {
            ptr_active_right[x_pix] = 0;
          }
        }
      }
    }
  }
}

void MergePatch::merge_patches_y(std::vector<MergePatch>& patches, const mat<std::vector<int>>& patch_indices, int y, int x_start, int x_end, int subpatch_size, int start_left, int start_right)
{
  /*
   * Build error matrix.
   */
  cv::Mat error_top(subpatch_size, (x_end - x_start) * subpatch_size, CV_32FC1, std::numeric_limits<float>::infinity());
  cv::Mat error_bottom(subpatch_size, (x_end - x_start) * subpatch_size, CV_32FC1, std::numeric_limits<float>::infinity());

  std::vector<std::vector<int>> indices_top_vec, indices_bottom_vec;
  std::vector<std::vector<cv::Mat>> active_top_vec, active_bottom_vec;
  std::vector<std::vector<cv::Mat>> error_top_vec, error_bottom_vec;

  for (int x = x_start; x < x_end; ++x)
  {
    const std::vector<int>& indices_coord = patch_indices(y, x);

    indices_top_vec.emplace_back();
    indices_bottom_vec.emplace_back();
    active_top_vec.emplace_back();
    active_bottom_vec.emplace_back();
    error_top_vec.emplace_back();
    error_bottom_vec.emplace_back();

    for (int i : indices_coord)
    {
      if (patches[i].region_subpatch.y < y)
      {
        indices_top_vec.back().push_back(i);
        active_top_vec.back().push_back(patches[i].get_active_pixel(y, x));
        error_top_vec.back().push_back(patches[i].get_error_mat(y, x));
      }
      else
      {
        indices_bottom_vec.back().push_back(i);
        active_bottom_vec.back().push_back(patches[i].get_active_pixel(y, x));
        error_bottom_vec.back().push_back(patches[i].get_error_mat(y, x));
      }
    }

    cv::Mat error_top_subpatch = error_top(cv::Rect((x - x_start) * subpatch_size, 0, subpatch_size, subpatch_size));
    cv::Mat error_bottom_subpatch = error_bottom(cv::Rect((x - x_start) * subpatch_size, 0, subpatch_size, subpatch_size));

    for (int y_pix = 0; y_pix < subpatch_size; ++y_pix)
    {
      float* ptr_error_top = reinterpret_cast<float*>(error_top_subpatch.ptr(y_pix));
      for (size_t i = 0; i < indices_top_vec.back().size(); ++i)
      {
        const uint8_t* ptr_active = reinterpret_cast<const uint8_t*>(active_top_vec.back()[i].ptr(y_pix));
        const float* ptr_error = reinterpret_cast<const float*>(error_top_vec.back()[i].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix])
          {
            ptr_error_top[x_pix] = ptr_error[x_pix];
          }
        }
      }

      float* ptr_error_bottom = reinterpret_cast<float*>(error_bottom_subpatch.ptr(y_pix));
      for (size_t i = 0; i < indices_bottom_vec.back().size(); ++i)
      {
        const uint8_t* ptr_active = reinterpret_cast<const uint8_t*>(active_bottom_vec.back()[i].ptr(y_pix));
        const float* ptr_error = reinterpret_cast<const float*>(error_bottom_vec.back()[i].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix])
          {
            ptr_error_bottom[x_pix] = ptr_error[x_pix];
          }
        }
      }
    }
  }

  cv::Mat active = merge_patches_y(error_top, error_bottom, start_left, start_right);
  for (int i = 0; i < static_cast<int>(indices_top_vec.size()); ++i)
  {
    cv::Mat active_subpatch = active(cv::Rect(i * subpatch_size, 0, subpatch_size, subpatch_size));
    for (int y_pix = 0; y_pix < subpatch_size; ++y_pix)
    {
      const uint8_t* ptr_active = reinterpret_cast<const uint8_t*>(active_subpatch.ptr(y_pix));
      for (size_t j = 0; j < indices_top_vec[i].size(); ++j)
      {
        uint8_t* ptr_active_top = reinterpret_cast<uint8_t*>(active_top_vec[i][j].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix] == bottom)
          {
            ptr_active_top[x_pix] = 0;
          }
        }
      }
      for (size_t j = 0; j < indices_bottom_vec[i].size(); ++j)
      {
        uint8_t* ptr_active_bottom = reinterpret_cast<uint8_t*>(active_bottom_vec[i][j].ptr(y_pix));
        for (int x_pix = 0; x_pix < subpatch_size; ++x_pix)
        {
          if (ptr_active[x_pix] == top)
          {
            ptr_active_bottom[x_pix] = 0;
          }
        }
      }
    }
  }
}

void MergePatch::merge_patches_x(MergePatch& patch_left, MergePatch& patch_right, int y_subpatch, int x_subpatch)
{
  cv::Mat error_left = patch_left.get_error_mat(y_subpatch, x_subpatch);
  cv::Mat error_right = patch_right.get_error_mat(y_subpatch, x_subpatch);
  cv::Mat active = merge_patches_x(error_left, error_right);

  cv::Mat active_left = patch_left.get_active_pixel(y_subpatch, x_subpatch);
  cv::Mat active_right = patch_right.get_active_pixel(y_subpatch, x_subpatch);

  for (int y = 0; y < active.rows; ++y)
  {
    uint8_t* ptr_left = active_left.ptr(y);
    uint8_t* ptr_right = active_right.ptr(y);
    const uint8_t* ptr = active.ptr(y);
    for (int x = 0; x < active.cols; ++x)
    {
      if (ptr[x] == left)
      {
        ptr_left[x] = 1;
        ptr_right[x] = 0;
      }
      else
      {
        ptr_left[x] = 0;
        ptr_right[x] = 1;
      }
    }
  }
}

void MergePatch::merge_patches_y(MergePatch& patch_top, MergePatch& patch_bottom, int y_subpatch, int x_subpatch)
{
  cv::Mat error_top = patch_top.get_error_mat(y_subpatch, x_subpatch);
  cv::Mat error_bottom = patch_bottom.get_error_mat(y_subpatch, x_subpatch);
  cv::Mat active = merge_patches_y(error_top, error_bottom);

  cv::Mat active_top = patch_top.get_active_pixel(y_subpatch, x_subpatch);
  cv::Mat active_bottom = patch_bottom.get_active_pixel(y_subpatch, x_subpatch);

  for (int y = 0; y < active.rows; ++y)
  {
    uint8_t* ptr_top = active_top.ptr(y);
    uint8_t* ptr_bottom = active_bottom.ptr(y);
    const uint8_t* ptr = active.ptr(y);
    for (int x = 0; x < active.cols; ++x)
    {
      if (ptr[x] == top)
      {
        ptr_top[x] = 1;
        ptr_bottom[x] = 0;
      }
      else
      {
        ptr_top[x] = 0;
        ptr_bottom[x] = 1;
      }
    }
  }
}

void MergePatch::merge_patches_cross(std::vector<MergePatch>& patches, const std::vector<int> cross_indices, int subpatch_size)
{
  if (cross_indices.size() != 4)
  {
    throw(std::invalid_argument("merge_patches_cross: cross_indices.size() != 4."));
  }

  /*
   * Get overlapping patch region.
   */
  cv::Rect overlap_region(0, 0, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
  for (int i : cross_indices)
  {
    overlap_region &= cv::Rect(patches[i].anchor_target, patches[i].size);
  }

  if (overlap_region.area() == 0)
  {
    throw(std::invalid_argument("merge_patches_cross regions do not overlap."));
  }

  /*
   * Find out which patch belongs to which corner of the overlapping region.
   */
  MergePatch* patch_tl = 0;
  MergePatch* patch_tr = 0;
  MergePatch* patch_bl = 0;
  MergePatch* patch_br = 0;
  cv::Point overlap_region_center = (overlap_region.tl() + overlap_region.br()) / 2;
  for (int i : cross_indices)
  {
    cv::Point p_center = patches[i].anchor_target + cv::Point(patches[i].size) / 2;
    if (p_center.x < overlap_region_center.x)
    {
      if (p_center.y < overlap_region_center.y)
      {
        patch_tl = &patches[i];
      }
      else
      {
        patch_bl = &patches[i];
      }
    }
    else
    {
      if (p_center.y < overlap_region_center.y)
      {
        patch_tr = &patches[i];
      }
      else
      {
        patch_br = &patches[i];
      }
    }
  }

  if (!patch_tl || !patch_tr || !patch_bl || !patch_br)
  {
    throw(std::runtime_error("merge_patches_cross failed."));
  }

  cv::Rect region_tl = overlap_region - patch_tl->anchor_target;
  cv::Rect region_tr = overlap_region - patch_tr->anchor_target;
  cv::Rect region_bl = overlap_region - patch_bl->anchor_target;
  cv::Rect region_br = overlap_region - patch_br->anchor_target;

  merge_patches_cross(*patch_tl, region_tl, *patch_bl, region_bl, *patch_tr, region_tr, *patch_br, region_br);
}

std::vector<int> MergePatch::get_cut_coordinates_horiz(int y_subpatch, int x_subpatch) const
{
  std::vector<int> coordinates;

  const cv::Point subpatch_center = (region_subpatch.tl() + region_subpatch.br()) / 2;
  const uint8_t val_top = subpatch_center.y < y_subpatch ? 1 : 0;
  std::cout << static_cast<int>(val_top) << std::endl;

  cv::Mat active_mat = get_active_pixel(y_subpatch, x_subpatch);



  for (int x = 0; x < active_mat.cols; ++x)
  {
    int y = 0;
    for (; y < active_mat.rows; ++y)
    {
      if (active_mat.at<uint8_t>(y, x) != val_top)
      {
        break;
      }
    }
    coordinates.push_back(y);
  }

  std::copy(coordinates.begin(), coordinates.end(), std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;

  cv::Mat active_mat_large;
  cv::resize(active_mat, active_mat_large, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);
  cv::imshow("active_mat", 255 * active_mat_large);
  cv::waitKey(0);

  return coordinates;
}

std::vector<int> MergePatch::get_cut_coordinates_vert(int y_subpatch, int x_subpatch) const
{
  std::vector<int> coordinates;

  const cv::Point subpatch_center = (region_subpatch.tl() + region_subpatch.br()) / 2;
  const uint val_left = subpatch_center.x < x_subpatch ? 1 : 0;

  cv::Mat active_mat = get_active_pixel(y_subpatch, x_subpatch);
  
  for (int y = 0; y < active_mat.rows; ++y)
  {
    int x = 0;
    for (; x < active_mat.cols; ++x)
    {
      if (active_mat.at<uint8_t>(y, x) != val_left)
      {
        break;
      }
    }
    coordinates.push_back(x);
  }

  return coordinates;
}

cv::Mat MergePatch::get_error_mat(int y_subpatch, int x_subpatch) const
{
  const int x_local = x_subpatch - region_subpatch.x;
  const int y_local = y_subpatch - region_subpatch.y;
  cv::Rect region(x_local * subpatch_size, y_local * subpatch_size, subpatch_size, subpatch_size);
  return error(region);
}

cv::Mat MergePatch::get_error_mat(cv::Rect region_global) const
{
  const cv::Rect region_local = region_global - anchor_target;
  return error(region_local);
}

cv::Mat MergePatch::get_active_pixel(int y_subpatch, int x_subpatch) const
{
  const int x_local = x_subpatch - region_subpatch.x;
  const int y_local = y_subpatch - region_subpatch.y;
  cv::Rect region(x_local * subpatch_size, y_local * subpatch_size, subpatch_size, subpatch_size);
  return active_pixel(region);
}

cv::Mat MergePatch::get_active_pixel(cv::Rect region_global) const
{
  cv::Rect region_local = region_global - anchor_target;
  return active_pixel(region_local);
}

void MergePatch::add_curve_horiz(int y_subpatch, int x_subpatch, const BezierCurve& curve)
{
  if (region_subpatch.y + region_subpatch.height / 2 < y_subpatch)
  {
    m_curves_bottom.push_back(curve);
  }
  else
  {
    m_curves_top.push_back(curve);
  }
}

void MergePatch::add_curve_vert(int y_subpatch, int x_subpatch, const BezierCurve& curve)
{
  if (region_subpatch.x + region_subpatch.width / 2 < x_subpatch)
  {
    m_curves_right.push_back(curve);
  }
  else
  {
    m_curves_left.push_back(curve);
  }
}

void MergePatch::draw_bezier_curves(cv::Mat image, const std::vector<MergePatch>& patches, double factor, const cv::Point2d& delta)
{
  for (const MergePatch& patch : patches)
  {
    for (const BezierCurve& curve : patch.m_curves_left)
    {
      (curve + delta).draw(image, patch.subpatch_size, factor);
    }
    for (const BezierCurve& curve : patch.m_curves_right)
    {
      (curve + delta).draw(image, patch.subpatch_size, factor);
    }
    for (const BezierCurve& curve : patch.m_curves_top)
    {
      (curve + delta).draw(image, patch.subpatch_size, factor);
    }
    for (const BezierCurve& curve : patch.m_curves_bottom)
    {
      (curve + delta).draw(image, patch.subpatch_size, factor);
    }
  }
}

void MergePatch::draw_bezier_curves_source(std::vector<cv::Mat>& images, const std::vector<MergePatch>& patches, double factor)
{
  for (const MergePatch& patch : patches)
  {
    cv::Mat image = images[patch.source_index];
    for (const BezierCurve& curve : patch.source_curves_left())
    {
      curve.draw(image, patch.subpatch_size, factor);
    }
    for (const BezierCurve& curve : patch.source_curves_right())
    {
      curve.draw(image, patch.subpatch_size, factor);
    }
    for (const BezierCurve& curve : patch.source_curves_top())
    {
      curve.draw(image, patch.subpatch_size, factor);
    }
    for (const BezierCurve& curve : patch.source_curves_bottom())
    {
      curve.draw(image, patch.subpatch_size, factor);
    }
  }
}

std::vector<int> MergePatch::get_cut_coordinates_horiz(const std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat, int y_subpatch, int x_subpatch, int subpatch_size)
{
  std::vector<int> coordinates;
  cv::Mat active_mat = cv::Mat::zeros(subpatch_size, subpatch_size, CV_8UC1);

  for (int i : subpatch_index_mat(y_subpatch, x_subpatch))
  {
    if (patches[i].center_subpatch().y < y_subpatch)
    {
      active_mat |= patches[i].get_active_pixel(y_subpatch, x_subpatch);
    }
  }

  for (int x = 0; x < active_mat.cols; ++x)
  {
    int y = 0;
    for (; y < active_mat.rows; ++y)
    {
      if (active_mat.at<uint8_t>(y, x) != 1)
      {
        break;
      }
    }
    coordinates.push_back(y);
  }

  return coordinates;
}

std::vector<int> MergePatch::get_cut_coordinates_vert(const std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat, int y_subpatch, int x_subpatch, int subpatch_size)
{
  std::vector<int> coordinates;
  cv::Mat active_mat = cv::Mat::zeros(subpatch_size, subpatch_size, CV_8UC1);

  for (int i : subpatch_index_mat(y_subpatch, x_subpatch))
  {
    if (patches[i].center_subpatch().x < x_subpatch)
    {
      active_mat |= patches[i].get_active_pixel(y_subpatch, x_subpatch);
    }
  }

  for (int y = 0; y < active_mat.rows; ++y)
  {
    int x = 0;
    for (; x < active_mat.cols; ++x)
    {
      if (active_mat.at<uint8_t>(y, x) != 1)
      {
        break;
      }
    }
    coordinates.push_back(x);
  }

  /*
  cv::Mat active_mat_debug;
  cv::resize(active_mat, active_mat_debug, cv::Size(), 4.0, 4.0, cv::INTER_NEAREST);
  cv::imshow("active_mat_debug", 255*active_mat_debug);
  cv::waitKey(0);
  */

  return coordinates;
}

void MergePatch::translate_all_curves(std::vector<MergePatch>& patches, const cv::Point2d delta)
{
  for (MergePatch& patch : patches)
  {
    for (BezierCurve& curve : patch.m_curves_left)
    {
      curve += delta;
    }
    for (BezierCurve& curve : patch.m_curves_right)
    {
      curve += delta;
    }
    for (BezierCurve& curve : patch.m_curves_top)
    {
      curve += delta;
    }
    for (BezierCurve& curve : patch.m_curves_bottom)
    {
      curve += delta;
    }
  }
}

static void put_front(std::vector<BezierCurve>& curves, const BezierCurve& b1, const BezierCurve& b2)
{
  curves.front() = b2;
  curves[1].set_control_point(0, b2.control_point(3));
  curves.insert(curves.begin(), b1);
}

static void put_back(std::vector<BezierCurve>& curves, const BezierCurve& b1, const BezierCurve& b2)
{
  curves.back() = b1;
  curves[curves.size() - 2].set_control_point(3, curves.back().control_point(0));
  curves.push_back(b2);
}

static void get_neighboring_patches(std::vector<MergePatch>& patches, const mat<std::vector<int>>& indices, int y_subpatch, int x_subpatch, MergePatch*& patch_tl, MergePatch*& patch_tr, MergePatch*& patch_bl, MergePatch*& patch_br)
{
  if (indices(y_subpatch, x_subpatch).size() != 4)
  {
    throw(std::invalid_argument("get_neighboring_patches only accepts corners with 4 patches."));
  }

  for (int i : indices(y_subpatch, x_subpatch))
  {
    const cv::Point p_center = patches[i].center_subpatch();
    if (p_center.x < x_subpatch)
    {
      if (p_center.y < y_subpatch)
      {
        patch_tl = &patches[i];
      }
      else
      {
        patch_bl = &patches[i];
      }
    }
    else
    {
      if (p_center.y < y_subpatch)
      {
        patch_tr = &patches[i];
      }
      else
      {
        patch_br = &patches[i];
      }
    }
  }
}

void MergePatch::fix_curve_corners(std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat)
{
  const double corner_threshold = 1.0;
  MergePatch* patch_tl;
  MergePatch* patch_tr;
  MergePatch* patch_bl;
  MergePatch* patch_br;

  const int rows = subpatch_index_mat.height();
  const int cols = subpatch_index_mat.width();
  
  mat<uint8_t> corner_handled(rows, cols, 0);

  for (MergePatch& patch : patches)
  {
    const int x_start = patch.region_subpatch.x;
    const int y_start = patch.region_subpatch.y;

    const int x_end = x_start + patch.region_subpatch.width - 1;
    const int y_end = y_start + patch.region_subpatch.height - 1;

    // Check top left corner
    if (x_start > 0 && y_start > 0 && !corner_handled(y_start, x_start))
    {
      BezierCurve b1 = patch.m_curves_top.front();
      BezierCurve b2 = patch.m_curves_left.front();

      const cv::Point2d p1 = b1.eval(1.0);
      const cv::Point2d p2 = b2.eval(1.0);

      if (cv::norm(p1 - p2) < corner_threshold)
      {
        std::pair<double, double> t_center = BezierCurve::intersect_curves_cubic(b1, b2);
        if (t_center.first >= 0.0)
        {
          get_neighboring_patches(patches, subpatch_index_mat, y_start, x_start, patch_tl, patch_tr, patch_bl, patch_br);

          std::pair<BezierCurve, BezierCurve> b1_split = b1.split_cubic(t_center.first);
          std::pair<BezierCurve, BezierCurve> b2_split = b2.split_cubic(t_center.second);
          BezierCurve b_average = BezierCurve::average(b1_split.second, b2_split.second);

          put_front(patch_tr->m_curves_bottom, b1_split.first, b_average);
          put_back(patch_tr->m_curves_left, b2_split.first, b_average);

          put_back(patch_bl->m_curves_top, b1_split.first, b_average);
          put_front(patch_bl->m_curves_right, b2_split.first, b_average);

          patch_tl->m_curves_bottom.back() = b1_split.first;
          patch_tl->m_curves_right.back() = b2_split.first;

          patch.m_curves_top.erase(patch.m_curves_top.begin());
          patch.m_curves_left.erase(patch.m_curves_left.begin());
          patch.m_curves_top.front().set_control_point(0, b_average.control_point(3));
          patch.m_curves_left.front().set_control_point(0, b_average.control_point(3));

          corner_handled(y_start, x_start) = 1;
        }
      }
    }

    // Check top right corner
    if (x_end < cols-1 && y_start > 0 && !corner_handled(y_start, x_end))
    {
      BezierCurve b1 = patch.m_curves_top.back();
      BezierCurve b2 = patch.m_curves_right.front();

      const cv::Point2d p1 = b1.eval(0.0);
      const cv::Point2d p2 = b2.eval(1.0);

      if (cv::norm(p1 - p2) < corner_threshold)
      {
        std::pair<double, double> t_center = BezierCurve::intersect_curves_cubic(b1, b2);
        if (t_center.first >= 0.0)
        {
          get_neighboring_patches(patches, subpatch_index_mat, y_start, x_end, patch_tl, patch_tr, patch_bl, patch_br);

          std::pair<BezierCurve, BezierCurve> b1_split = b1.split_cubic(t_center.first);
          std::pair<BezierCurve, BezierCurve> b2_split = b2.split_cubic(t_center.second);
          BezierCurve b_average = BezierCurve::average(b1_split.first, b2_split.second.reversed());

          put_back(patch_tl->m_curves_bottom, b_average, b1_split.second);
          put_back(patch_tl->m_curves_right, b2_split.first, b_average.reversed());

          put_front(patch_br->m_curves_top, b_average, b1_split.second);
          put_front(patch_br->m_curves_left, b2_split.first, b_average.reversed());

          patch_tr->m_curves_bottom.front() = b1_split.second;
          patch_tr->m_curves_left.back() = b2_split.first;

          patch.m_curves_top.erase(patch.m_curves_top.end() - 1);
          patch.m_curves_right.erase(patch.m_curves_right.begin());
          patch.m_curves_top.back().set_control_point(3, b_average.control_point(0));
          patch.m_curves_right.front().set_control_point(0, b_average.control_point(0));

          corner_handled(y_start, x_end) = 1;
        }
      }
    }

    // Check bottom left corner
    if (x_start > 0 && y_end < rows-1 && !corner_handled(y_end, x_start))
    {
      BezierCurve b1 = patch.m_curves_bottom.front();
      BezierCurve b2 = patch.m_curves_left.back();

      const cv::Point2d p1 = b1.eval(1.0);
      const cv::Point2d p2 = b2.eval(0.0);

      if (cv::norm(p1 - p2) < corner_threshold)
      {
        std::pair<double, double> t_center = BezierCurve::intersect_curves_cubic(b1, b2);
        if (t_center.first >= 0.0)
        {
          get_neighboring_patches(patches, subpatch_index_mat, y_end, x_start, patch_tl, patch_tr, patch_bl, patch_br);

          std::pair<BezierCurve, BezierCurve> b1_split = b1.split_cubic(t_center.first);
          std::pair<BezierCurve, BezierCurve> b2_split = b2.split_cubic(t_center.second);
          BezierCurve b_average = BezierCurve::average(b1_split.second, b2_split.first.reversed());

          put_front(patch_br->m_curves_top, b1_split.first, b_average);
          put_front(patch_br->m_curves_left, b_average.reversed(), b2_split.second);

          put_back(patch_tl->m_curves_bottom, b1_split.first, b_average);
          put_back(patch_tl->m_curves_right, b_average.reversed(), b2_split.second);

          patch_bl->m_curves_top.back() = b1_split.first;
          patch_bl->m_curves_right.front() = b2_split.second;

          patch.m_curves_bottom.erase(patch.m_curves_bottom.begin());
          patch.m_curves_left.erase(patch.m_curves_left.end() - 1);
          patch.m_curves_bottom.front().set_control_point(0, b_average.control_point(3));
          patch.m_curves_left.back().set_control_point(3, b_average.control_point(3));

          corner_handled(y_end, x_start) = 1;
        }
      }
    }

    // Check bottom right corner
    if (x_end < cols-1 && y_end < rows-1 && !corner_handled(y_end, x_end))
    {
      BezierCurve b1 = patch.m_curves_bottom.back();
      BezierCurve b2 = patch.m_curves_right.back();

      const cv::Point2d p1 = b1.eval(0.0);
      const cv::Point2d p2 = b2.eval(0.0);

      if (cv::norm(p1 - p2) < corner_threshold)
      {
        std::pair<double, double> t_center = BezierCurve::intersect_curves_cubic(b1, b2);
        if (t_center.first >= 0.0)
        {
          get_neighboring_patches(patches, subpatch_index_mat, y_end, x_end, patch_tl, patch_tr, patch_bl, patch_br);

          std::pair<BezierCurve, BezierCurve> b1_split = b1.split_cubic(t_center.first);
          std::pair<BezierCurve, BezierCurve> b2_split = b2.split_cubic(t_center.second);
          BezierCurve b_average = BezierCurve::average(b1_split.first, b2_split.first);

          put_back(patch_bl->m_curves_top, b_average, b1_split.second);
          put_front(patch_bl->m_curves_right, b_average, b2_split.second);

          put_front(patch_tr->m_curves_bottom, b_average, b1_split.second);
          put_back(patch_tr->m_curves_left, b_average, b2_split.second);

          patch_br->m_curves_top.front() = b1_split.second;
          patch_br->m_curves_left.front() = b2_split.second;

          patch.m_curves_bottom.erase(patch.m_curves_bottom.end() - 1);
          patch.m_curves_right.erase(patch.m_curves_right.end() - 1);
          patch.m_curves_bottom.back().set_control_point(3, b_average.control_point(0));
          patch.m_curves_right.back().set_control_point(3, b_average.control_point(0));

          corner_handled(y_end, x_end) = 1;
        }
      }
    }
  }
}

cv::Mat MergePatch::draw_rect(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, double scale, int subpatch_size, bool draw_boundaries)
{
  if (patches.empty())
  {
    return cv::Mat();
  }

  cv::Rect bbox(patches[0].anchor_target, patches[0].size);
  for (const MergePatch& patch : patches)
  {
    bbox = cv::boundingRect(std::vector<cv::Point>({bbox.tl(), bbox.br() - cv::Point(1, 1), patch.anchor_target, patch.anchor_target + cv::Point(patch.size)}));
  }

  const int boundary_size = subpatch_size / 2;
  const float dx = static_cast<float>(1.0f / scale);
  const float dy = static_cast<float>(1.0f / scale);

  cv::Mat image(static_cast<int>((bbox.height-1) * scale), static_cast<int>((bbox.width-1) * scale), CV_8UC3, cv::Scalar(255, 255, 255));
  cv::Mat boundary_mask;
  if (draw_boundaries)
  {
    boundary_mask = cv::Mat::zeros(image.size(), CV_32FC1);
  }

  for (const MergePatch& patch : patches)
  {
    float x_start = static_cast<float>(patch.anchor_target.x);
    float x_end = static_cast<float>(patch.anchor_target.x + patch.size.width);
    float y_start = static_cast<float>(patch.anchor_target.y);
    float y_end = static_cast<float>(patch.anchor_target.y + patch.size.height);

    if (patch.anchor_target.x > bbox.x)
    {
      x_start += static_cast<float>(boundary_size);
    }

    if (patch.anchor_target.x + patch.size.width < bbox.x + bbox.width - 1)
    {
      x_end -= static_cast<float>(boundary_size);
    }

    if (patch.anchor_target.y > bbox.y)
    {
      y_start += static_cast<float>(boundary_size);
    }

    if (patch.anchor_target.y + patch.size.height < bbox.y + bbox.height - 1)
    {
      y_end -= static_cast<float>(boundary_size);
    }

    cv::Point2f p_target;
    for (p_target.y = y_start; p_target.y < y_end; p_target.y += dy)
    {
      for (p_target.x = x_start; p_target.x < x_end; p_target.x += dx)
      {
        cv::Point p_target_scale((p_target - cv::Point2f(bbox.tl())) * scale);
        cv::Point2f p_source = p_target - cv::Point2f(patch.anchor_target) + cv::Point2f(patch.anchor_source);
        p_source = AffineTransformation::transform(patch.transformation_source_inv, p_source);
        image.at<cv::Vec3b>(p_target_scale) = textures[patch.source_index].interpolate_texture(p_source);
      }
    }

    if (draw_boundaries)
    {
      cv::Point p_target_scale(static_cast<int>((x_start - bbox.x)*scale), static_cast<int>((y_start - bbox.y) * scale));
      cv::Size size_scale(static_cast<int>((x_end-x_start)*scale+1), static_cast<int>((y_end-y_start)*scale+1));
      cv::rectangle(boundary_mask, cv::Rect(p_target_scale, size_scale), cv::Scalar(1.0f));
    }
  }

  if (draw_boundaries)
  {
    const cv::Vec3b dark_brown(14, 29, 43);
    cv::GaussianBlur(boundary_mask, boundary_mask, cv::Size(5, 5), 0.0);

    for (int y = 0; y < image.rows; ++y)
    {
      cv::Vec3b* ptr_image = reinterpret_cast<cv::Vec3b*>(image.ptr(y));
      const float* ptr_mask = reinterpret_cast<const float*>(boundary_mask.ptr(y));
      for (int x = 0; x < image.cols; ++x)
      {
        ptr_image[x] = ptr_mask[x] * dark_brown + (1.0f - ptr_mask[x]) * ptr_image[x];
      }
    }
  }

  return image;
}

cv::Mat MergePatch::draw_rect_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale, int subpatch_size, bool draw_boundaries)
{
  std::vector<Texture> textures_fullres;
  for (const Texture& t : textures)
  {
    textures_fullres.emplace_back(base_path / t.filename, t.dpi, 1.0);
  }
  return draw_rect_fullres(patches, textures, textures_fullres, scale, subpatch_size, draw_boundaries);
}

cv::Mat MergePatch::draw_rect_fullres(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale, int subpatch_size, bool draw_boundaries)
{
  if (patches.empty())
  {
    return cv::Mat();
  }

  cv::Rect bbox(patches[0].anchor_target, patches[0].size);
  for (const MergePatch& patch : patches)
  {
    bbox = cv::boundingRect(std::vector<cv::Point>({bbox.tl(), bbox.br() - cv::Point(1, 1), patch.anchor_target, patch.anchor_target + cv::Point(patch.size)}));
  }

  const int boundary_size = subpatch_size / 2;
  const float dx = static_cast<float>(1.0f / scale);
  const float dy = static_cast<float>(1.0f / scale);

  cv::Mat image(static_cast<int>((bbox.height-1) * scale), static_cast<int>((bbox.width-1) * scale), CV_8UC3, cv::Scalar(255, 255, 255));
  cv::Mat boundary_mask;
  if (draw_boundaries)
  {
    boundary_mask = cv::Mat::zeros(image.size(), CV_32FC1);
  }

  for (const MergePatch& patch : patches)
  {
    const double scale_source = 1.0 / textures[patch.source_index].scale;
    const cv::Mat T_source = AffineTransformation::concat(AffineTransformation::T_scale(scale_source, scale_source), patch.transformation_source_inv);

    float x_start = static_cast<float>(patch.anchor_target.x);
    float x_end = static_cast<float>(patch.anchor_target.x + patch.size.width);
    float y_start = static_cast<float>(patch.anchor_target.y);
    float y_end = static_cast<float>(patch.anchor_target.y + patch.size.height);

    if (patch.anchor_target.x > bbox.x)
    {
      x_start += static_cast<float>(boundary_size);
    }

    if (patch.anchor_target.x + patch.size.width < bbox.x + bbox.width - 1)
    {
      x_end -= static_cast<float>(boundary_size);
    }

    if (patch.anchor_target.y > bbox.y)
    {
      y_start += static_cast<float>(boundary_size);
    }

    if (patch.anchor_target.y + patch.size.height < bbox.y + bbox.height - 1)
    {
      y_end -= static_cast<float>(boundary_size);
    }

    cv::Point2f p_target;
    for (p_target.y = y_start; p_target.y < y_end; p_target.y += dy)
    {
      for (p_target.x = x_start; p_target.x < x_end; p_target.x += dx)
      {
        cv::Point p_target_scale((p_target - cv::Point2f(bbox.tl())) * scale);
        cv::Point2f p_source = p_target - cv::Point2f(patch.anchor_target) + cv::Point2f(patch.anchor_source);
        p_source = AffineTransformation::transform(T_source, p_source);
        image.at<cv::Vec3b>(p_target_scale) = textures_fullres[patch.source_index].interpolate_texture(p_source);
      }
    }

    if (draw_boundaries)
    {
      cv::Point p_target_scale(static_cast<int>((x_start - bbox.x)*scale), static_cast<int>((y_start - bbox.y) * scale));
      cv::Size size_scale(static_cast<int>((x_end-x_start)*scale+1), static_cast<int>((y_end-y_start)*scale+1));
      cv::rectangle(boundary_mask, cv::Rect(p_target_scale, size_scale), cv::Scalar(1.0f));
    }
  }

  if (draw_boundaries)
  {
    const cv::Vec3b dark_brown(14, 29, 43);
    cv::GaussianBlur(boundary_mask, boundary_mask, cv::Size(5, 5), 0.0);

    for (int y = 0; y < image.rows; ++y)
    {
      cv::Vec3b* ptr_image = reinterpret_cast<cv::Vec3b*>(image.ptr(y));
      const float* ptr_mask = reinterpret_cast<const float*>(boundary_mask.ptr(y));
      for (int x = 0; x < image.cols; ++x)
      {
        ptr_image[x] = ptr_mask[x] * dark_brown + (1.0f - ptr_mask[x]) * ptr_image[x];
      }
    }
  }

  return image;
}

MergePatch::MergePatch(cv::Rect region_subpatch, cv::Point anchor_target, int subpatch_size) :
  region_subpatch(region_subpatch),
  subpatch_size(subpatch_size)
{
  size.height = region_subpatch.height * subpatch_size;
  size.width = region_subpatch.width * subpatch_size;
  error = cv::Mat(size, CV_32FC1, std::numeric_limits<float>::infinity());
  this->anchor_target = anchor_target;
  pos_patch.x = -1;
  pos_patch.y = -1;
  active_pixel =  cv::Mat::zeros(size, CV_8UC1);
  target_id = -1;
}

void MergePatch::trim_bezier_curves()
{
  try
  {
    m_curves_top.insert(m_curves_top.begin(), m_curves_top.front().extend_curve_left());
    m_curves_top.push_back(m_curves_top.back().extend_curve_right());

    m_curves_left.insert(m_curves_left.begin(), m_curves_left.front().extend_curve_left());
    m_curves_left.push_back(m_curves_left.back().extend_curve_right());

    m_curves_bottom.insert(m_curves_bottom.begin(), m_curves_bottom.front().extend_curve_left());
    m_curves_bottom.push_back(m_curves_bottom.back().extend_curve_right());

    m_curves_right.insert(m_curves_right.begin(), m_curves_right.front().extend_curve_left());
    m_curves_right.push_back(m_curves_right.back().extend_curve_right());

    BezierCurve::trim_curves_cubic(m_curves_top, true, m_curves_left, true);
    BezierCurve::trim_curves_cubic(m_curves_top, false, m_curves_right, true);
    BezierCurve::trim_curves_cubic(m_curves_bottom, true, m_curves_left, false);
    BezierCurve::trim_curves_cubic(m_curves_bottom, false, m_curves_right, false);

    m_curves_top.erase(std::remove_if(m_curves_top.begin(), m_curves_top.end(), std::mem_fn(&BezierCurve::is_empty_curve)), m_curves_top.end());
    m_curves_bottom.erase(std::remove_if(m_curves_bottom.begin(), m_curves_bottom.end(), std::mem_fn(&BezierCurve::is_empty_curve)), m_curves_bottom.end());
    m_curves_left.erase(std::remove_if(m_curves_left.begin(), m_curves_left.end(), std::mem_fn(&BezierCurve::is_empty_curve)), m_curves_left.end());
    m_curves_right.erase(std::remove_if(m_curves_right.begin(), m_curves_right.end(), std::mem_fn(&BezierCurve::is_empty_curve)), m_curves_right.end());
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl
      << pos_patch << std::endl;
    cv::Mat image = draw(16.0);
    cv::imshow("Image", image);
    cv::waitKey(0);
    std::exit(-1);
  }
}

cv::Mat MergePatch::draw(double scale) const
{
  const float dx = static_cast<float>(1.0f / scale);
  const float dy = static_cast<float>(1.0f / scale);

  cv::Mat image = cv::Mat::zeros(static_cast<int>((size.height+8) * scale), static_cast<int>((size.width+8) * scale), CV_8UC3);

  cv::line(image, cv::Point(static_cast<int>(4.0*scale), 0), cv::Point(static_cast<int>(4.0*scale), image.rows - 1), cv::Scalar(255));
  cv::line(image, cv::Point(image.cols - static_cast<int>(4.0*scale), 0), cv::Point(image.cols - static_cast<int>(4.0*scale), image.rows - 1), cv::Scalar(255));
  cv::line(image, cv::Point(0, static_cast<int>(4.0*scale)), cv::Point(image.cols - 1, static_cast<int>(4.0*scale)), cv::Scalar(255));
  cv::line(image, cv::Point(0, image.rows - static_cast<int>(4.0*scale)), cv::Point(image.cols - 1, image.rows - static_cast<int>(4.0*scale)), cv::Scalar(255));

  for (const BezierCurve& curve : m_curves_left)
  {
    (curve - anchor_target + cv::Point2d(4.0, 4.0)).draw(image, subpatch_size, scale);
  }
  for (const BezierCurve& curve : m_curves_right)
  {
    (curve - anchor_target + cv::Point2d(4.0, 4.0)).draw(image, subpatch_size, scale);
  }
  for (const BezierCurve& curve : m_curves_top)
  {
    (curve - anchor_target + cv::Point2d(4.0, 4.0)).draw(image, subpatch_size, scale);
  }
  for (const BezierCurve& curve : m_curves_bottom)
  {
    (curve - anchor_target + cv::Point2d(4.0, 4.0)).draw(image, subpatch_size, scale);
  }

  /*
  std::vector<float> x_start, x_end, y_start, y_end, x_vals, y_vals;
  get_patch_boundaries(*this, dx, dy, x_vals, y_vals, x_start, x_end, y_start, y_end);
  for (int i = 0; i < static_cast<int>(y_vals.size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(x_vals.size()); ++j)
    {
      if (x_vals[j] >= x_start[i] - 1 && x_vals[j] <= x_end[i] + 1 &&
        y_vals[i] >= y_start[j] - 1 && y_vals[i] <= y_end[j] + 1)
      {
        float x = x_vals[j];
        float y = y_vals[i];
        int x_target = static_cast<int>((x - anchor_target.x) * scale);
        int y_target = static_cast<int>((y - anchor_target.y) * scale);
        image.at<cv::Vec3b>(y_target, x_target) = cv::Vec3b::all(255);
      }
    }
  }
  */

  return image;
}

cv::Mat MergePatch::draw_half_rect_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale, int subpatch_size, bool draw_boundaries)
{
  std::vector<Texture> textures_fullres;
  for (const Texture& t : textures)
  {
    textures_fullres.emplace_back(base_path / t.filename, t.dpi, 1.0);
  }
  return draw_half_rect_fullres(patches, textures, textures_fullres, scale, subpatch_size, draw_boundaries);
}

cv::Mat MergePatch::draw_half_rect_fullres(std::vector<MergePatch> patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale, int subpatch_size, bool draw_boundaries)
{
  for (MergePatch& p : patches)
  {
    const float x_left = static_cast<float>(p.anchor_target.x + subpatch_size / 2);
    const float x_right = static_cast<float>(p.anchor_target.x + p.size.width - subpatch_size / 2);
    const float y_top = static_cast<float>(p.anchor_target.y);
    const float y_bottom = static_cast<float>(p.anchor_target.y + p.size.width);

    const float y_delta = 0.25f * (y_bottom - y_top);

    std::vector<BezierCurve> curves_left, curves_right;
    for (int i = 0; i < 4; ++i)
    {
      curves_left.emplace_back(cv::Point2f(x_left, y_top + i * y_delta),
                               cv::Point2f(x_left, y_top + (i+1) * y_delta), 3);

      curves_right.emplace_back(cv::Point2f(x_right, y_top + i * y_delta),
                                cv::Point2f(x_right, y_top + (i+1) * y_delta), 3);
    }

    p.m_curves_left.swap(curves_left);
    p.m_curves_right.swap(curves_right);
    p.trim_bezier_curves();    
  }

  return draw_fullres(patches, textures, textures_fullres, scale, draw_boundaries);
}
