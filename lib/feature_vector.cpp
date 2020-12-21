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

#include "feature_vector.hpp"

#include <boost/format.hpp>

#include "print.hpp"
#include "timer.hpp"

namespace fs = boost::filesystem;

void FeatureVector::render() const
{
  const cv::Size feature_size = feature_vector.front().size();
  const int num_channels = static_cast<int>(feature_vector.size());
  cv::Mat image_out = cv::Mat::zeros(feature_size, CV_8UC1);

  enum key_press
  {
    left = 2424832,
    right = 2555904,
    up = 2490368,
    down = 2621440,
    space = 32
  };
  
  int cur_index = 0;  
  bool run_visualization = true;
  while (run_visualization)
  {
    cv::Mat im_out;
    feature_vector[cur_index].convertTo(im_out, CV_8UC1, 1.0/255);
    cv::imshow("Response", im_out);
    int key = cv::waitKeyEx();

    switch (key)
    {
    case left:
    case down:
      --cur_index;
      if (cur_index < 0)
      {
        cur_index = num_channels - 1;
      }
      break;
    case right:
    case up:
      ++cur_index;
      if (cur_index >= num_channels)
      {
        cur_index  = 0;
      }
      break;
    case space:
      run_visualization = false;
    }
  }
}

void FeatureVector::save(const boost::filesystem::path& path) const
{
  const int num_channels = static_cast<int>(feature_vector.size());
  cv::Mat image_out;

  if (!fs::exists(path))
  {
    fs::create_directories(path);
  }

  if (!fs::exists(path) || !fs::is_directory(path))
  {
    throw(std::runtime_error("Unable to create output path: " + path.string()));
  }

  for (int i = 0; i < num_channels; ++i)
  {
    feature_vector[i].convertTo(image_out, CV_8UC1, 1.0 / 255.0);
    cv::imwrite((path / (boost::format("%04d.png") % i).str()).string(), image_out);
  }
}

void FeatureVector::render(const std::vector<FeatureVector>& features)
{
  const int num_channels = static_cast<int>(features.front().num_channels());

  int cols_out = 0;
  int rows_out = 0;
  double global_min = std::numeric_limits<double>::max();
  double global_max = -std::numeric_limits<double>::max();
  for (const FeatureVector& f : features)
  {
    cols_out += f[0].cols;
    rows_out = std::max(rows_out, f[0].rows);
    for (int i = 0; i < num_channels; ++i)
    {
      double local_min, local_max;
      cv::minMaxLoc(f[i], &local_min, &local_max);
      global_min = std::min(local_min, global_min);
      global_max = std::max(local_max, global_max);
    }
  }
  //global_min = 0.0;

  std::cout << global_max << std::endl;

  enum key_press
  {
    left = 2424832,
    right = 2555904,
    up = 2490368,
    down = 2621440,
    space = 32
  };

  cv::Mat image_out = cv::Mat::zeros(rows_out, cols_out, CV_16UC1);

  int cur_index = 0;
  bool run_visualization = true;
  while (run_visualization)
  {
    std::vector<double> response_min, response_max;
    for (const FeatureVector& f : features)
    {
      double local_min, local_max;
      cv::minMaxLoc(f[cur_index], &local_min, &local_max);
      response_min.push_back(local_min);
      response_max.push_back(local_max);
    }

    //const double global_min = *std::min_element(response_min.begin(), response_min.end());
    //const double global_max = *std::max_element(response_max.begin(), response_max.end());

    int current_col = 0;
    for (size_t i = 0; i < features.size(); ++i)
    {
      const double scale_min = 65535.0 * (response_min[i] - global_min) / (global_max - global_min);
      const double scale_max = 65535.0 * (response_max[i] - global_min) / (global_max - global_min);

      cv::Mat response;
      cv::normalize(features[i][cur_index], response, scale_max, scale_min, cv::NORM_MINMAX);

      response.copyTo(image_out(cv::Rect(current_col, 0, features[i][cur_index].cols, features[i][cur_index].rows)));
      current_col += features[i][cur_index].cols;
    }

    cv::imshow("Response", image_out);
    int key = cv::waitKeyEx();

    switch (key)
    {
    case left:
    case down:
      --cur_index;
      if (cur_index < 0)
      {
        cur_index = num_channels - 1;
      }
      break;
    case right:
    case up:
      ++cur_index;
      if (cur_index >= num_channels)
      {
        cur_index = 0;
      }
      break;
    case space:
      run_visualization = false;
    }
  }
}

float FeatureVector::dist(cv::Point p_target, const FeatureVector& rhs, cv::Point p_rhs) const
{
  float dist = 0.0f;
  for (size_t i = 0; i < feature_vector.size(); ++i)
  {
    const float diff = (static_cast<float>(feature_vector[i].at<uint16_t>(p_target)) - static_cast<float>(rhs.feature_vector[i].at<uint16_t>(p_rhs))) / 65535.0f;
    dist += diff * diff;
  }
  return std::sqrt(dist);
}

cv::Mat FeatureVector::dist_sqr_mat(const FeatureVector& rhs) const
{
  if (size() != rhs.size())
  {
    throw(std::invalid_argument("FeatureVector sizes differ."));
  }

  cv::Mat channel_lhs, channel_rhs, channel_diff;

  cv::Mat dist = cv::Mat::zeros(size(), CV_32FC1);

  for (int i = 0; i < num_channels(); ++i)
  {
    feature_vector[i].convertTo(channel_lhs, CV_32FC1, 1.0 / 65535.0);
    rhs.feature_vector[i].convertTo(channel_rhs, CV_32FC1, 1.0 / 65535.0);
    channel_diff = channel_lhs - channel_rhs;
    dist += channel_diff.mul(channel_diff);
  }

  return dist;
}

FeatureVector FeatureVector::clone() const
{
  FeatureVector rhs;
  for (const cv::Mat& feature : feature_vector)
  {
    rhs.feature_vector.push_back(feature.clone());
  }
  return rhs;
}

void FeatureVector::downsample_nn(int factor)
{
  for (cv::Mat& feature : feature_vector)
  {
    cv::resize(feature, feature, cv::Size(), 1.0/factor, 1.0/factor, cv::INTER_NEAREST);
  }
}