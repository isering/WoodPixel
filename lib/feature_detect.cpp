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

#include "feature_detect.hpp"

#include <stdexcept>

std::vector<cv::Point2f> feature_detect(cv::Mat texture, int max_corners, double quality_level, double min_distance)
{
  std::vector<cv::Point2f> corners;
  cv::Mat texture_gray;

  if (texture.channels() == 3)
  {
    cv::cvtColor(texture, texture_gray, cv::COLOR_BGR2GRAY);
  }
  else if (texture.channels() == 1)
  {
    texture_gray = texture;
  }
  else
  {
    throw std::invalid_argument("Input texture must be either 1 or 3 channel.");
  }

  cv::goodFeaturesToTrack(texture_gray, corners, max_corners, quality_level, min_distance);

  return corners;
}

void draw_features(const std::string& window_name, cv::Mat texture, const std::vector<cv::Point2f>& features)
{
  std::vector<cv::KeyPoint> key_points;
  for (cv::Point2f p : features)
  {
    key_points.emplace_back(p, 1.0f);
  }

  cv::Mat texture_uchar;
  texture.convertTo(texture_uchar, CV_8UC3, 255.0);
  cv::drawKeypoints(texture_uchar, key_points, texture_uchar);

  cv::imshow("Keypoints", texture_uchar);
  cv::waitKey(1);
}

cv::Mat features_to_binary_image(const std::vector<cv::Point2f>& features, int height, int width)
{
  cv::Mat binary_image = cv::Mat::zeros(height, width, CV_32FC1);

  for (const cv::Point2f& p : features)
  {
    cv::circle(binary_image, static_cast<cv::Point>(p), 5, 1.0, -1);
  }

  return binary_image;
}