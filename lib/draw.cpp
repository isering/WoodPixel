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

#include "draw.hpp"

static cv::Scalar hsv_to_bgr(cv::Scalar hsv)
{
  cv::Mat hsv_mat(1, 1, CV_8UC3, hsv);
  cv::Mat rgb_mat;
  cv::cvtColor(hsv_mat, rgb_mat, cv::COLOR_HSV2BGR);
  return cv::Scalar(rgb_mat.data[0], rgb_mat.data[1], rgb_mat.data[2]);
}

static std::vector<cv::Scalar> get_colors(int num_colors)
{
  std::vector<cv::Scalar> colors;
  switch (num_colors)
  {
  case 1:
    colors.emplace_back(0, 0, 0);
    break;
  case 3:
    colors.emplace_back(255, 0, 0);
    colors.emplace_back(0, 255, 0);
    colors.emplace_back(0, 0, 255);
    break;
  default:
    for (int i = 0; i < num_colors; ++i)
    {
      cv::Scalar color_hsv(180.0 * i / num_colors, 255, 255);
      colors.emplace_back(hsv_to_bgr(color_hsv));
    }
  }
  return colors;
}

cv::Mat draw_bars(const std::vector<cv::Mat>& hist)
{
  const int num_colors = static_cast<int>(hist.size());
  const int num_bins = hist.front().rows;
  const int bin_w = 3;
  const int image_h = 512;

  double local_max;
  double global_max = 0.0;
  for (const cv::Mat& mat : hist)
  {
    cv::minMaxLoc(mat, 0, &local_max);
    global_max = std::max(global_max, local_max);
  }

  std::cout << "global_max: " << global_max << std::endl;

  std::vector<cv::Scalar> colors = get_colors(num_colors);

  cv::Mat hist_image(image_h, bin_w * num_colors * num_bins, CV_8UC3, cv::Scalar(255, 255, 255));
  for (int color = 0; color < num_colors; ++color)
  {
    for (int x = 0; x < num_bins; ++x)
    {
      const int x1 = (num_colors * x + color) * bin_w;
      const int x2 = x1 + bin_w - 1;
      const int y = static_cast<int>((hist_image.rows-1) * (1.0 - (hist[color].at<float>(x, 0) / global_max)));

      cv::rectangle(hist_image, cv::Point(x1, image_h-1), cv::Point(x2, y), colors[color], -1);
    }
  }

  return hist_image;
}

cv::Mat draw_lines(const std::vector<cv::Mat>& hist)
{
  const int num_colors = static_cast<int>(hist.size());
  const int num_bins = hist.front().rows;
  const int bin_w = 4;
  const int image_h = 512;

  double local_max;
  double global_max = 0.0;
  for (const cv::Mat& mat : hist)
  {
    cv::minMaxLoc(mat, 0, &local_max);
    global_max = std::max(global_max, local_max);
  }

  std::vector<cv::Scalar> colors = get_colors(num_colors);

  cv::Mat hist_image(image_h, bin_w * num_bins, CV_8UC3, cv::Scalar(255, 255, 255));
  for (int color = 0; color < num_colors; ++color)
  {
    cv::Point p_last(0, static_cast<int>((hist_image.rows-1) * (1.0 - (hist[color].at<float>(0, 0) / global_max))));
    for (int x = 1; x < num_bins; ++x)
    {
      cv::Point p(x * bin_w, static_cast<int>((hist_image.rows-1) * (1.0 - (hist[color].at<float>(x, 0) / global_max))));
      cv::line(hist_image, p_last, p, colors[color], 2);
      p_last = p;
    }
  }

  return hist_image;
}