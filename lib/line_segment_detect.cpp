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

#include "line_segment_detect.hpp"

std::vector<cv::Vec4f> line_segment_detect(cv::Mat texture)
{
  cv::Mat texture_gray;
  std::vector<cv::Vec4f> lines;

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

  texture_gray.convertTo(texture_gray, CV_8UC1, 255.0);

  cv::Ptr<cv::LineSegmentDetector> detector = cv::createLineSegmentDetector(cv::LSD_REFINE_STD, 0.9, 0.6, 1.25, 22.5);
  detector->detect(texture_gray, lines);

  return lines;
}

void draw_line_segments(const std::string& window_name, cv::Mat texture, const std::vector<cv::Vec4f>& line_segments)
{
  cv::Mat texture_uchar;
  texture.convertTo(texture_uchar, CV_8UC3, 255.0);

  cv::Ptr<cv::LineSegmentDetector> detector = cv::createLineSegmentDetector(cv::LSD_REFINE_STD);
  detector->drawSegments(texture_uchar, line_segments);

  cv::imshow(window_name, texture_uchar);
  cv::waitKey(1);
}

cv::Mat line_segments_to_binary_image(const std::vector<cv::Vec4f>& line_segments, int height, int width)
{
  cv::Mat binary_image = cv::Mat::zeros(height, width, CV_32FC1);

  for (const cv::Vec4f& line : line_segments)
  {
    const cv::Point p1(static_cast<int>(line[0]), static_cast<int>(line[1]));
    const cv::Point p2(static_cast<int>(line[2]), static_cast<int>(line[3]));

    cv::line(binary_image, p1, p2, 1.0f, 1, 8);
  }

  return binary_image;
}