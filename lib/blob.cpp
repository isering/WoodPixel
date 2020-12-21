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

#include "blob.hpp"

Blob Blob::detect(cv::Mat image, const cv::Point &p_start, unsigned char color_fill)
{
  Blob blob;
  unsigned char color_start = image.at<unsigned char>(p_start);

  std::deque<cv::Point> queue({p_start});
  while (!queue.empty())
  {
    const cv::Point p = queue.front();
    queue.pop_front();

    unsigned char &color = image.at<unsigned char>(p);

    if (color == color_start)
    {
      blob.m_points.push_back(p);

      if (p.x > 0)
      {
        queue.emplace_back(p.x - 1, p.y);
      }

      if (p.x < image.cols - 1)
      {
        queue.emplace_back(p.x + 1, p.y);
      }

      if (p.y > 0)
      {
        queue.emplace_back(p.x, p.y - 1);
      }

      if (p.y < image.rows - 1)
      {
        queue.emplace_back(p.x, p.y + 1);
      }

      color = color_fill;
    }
  }

  return blob;
}

std::vector<Blob> Blob::detect(cv::Mat image, unsigned char fg)
{
  cv::Mat blob_image = image.clone();

  std::vector<Blob> blobs;
  for (int y = 0; y < blob_image.rows; ++y)
  {
    const unsigned char *ptr = blob_image.ptr(y);
    for (int x = 0; x < blob_image.cols; ++x)
    {
      if (ptr[x] != fg)
      {
        blobs.push_back(Blob::detect(blob_image, cv::Point(x, y), fg));
      }
    }
  }

  return blobs;
}

std::vector<BezierCurve> Blob::contours(cv::Size size) const
{
  cv::Mat mask = cv::Mat::zeros(size, CV_8UC1);

  for (const cv::Point &p : m_points)
  {
    mask.at<unsigned char>(p) = 255;
  }

  std::vector<std::vector<cv::Point>> contours;
  cv::findContours(mask, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

  std::vector<BezierCurve> curves;
  for (const std::vector<cv::Point> &contour : contours)
  {
    for (size_t i = 1; i < contours.size(); ++i)
    {
      curves.emplace_back(contour[i - 1], contour[i], 3);
    }
    if (contour.size() > 1)
    {
      curves.emplace_back(contour.front(), contour.back(), 3);
    }
  }

  return curves;
}

void Blob::draw_contour(cv::Mat mask) const
{
  cv::Mat mask_temp = cv::Mat::zeros(mask.size(), CV_8UC1);
  for (const cv::Point &p : m_points)
  {
    mask_temp.at<unsigned char>(p) = 255;
  }

  std::vector<std::vector<cv::Point>> contours;
  cv::findContours(mask_temp, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

  for (const std::vector<cv::Point> &contour : contours)
  {
    for (const cv::Point &p : contour)
    {
      mask.at<unsigned char>(p) = 255;
    }
  }
}