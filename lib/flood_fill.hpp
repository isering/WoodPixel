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

#ifndef TRLIB_FLOOD_FILL_HPP_
#define TRLIB_FLOOD_FILL_HPP_

#include <stack>

#include <opencv2/opencv.hpp>

#include "bezier_curve.hpp"

template <typename T>
void flood_fill_patch(cv::Mat mask, T mask_val)
{
  std::stack<cv::Point> fill_stack;

  cv::Mat mask_border;
  cv::copyMakeBorder(mask, mask_border, 1, 1, 1, 1, cv::BORDER_CONSTANT, 0);
  cv::floodFill(mask_border, cv::Point(0, 0), mask_val, 0, 0, 0, 4 | cv::FLOODFILL_FIXED_RANGE);

  mask += mask_val - mask - mask_border(cv::Rect(1, 1, mask.cols, mask.rows));
}

template <typename T>
void flood_fill_separating_curve(cv::Mat mask, BezierCurve curve, T val_curve, T val_left, T val_right, T val_mask)
{
  std::stack<cv::Point> fill_stack_left, fill_stack_right;

  curve.draw(mask, val_curve);

  const std::vector<CurveDrawPoint> points = curve.get_draw_points();
  const BezierCurve curve_deriv = curve.deriv(1);

  std::vector<cv::Point2d> normals(points.size());
  for (size_t i = 0; i < points.size(); ++i)
  {
    const cv::Vec2d tangent = cv::normalize(cv::Vec2d(curve_deriv.eval(points[i].t)));
    const cv::Vec2d normal(tangent[1], -tangent[0]);

    const int delta_x = normal[0] >= 0.383 ? 1 : normal[0] >= -0.383 ? 0 : -1;
    const int delta_y = normal[1] >= 0.383 ? 1 : normal[1] >= -0.383 ? 0 : -1;

    fill_stack_left.emplace(points[i].p.x - delta_x, points[i].p.y - delta_y);
    fill_stack_right.emplace(points[i].p.x + delta_x, points[i].p.y + delta_y);
  }

  while (!fill_stack_left.empty())
  {
    const cv::Point p = fill_stack_left.top();
    fill_stack_left.pop();

    if (p.x >= 0 && p.x < mask.cols && p.y >= 0 && p.y < mask.rows)
    {
      T &val = mask.at<T>(p);
      if (val == val_mask)
      {
        val = val_left;
        fill_stack_left.emplace(p.x - 1, p.y);
        fill_stack_left.emplace(p.x + 1, p.y);
        fill_stack_left.emplace(p.x, p.y - 1);
        fill_stack_left.emplace(p.x, p.y + 1);
      }
    }
  }

  while (!fill_stack_right.empty())
  {
    const cv::Point p = fill_stack_right.top();
    fill_stack_right.pop();

    if (p.x >= 0 && p.x < mask.cols && p.y >= 0 && p.y < mask.rows)
    {
      T &val = mask.at<T>(p);
      if (val == val_mask)
      {
        val = val_right;
        fill_stack_right.emplace(p.x - 1, p.y);
        fill_stack_right.emplace(p.x + 1, p.y);
        fill_stack_right.emplace(p.x, p.y - 1);
        fill_stack_right.emplace(p.x, p.y + 1);
      }
    }
  }
}

#endif /* TRLIB_FLOOD_FILL_HPP_ */