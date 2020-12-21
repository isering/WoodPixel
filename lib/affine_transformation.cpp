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

cv::Mat AffineTransformation::concat(const cv::Mat& lhs, const cv::Mat& rhs)
{
  cv::Mat result(2, 3, lhs.type());

  const cv::Mat rot_lhs = lhs(cv::Range(0, 2), cv::Range(0, 2));
  const cv::Mat rot_rhs = rhs(cv::Range(0, 2), cv::Range(0, 2));
  cv::Mat rot_result = result(cv::Range(0, 2), cv::Range(0, 2));

  const cv::Mat t_lhs = lhs(cv::Range(0, 2), cv::Range(2, 3));
  const cv::Mat t_rhs = rhs(cv::Range(0, 2), cv::Range(2, 3));
  cv::Mat t_result = result(cv::Range(0, 2), cv::Range(2, 3));

  rot_result = rot_lhs * rot_rhs;
  t_result = rot_lhs * t_rhs + t_lhs;

  return result;
}

cv::Mat AffineTransformation::concat(std::initializer_list<cv::Mat> mat_list)
{
  cv::Mat M = cv::Mat::eye(2, 3, CV_64FC1);
  for (auto iter = mat_list.begin(); iter != mat_list.end(); ++iter)
  {
    M = concat(M, *iter);
  }
  return M;
}

cv::Mat AffineTransformation::fit(const std::vector<cv::Point2d>& points_from, const std::vector<cv::Point2d>& points_to)
{
  const std::vector<cv::Point2f> points_from_float(points_from.begin(), points_from.end());
  const std::vector<cv::Point2f> points_to_float(points_to.begin(), points_to.end());

  if (points_from.size() < 3 || points_to.size() < 3)
  {
    throw(std::invalid_argument("AffineTransformation::fit expects at least three points in points_from, points_to."));
  }

  return cv::getAffineTransform(points_from_float, points_to_float);
}