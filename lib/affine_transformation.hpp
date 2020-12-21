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

#ifndef TRLIB_AFFINE_TRANSFORMATION_HPP_
#define TRLIB_AFFINE_TRANSFORMATION_HPP_

#include <boost/math/constants/constants.hpp>
#include <opencv2/opencv.hpp>

namespace AffineTransformation
{

template <typename T, typename Tout = T>
cv::Point_<Tout> transform(cv::Mat A, const cv::Point_<T>& p)
{
  cv::Mat p_transformed(A * cv::Mat(cv::Vec3d(static_cast<double>(p.x), static_cast<double>(p.y), 1.0)));
  return cv::Point_<Tout>(static_cast<Tout>(p_transformed.at<double>(0)), static_cast<Tout>(p_transformed.at<double>(1)));
}

template <typename T, typename Tout = T>
std::vector<cv::Point_<Tout>> transform(cv::Mat A, const std::vector<cv::Point_<T>>& points)
{
  std::vector<cv::Point_<Tout>> points_transformed(points.size());
  for (size_t i = 0; i < points.size(); ++i)
  {
    points_transformed[i] = transform<T, Tout>(A, points[i]);
  }
  return points_transformed;
}

inline double rotation_rad(cv::Mat T)
{
  return std::atan2(T.at<double>(1, 0), T.at<double>(0, 0));
}

inline double rotation_deg(cv::Mat T)
{
  double rotation_rad = AffineTransformation::rotation_rad(T);
  return 180.0 * rotation_rad / boost::math::double_constants::pi;
}

cv::Mat concat(const cv::Mat& lhs, const cv::Mat& rhs);
cv::Mat concat(std::initializer_list<cv::Mat> mat_list);

inline cv::Mat T_scale(double scale_x, double scale_y)
{
  cv::Mat T = cv::Mat::zeros(2, 3, CV_64F);
  T.at<double>(0, 0) = scale_x;
  T.at<double>(1, 1) = scale_y;
  return T;
}

inline cv::Mat T_translate(double delta_x, double delta_y)
{
  cv::Mat T = cv::Mat::eye(2, 3, CV_64F);
  T.at<double>(0, 2) = delta_x;
  T.at<double>(1, 2) = delta_y;
  return T;
}

inline cv::Mat T_rotate_rad(double angle_rad)
{
  cv::Mat T = cv::Mat::eye(2, 3, CV_64F);
  T.at<double>(0, 0) = std::cos(angle_rad);
  T.at<double>(0, 1) = -std::sin(angle_rad);
  T.at<double>(1, 0) = std::sin(angle_rad);
  T.at<double>(1, 1) = std::cos(angle_rad);
  return T;
}

inline cv::Mat T_rotate_deg(double angle_deg)
{
  return T_rotate_rad(angle_deg / 180.0 * boost::math::double_constants::pi);
}

cv::Mat fit(const std::vector<cv::Point2d>& points_from, const std::vector<cv::Point2d>& points_to);

}

#endif /* TRLIB_AFFINE_TRANSFORMATION_HPP_ */
