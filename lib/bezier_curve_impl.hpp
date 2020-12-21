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

#ifndef TRLIB_BEZIER_CURVE_IMPL_HPP_
#define TRLIB_BEZIER_CURVE_IMPL_HPP_

#include "bezier_curve.hpp"

template <typename T>
void BezierCurve::draw(cv::Mat image, const T& color, const cv::Point& offset) const
{
  if (m_control_points.empty())
  {
    return;
  }

  if (m_draw_points.empty())
  {
    generate_draw_points();
  }

  // Draw points onto image.
  for (const CurveDrawPoint& d : m_draw_points)
  {
    const cv::Point p = d.p - offset;
    if (p.x >= 0 && p.x < image.cols &&
        p.y >= 0 && p.y < image.rows)
    {
      image.at<T>(p) = color;
    }
  }
}

template <typename T>
void BezierCurve::operator+=(const cv::Point_<T>& p)
{
  for (cv::Point2d& c : m_control_points)
  {
    c += static_cast<cv::Point2d>(p);
  }
  clear_draw_points();
}

template <typename T>
void BezierCurve::operator-=(const cv::Point_<T>& p)
{
  for (cv::Point2d& c : m_control_points)
  {
    c -= static_cast<cv::Point2d>(p);
  }
  clear_draw_points();
}

template <typename T>
void BezierCurve::operator*=(const T& scalar)
{
  for (cv::Point2d& c : m_control_points)
  {
    c *= static_cast<double>(scalar);
  }
  clear_draw_points();
}

template <typename T>
void BezierCurve::operator/=(const T& scalar)
{
  for (cv::Point2d& c : m_control_points)
  {
    c /= static_cast<double>(scalar);
  }
  clear_draw_points();
}

template <typename T>
BezierCurve operator+(const BezierCurve& curve, const cv::Point_<T>& p)
{
  BezierCurve curve_out(curve);
  curve_out += p;
  return curve_out;
}

template <typename T>
BezierCurve operator+(const cv::Point_<T>& p, const BezierCurve& curve)
{
  BezierCurve curve_out(curve);
  curve_out += p;
  return curve_out;
}

template <typename T>
BezierCurve operator-(const BezierCurve& curve, const cv::Point_<T>& p)
{
  BezierCurve curve_out(curve);
  curve_out -= p;
  return curve_out;
}

template <typename T>
BezierCurve operator-(const cv::Point_<T>& p, const BezierCurve& curve)
{
  BezierCurve curve_out(curve);
  curve_out *= -1.0;
  curve_out += p;
  return curve_out;
}

template <typename T>
BezierCurve operator*(const BezierCurve& curve, const T& scalar)
{
  BezierCurve curve_out(curve);
  curve_out *= scalar;
  return curve_out;
}

template <typename T>
BezierCurve operator*(const T& scalar, const BezierCurve& curve)
{
  BezierCurve curve_out(curve);
  curve_out *= scalar;
  return curve_out;
}

template <typename T>
BezierCurve operator/(const BezierCurve& curve, const T& scalar)
{
  BezierCurve curve_out(curve);
  curve_out /= scalar;
  return curve_out;
}

#endif /* TRLIB_BEZIER_CURVE_IMPL_HPP_ */