/*
WoodPixel - Supplementary code for Computational Parquetry:
            Fabricated Style Transfer with Wood Pixels
            ACM Transactions on Graphics 39(2), 2020

Copyright (C) 2020  Julian Iseringhausen, University of Bonn, <iseringhausen@cs.uni-bonn.de>
Copyright (C) 2020  Matthias Hullin, University of Bonn, <hullin@cs.uni-bonn.de>

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

#ifndef TRLIB_PATCH_REGION_IMPL_HPP_
#define TRLIB_PATCH_REGION_IMPL_HPP_

#include "patch_region.hpp"

template<typename T>
void PatchRegion::draw(cv::Mat image, const T& color, double scale) const
{
  if (has_sub_regions())
  {
    for (const PatchRegion& sub_region : m_sub_regions)
    {
      sub_region.draw(image, color, scale);
    }
  }
  else
  {
    for (const BezierCurve& c : m_curves_top)
    {
      c.scaled(scale).draw<T>(image, color);
    }
    for (const BezierCurve& c : m_curves_bot)
    {
      c.scaled(scale).draw<T>(image, color);
    }
    for (const BezierCurve& c : m_curves_left)
    {
      c.scaled(scale).draw<T>(image, color);
    }
    for (const BezierCurve& c : m_curves_right)
    {
      c.scaled(scale).draw<T>(image, color);
    }
    for (const BezierCurve& c : m_curves_diag)
    {
      c.scaled(scale).draw<T>(image, color);
    }
  }
}

template <typename T>
void PatchRegion::draw_edge_mask(cv::Mat image, const T& color, int dilate, cv::Point offset) const
{
  cv::Mat mask = cv::Mat::zeros(image.size(), CV_8UC1);

  for (const BezierCurve& curve : m_curves_top)
  {
    const std::vector<CurveDrawPoint> curve_points = (curve-offset).get_draw_points();
    for (const CurveDrawPoint& c : curve_points)
    {
      if (c.p.x >= 0 && c.p.x < image.cols && c.p.y >= 0 && c.p.y < image.rows)
      {
        mask.at<unsigned char>(c.p) = 255;
      }
    }
  }
  for (const BezierCurve& curve : m_curves_bot)
  {
    const std::vector<CurveDrawPoint> curve_points = (curve-offset).get_draw_points();
    for (const CurveDrawPoint& c : curve_points)
    {
      if (c.p.x >= 0 && c.p.x < image.cols && c.p.y >= 0 && c.p.y < image.rows)
      {
        mask.at<unsigned char>(c.p) = 255;
      }
    }
  }
  for (const BezierCurve& curve : m_curves_left)
  {
    const std::vector<CurveDrawPoint> curve_points = (curve-offset).get_draw_points();
    for (const CurveDrawPoint& c : curve_points)
    {
      if (c.p.x >= 0 && c.p.x < image.cols && c.p.y >= 0 && c.p.y < image.rows)
      {
        mask.at<unsigned char>(c.p) = 255;
      }
    }
  }
  for (const BezierCurve& curve : m_curves_right)
  {
    const std::vector<CurveDrawPoint> curve_points = (curve-offset).get_draw_points();
    for (const CurveDrawPoint& c : curve_points)
    {
      if (c.p.x >= 0 && c.p.x < image.cols && c.p.y >= 0 && c.p.y < image.rows)
      {
        mask.at<unsigned char>(c.p) = 255;
      }
    }
  }
  for (const BezierCurve& curve : m_curves_diag)
  {
    const std::vector<CurveDrawPoint> curve_points = (curve-offset).get_draw_points();
    for (const CurveDrawPoint& c : curve_points)
    {
      if (c.p.x >= 0 && c.p.x < image.cols && c.p.y >= 0 && c.p.y < image.rows)
      {
        mask.at<unsigned char>(c.p) = 255;
      }
    }
  }

  if (dilate > 0)
  {
    cv::dilate(mask, mask, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)), cv::Point(-1, -1), dilate);
  }

  image.setTo(color, mask);
}

template <typename T, typename Iter>
void PatchRegion::draw(cv::Mat image, const T& color, double scale, Iter begin, Iter end)
{
  for (Iter iter = begin; iter != end; ++iter)
  {
    if (iter->has_sub_regions())
    {
      for (const PatchRegion& p : iter->sub_regions())
      {        
        p.draw(image, color, scale);
      }
    }
    else
    {
      iter->draw(image, color, scale);
    }
  }
}

#endif /* TRLIB_PATCH_REGION_IMPL_HPP_ */