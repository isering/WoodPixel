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

#include "convex_hull.hpp"

ConvexHull::ConvexHull(const std::vector<cv::Point> &points)
{
  from_points(points);
}

void ConvexHull::add_point(cv::Point p)
{
  std::vector<cv::Point> points;
  points.swap(m_convex_hull);
  points.push_back(p);
  cv::convexHull(points, m_convex_hull);
}

void ConvexHull::from_points(const std::vector<cv::Point> &points)
{
  m_convex_hull.clear();
  cv::convexHull(points, m_convex_hull);
}

void ConvexHull::fill(cv::Mat &mask) const
{
  cv::fillConvexPoly(mask, m_convex_hull, cv::Scalar(1));
}
