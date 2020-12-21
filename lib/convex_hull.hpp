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

#ifndef TRLIB_CONVEX_HULL_HPP_
#define TRLIB_CONVEX_HULL_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

class ConvexHull
{
public:
  ConvexHull() = default;
  ConvexHull(const ConvexHull&) = default;
  ConvexHull(const std::vector<cv::Point>& points);

  const std::vector<cv::Point>& convex_hull() const
  {
    return m_convex_hull;
  }

  void fill(cv::Mat& mask) const;
  std::vector<cv::Point> fill_convex_hull() const;

  void add_point(cv::Point p);
  void from_points(const std::vector<cv::Point>& points);

private:
  std::vector<cv::Point> m_convex_hull;
};

#endif /* TRLIB_CONVEX_POLYGON_HPP_ */