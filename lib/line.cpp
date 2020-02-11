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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "line.hpp"

double project_point_line(const cv::Point2d& p, const cv::Point2d& line_1, const cv::Point2d& line_2)
{
  const cv::Point2d s = line_2 - line_1;
  return s.dot(p - line_1) / s.dot(s);
}

double dist_point_line(const cv::Point2d& p, const cv::Point2d& line_1, const cv::Point2d& line_2)
{
  const cv::Point2d d = line_2 - line_1;
  return std::abs(d.y*p.x - d.x*p.y + line_2.x*line_1.y - line_2.y*line_1.x) / cv::norm(d);
}