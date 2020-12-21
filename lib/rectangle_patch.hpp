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

#ifndef TRLIB_RECTANGLE_PATCH_HPP_
#define TRLIB_RECTANGLE_PATCH_HPP_

#include <opencv2/opencv.hpp>

struct RectanglePatch
{
  RectanglePatch(int index, cv::Rect region_source, cv::Rect region_target) : index(index),
                                                                              region_source(region_source),
                                                                              region_target(region_target)
  {
  }

  int index;
  cv::Rect region_source;
  cv::Rect region_target;
};

#endif /* TRLIB_RECTANGLE_PATCH_HPP_ */
