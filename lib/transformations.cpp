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

#include "transformations.hpp"

cv::Mat transformations::scale(double scale_x, double scale_y)
{
  return (cv::Mat_<double>(3, 3) << scale_x, 0.0, 0.0, 0.0, scale_y, 0.0, 0.0, 0.0, 1.0);
}

cv::Mat transformations::shear_x(double shear)
{
  return (cv::Mat_<double>(3, 3) << 1.0, shear, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
}

cv::Mat transformations::shear_y(double shear)
{
  return (cv::Mat_<double>(3, 3) << 1.0, 0.0, 0.0, shear, 1.0, 0.0, 0.0, 0.0, 1.0);
}
