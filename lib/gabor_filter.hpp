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

#ifndef TRLIB_GABOR_FILTER_HPP_
#define TRLIB_GABOR_FILTER_HPP_

#include <boost/math/constants/constants.hpp>
#include <opencv2/opencv.hpp>

struct GaborFilter
{
  GaborFilter() = default;
  GaborFilter(double frequency, double theta, double sigma_x, double sigma_y);

  cv::Mat apply(cv::Mat texture) const;
  cv::Mat mask(cv::Mat mask_texture) const;

  cv::Mat kernel_real;
  cv::Mat kernel_imag;

  double frequency;
  double theta;
  double sigma_x;
  double sigma_y;

  int kernel_size_x;
  int kernel_size_y;
};

double compute_sigma_x(double frequency, double bandwidth_frequency_octaves);
double compute_sigma_y(double frequency, double bandwidth_angular);

#endif /* TRLIB_GABOR_FILTER_HPP_ */