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

#ifndef TRLIB_GABOR_FILTER_BANK_HPP_
#define TRLIB_GABOR_FILTER_BANK_HPP_

#include <opencv2/opencv.hpp>

#include "gabor_filter.hpp"
//#include "gabor_filter_response.hpp"
#include "mat.hpp"

class GaborFilterBank
{
public:
  GaborFilterBank(int filter_resolution, double frequency_octaves, int num_directions);

  mat<cv::Mat> compute_response(cv::Mat texture) const;

  cv::Mat draw() const;

  void print_filter_sums() const;

  int num_frequencies() const
  {
    return m_gabor_filters.height();
  }

  int num_directions() const
  {
    return m_gabor_filters.width();
  }

  const mat<GaborFilter>& filters() const
  {
    return m_gabor_filters;
  }

  cv::Size max_filter_size() const;

private:
  mat<GaborFilter> m_gabor_filters;
};

#endif /* TRLIB_GABOR_FILTER_BANK_HPP_ */