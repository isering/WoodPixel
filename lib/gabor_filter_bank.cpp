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

#include "gabor_filter_bank.hpp"

#include <boost/math/constants/constants.hpp>

#include "linspace.hpp"

GaborFilterBank::GaborFilterBank(int filter_resolution, double frequency_octaves, int num_directions)
{
  const double bandwidth_angular = boost::math::constants::pi<double>() / num_directions;

  std::vector<double> frequencies;
  for (int freq = 4; freq <= filter_resolution / 4; freq *= 2)
  {
    frequencies.push_back((sqrt(2.0) * freq) / filter_resolution);
  }
  const int num_frequencies = static_cast<int>(frequencies.size());

  std::vector<double> directions = linspace(0.0, boost::math::constants::pi<double>(), num_directions, false);

  m_gabor_filters.resize(num_frequencies, num_directions);

  for (int i = 0; i < num_frequencies; ++i)
  {
    const double sigma_x = compute_sigma_x(frequencies[i], frequency_octaves);
    const double sigma_y = compute_sigma_y(frequencies[i], bandwidth_angular);

    for (int j = 0; j < num_directions; ++j)
    {
      m_gabor_filters(i, j) = GaborFilter(frequencies[i], directions[j], sigma_x, sigma_y);
    }
  }
}

mat<cv::Mat> GaborFilterBank::compute_response(cv::Mat texture) const
{
  const double gabor_max = 0.2;
  mat<cv::Mat> response(m_gabor_filters.height(), m_gabor_filters.width());
  for (int i = 0; i < m_gabor_filters.size(); ++i)
  {
    response[i] = m_gabor_filters[i].apply(texture);
    response[i] = cv::min(response[i], gabor_max) / gabor_max;
    response[i].convertTo(response[i], CV_16UC1, 65535.0);
  }
  return response;
}

static void minmax(cv::Mat mat)
{
  double min_val, max_val;
  cv::minMaxLoc(mat, &min_val, &max_val);
  std::cout << min_val << " / " << max_val << std::endl;
}

cv::Mat GaborFilterBank::draw() const
{
  int max_rows = 0;
  int max_cols = 0;

  for (int i = 0; i < m_gabor_filters.size(); ++i)
  {
    max_rows = std::max(max_rows, m_gabor_filters[i].kernel_real.rows);
    max_cols = std::max(max_cols, m_gabor_filters[i].kernel_real.cols);
  }

  cv::Mat image = cv::Mat::zeros(m_gabor_filters.height() * max_rows, m_gabor_filters.width() * max_cols, CV_32FC1);

  for (int y = 0; y < m_gabor_filters.height(); ++y)
  {
    for (int x = 0; x < m_gabor_filters.width(); ++x)
    {
      cv::Mat kernel;
      cv::normalize(m_gabor_filters(y, x).kernel_real, kernel, 1.0, 0.0, cv::NORM_MINMAX);
      cv::Rect rect_target(x * max_cols, y * max_rows, kernel.cols, kernel.rows);
      kernel.copyTo(image(rect_target));
    }
  }

  return image;
}

void GaborFilterBank::print_filter_sums() const
{
  for (int y = 0; y < num_frequencies(); ++y)
  {
    for (int x = 0; x < num_directions(); ++x)
    {
      cv::Scalar sum_real = cv::sum(m_gabor_filters(y, x).kernel_real);
      cv::Scalar sum_imag = cv::sum(m_gabor_filters(y, x).kernel_imag);
      std::cout << "Filter " << y << " " << x << ": " << sum_real[0] << " / " << sum_imag[0] << std::endl;
    }
  }
}

cv::Size GaborFilterBank::max_filter_size() const
{
  int kernel_rows_max = 0;
  int kernel_cols_max = 0;
  for (int i = 0; i < m_gabor_filters.size(); ++i)
  {
    kernel_rows_max = std::max(kernel_rows_max, m_gabor_filters[i].kernel_size_y);
    kernel_cols_max = std::max(kernel_cols_max, m_gabor_filters[i].kernel_size_x);
  }
  return cv::Size(kernel_cols_max, kernel_rows_max);
}