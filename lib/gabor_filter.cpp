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

#include "gabor_filter.hpp"

#include <cmath>

#include <boost/math/constants/constants.hpp>

cv::Mat GaborFilter::apply(cv::Mat texture) const
{
  cv::Mat response;
  cv::Mat texture_gray;

  if (texture.channels() == 3)
  {
    cv::cvtColor(texture, texture_gray, cv::COLOR_BGR2GRAY);
  }
  else
  {
    texture_gray = texture.clone();
  }

  if (texture_gray.depth() == CV_8U)
  { 
    texture_gray.convertTo(texture_gray, CV_32FC1, 1.0 / 255.0);
  }
  else if (texture_gray.depth() == CV_16U)
  {
    texture_gray.convertTo(texture_gray, CV_32FC1, 1.0 / 65535.0);
  }

  cv::Mat response_real, response_imag;
  cv::filter2D(texture_gray, response_real, -1, kernel_real, cv::Point(-1, -1), 0.0, cv::BORDER_REFLECT_101);
  cv::filter2D(texture_gray, response_imag, -1, kernel_imag, cv::Point(-1, -1), 0.0, cv::BORDER_REFLECT_101);
  cv::magnitude(response_real, response_imag, response);
  
  return response;
}

cv::Mat GaborFilter::mask(cv::Mat mask_texture) const
{
  cv::Mat mask_filter;
  cv::Mat kernel = cv::Mat::ones(kernel_size_y, kernel_size_x, CV_8UC1);
  cv::erode(mask_texture, mask_filter, kernel, cv::Point(-1, -1), 1, cv::BORDER_CONSTANT, cv::Scalar(0));
  return mask_filter;
}

GaborFilter::GaborFilter(double frequency, double theta, double sigma_x, double sigma_y) :
  frequency(frequency),
  theta(theta),
  sigma_x(sigma_x),
  sigma_y(sigma_y)
{
  const double cos_theta = std::cos(theta);
  const double sin_theta = std::sin(theta);

  const double sigma_x_squared = sigma_x * sigma_x;
  const double sigma_y_squared = sigma_y * sigma_y;

  const double pi = boost::math::constants::pi<double>();
  const double factor = 1.0 / (2.0 * pi * sigma_x * sigma_y);

  int kernel_size = static_cast<int>(pi * std::max(sigma_x, sigma_y));
  kernel_size_x = kernel_size;
  kernel_size_y = kernel_size;

  /*
  kernel_size_x = static_cast<int>(pi * std::max(sigma_x * std::abs(std::cos(theta)), sigma_y * std::abs(std::sin(theta))));
  kernel_size_y = static_cast<int>(pi * std::max(sigma_x * std::abs(std::sin(theta)), sigma_y * std::abs(std::cos(theta))));
  */

  kernel_real = cv::Mat(2*kernel_size_y+1, 2*kernel_size_x+1, CV_32FC1);
  kernel_imag = cv::Mat(2*kernel_size_y+1, 2*kernel_size_x+1, CV_32FC1);

  for (int y = -kernel_size_y; y <= kernel_size_y; ++y)
  {
    for (int x = -kernel_size_x; x <= kernel_size_x; ++x)
    {
      const double x_theta = x * cos_theta + y * sin_theta;
      const double y_theta = -x * sin_theta + y * cos_theta;
      const double gaussian = std::exp(-0.5 * (x_theta*x_theta/sigma_x_squared + y_theta*y_theta/sigma_y_squared));

      kernel_real.at<float>(y+kernel_size_y, x+kernel_size_x) = static_cast<float>(
          factor * gaussian * std::cos(2.0 * pi * frequency * x_theta));
      kernel_imag.at<float>(y+kernel_size_y, x+kernel_size_x) = static_cast<float>(
          factor * gaussian * std::sin(2.0 * pi * frequency * x_theta));
    }
  }
}

double compute_sigma_x(double frequency, double bandwidth_frequency_octaves)
{
  const double bandwidth_frequency = std::pow(2.0, bandwidth_frequency_octaves);
  const double pi = boost::math::constants::pi<double>();

  return std::sqrt(std::log(2.0)) * (bandwidth_frequency + 1.0) / (std::sqrt(2.0) * pi * frequency * (bandwidth_frequency - 1.0));
}

double compute_sigma_y(double frequency, double bandwidth_angular)
{
  const double pi = boost::math::constants::pi<double>();

  return std::sqrt(std::log(2.0)) / (std::sqrt(2.0) * pi * frequency * std::tan(0.5 * bandwidth_angular));
}
