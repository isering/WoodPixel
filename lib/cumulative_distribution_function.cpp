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

#include "cumulative_distribution_function.hpp"

CumulativeDistributionFunction::CumulativeDistributionFunction(const Histogram& histogram)
{
  from_histogram(histogram);
}

CumulativeDistributionFunction::CumulativeDistributionFunction(const std::vector<cv::Mat>& cdf, float range_min, float range_max, int num_bins) :
  cdf(cdf),
  range_min(range_min),
  range_max(range_max),
  num_bins(num_bins)
{}

CumulativeDistributionFunction::CumulativeDistributionFunction(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, int num_bins)
{
  Histogram hist(data, num_bins);
  from_histogram(hist);
}

CumulativeDistributionFunction::CumulativeDistributionFunction(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, float range_min, float range_max, int num_bins)
{
  Histogram hist(data, range_min, range_max, num_bins);
  from_histogram(hist);
}

void CumulativeDistributionFunction::from_histogram(const Histogram& histogram)
{
  this->cdf = histogram.histogram;
  this->range_min = histogram.range_min;
  this->range_max = histogram.range_max;
  this->num_bins = histogram.num_bins;

  // Compute CDF for  histogram equalization.
  for (cv::Mat cdf_channel : cdf)
  {
    float* ptr = reinterpret_cast<float*>(cdf_channel.data);
    for (int i = 1; i < cdf_channel.rows; ++i)
    {
      ptr[i] += ptr[i-1];
    }

    // Normalize CDF.
    cdf_channel /= ptr[cdf_channel.rows-1];
  }
}

CumulativeDistributionFunction CumulativeDistributionFunction::inverse() const
{
  std::vector<cv::Mat> cdf_inverse;

  for (cv::Mat cdf_channel : cdf)
  {
    cv::Mat cdf_inverse_channel(cdf_channel.size(), CV_32FC1);
    const float* ptr_cdf = reinterpret_cast<const float*>(cdf_channel.ptr());
    float* ptr_cdf_inverse = reinterpret_cast<float*>(cdf_inverse_channel.ptr());

    int c = 0;
    for (int i = 0; i < num_bins; ++i)
    {
      float x = static_cast<float>(i+1) / num_bins;
      while (c < num_bins && ptr_cdf[c] < x)
      {
        ++c;
      }
      ptr_cdf_inverse[i] = static_cast<float>(c+1) / num_bins;
    }
    cdf_inverse.push_back(cdf_inverse_channel);
  }
  
  return CumulativeDistributionFunction(cdf_inverse, 0.0f, 1.0f, num_bins);
}