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

#ifndef TRLIB_CUMULATIVE_DISTRIBUTION_FUNCTION_HPP_
#define TRLIB_CUMULATIVE_DISTRIBUTION_FUNCTION_HPP_

#include <opencv2/opencv.hpp>

#include "draw.hpp"
#include "histogram.hpp"

struct CumulativeDistributionFunction
{
  CumulativeDistributionFunction(const Histogram& histogram);
  CumulativeDistributionFunction(const std::vector<cv::Mat>& cdf, float range_min, float range_max, int num_bins);
  CumulativeDistributionFunction(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, int num_bins);
  CumulativeDistributionFunction(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, float range_min, float range_max, int num_bins);

  void from_histogram(const Histogram& histogram);

  std::vector<cv::Mat> cdf;
  float range_min, range_max;
  int num_bins;

  float operator()(int c, float x) const
  {
    return cdf[c].template at<float>(bin(x), 0);
  }

  int bin(float val) const
  {
    return static_cast<int>((num_bins - 1) * (static_cast<double>(val - range_min) / static_cast<double>(range_max - range_min)));
  }

  template <typename T, typename U>
  T equalize(U val, int channel, unsigned char mask, T range_out_min, T range_out_max)
  {
    if (mask)
    {
      float val_normalized = cdf[channel].template at<float>(bin(val), 0);
      return range_out_min + static_cast<T>(val_normalized * (range_out_max - range_out_min));
    }
    else
    {
      return T();
    }
  }
  
  double get_percentile(int channel, double perc)
  {
    int index = 0;
    float* ptr = reinterpret_cast<float*>(cdf[channel].data);
    while (index < num_bins && ptr[index] < perc)
    {
      ++index;
    }
    return range_min + (static_cast<double>(index) / num_bins) * (range_max - range_min);
  }

  CumulativeDistributionFunction inverse() const;

  cv::Mat draw() const
  {
    return draw_lines(cdf);
  }
};

#endif /* TRLIB_CUMULATIVE_DISTRIBUTION_FUNCTION_HPP_ */