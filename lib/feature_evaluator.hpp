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

#ifndef TRLIB_FEATURE_EVALUATOR_HPP_
#define TRLIB_FEATURE_EVALUATOR_HPP_

#include <opencv2/opencv.hpp>

#include "feature_vector.hpp"
#include "gabor_filter_bank.hpp"
#include "histogram_vector.hpp"
#include "texture.hpp"

class FeatureEvaluator
{
public:
  FeatureEvaluator(double weight_intensity, double weight_sobel, double weight_gabor, const GaborFilterBank& filter_bank) :
    m_weight_intensity(weight_intensity),
    m_weight_sobel(weight_sobel),
    m_weight_gabor(weight_gabor),
    m_filter_bank(filter_bank)
  {
    m_num_channels = filter_bank.num_directions() * filter_bank.num_frequencies() + 3;
  }

  FeatureVector evaluate(cv::Mat texture, cv::Mat mask=cv::Mat()) const;
  FeatureVector evaluate_with_histogram_matching(cv::Mat texture, const std::vector<Texture>& texture_target, cv::Mat mask, double dampening_factor) const;

  HistogramVector compute_feature_histogram(const FeatureVector& feature_vec, const cv::Size& patch_size) const;

  cv::Size max_filter_size() const
  {
    cv::Size filter_size(1, 1);
    if (m_weight_sobel > 0.0)
    {
      filter_size = cv::Size(3, 3);
    }
    if (m_weight_gabor > 0.0)
    {
      cv::Size gabor_size = m_filter_bank.max_filter_size();
      filter_size.width = std::max(filter_size.width, gabor_size.width);
      filter_size.height = std::max(filter_size.height, gabor_size.height);
    }
    return filter_size;
  }

private:
  double m_weight_intensity, m_weight_sobel, m_weight_gabor;
  const GaborFilterBank& m_filter_bank;
  int m_num_channels;
};

#endif /* TRLIB_FEATURE_EVALUATOR_HPP_ */