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

#include "feature_evaluator.hpp"

#include "histogram.hpp"

FeatureVector FeatureEvaluator::evaluate(cv::Mat texture, cv::Mat mask) const
{
  mat<cv::Mat> gabor_response;
  if (m_weight_gabor > 0.0)
  {
    gabor_response = m_filter_bank.compute_response(texture);
  }

  cv::Mat texture_conv;
  cv::cvtColor(texture, texture_conv, cv::COLOR_BGR2GRAY);

  std::vector<cv::Mat> response_vec;
  if (m_weight_intensity > 0.0)
  {
    response_vec.push_back(m_weight_intensity * texture_conv);
  }
  cv::Mat texture_float;
  texture_conv.convertTo(texture_float, CV_32FC1, 1.0 / 65535.0);

  cv::Mat texture_sobel_x, texture_sobel_y, texture_sobel_mag;
  cv::Sobel(texture_float, texture_sobel_x, CV_32F, 1, 0);
  cv::Sobel(texture_float, texture_sobel_y, CV_32F, 0, 1);
  cv::magnitude(texture_sobel_x, texture_sobel_y, texture_sobel_mag);
  texture_sobel_mag *= 0.25;
  texture_sobel_mag.convertTo(texture_sobel_mag, CV_16UC1, 65535.0);

  if (m_weight_sobel > 0.0)
  {
    response_vec.push_back(m_weight_sobel * texture_sobel_mag);
  }

  if (m_weight_gabor > 0.0)
    for (int i = 0; i < gabor_response.size(); ++i)
    {
      response_vec.push_back(m_weight_gabor * gabor_response[i]);
    }

  return FeatureVector(response_vec);
}

FeatureVector FeatureEvaluator::evaluate_with_histogram_matching(cv::Mat texture, const std::vector<Texture> &texture_target, cv::Mat mask, double dampening_factor) const
{
  mat<cv::Mat> gabor_response;
  if (m_weight_gabor > 0.0)
  {
    gabor_response = m_filter_bank.compute_response(texture);
  }

  cv::Mat texture_gray;
  cv::cvtColor(texture, texture_gray, cv::COLOR_BGR2GRAY);

  cv::Mat texture_float;
  texture_gray.convertTo(texture_float, CV_32FC1, 1.0 / 65535.0);

  cv::Mat texture_sobel_x, texture_sobel_y, texture_sobel_mag;
  cv::Sobel(texture_float, texture_sobel_x, CV_32F, 1, 0);
  cv::Sobel(texture_float, texture_sobel_y, CV_32F, 0, 1);
  cv::magnitude(texture_sobel_x, texture_sobel_y, texture_sobel_mag);
  texture_sobel_mag *= 0.25;
  texture_sobel_mag.convertTo(texture_sobel_mag, CV_16UC1, 65535.0);

  std::vector<cv::Mat> textures_target_gray(texture_target.size());
  for (size_t i = 0; i < texture_target.size(); ++i)
  {
    if (texture_target[i].texture.channels() == 3)
    {
      cv::cvtColor(texture_target[i].texture, textures_target_gray[i], cv::COLOR_BGR2GRAY);
    }
    else
    {
      textures_target_gray[i] = texture_target[i].texture.clone();
    }
  }

  cv::Mat texture_matched = Histogram::histogram_matching(texture_gray, mask, textures_target_gray);
  cv::Mat texture_matched_float;
  texture_matched.convertTo(texture_matched_float, CV_32FC1, 1.0 / 65535.0);
  texture_matched_float = dampening_factor * texture_matched_float + (1.0 - dampening_factor) * texture_float;
  texture_matched_float = cv::min(texture_matched_float, 1.0f);
  texture_matched_float = cv::max(texture_matched_float, 0.0f);
  texture_matched_float.convertTo(texture_matched, CV_16UC1, 65535.0);

  std::vector<cv::Mat> response_vec;
  if (m_weight_intensity > 0.0)
  {
    response_vec.push_back(m_weight_intensity * texture_matched);
  }

  if (m_weight_sobel > 0.0)
  {
    response_vec.push_back(m_weight_sobel * texture_sobel_mag);
  }

  if (m_weight_gabor > 0.0)
  {
    for (int i = 0; i < gabor_response.size(); ++i)
    {
      response_vec.push_back(m_weight_gabor * gabor_response[i]);
    }
  }
  return FeatureVector(response_vec);
}

HistogramVector FeatureEvaluator::compute_feature_histogram(const FeatureVector &feature_vec, const cv::Size &patch_size) const
{
  cv::Mat feature_float;
  std::vector<cv::Mat> histogram_vec(feature_vec.num_channels());

  for (int i = 0; i < feature_vec.num_channels(); ++i)
  {
    feature_vec[i].convertTo(feature_float, CV_32FC1, 1.0 / 65535.0);
    cv::blur(feature_float, histogram_vec[i], patch_size, cv::Point(0, 0));
    //cv::dilate(feature_vec[i], histogram_vec[i], cv::Mat::ones(patch_size, CV_8UC1), cv::Point(0, 0));

    //cv::imshow("Before", feature_vec[i]);
    //cv::imshow("After", histogram_vec[i]);
    //cv::waitKey();

    histogram_vec[i] = histogram_vec[i](cv::Range(0, histogram_vec[i].rows - patch_size.height + 1),
                                        cv::Range(0, histogram_vec[i].cols - patch_size.width + 1));
  }

  return HistogramVector(histogram_vec, std::vector<cv::Mat>());
}