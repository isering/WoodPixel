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

#ifndef TRLIB_HISTOGRAM_HPP_
#define TRLIB_HISTOGRAM_HPP_

#include <opencv2/opencv.hpp>

#include "mat.hpp"
#include "texture.hpp"

class Histogram
{ 
public:
  Histogram(cv::Mat texture, cv::Mat mask, int num_bins);
  Histogram(cv::Mat texture, cv::Mat mask, float range_min, float range_max, int num_bins);
  Histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, int num_bins);
  Histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, float range_min, float range_max, int num_bins);

  cv::Mat draw(bool lines) const;

  static void linear_normalization(cv::Mat image, cv::Mat mask, double perc_min=0.05, double perc_max=0.95);
  static void linear_normalization_rows(mat<cv::Mat>& matrix, cv::Mat mask);
  static cv::Mat histogram_matching(cv::Mat texture_source, cv::Mat mask_source, const std::vector<cv::Mat>& target);
  static cv::Mat histogram_matching_hsv(cv::Mat texture_source, cv::Mat mask_source, const std::vector<Texture>& target);
  static cv::Mat histogram_matching_hsv_nn(cv::Mat texture_source, cv::Mat mask_source, const std::vector<Texture>& target);

  static void histogram_equalization_rows(mat<cv::Mat>& matrix, cv::Mat mask);

  std::vector<cv::Mat> histogram;
  float range_min, range_max;
  int num_bins;
  int num_channels;

private:
  void build_histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, int num_bins);
  void build_histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, float range_min, float range_max, int num_bins);
  void build_histogram(const std::vector<std::vector<cv::Mat>>& texture_planes, const std::vector<cv::Mat>& masks, float range_min, float range_max, int num_bins);
};

#endif /* TRLIB_HISTOGRAM_HPP_ */