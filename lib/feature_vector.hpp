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

#ifndef TRLIB_FEATURE_VECTOR_HPP_
#define TRLIB_FEATURE_VECTOR_HPP_

#include <vector>

#include <boost/filesystem.hpp>
#include <opencv2/opencv.hpp>

class FeatureVector
{
public:
  FeatureVector() = default;
  FeatureVector(const std::vector<cv::Mat>& feature_vector) :
    feature_vector(feature_vector)
  {}

  FeatureVector clone() const;

  void render() const;
  static void render(const std::vector<FeatureVector>& features);

  void save(const boost::filesystem::path& path) const;

  float dist(cv::Point p, const FeatureVector& rhs, cv::Point p_rhs) const;

  int num_channels() const
  {
    return static_cast<int>(feature_vector.size());
  }

  int rows() const
  {
    return num_channels() > 0 ? feature_vector.front().rows : 0;
  }

  int cols() const
  {
    return num_channels() > 0 ? feature_vector.front().cols : 0;
  }

  int depth() const
  {
    return num_channels() > 0 ? feature_vector.front().depth() : -1;
  }

  cv::Size size() const
  {
    return num_channels() > 0 ? feature_vector.front().size() : cv::Size();
  }

  cv::Mat operator[](int i) const
  {
    return feature_vector[i];
  }

  FeatureVector operator()(const cv::Rect& region) const
  {
    FeatureVector result;
    for (cv::Mat feature : feature_vector)
    {
      result.feature_vector.push_back(feature(region));
    }
    return result;
  }

  FeatureVector operator()(const cv::Range& row_range, const cv::Range& col_range) const
  {
    FeatureVector result;
    for (cv::Mat feature : feature_vector)
    {
      result.feature_vector.push_back(feature(row_range, col_range));
    }
    return result;
  }

  void update_weights(float weight_intensity, float weight_gabor, float weight_laplacian)
  {
    feature_vector[0] *= weight_intensity;
    feature_vector[1] *= weight_laplacian;
    for (int i = 2; i < num_channels(); ++i)
    {
      feature_vector[i] *= weight_gabor;
    }
  }

  cv::Mat dist_sqr_mat(const FeatureVector& rhs) const;

  void downsample_nn(int factor);

  std::vector<cv::Mat> feature_vector;
};

#endif /* TRLIB_FEATURE_VECTOR_HPP_ */