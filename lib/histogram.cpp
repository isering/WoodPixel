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

#include "histogram.hpp"

#include <numeric>

#include "cumulative_distribution_function.hpp"
#include "draw.hpp"
#include "linspace.hpp"

Histogram::Histogram(cv::Mat texture, cv::Mat mask, float range_min, float range_max, int num_bins)
{
  std::vector<std::pair<cv::Mat, cv::Mat>> data(1, std::make_pair(texture, mask));
  build_histogram(data, range_min, range_max, num_bins);
}

Histogram::Histogram(cv::Mat texture, cv::Mat mask, int num_bins)
{
  std::vector<std::pair<cv::Mat, cv::Mat>> data(1, std::make_pair(texture, mask));
  build_histogram(data, num_bins);
}

Histogram::Histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, float range_min, float range_max, int num_bins)
{
  build_histogram(data, range_min, range_max, num_bins);
}

Histogram::Histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, int num_bins)
{
  build_histogram(data, num_bins);
}

void Histogram::build_histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, float range_min, float range_max, int num_bins)
{
  std::vector<std::vector<cv::Mat>> texture_planes(data.size());
  std::vector<cv::Mat> masks(data.size());
  for (size_t i = 0; i < data.size(); ++i)
  {
    cv::split(data[i].first, texture_planes[i]);
    masks[i] = data[i].second;
  }

  build_histogram(texture_planes, masks, range_min, range_max, num_bins);
}

void Histogram::build_histogram(const std::vector<std::pair<cv::Mat, cv::Mat>>& data, int num_bins)
{
  std::vector<std::vector<cv::Mat>> texture_planes(data.size());
  std::vector<cv::Mat> masks(data.size());
  for (size_t i = 0; i < data.size(); ++i)
  {
    cv::split(data[i].first, texture_planes[i]);
    masks[i] = data[i].second;
  }

  double global_min = std::numeric_limits<double>::max();
  double global_max = -std::numeric_limits<double>::max();
  for (size_t i = 0; i < data.size(); ++i)
  {
    for (cv::Mat channel : texture_planes[i])
    {
      double local_min, local_max;
      cv::minMaxLoc(channel, &local_min, &local_max, 0, 0, masks[i]);
      global_min = std::min(global_min, local_min);
      global_max = std::max(global_max, local_max);
    }
  }

  build_histogram(texture_planes, masks, static_cast<float>(global_min), static_cast<float>(global_max), num_bins);
}

void Histogram::build_histogram(const std::vector<std::vector<cv::Mat>>& texture_planes, const std::vector<cv::Mat>& masks, float range_min, float range_max, int num_bins)
{
  const int num_textures = static_cast<int>(texture_planes.size());
  this->range_min = range_min;
  this->range_max = range_max;
  this->num_bins = num_bins;
  this->num_channels = static_cast<int>(texture_planes.front().size());
  this->histogram.resize(this->num_channels);
  
  for (int c = 0; c < num_channels; ++c)
  {
    cv::Mat hist = cv::Mat();
    std::vector<int> channels(1, 0);
    std::vector<int> hist_size(1, num_bins);
    std::vector<float> ranges({static_cast<float>(range_min), static_cast<float>(range_max)});
    for (int i = 0; i < num_textures; ++i)
    {
      std::vector<cv::Mat> images(1);
      texture_planes[i][c].convertTo(images[0], CV_32FC1);
      cv::calcHist(images, channels, masks[i], hist, hist_size, ranges, true);
    }
    this->histogram[c] = hist;
  }
}

cv::Mat Histogram::draw(bool lines) const
{
  cv::Mat hist_image;
  if (lines)
  {
    hist_image = draw_lines(histogram);
  }
  else
  {
    hist_image = draw_bars(histogram);
  }
  return hist_image;
}

void Histogram::histogram_equalization_rows(mat<cv::Mat>& matrix, cv::Mat mask)
{
  // 16 bit histograms.
  const int num_bins = 1 << 16;

  // Normalize rows independently.
  for (int i = 0; i < matrix.height(); ++i)
  {
    // Accumulate data to normalize.
    std::vector<std::pair<cv::Mat, cv::Mat>> data;
    for (int j = 0; j < matrix.width(); ++j)
    {
      data.emplace_back(matrix(i, j), mask);
    }

    // Get CDF for normalization.
    CumulativeDistributionFunction cdf(data, num_bins);
    auto cdf_equalize = std::bind(&CumulativeDistributionFunction::equalize<uint16_t, uint16_t>, cdf, std::placeholders::_1, 0, std::placeholders::_2, 0, 65535);

    // Equalize feature vectors.
    for (int j = 0; j < matrix.width(); ++j)
    {
      cv::Mat response = matrix(i, j);
      for (int row = 0; row < response.rows; ++row)
      {
        uint16_t* ptr = reinterpret_cast<uint16_t*>(response.ptr(row));
        unsigned char* ptr_mask = reinterpret_cast<unsigned char*>(mask.ptr(row));
        std::transform(ptr, ptr+response.cols, ptr_mask, ptr, cdf_equalize);
      }
    }
  }
}

void Histogram::linear_normalization(cv::Mat image, cv::Mat mask, double perc_min, double perc_max)
{
  // 16 bit histograms.
  const int num_bins = 1 << 16;

  std::vector<std::pair<cv::Mat, cv::Mat>> data(1, std::make_pair(image, mask));

  // Get CDF for normalization.
  CumulativeDistributionFunction cdf(data, num_bins);
  double range_min = cdf.get_percentile(0, perc_min);
  double range_max = cdf.get_percentile(0, perc_max);

  struct LinearNormalizeOperator
  {
    LinearNormalizeOperator(double range_min_in, double range_max_in, double range_min_out, double range_max_out) :
      range_min_in(range_min_in),
      range_max_in(range_max_in),
      range_min_out(range_min_out),
      range_max_out(range_max_out)
    {
    }

    uint16_t operator()(uint16_t val_in) const
    {
      double val = static_cast<double>(val_in - range_min_in) / static_cast<double>(range_max_in - range_min_in);
      val = range_min_out + val * (range_max_out - range_min_out);
      val = std::max(val, range_min_out);
      val = std::min(val, range_max_out);
      return static_cast<uint16_t>(val);
    }

    double range_min_in, range_max_in;
    double range_min_out, range_max_out;
  };

  LinearNormalizeOperator normalize_op(range_min, range_max, 0.0, 65535.0);

  for (int row = 0; row < image.rows; ++row)
  {
    uint16_t* ptr = reinterpret_cast<uint16_t*>(image.ptr(row));
    std::transform(ptr, ptr+image.cols, ptr, normalize_op);
  }
}

void Histogram::linear_normalization_rows(mat<cv::Mat>& matrix, cv::Mat mask)
{
  // Set percentile for linear transformation
  const double percentile_min = 0.05;
  const double percentile_max = 0.95;

  // 16 bit histograms.
  const int num_bins = 1 << 16;

  // Normalize rows independently.
  for (int i = 0; i < matrix.height(); ++i)
  {
    // Accumulate data to normalize.
    std::vector<std::pair<cv::Mat, cv::Mat>> data;
    for (int j = 0; j < matrix.width(); ++j)
    {
      data.emplace_back(matrix(i, j), mask);
    }

    // Get CDF for normalization.
    CumulativeDistributionFunction cdf(data, num_bins);
    double range_min = cdf.get_percentile(0, percentile_min);
    double range_max = cdf.get_percentile(0, percentile_max);

    struct LinearNormalizeOperator
    {
      LinearNormalizeOperator(double range_min_in, double range_max_in, double range_min_out, double range_max_out) :
        range_min_in(range_min_in),
        range_max_in(range_max_in),
        range_min_out(range_min_out),
        range_max_out(range_max_out)
      {}

      uint16_t operator()(uint16_t val_in) const
      {
        double val = static_cast<double>(val_in - range_min_in) / static_cast<double>(range_max_in - range_min_in);
        val = range_min_out + val * (range_max_out - range_min_out);
        val = std::max(val, range_min_out);
        val = std::min(val, range_max_out);
        return static_cast<uint16_t>(val);
      }

      double range_min_in, range_max_in;
      double range_min_out, range_max_out;
    };

    LinearNormalizeOperator normalize_op(range_min, range_max, 0.0, 65535.0);

    // Equalize feature vectors.
    for (int j = 0; j < matrix.width(); ++j)
    {
      cv::Mat response = matrix(i, j);
      for (int row = 0; row < response.rows; ++row)
      {
        uint16_t* ptr = reinterpret_cast<uint16_t*>(response.ptr(row));
        std::transform(ptr, ptr+response.cols, ptr, normalize_op);
      }
    }
  }
}

cv::Mat Histogram::histogram_matching(cv::Mat texture_source, cv::Mat mask_source, const std::vector<cv::Mat>& target)
{
  const int num_bins = 1 << 8;

  std::vector<std::pair<cv::Mat, cv::Mat>> data_source(1, std::make_pair(texture_source, mask_source));
  std::vector<std::pair<cv::Mat, cv::Mat>> data_target;
  for (const cv::Mat& t : target)
  {
    data_target.emplace_back(t, cv::Mat());
  }

  CumulativeDistributionFunction cdf_source(data_source, 0, 65535, num_bins);
  CumulativeDistributionFunction cdf_target(data_target, 0, 65535, num_bins);
  CumulativeDistributionFunction cdf_target_inverse = cdf_target.inverse();

  std::vector<cv::Mat> texture_out_planes;
  std::vector<cv::Mat> texture_in_planes;
  cv::split(texture_source, texture_in_planes);

  for (int channel = 0; channel < static_cast<int>(texture_in_planes.size()); ++channel)
  {
    cv::Mat texture_in_channel = texture_in_planes[channel];
    cv::Mat texture_out_channel = cv::Mat::zeros(texture_in_channel.size(), CV_16UC1);

    for (int y = 0; y < texture_in_channel.rows; ++y)
    {
      const uint16_t* ptr_texture_in = reinterpret_cast<const uint16_t*>(texture_in_channel.ptr(y));
      uint16_t* ptr_texture_out = reinterpret_cast<uint16_t*>(texture_out_channel.ptr(y));
      const unsigned char* ptr_mask = mask_source.ptr(y);

      for (int x = 0; x < texture_in_channel.cols; ++x)
      {
        if (ptr_mask[x])
        {
          const float cdf_val_in = cdf_source(channel, ptr_texture_in[x]);
          const float cdf_val_out = cdf_target_inverse(channel, cdf_val_in);
          ptr_texture_out[x] = static_cast<uint16_t>(65535.0f * cdf_val_out);
        }
      }
    }
    texture_out_planes.push_back(texture_out_channel);
  }

  cv::Mat texture_out;
  cv::merge(texture_out_planes, texture_out);

  return texture_out;
}

void temp_histogram_matching(std::vector<cv::Mat>& source_planes, cv::Mat mask_source, const std::vector<std::vector<cv::Mat>>& target_planes, int channel, float range_min, float range_max, int num_bins)
{

  double min_val, max_val;
  cv::minMaxLoc(source_planes[channel], &min_val, &max_val);
  std::cout << "before: " << min_val << " " << max_val << std::endl;

  std::vector<std::pair<cv::Mat, cv::Mat>> data_source(1, std::make_pair(source_planes[channel], mask_source));
  std::vector<std::pair<cv::Mat, cv::Mat>> data_target;
  for (const std::vector<cv::Mat>& plane : target_planes)
  {
    data_target.emplace_back(plane[channel], cv::Mat());
  }

  CumulativeDistributionFunction cdf_source(data_source, range_min, range_max, num_bins);
  CumulativeDistributionFunction cdf_target(data_target, range_min, range_max, num_bins);
  CumulativeDistributionFunction cdf_target_inverse = cdf_target.inverse();

  for (int y = 0; y < source_planes[channel].rows; ++y)
  {
    uint16_t* ptr_texture = reinterpret_cast<uint16_t*>(source_planes[channel].ptr(y));
    const unsigned char* ptr_mask = mask_source.ptr(y);
    for (int x = 0; x < source_planes[channel].cols; ++x)
    {
      if (ptr_mask[x])
      {
        const float cdf_val_in = cdf_source(0, ptr_texture[x]);
        const float cdf_val_out = cdf_target_inverse(0, cdf_val_in);
        ptr_texture[x] = static_cast<uint16_t>(65535.0f * cdf_val_out);
      }
    }
  }

  cv::minMaxLoc(source_planes[channel], &min_val, &max_val);
  std::cout << "after: " << min_val << " " << max_val << std::endl;
}

cv::Mat Histogram::histogram_matching_hsv(cv::Mat texture_source, cv::Mat mask_source, const std::vector<Texture>& target)
{  
  cv::Mat texture_source_hsv;
  cv::cvtColor(texture_source, texture_source_hsv, cv::COLOR_BGR2HSV);

  std::vector<cv::Mat> source_hsv_planes;
  cv::split(texture_source_hsv, source_hsv_planes);

  std::vector<std::vector<cv::Mat>> target_hsv_planes;
  for (const Texture& t : target)
  {
    cv::Mat texture_target_hsv;
    cv::cvtColor(t.texture, texture_target_hsv, cv::COLOR_BGR2HSV);
    std::vector<cv::Mat> temp;
    cv::split(texture_target_hsv, temp);
    target_hsv_planes.push_back(temp);
  }

  //temp_histogram_matching(source_hsv_planes, mask_source, target_hsv_planes, 0, 0.0f, 360.0f, 360);
  //temp_histogram_matching(source_hsv_planes, mask_source, target_hsv_planes, 1, 0.0f, 1.0f, 255);
  temp_histogram_matching(source_hsv_planes, mask_source, target_hsv_planes, 2, 0.0f, 1.0f, 255);
  
  //source_hsv_planes[0] *= 360.0f;

  cv::Mat source_hsv_out;
  cv::merge(source_hsv_planes, source_hsv_out);

  cv::Mat source_bgr_out;
  cv::cvtColor(source_hsv_out, source_bgr_out, cv::COLOR_HSV2BGR);

  return source_bgr_out;
}

cv::Mat Histogram::histogram_matching_hsv_nn(cv::Mat texture_source, cv::Mat mask_source, const std::vector<Texture>& target)
{
  std::vector<cv::Mat> nn_vec;
  for (const Texture& t : target)
  {
    nn_vec.push_back(histogram_matching_hsv(texture_source, mask_source, std::vector<Texture>(1, t)));
  }

  cv::Mat texture_out = cv::Mat::zeros(texture_source.size(), CV_32FC3);
  cv::Mat dist_mat(texture_source.size(), CV_32FC1, 1000.0f);

  for (const cv::Mat& nn : nn_vec)
  {
    for (int y = 0; y < texture_out.rows; ++y)
    {
      const cv::Vec3f* ptr_in = reinterpret_cast<const cv::Vec3f*>(texture_source.ptr(y));
      const cv::Vec3f* ptr_hist_in = reinterpret_cast<const cv::Vec3f*>(nn.ptr(y));
      const unsigned char* ptr_mask = reinterpret_cast<const unsigned char*>(mask_source.ptr(y));
      float* ptr_dist = reinterpret_cast<float*>(dist_mat.ptr(y));
      cv::Vec3f* ptr_out = reinterpret_cast<cv::Vec3f*>(texture_out.ptr(y));
      for (int x = 0; x < texture_out.cols; ++x)
      {
        if (ptr_mask[x])
        {
          const float dist = static_cast<float>(cv::norm(ptr_in[x] - ptr_hist_in[x]));
          if (dist < ptr_dist[x])
          {
            ptr_out[x] = ptr_hist_in[x];
            ptr_dist[x] = dist;
          }
        }
      }
    }
  }

  return texture_out;
}
