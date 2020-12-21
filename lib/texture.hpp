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

#ifndef TRLIB_TEXTURE_HPP_
#define TRLIB_TEXTURE_HPP_

#include <boost/filesystem.hpp>
#include <opencv2/opencv.hpp>

#include "bezier_curve.hpp"
#include "feature_vector.hpp"
#include "texture_marker.hpp"
#include "serializable.hpp"

struct TextureRegion
{
  cv::Rect region;
  cv::Mat mask_1, mask_2;
  BezierCurve separating_curve;

  bool empty() const
  {
    return (region.area() == 0) || (!mask_1.data) || (!mask_2.data);
  }
};

class Texture : public Serializable
{
public:
  Texture() = default;
  Texture(const boost::filesystem::path& filename, double dpi, double scale, const TextureMarker marker=TextureMarker(), const std::string id=std::string()) :
    dpi(scale * dpi),
    scale(scale),
    id(id)
  {
    load_texture(filename, scale);

    this->marker = marker.scaled(scale);
    mask_done = cv::Mat(texture.size(), CV_8UC1, cv::Scalar(255));
    mask_rotation = cv::Mat(texture.size(), CV_8UC1, cv::Scalar(255));
    transformation_matrix = cv::Mat::eye(2, 3, CV_64FC1);
    transformation_matrix_inv = cv::Mat::eye(2, 3, CV_64FC1);
    angle_rad = 0.0;
    this->filename = filename;
  }

  Texture(const boost::filesystem::path& filename, const boost::filesystem::path& filename_mask, double dpi, double scale, TextureMarker marker=TextureMarker(), const std::string id=std::string()) :
    dpi(scale * dpi),
    scale(scale),
    id(id)
  {
	  load_texture(filename, scale);
	  load_mask(filename_mask, scale);

    this->marker = marker.scaled(scale);
    mask_rotation = cv::Mat(texture.size(), CV_8UC1, cv::Scalar(255));
    transformation_matrix = cv::Mat::eye(2, 3, CV_64FC1);
    transformation_matrix_inv = cv::Mat::eye(2, 3, CV_64FC1);
    angle_rad = 0.0;
    this->filename = filename;
  }

  virtual void load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const override;

  Texture clone() const;

  Texture operator()(const cv::Range& row_range, const cv::Range& col_range);
  Texture operator()(const cv::Rect& rect);

  TextureRegion get_regions(cv::Rect region, cv::Mat edge_image);

  cv::Mat mask() const
  {
    return cv::min(mask_done, mask_rotation);
  }

  Texture rotate(double angle_rad) const;

  static void mask_patch(std::vector<Texture>& rotated_textures, int index, cv::Point anchor, const std::vector<cv::Point>& patch);

  cv::Mat mask_rotation_inv() const
  {
    return cv::Scalar(255) - mask_rotation;
  }

  void downsample_nn(int factor);

  static cv::Mat compute_transformation_matrix(const cv::Mat& texture, double angle_rad);

  cv::Mat template_match(const Texture& kernel) const;
  cv::Mat template_match(const Texture& kernel, cv::Mat mask) const;

  std::vector<cv::Vec3f> find_markers(double marker_size_mm, int num_markers);

  cv::Vec3b interpolate_texture(const cv::Point2f& p) const;
  
  cv::Mat texture;
  cv::Mat mask_done;
  cv::Mat mask_rotation;
  cv::Mat transformation_matrix;
  cv::Mat transformation_matrix_inv;

  TextureMarker marker;
  std::string id;

  double angle_rad;
  double scale;
  double dpi;

  FeatureVector response;

  boost::filesystem::path filename;

private:
  void load_texture(const boost::filesystem::path& filename, double scale);
  void load_mask(const boost::filesystem::path& filename, double scale);
};


#endif /* TRLIB_TEXTURE_HPP_ */