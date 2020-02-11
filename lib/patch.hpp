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

#ifndef TRLIB_PATCH_HPP_
#define TRLIB_PATCH_HPP_

#include <vector>
#include <string>

#include <opencv2/opencv.hpp>

#include "patch_region.hpp"
#include "serializable.hpp"
#include "texture.hpp"

class Patch : public Serializable
{
public:
  Patch() = default;
  Patch(const PatchRegion& region_target, const cv::Point& anchor_source, const cv::Mat& transformation_source, int source_index, int source_rot, const cv::Mat& error, double cost) :
    region_target(region_target),
    anchor_source(anchor_source),
    transformation_source(transformation_source.clone()),
    source_index(source_index),
    source_rot(source_rot),
    error(error.clone()),
    m_cost(cost)
  {
    cv::invertAffineTransform(transformation_source, transformation_source_inv);
  }

  virtual cv::Size size() const
  {
    return region_target.size();
  }

  cv::Point anchor_target() const
  {
    return region_target.anchor();
  }

  cv::Rect bounding_box() const
  {
    return region_target.bounding_box();
  }

  Patch scaled(double scale) const;

  void draw(cv::Mat image, double image_scale, const Texture& texture, double texture_scale) const;
  void draw_edge_mask(cv::Mat image, double image_scale) const;

  cv::Mat mask() const
  {
    return region_target.mask();
  }

  static cv::Rect bounding_box(const std::vector<Patch>& patches, double scale);

  PatchRegion region_target;

  cv::Point anchor_source;
  cv::Mat transformation_source;
  cv::Mat transformation_source_inv;
  int source_index;
  int source_rot;

  cv::Mat error;
  
  virtual double cost() const
  {
    return m_cost;
  }

  virtual void load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const override;

protected:
  double m_cost;
};


#endif /* TRLIB_PATCH_HPP_ */