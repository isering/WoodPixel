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

#include "affine_transformation.hpp"
#include "patch.hpp"

cv::Rect Patch::bounding_box(const std::vector<Patch>& patches, double scale)
{
  std::vector<cv::Point> points;
  for (const Patch& patch : patches)
  {
    const std::vector<cv::Point> points_patch = patch.scaled(scale).region_target.get_draw_points();
    points.insert(points.end(), points_patch.begin(), points_patch.end());
  }
  return cv::boundingRect(points);
}

static cv::Vec3b interpolate_bilinear(cv::Mat image, const cv::Point2d& p)
{
  const int x = static_cast<int>(p.x);
  const int y = static_cast<int>(p.y);
  const double dx = p.x - x;
  const double dy = p.y - y;

  const cv::Vec3f c11 = cv::Vec3f(image.at<cv::Vec3b>(y, x));
  const cv::Vec3f c12 = cv::Vec3f(image.at<cv::Vec3b>(y, x + 1));
  const cv::Vec3f c21 = cv::Vec3f(image.at<cv::Vec3b>(y + 1, x));
  const cv::Vec3f c22 = cv::Vec3f(image.at<cv::Vec3b>(y + 1, x + 1));
  const cv::Vec3f c1 = (1.0f - dx) * c11 + dx * c12;
  const cv::Vec3f c2 = (1.0f - dx) * c21 + dx * c22;
  const cv::Vec3f c = (1.0f - dy) * c1 + dy * c2;

  return cv::Vec3b(c);
}

void Patch::draw(cv::Mat image, double image_scale, const Texture& texture, double texture_scale) const
{
  const cv::Rect bbox = bounding_box();

  const cv::Mat T_1 = AffineTransformation::T_scale(1.0 / texture_scale, 1.0 / texture_scale);
  const cv::Mat T_2 = transformation_source_inv;
  const cv::Mat T_3 = AffineTransformation::T_scale(1.0 / image_scale, 1.0 / image_scale);
  const cv::Mat T_source = AffineTransformation::concat({T_1, T_2});
  
  const PatchRegion region_scaled = region_target.scaled(image_scale);
  const cv::Rect bbox_scaled = region_scaled.bounding_box();
  cv::Mat mask_scaled = region_scaled.mask();

  cv::Point anchor_target_scaled = image_scale * anchor_target();
  
  for (int y = 0; y < bbox_scaled.height; ++y)
  {
    cv::Vec3b* image_ptr = reinterpret_cast<cv::Vec3b*>(image.ptr(bbox_scaled.y+y));
    const unsigned char* mask_ptr = reinterpret_cast<const unsigned char*>(mask_scaled.ptr(y));
    for (int x = 0; x < bbox_scaled.width; ++x)
    {
      if (mask_ptr[x])
      {
        cv::Point2f p_source = AffineTransformation::transform(
          T_source,
          cv::Point2f(
            static_cast<float>(anchor_source.x + (anchor_target_scaled.x - bbox_scaled.x + x) / image_scale),                       
            static_cast<float>(anchor_source.y + (anchor_target_scaled.y - bbox_scaled.y + y) / image_scale)));
        image_ptr[bbox_scaled.x+x] = texture.interpolate_texture(p_source);
      }
    }
  }
}

void Patch::draw_edge_mask(cv::Mat image, double image_scale) const
{
  region_target.scaled(image_scale).draw_edge_mask<unsigned char>(image, 255, 0);
}

Patch Patch::scaled(double scale) const
{
  cv::Mat transformation_scaled = AffineTransformation::concat(transformation_source, AffineTransformation::T_scale(1.0/scale, 1.0/scale));
  return Patch(region_target.scaled(scale), scale*anchor_source, transformation_scaled, source_index, source_rot, cv::Mat(), -1.0);
}

boost::property_tree::ptree Patch::save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const
{
  boost::property_tree::ptree tree;
  serialize(tree, "region_target", region_target, base_path, path);
  serialize(tree, "anchor_source", anchor_source, base_path, path);
  serialize_mat<double>(tree, "transformation_source", transformation_source, base_path, path);
  serialize_mat<double>(tree, "transformation_source_inv", transformation_source_inv, base_path, path);
  serialize(tree, "source_index", source_index, base_path, path);
  serialize(tree, "source_rot", source_rot, base_path, path);
  serialize_image(tree, "error", error, base_path, path);
  return tree;
}

void Patch::load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree)
{
  deserialize(tree, "region_target", region_target, base_path);
  deserialize(tree, "anchor_source", anchor_source, base_path);
  deserialize_mat<double>(tree, "transformation_source", transformation_source, base_path);
  deserialize_mat<double>(tree, "transformation_source_inv", transformation_source_inv, base_path);
  deserialize(tree, "source_index", source_index, base_path);
  deserialize(tree, "source_rot", source_rot, base_path);
  deserialize_image(tree, "error", error, base_path);
  if (!error.data || error.cols == 0 || error.rows == 0)
  {
    throw(std::runtime_error("Unable to read error image."));
  }
}
