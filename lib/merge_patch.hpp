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

#ifndef TRLIB_MERGE_PATCH_HPP_
#define TRLIB_MERGE_PATCH_HPP_

#include <vector>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "bezier_curve.hpp"
#include "mat.hpp"
#include "patch.hpp"
#include "texture.hpp"

class MergePatch
{
public:
  MergePatch(const Patch& patch, cv::Rect region_subpatch, cv::Point pos_patch, int subpatch_size, int target_id) :
    anchor_target(patch.region_target.anchor()),
    anchor_source(patch.anchor_source),
    source_index(patch.source_index),
    source_rot(patch.source_rot),
    transformation_source(patch.transformation_source),
    transformation_source_inv(patch.transformation_source_inv),
    error(patch.error),
    region_subpatch(region_subpatch),
    pos_patch(pos_patch),
    subpatch_size(subpatch_size),
    active_pixel(cv::Mat::ones(patch.region_target.mask().size(), CV_8UC1)),
    size(patch.region_target.mask().size()),
    target_id(target_id)
  {}

  MergePatch(cv::Rect region_subpatch, cv::Point anchor_target, int subpatch_size);

  static std::vector<MergePatch> from_patches(std::vector<Patch>& patches_in, int subpatch_size, int target_id);
  static mat<std::vector<int>> get_subpatch_index_mat(std::vector<MergePatch>& patches, int subpatch_size);

  cv::Mat draw(double scale) const;

  static cv::Mat draw(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, double scale, bool draw_boundaries);
  static cv::Mat draw_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale, bool draw_boundaries);
  static cv::Mat draw_fullres(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale, bool draw_boundaries);

  static std::vector<cv::Mat> draw_source_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale);
  static std::vector<cv::Mat> draw_source_fullres(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale);

  static cv::Mat draw_rect(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, double scale, int subpatch_size, bool draw_boundaries);
  static cv::Mat draw_rect_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale, int subpatch_size, bool draw_boundaries);
  static cv::Mat draw_rect_fullres(const std::vector<MergePatch>& patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale, int subpatch_size, bool draw_boundaries);

  static cv::Mat draw_half_rect_fullres(const std::vector<MergePatch>& patches, const boost::filesystem::path& base_path, const std::vector<Texture>& textures, double scale, int subpatch_size, bool draw_boundaries);
  static cv::Mat draw_half_rect_fullres(std::vector<MergePatch> patches, const std::vector<Texture>& textures, const std::vector<Texture>& textures_fullres, double scale, int subpatch_size, bool draw_boundaries);

  cv::Mat get_error_mat(int y_subpatch, int x_subpatch) const;
  cv::Mat get_error_mat(cv::Rect region_global) const;

  cv::Mat get_active_pixel(int y_subpatch, int x_subpatch) const;
  cv::Mat get_active_pixel(cv::Rect region_global) const;

  std::string get_id_string() const
  {
    std::string id_string;
    id_string += output_characters[target_id % output_characters.size()];
    id_string += output_characters[pos_patch.y % output_characters.size()];
    id_string += output_characters[pos_patch.x % output_characters.size()];
    return id_string;
  }

  void trim_bezier_curves();

  static void merge_patches_x(std::vector<MergePatch>& patches, const mat<std::vector<int>>& patch_indices, int x, int y_start, int y_end, int subpatch_size, int start_top=-1, int start_bottom=-1);
  static void merge_patches_y(std::vector<MergePatch>& patches, const mat<std::vector<int>>& patch_indices, int y, int x_start, int x_end, int subpatch_size, int start_left=-1, int start_right=-1);
  static void merge_patches_cross(std::vector<MergePatch>& patches, const std::vector<int> cross_indices, int subpatch_size);

  static void merge_patches_x(MergePatch& patch_left, cv::Rect region_left, MergePatch& patch_right, cv::Rect region_right, int start_top = -1, int start_bottom = -1);
  static void merge_patches_y(MergePatch& patch_top, cv::Rect region_top, MergePatch& patch_bottom, cv::Rect region_bottom, int start_left = -1, int start_right = -1);
  static void merge_patches_cross(MergePatch& patch_tl, cv::Rect region_tl, MergePatch& patch_bl, cv::Rect region_bl, MergePatch& patch_tr, cv::Rect region_tr, MergePatch& patch_br, cv::Rect region_br);

  static void merge_patches_x(MergePatch& patch_left, MergePatch& patch_right, int y_subpatch, int x_subpatch);
  static void merge_patches_y(MergePatch& patch_top, MergePatch& patch_bottom, int y_subpatch, int x_subpatch);

  static cv::Mat merge_patches_x(cv::Mat error_left, cv::Mat error_right, int start_top = -1, int start_bottom = -1);
  static cv::Mat merge_patches_y(cv::Mat error_top, cv::Mat error_bottom, int start_left = -1, int start_right = -1);

  static std::vector<int> get_cut_coordinates_horiz(const std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat, int y_subpatch, int x_subpatch, int subpatch_size);

  static std::vector<int> get_cut_coordinates_vert(const std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat, int y_subpatch, int x_subpatch, int subpatch_size);

  std::vector<int> get_cut_coordinates_horiz(int y_subpatch, int x_subpatch) const;
  std::vector<int> get_cut_coordinates_vert(int y_subpatch, int x_subpatch) const;

  void add_curve_horiz(int y_subpatch, int x_subpatch, const BezierCurve& curve);
  void add_curve_vert(int y_subpatch, int x_subpatch, const BezierCurve& curve);

  static void fix_curve_corners(std::vector<MergePatch>& patches, const mat<std::vector<int>>& subpatch_index_mat);
  static void translate_all_curves(std::vector<MergePatch>& patches, const cv::Point2d delta);

  static void draw_bezier_curves(cv::Mat image, const std::vector<MergePatch>& patches, double factor, const cv::Point2d& delta);
  static void draw_bezier_curves_source(std::vector<cv::Mat>& image, const std::vector<MergePatch>& patches, double factor);

  const std::vector<BezierCurve>& target_curves_top() const
  {
    return m_curves_top;
  }

  const std::vector<BezierCurve>& target_curves_bottom() const
  {
    return m_curves_bottom;
  }

  const std::vector<BezierCurve>& target_curves_left() const
  {
    return m_curves_left;
  }

  const std::vector<BezierCurve>& target_curves_right() const
  {
    return m_curves_right;
  }

  std::vector<BezierCurve> source_curves_top() const
  {
    std::vector<BezierCurve> curves = BezierCurve::translated(m_curves_top, anchor_source - anchor_target);
    return BezierCurve::transformed(curves, transformation_source_inv);
  }

  std::vector<BezierCurve> source_curves_bottom() const
  {
    std::vector<BezierCurve> curves = BezierCurve::translated(m_curves_bottom, anchor_source - anchor_target);
    return BezierCurve::transformed(curves, transformation_source_inv);
  }

  std::vector<BezierCurve> source_curves_left() const
  {
    std::vector<BezierCurve> curves = BezierCurve::translated(m_curves_left, anchor_source - anchor_target);
    return BezierCurve::transformed(curves, transformation_source_inv);
  }

  std::vector<BezierCurve> source_curves_right() const
  {
    std::vector<BezierCurve> curves = BezierCurve::translated(m_curves_right, anchor_source - anchor_target);
    return BezierCurve::transformed(curves, transformation_source_inv);
  }

  cv::Point center_subpatch() const
  {
    return (region_subpatch.tl() + region_subpatch.br()) / 2;
  }

  cv::Mat active_pixel;
  cv::Rect region_subpatch;
  int subpatch_size;

  cv::Point  pos_patch;

  cv::Size size;

  cv::Point anchor_target;
  cv::Point anchor_source;

  int source_index, source_rot;

  cv::Mat transformation_source, transformation_source_inv;

  cv::Mat error;

private:
  std::vector<BezierCurve> m_curves_top;
  std::vector<BezierCurve> m_curves_bottom;
  std::vector<BezierCurve> m_curves_left;
  std::vector<BezierCurve> m_curves_right;

  int target_id;

  static const std::string output_characters;
};

#endif /* TRLIB_MERGE_PATCH_HPP_ */