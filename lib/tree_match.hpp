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

#ifndef TRLIB_TREE_MATCH_HPP_
#define TRLIB_TREE_MATCH_HPP_

#include <deque>
#include <vector>

#include <boost/filesystem.hpp>
#include <opencv2/opencv.hpp>

#include "adaptive_patch.hpp"
#include "feature_evaluator.hpp"
#include "gabor_filter_bank.hpp"
#include "grid.hpp"
#include "patch.hpp"
#include "texture.hpp"

class TreeMatch
{
public:
  TreeMatch(int min_patch_size, int patch_levels, double patch_quality_factor, int filter_resolution, double frequency_octaves, int num_filter_directions);

  static TreeMatch load(const boost::filesystem::path& path, bool load_textures);

  void add_target(const boost::filesystem::path& path, double dpi, double scale);
  void add_texture(const boost::filesystem::path& path, double dpi, double scale, int num_rotations, const TextureMarker& marker=TextureMarker(), const std::string id=std::string());
  void add_texture(const boost::filesystem::path& path, const boost::filesystem::path& mask, double dpi, double scale, int num_rotations, const TextureMarker& marker=TextureMarker(), const std::string id=std::string());

  void generate_patches(int target_index, const Grid& morphed_grid, cv::Mat edge_image);
  void generate_patches_square(int target_index);
  void add_patches(int target_index, const std::vector<PatchRegion>& patches, double scale);

  void compute_responses(double weight_intensity, double weight_sobel, double weight_gabor, double histogram_matching_factor);

  bool find_next_patch();
  bool find_next_patch_adaptive();

  cv::Mat fit_single_patch(const boost::filesystem::path& path);

  cv::Mat draw(int target_index, bool draw_target) const;
  cv::Mat draw_matched_target(int target_index, double histogram_matching_factor) const;
  cv::Mat draw_masked_target(int target_index) const;
  cv::Mat draw_patch(const Patch& patch) const;
  std::pair<cv::Mat, cv::Mat> draw_saliency(int target_index) const;

  std::vector<cv::Mat> draw_masked_textures() const;
  std::vector<cv::Mat> draw_masked_textures_patch(const std::vector<Patch>& patch, cv::Scalar color, double alpha) const;
  std::vector<cv::Mat> draw_masked_textures_patch(const std::vector<Patch>& patch, const std::vector<cv::Scalar>& color, const std::vector<double>& alpha) const;

  std::vector<cv::Mat> draw_masked_textures_patch(const Patch& patches, cv::Scalar color, double alpha) const;
  std::vector<cv::Mat> draw_masked_textures_patch_last(const std::vector<Patch>& patches, cv::Scalar color_1, double alpha_1, cv::Scalar color_2, double alpha_2, double scale) const;

  const std::vector<Texture>& targets() const
  {
    return m_targets;
  }

  const std::vector<std::vector<Texture>>& textures() const
  {
    return m_textures;
  }

  const GaborFilterBank& filter_bank() const
  {
    return m_filter_bank;
  }

  int num_targets() const
  {
    return static_cast<int>(m_targets.size());
  }

  int num_textures() const
  {
    return static_cast<int>(m_textures.size());
  }

  void downsample(int factor);

  void save(int target_index, boost::filesystem::path path) const;
  void find_markers(double marker_size_mm, int num_marker);

  cv::Size max_filter_size(double weight_intensity, double weight_sobel, double weight_gabor) const;
  
  const std::vector<Patch>& patches() const
  {
    return m_patches;
  }

  const std::deque<PatchRegion>& reconstruction_regions() const
  {
    return m_reconstruction_regions;
  }

  void sort_patches_by_saliency();
  void sort_patches_by_center_distance();

private:
  void mask_patch_resources(const Patch& patch);
  void mask_patch_resources(const AdaptivePatch& adaptive_patch);
  void mask_patch_resources(const Patch& patch, const cv::Mat& mask);
  
  void unmask_patch_resources(const Patch& patch);
  void unmask_patch_resources(const AdaptivePatch& adaptive_patch);
  void unmask_patch_resources(const Patch& patch, const cv::Mat& mask);

  void add_patch(const Patch& match);
 
  std::vector<Patch> match_patch(const PatchRegion& region);
  Patch match_patch_impl(const PatchRegion& region, cv::Mat mask);

  static cv::Mat compute_priority_map(const cv::Mat& texture);

  std::vector<cv::Size> m_patch_sizes;
  double m_patch_quality_factor;
  cv::Size m_subpatch_size;

  std::vector<Texture> m_targets;
  std::vector<std::vector<Texture>> m_textures;

  std::vector<cv::Mat> m_target_images;
  
  std::deque<PatchRegion> m_reconstruction_regions;

  GaborFilterBank m_filter_bank;
  std::vector<Patch> m_patches;
};

#endif /* TRLIB_TREE_MATCH_HPP_ */