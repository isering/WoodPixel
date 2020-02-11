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

#ifndef TRLIB_CUT_SAVER_HPP_
#define TRLIB_CUT_SAVER_HPP_

#include <ostream>
#include <vector>

#include <boost/filesystem.hpp>

#include "bezier_curve.hpp"
#include "eps_saver.hpp"
#include "material_panel.hpp"
#include "merge_patch.hpp"
#include "svg_saver.hpp"
#include "texture.hpp"

namespace CutSaver
{
  template <typename SaverT>
  void save(const std::vector<SaverT>& saver, const boost::filesystem::path& base_path, const std::string& base_filename, const cv::Size2d& table_dimensions_mm);

  template <typename SaverT>
  void add_sources_bezier(std::vector<SaverT>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel);
  void add_sources_bezier_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel);
  void add_sources_bezier_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel);

  template <typename SaverT>
  void add_sources_bezier(std::vector<SaverT>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<Patch>& patches, const std::vector<MaterialPanel>& material_panel);
  void add_sources_bezier_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<Patch>& patches, const std::vector<MaterialPanel>& material_panel);
  void add_sources_bezier_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<Patch>& patches, const std::vector<MaterialPanel>& material_panel);

  template <typename SaverT>
  void save_sources_bezier(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel, const cv::Size2d& table_dimensions_mm);
  void save_sources_bezier_svg(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel, const cv::Size2d& table_dimensions_mm);
  void save_sources_bezier_eps(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel, const cv::Size2d& table_dimensions_mm);

  template <typename SaverT>
  void add_target_bezier(std::vector<SaverT>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches);
  void add_target_bezier_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches);
  void add_target_bezier_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches);

  template <typename SaverT>
  void save_target_bezier(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);
  void save_target_bezier_svg(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);
  void save_target_bezier_eps(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);

  template <typename SaverT>
  void add_sources_rect(std::vector<SaverT>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches);
  void add_sources_rect_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches);
  void add_sources_rect_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches);

  template <typename SaverT>
  void save_sources_rect(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);
  void save_sources_rect_svg(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);
  void save_sources_rect_eps(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);

  template <typename SaverT>
  void add_target_rect(std::vector<SaverT>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches);
  void add_target_rect_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches);
  void add_target_rect_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches);

  template <typename SaverT>
  void save_target_rect(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);
  void save_target_rect_svg(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);
  void save_target_rect_eps(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm);
};

#include "cut_saver_impl.hpp"

#endif /* TRLIB_CUT_SAVER_HPP_ */