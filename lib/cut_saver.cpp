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

#include "cut_saver.hpp"

void CutSaver::add_sources_bezier_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel)
{
  add_sources_bezier<SVGSaver>(saver, path, textures, patches, material_panel);
}

void CutSaver::add_sources_bezier_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel)
{
  add_sources_bezier<EPSSaver>(saver, path, textures, patches, material_panel);
}

void CutSaver::add_sources_bezier_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<Patch>& patches, const std::vector<MaterialPanel>& material_panel)
{
  add_sources_bezier<SVGSaver>(saver, path, textures, patches, material_panel);
}

void CutSaver::add_sources_bezier_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<Patch>& patches, const std::vector<MaterialPanel>& material_panel)
{
  add_sources_bezier<EPSSaver>(saver, path, textures, patches, material_panel);
}

void CutSaver::save_sources_bezier_svg(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel, const cv::Size2d& table_dimensions_mm)
{
  save_sources_bezier<SVGSaver>(path, textures, patches, material_panel, table_dimensions_mm);
}

void CutSaver::save_sources_bezier_eps(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel, const cv::Size2d& table_dimensions_mm)
{
  save_sources_bezier<EPSSaver>(path, textures, patches, material_panel, table_dimensions_mm);
}

void CutSaver::add_target_bezier_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches)
{
  add_target_bezier<SVGSaver>(saver, path, texture, patches);
}

void CutSaver::add_target_bezier_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches)
{
  add_target_bezier<EPSSaver>(saver, path, texture, patches);
}

void CutSaver::save_target_bezier_svg(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  save_target_bezier<SVGSaver>(path, texture, patches, table_dimensions_mm);
}

void CutSaver::save_target_bezier_eps(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  save_target_bezier<EPSSaver>(path, texture, patches, table_dimensions_mm);
}

void CutSaver::add_sources_rect_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches)
{
  add_sources_rect<SVGSaver>(saver, path, textures, patches);
}

void CutSaver::add_sources_rect_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches)
{
  add_sources_rect<EPSSaver>(saver, path, textures, patches);
}

void CutSaver::save_sources_rect_svg(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  save_sources_rect<SVGSaver>(path, textures, patches, table_dimensions_mm);
}

void CutSaver::save_sources_rect_eps(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  save_sources_rect<EPSSaver>(path, textures, patches, table_dimensions_mm);
}

void CutSaver::add_target_rect_svg(std::vector<SVGSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches)
{
  add_target_rect<SVGSaver>(saver, path, texture, patches);
}

void CutSaver::add_target_rect_eps(std::vector<EPSSaver>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches)
{
  add_target_rect<EPSSaver>(saver, path, texture, patches);
}

void CutSaver::save_target_rect_svg(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  save_target_rect<SVGSaver>(path, texture, patches, table_dimensions_mm);
}

void CutSaver::save_target_rect_eps(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  save_target_rect<EPSSaver>(path, texture, patches, table_dimensions_mm);
}
