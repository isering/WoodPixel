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

#ifndef TRLIB_CUT_SAVER_IMPL_HPP_
#define TRLIB_CUT_SAVER_IMPL_HPP_

#include "cut_saver.hpp"

#include "affine_transformation.hpp"

template <typename SaverT>
void CutSaver::save(const std::vector<SaverT>& saver, const boost::filesystem::path& base_path, const std::string& base_filename, const cv::Size2d& table_dimensions_mm)
{
  for (size_t i = 0; i < saver.size(); ++i)
  {
    if (saver[i].get_output_name() == "")
    {
      saver[i].save(base_path / (boost::format("%s_saver%04d%s") % base_filename % i % saver[i].extension()).str(), table_dimensions_mm);
    }
    else
    {
      saver[i].save(base_path / (boost::format("%s_%s%s") % base_filename % saver[i].get_output_name() % saver[i].extension()).str(), table_dimensions_mm);
    }
  }
}

template <typename SaverT>
void CutSaver::add_sources_bezier(std::vector<SaverT>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel)
{
  if (saver.empty())
  {
    for (const Texture& texture : textures)
    {
      saver.emplace_back(texture);
    }
  }

  for (const MergePatch& p : patches)
  {
    SaverT& s = saver[p.source_index];
    const Texture& t = textures[p.source_index];

    MergePatch p_trimmed = p;
    p_trimmed.trim_bezier_curves();

    /*
     * Add text.
     */
    cv::Point2f source_tl = AffineTransformation::transform<int, float>(p.transformation_source_inv, p.anchor_source);
    cv::Point2f source_br = AffineTransformation::transform<int, float>(p.transformation_source_inv, p.anchor_source + cv::Point(p.size));
    cv::Point2f pos_text(0.5f * (source_tl + source_br));
    //pos_text.x = static_cast<float>(t.texture.cols - pos_text.x);
    s.add_text(pos_text, p.get_id_string(), AffineTransformation::rotation_rad(p.transformation_source_inv));

    /*
     * Get the four Bezier contours, bring to right orientation, trim and add to saver.
     */
    std::vector<BezierCurve> curves_top = p_trimmed.source_curves_top();
    std::vector<BezierCurve> curves_right = p_trimmed.source_curves_right();

    std::vector<BezierCurve> curves_bottom = p_trimmed.source_curves_bottom();
    std::reverse(curves_bottom.begin(), curves_bottom.end());
    for (BezierCurve& c : curves_bottom)
    {
      c = c.reversed();
    }

    std::vector<BezierCurve> curves_left = p_trimmed.source_curves_left();
    std::reverse(curves_left.begin(), curves_left.end());
    for (BezierCurve& c : curves_left)
    {
      c = c.reversed();
    }

    std::vector<BezierCurve> curves_complete = curves_top;
    curves_complete.insert(curves_complete.end(), curves_right.begin(), curves_right.end());
    curves_complete.insert(curves_complete.end(), curves_bottom.begin(), curves_bottom.end());
    curves_complete.insert(curves_complete.end(), curves_left.begin(), curves_left.end());

    s.add_bezier(curves_complete);
  }
}

static std::vector<BezierCurve> get_bezier_curves(const PatchRegion& patch)
{
  std::vector<std::vector<BezierCurve>> curves;
  if (patch.has_diag_curve())
  {
    if (patch.diag_origin_tl())
    {
      if (!patch.curves_top().empty())
      {
        curves.push_back(patch.curves_diag_reversed());
        curves.push_back(patch.curves_top());
        curves.push_back(patch.curves_right());
      }
      else
      {
        curves.push_back(patch.curves_diag_reversed());
        curves.push_back(patch.curves_left());
        curves.push_back(patch.curves_bot());
      }
    }
    else
    {
      if (!patch.curves_top().empty())
      {
        curves.push_back(patch.curves_diag());
        curves.push_back(patch.curves_top_reversed());
        curves.push_back(patch.curves_left());
      }
      else
      {
        curves.push_back(patch.curves_diag());
        curves.push_back(patch.curves_right());
        curves.push_back(patch.curves_bot_reversed());
      }
    }
  }
  else
  {
    curves.push_back(patch.curves_top());
    curves.push_back(patch.curves_right());
    curves.push_back(patch.curves_bot_reversed());
    curves.push_back(patch.curves_left_reversed());
  }

  std::vector<BezierCurve> curves_out;
  for (const std::vector<BezierCurve>& curve : curves)
  {
    if (curve.empty())
    {
      throw(std::runtime_error("get_bezier_curves: Got empty curve where non-empty curve was expected."));
    }
    curves_out.insert(curves_out.end(), curve.begin(), curve.end());
  }

  return curves_out;
}

static std::vector<MaterialPanel>::const_iterator find_panel(int source_index, const std::vector<Texture>& textures, const std::vector<MaterialPanel>& material_panels)
{
  const std::string id = textures[source_index].id;

  for (std::vector<MaterialPanel>::const_iterator iter = material_panels.begin(); iter != material_panels.end(); ++iter)
  {
    if (id == iter->id)
    {
      return iter;
    }
  }
  return material_panels.end();
}

template <typename SaverT>
void CutSaver::add_sources_bezier(std::vector<SaverT>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<Patch>& patches, const std::vector<MaterialPanel>& material_panels)
{
  // First, determine which textures are active in current config.
  std::vector<int> source_ids;
  for (int i = 0; i < static_cast<int>(textures.size()); ++i)
  {
    // Check if texture is in current list of material panels.
    if (find_panel(i, textures, material_panels) == material_panels.end())
    {
      continue;
    }

    if (std::find(source_ids.begin(), source_ids.end(), i) == source_ids.end())
    {
      source_ids.push_back(i);
    }
  }

  // Then, determine which target_ids are active in current config.
  const size_t num_textures_used = source_ids.size();
  std::vector<int> target_ids;
  for (const Patch& p : patches)
  {
    if (std::find(target_ids.begin(), target_ids.end(), p.region_target.target_index()) == target_ids.end())
    {
      target_ids.push_back(p.region_target.target_index());
    }
  }

  if (saver.empty())
  {
    for (int i = 0; i < num_textures_used; ++i)
    {
      saver.emplace_back(textures[source_ids[i]], "panel_" + textures[source_ids[i]].id);
      std::cout << "Adding texture " << textures[source_ids[i]].id << std::endl;
    }
  }

  for (int i : target_ids)
  {
    saver.emplace_back(SaverT("target_" + std::to_string(i)));
    std::cout << "Adding target " << i << std::endl;
  }

  for (const Patch& p : patches)
  {
    std::vector<MaterialPanel>::const_iterator panel = find_panel(p.source_index, textures, material_panels);
    if (panel == material_panels.end())
    {
      continue;
    }
 
    const size_t source_index = std::distance(source_ids.begin(), std::find(source_ids.begin(), source_ids.end(), p.source_index));
    if (source_index >= source_ids.size())
    {
      continue;
    }
    SaverT& s = saver[source_index];
 
    const size_t target_index = std::distance(target_ids.begin(), std::find(target_ids.begin(), target_ids.end(), p.region_target.target_index()));
    if (target_index < 0)
    {
      continue;
    }    
    SaverT& t = saver[num_textures_used + target_index];

    const Texture& tex = textures[p.source_index];
    const cv::Mat T_texture_table = AffineTransformation::fit(tex.marker.markers_pix, panel->markers_table_mm);

    if (p.region_target.has_sub_regions())
    {
      throw(std::invalid_argument("CutSaver::add_sources_bezier: only patches without sub-regions expected here."));
    }

    std::vector<BezierCurve> curves_target = get_bezier_curves(p.region_target);
    cv::Mat T = AffineTransformation::concat(
      p.transformation_source_inv,
      AffineTransformation::T_translate(p.anchor_source.x - p.anchor_target().x, p.anchor_source.y - p.anchor_target().y));
    T = AffineTransformation::concat(T_texture_table, T);

    std::vector<BezierCurve> curves_source = BezierCurve::transformed(curves_target, T);

    s.add_bezier(curves_source);
    t.add_bezier(curves_source);
    
    /*
     * Add text.
     * Get barycenter of curve midpoints as text anchor.
     */
    cv::Point2f text_pos(0.0f, 0.0f);
    for (int j = 0; j < curves_source.size(); ++j)
    {
      text_pos += cv::Point2f(curves_source[j].eval(0.5));
    }
    text_pos /= static_cast<float>(curves_source.size());

    cv::Mat T_text = AffineTransformation::concat(T_texture_table, p.transformation_source_inv);
    cv::Point2f source_tl = AffineTransformation::transform<int, float>(T_text, p.anchor_source);
    cv::Point2f source_br = AffineTransformation::transform<int, float>(T_text, p.anchor_source + cv::Point(p.size()));

    s.add_text(text_pos, p.region_target.get_short_id_string(), AffineTransformation::rotation_rad(T_text));
    t.add_text(text_pos, p.region_target.get_short_id_string(), AffineTransformation::rotation_rad(T_text));
  }

  for (SaverT& s : saver)
  {
    for (const MaterialPanel& panel : material_panels)
    {
      for (const cv::Point2d& p : panel.markers_table_mm)
      {
        s.add_debug_circle(cv::Point2f(p), 5.0);
      }
    }
  }
}

template <typename SaverT>
void CutSaver::save_sources_bezier(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const std::vector<MaterialPanel>& material_panel, const cv::Size2d& table_dimensions_mm)
{
  std::vector<SaverT> saver;
  add_sources_bezier(saver, path, textures, patches, material_panel);
  save(saver, path, "source_bezier", table_dimensions_mm);
}

template <typename SaverT>
void CutSaver::add_target_bezier(std::vector<SaverT>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches)
{
  saver.emplace_back(texture);

  for (const MergePatch& p : patches)
  {
    MergePatch p_trimmed = p;
    p_trimmed.trim_bezier_curves();

    /*
    * Add text.
    */
    cv::Point2f pos_text(cv::Point2f(p.anchor_target) + 0.5f * cv::Point2f(p.size));
    //pos_text.x = static_cast<float>(texture.texture.cols / texture.scale - pos_text.x);
    saver.back().add_text(pos_text, p.get_id_string(), 0.0);

    /*
    * Get the four Bezier contours, bring to right orientation, trim and add to saver.
    */
    const std::vector<BezierCurve> curves_top = p_trimmed.target_curves_top();
    const std::vector<BezierCurve> curves_right = p_trimmed.target_curves_right();

    std::vector<BezierCurve> curves_bottom = p_trimmed.target_curves_bottom();
    std::reverse(curves_bottom.begin(), curves_bottom.end());
    for (BezierCurve& c : curves_bottom)
    {
      c = c.reversed();
    }

    std::vector<BezierCurve> curves_left = p_trimmed.target_curves_left();
    std::reverse(curves_left.begin(), curves_left.end());
    for (BezierCurve& c : curves_left)
    {
      c = c.reversed();
    }

    std::vector<BezierCurve> curves_complete = curves_top;
    curves_complete.insert(curves_complete.end(), curves_right.begin(), curves_right.end());
    curves_complete.insert(curves_complete.end(), curves_bottom.begin(), curves_bottom.end());
    curves_complete.insert(curves_complete.end(), curves_left.begin(), curves_left.end());

    saver.back().add_bezier(curves_complete);
  }
}

template <typename SaverT>
void CutSaver::save_target_bezier(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  std::vector<SaverT> saver;
  add_target_bezier(saver, path, texture, patches);
  save(saver, path, "target_bezier", table_dimensions_mm);
}

template <typename SaverT>
void CutSaver::add_sources_rect(std::vector<SaverT>& saver, const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches)
{
  if (saver.empty())
  {
    for (const Texture& texture : textures)
    {
      saver.emplace_back(texture);
    }
  }

  for (const MergePatch& p : patches)
  {
    SaverT& s = saver[p.source_index];
    const Texture& t = textures[p.source_index];

    /*
    * Add text.
    */
    cv::Point2f source_tl = AffineTransformation::transform<int, float>(p.transformation_source_inv, p.anchor_source);
    cv::Point2f source_br = AffineTransformation::transform<int, float>(p.transformation_source_inv, p.anchor_source + cv::Point(p.size));
    cv::Point2f pos_text(0.5f * (source_tl + source_br));
    //pos_text.x = static_cast<float>(t.texture.cols - pos_text.x);
    s.add_text(pos_text, p.get_id_string(), AffineTransformation::rotation_rad(p.transformation_source_inv));

    /*
    * Add rectangle patch.
    */
    // FIXME
    const float tl_x = p.anchor_source.x + 0.125f * p.size.width;
    const float tl_y = p.anchor_source.y + 0.125f * p.size.height;
    const float br_x = p.anchor_source.x + 0.875f * p.size.width;
    const float br_y = p.anchor_source.y + 0.875f * p.size.height;

    std::vector<cv::Point2f> polygon({ { tl_x, tl_y }, { br_x, tl_y }, { br_x, br_y }, { tl_x, br_y } });
    s.add_polygon(AffineTransformation::transform(p.transformation_source_inv, polygon));
  }
}

template <typename SaverT>
static void CutSaver::save_sources_rect(const boost::filesystem::path& path, const std::vector<Texture>& textures, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  std::vector<SaverT> saver;
  add_sources_rect(saver, path, textures, patches);
  save(saver, path, "source_rect", table_dimensions_mm);
}

template <typename SaverT>
void CutSaver::add_target_rect(std::vector<SaverT>& saver, const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches)
{
  saver.emplace_back(texture);

  for (const MergePatch& p : patches)
  {
    /*
    * Add text.
    */
    cv::Point2f pos_text(cv::Point2f(p.anchor_target) + 0.5f * cv::Point2f(p.size));
    //pos_text.x = static_cast<float>(texture.texture.cols - pos_text.x);
    saver.back().add_text(pos_text, p.get_id_string(), 0.0);

    /*
    * Add rectangle polygon.
    */
    // FIXME
    float tl_x = p.anchor_target.x + 0.125f * p.size.width;
    float tl_y = p.anchor_target.y + 0.125f * p.size.height;
    float br_x = p.anchor_target.x + 0.875f * p.size.width;
    float br_y = p.anchor_target.y + 0.875f * p.size.height;

    std::vector<cv::Point2f> polygon({ { tl_x, tl_y }, { br_x, tl_y }, { br_x, br_y }, { tl_x, br_y } });
    saver.back().add_polygon(polygon);
  }
}

template <typename SaverT>
void CutSaver::save_target_rect(const boost::filesystem::path& path, const Texture& texture, const std::vector<MergePatch>& patches, const cv::Size2d& table_dimensions_mm)
{
  std::vector<SaverT> saver;
  add_target_rect(saver, path, texture, patches);
  save(saver, path, "target_rect", table_dimensions_mm);
}

#endif /* TRLIB_CUT_SAVER_IMPL_HPP_ */