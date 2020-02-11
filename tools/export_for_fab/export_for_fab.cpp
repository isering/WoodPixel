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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "cut_saver.hpp"
#include "patch.hpp"
#include "material_panel.hpp"
#include "svg_saver.hpp"


namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

cv::Rect get_roi(const cv::Point& p, const cv::Size& image_size, const cv::Size& roi_size)
{
  cv::Point p_tl = p - cv::Point(roi_size);
  cv::Point p_br = p + cv::Point(roi_size);

  p_tl.x = std::max(p_tl.x, 0);
  p_tl.x = std::min(p_tl.x, image_size.width-1);
  p_tl.y = std::max(p_tl.y, 0);
  p_tl.y = std::min(p_tl.y, image_size.height-1);

  p_br.x = std::max(p_br.x, 0);
  p_br.x = std::min(p_br.x, image_size.width-1);
  p_br.y = std::max(p_br.y, 0);
  p_br.y = std::min(p_br.y, image_size.height-1);

  return cv::Rect(p_tl, p_br);
}

int main(int argc, char* argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show this help message")
    ("table_layout,t", po::value<fs::path>(), "Input JSON file with material panels on table (table.json)")
    ("verify_markers,v", "Visually verify marker positions in source textures")
    ("out,o", po::value<fs::path>(), "Output path")
    ("cut_pattern,c", po::value<std::vector<fs::path> >()->composing(), "Input JSON file(s) with cut pattern(s) (result.json)");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
  po::notify(vm);

  bool verify_markers = false;
  if (vm.count("verify_markers"))
  {
    verify_markers = true;
  }

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 0;
  }

  std::vector<fs::path> paths_cut_pattern;
  fs::path path_table_layout;
  fs::path path_out;

  std::vector<Patch> patches;
  std::vector<Texture> textures;
  std::vector<MaterialPanel> material_panels;
  cv::Size2d table_dimensions_mm;

  try
  {
    if (vm.count("cut_pattern"))
    {
      paths_cut_pattern = vm["cut_pattern"].as<std::vector<fs::path>>();

      for (const fs::path path_cut_pattern : paths_cut_pattern)
      {
        if (!fs::exists(path_cut_pattern) || !fs::is_regular_file(path_cut_pattern))
        {
          std::cerr << "Input JSON path " << path_cut_pattern << "does not exist or is no regular file." << std::endl
            << desc << std::endl;
          return -1;
        }
      }
    }
    else
    {
      std::cerr << "No cut_pattern JSON specified." << std::endl
        << desc << std::endl;
      return -1;
    }

    if (vm.count("table_layout"))
    {
      path_table_layout = vm["table_layout"].as<fs::path>();

      if (!fs::exists(path_table_layout) || !fs::is_regular_file(path_table_layout))
      {
        std::cerr << "Input JSON path does not exist or is no regular file." << std::endl
          << desc << std::endl;
        return -1;
      }
    }
    else
    {
      std::cerr << "No table_layout JSON specified." << std::endl
        << desc << std::endl;
      return -1;
    }

    if (vm.count("out"))
    {
      path_out = vm["out"].as<fs::path>();
      if (fs::exists(path_out) && fs::is_directory(path_out))
      {
        std::cerr << "Output path exists but is a directory." << std::endl;
        return -1;
      }
    }
    else
    {
      std::cerr << "No output path specified." << std::endl
        << desc << std::endl;
    }

    pt::ptree tree_cut_pattern;
    for (const fs::path path_cut_pattern : paths_cut_pattern)
    {
      pt::read_json(path_cut_pattern.string(), tree_cut_pattern);
      std::vector<Patch> newpatches = Serializable::deserialize_vec<Patch>(tree_cut_pattern, "patches", path_cut_pattern.parent_path());
      patches.insert(patches.end(), newpatches.begin(), newpatches.end());
    }

    const fs::path path_cut_pattern = paths_cut_pattern.front();
    textures = Serializable::deserialize_vec<Texture>(tree_cut_pattern, "textures_source", path_cut_pattern.parent_path());

    pt::ptree tree_table_layout;
    pt::read_json(path_table_layout.string(), tree_table_layout);
    material_panels = Serializable::deserialize_vec<MaterialPanel>(tree_table_layout, "panels", path_table_layout.parent_path());
    Serializable::deserialize(tree_table_layout, "table_dimensions_mm", table_dimensions_mm, path_table_layout.parent_path());

    /*
     * Check that all patches for validity.
     */
    for (const Patch& patch : patches)
    {
      if (!patch.region_target.valid())
      {
        std::cerr << "Loaded invalid patch." << std::endl
          << "target_index: " << patch.region_target.target_index() << std::endl
          << "coordinate: " << patch.region_target.coordinate() << std::endl
          << "Quitting..." << std::endl;
        return -1;
      }
    }

    std::cout << "Patches loaded: " << patches.size() << std::endl
      << "Textures loaded: " << textures.size() << std::endl;

    if (verify_markers)
    {
      for (const Texture& texture : textures)
      {
        fs::path path_texture = texture.filename;
        if (!fs::exists(path_texture))
        {
          path_texture = path_cut_pattern.parent_path() / path_texture;
        }
        cv::Mat image = cv::imread(path_texture.string());

        for (cv::Point2d p : texture.marker.markers_pix)
        {
          p /= texture.scale;

          cv::circle(image, cv::Point(p), 10, cv::Scalar(255, 0, 255), 1);

          cv::Mat image_roi;
          cv::Rect roi = get_roi(cv::Point(p), image.size(), cv::Size(64, 64));
          cv::resize(image(roi), image_roi, cv::Size(), 10.0, 10.0, cv::INTER_NEAREST);
          cv::imshow(texture.id, image_roi);
          cv::waitKey();
        }
        cv::destroyWindow(texture.id);
      }
    }

    if (!fs::exists(path_out.parent_path()))
    {
      fs::create_directories(path_out.parent_path());
    }

    std::vector<SVGSaver> svg_saver;
    CutSaver::add_sources_bezier_svg(svg_saver, path_out.parent_path(), textures, patches, material_panels);
    CutSaver::save(svg_saver, path_out.parent_path(), path_out.filename().string(), table_dimensions_mm);
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  return 0;
}
