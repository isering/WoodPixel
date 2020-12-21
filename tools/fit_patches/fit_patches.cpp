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

#include <iomanip>
#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "tree_match.hpp"

//const float pi = boost::math::constants::pi<float>();

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

fs::path get_unique_path(const fs::path &base_path)
{
  fs::path path_unique;
  int i = 0;
  do
  {
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << i++;
    path_unique = base_path / ss.str();
  } while (fs::exists(path_unique));
  return path_unique;
}

std::vector<Patch> load_patches(const std::vector<fs::path> &paths_patches)
{
  std::vector<Patch> patches_all;

  for (const fs::path &p : paths_patches)
  {
    pt::ptree root;
    pt::read_json(p.string(), root);
    const fs::path base_path = p.parent_path();
    const std::vector<Patch> patches = Serializable::deserialize_vec<Patch>(root, "patches", base_path);
    patches_all.insert(patches_all.end(), patches.begin(), patches.end());
  }

  return patches_all;
}

int main(int argc, char *argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "Show this help message")("in,i", po::value<fs::path>(), "Input JSON file")("out,o", po::value<fs::path>(), "Output directory")("vis,v", "Visualization")("steps,s", po::value<fs::path>(), "Intermediate output directory")("patches,p", po::value<std::vector<fs::path>>(), "Old patches for visualization");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 0;
  }

  fs::path path_in;
  if (vm.count("in"))
  {
    path_in = vm["in"].as<fs::path>();
    if (!fs::exists(path_in) || !fs::is_regular_file(path_in))
    {
      std::cerr << "Input JSON file does not exist or is no regular file." << std::endl;
      return -1;
    }
  }
  else
  {
    std::cerr << "No input JSON file specified." << std::endl
              << desc << std::endl;
    return -1;
  }

  fs::path path_out;
  if (vm.count("out"))
  {
    path_out = vm["out"].as<fs::path>();

    if (!fs::exists(path_out))
    {
      fs::create_directories(path_out);
    }

    if (!fs::exists(path_out) || !fs::is_directory(path_out))
    {
      std::cerr << "Unable to create output directory." << std::endl
                << desc << std::endl;
      return -1;
    }

    path_out = get_unique_path(path_out);
  }

  bool visualization = false;
  if (vm.count("vis"))
  {
    visualization = true;
  }

  bool steps = false;
  fs::path path_steps;
  if (vm.count("steps"))
  {
    steps = true;
    path_steps = vm["steps"].as<fs::path>();

    if (!fs::exists(path_steps / "image"))
    {
      fs::create_directories(path_steps / "image");
    }

    if (!fs::exists(path_steps / "image_target"))
    {
      fs::create_directories(path_steps / "image_target");
    }

    if (!fs::exists(path_steps / "patch"))
    {
      fs::create_directories(path_steps / "patch");
    }

    if (!fs::exists(path_steps / "textures"))
    {
      fs::create_directories(path_steps / "textures");
    }
  }

  std::vector<fs::path> paths_patches;
  if (vm.count("patches"))
  {
    paths_patches = vm["patches"].as<std::vector<fs::path>>();
    for (const fs::path &p : paths_patches)
    {
      if (!fs::exists(p) || !fs::is_regular_file(p))
      {
        std::cerr << "Input JSON file does not exist or is no regular file: " << p << std::endl;
        return -1;
      }
    }
  }

  std::vector<Patch> patches_old = load_patches(paths_patches);

  TreeMatch matcher = TreeMatch::load(path_in, true);

  for (int i = 0; i < matcher.num_targets(); ++i)
  {
    std::cout << "Target " << i << " approx. size in mm: "
              << matcher.targets()[i].texture.size().width * 25.4 / matcher.textures()[0][0].dpi << " x "
              << matcher.targets()[i].texture.size().height * 25.4 / matcher.textures()[0][0].dpi << std::endl;
  }

  /*
  cv::Mat image_debug;
  matcher.targets()[0].texture.convertTo(image_debug, CV_8UC3, 1.0/255.0);
  PatchRegion::draw(image_debug, cv::Vec3b(255, 0, 255), 1.0, matcher.reconstruction_regions().begin(), matcher.reconstruction_regions().end());
  cv::imshow("image_debug", image_debug);
  cv::waitKey();
  std::exit(EXIT_SUCCESS);
  */

  /*
  if (!fs::exists(path_out))
  {
  fs::create_directories(path_out);
  }
  cv::Mat saliency_image, saliency_color_map;
  std::tie(saliency_image, saliency_color_map) = matcher.draw_saliency(0);
  cv::imwrite((path_out / "saliency_map.png").string(), saliency_image);
  cv::imwrite((path_out / "saliency_map_colormap.png").string(), saliency_color_map);
  std::exit(EXIT_SUCCESS);
  */

  /*
  matcher.target().response.render();
  matcher.textures()[0][0].response.render();
  */
  try
  {
    while (matcher.find_next_patch_adaptive())
    {
      if (steps)
      {
        static int iteration;
        const std::string filename = (boost::format("%08d.png") % iteration++).str();

        for (int i = 0; i < matcher.num_targets(); ++i)
        {
          cv::Mat image_patches = matcher.draw(i, false);
          cv::imwrite((path_steps / (boost::format("image_%04d_%s") % i % filename).str()).string(), image_patches);

          cv::Mat image_patches_target = matcher.draw(i, true);
          cv::imwrite((path_steps / (boost::format("image_target_%04d_%s") % i % filename).str()).string(), image_patches_target);
        }

        cv::Mat patch = matcher.draw_patch(matcher.patches().back());
        cv::imwrite((path_steps / "patch" / filename).string(), patch);

        std::vector<Patch> patches_temp = patches_old;
        patches_temp.insert(patches_temp.end(), matcher.patches().begin(), matcher.patches().end());

        std::vector<cv::Mat> masked_textures = matcher.draw_masked_textures_patch_last(
            patches_temp,
            cv::Scalar::all(0.0), 0.5,
            cv::Scalar(0.0, 1.0, 0.0), 0.5,
            1.0);
        for (size_t i = 0; i < masked_textures.size(); ++i)
        {
          cv::Mat masked_texture_scaled;
          cv::resize(masked_textures[i], masked_texture_scaled, cv::Size(), 0.5, 0.5);
          cv::imwrite((path_steps / "textures" / (std::to_string(i) + "_" + filename)).string(), masked_texture_scaled);
        }
      }

      if (visualization)
      {
        for (int i = 0; i < matcher.num_targets(); ++i)
        {
          cv::Mat image = matcher.draw(i, true);
          cv::imshow((boost::format("Image %04d") % i).str(), image);
        }

        /*
        const Patch& last_patch = matcher.patches().back();
        cv::Mat mask = cv::Mat::zeros(last_patch.region_target.size(), CV_8UC1);
        last_patch.region_target.draw_edge_mask(mask, last_patch.region_target.bbox_tl());
        cv::resize(mask, mask, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);
        cv::imshow("Mask", mask);
        */

        // If Ctrl+Escape, stop the matching and terminate
        const int key = cv::waitKey(1);
        if (((key & 0xff) == 0x1b))
        {
          break;
        }
      }
    }

    for (int i = 0; i < matcher.num_targets(); ++i)
    {
      matcher.save(i, path_out);
    }

    fs::copy_file(path_in, path_out / path_in.filename());
  }
  catch (std::exception &e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  return 0;
}