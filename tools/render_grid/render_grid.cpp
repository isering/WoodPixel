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
#include <opencv2/opencv.hpp>
#include <opencv2/ximgproc.hpp>

#include "blob.hpp"
#include "patch.hpp"
#include "texture.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show this help message")
    ("grid,g", po::value<fs::path>(), "Input JSON containing grid")
    ("image,i", po::value<fs::path>(), "Input image")
    ("out,o", po::value<fs::path>(), "Output image");

  cv::Mat image;
  std::vector<PatchRegion> patches;

  bool has_path_out =  false;
  fs::path path_out;

  try
  {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
      std::cout << desc << std::endl;
      return 0;
    }

    if (vm.count("grid"))
    {
      fs::path path_grid = vm["grid"].as<fs::path>();

      if (!fs::exists(path_grid) || !fs::is_regular_file(path_grid))
      {
        std::cerr << "Input JSON does not exist or is no regular file." << std::endl;
        return -1;
      }

      pt::ptree tree_json;
      pt::read_json(path_grid.string(), tree_json);
      patches = Serializable::deserialize_vec<PatchRegion>(tree_json, "patches", path_grid.parent_path());
    }
    else
    {
      std::cerr << "No input image specified." << std::endl
        << desc << std::endl;
      return -1;
    }

    if (vm.count("image"))
    {
      fs::path path_image = vm["image"].as<fs::path>();

      if (!fs::exists(path_image) || !fs::is_regular_file(path_image))
      {
        std::cerr << "Input image does not exist or is no regular file." << std::endl;
        return -1;
      }

      image = cv::imread(path_image.string());

      if (image.empty())
      {
        std::cerr << "Unable to load image." << std::endl;
        return -1;
      }
    }

    if (vm.count("out"))
    {
      path_out = vm["out"].as<fs::path>();
      has_path_out = true;
      if (fs::exists(path_out) && !fs::is_regular_file(path_out))
      {
        std::cerr << "Output image exists but is no regular file." << std::endl;
        return -1;
      }

      if (path_out.has_parent_path() && !fs::exists(path_out.parent_path()))
      {
        fs::create_directories(path_out.parent_path());
      }
    }
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  cv::Mat mask = cv::Mat::zeros(image.size(), CV_8UC1);
  for (const PatchRegion& p : patches)
  {
    if (p.has_sub_regions())
    {
      for (const PatchRegion& sub_region : p.sub_regions())
      {
        sub_region.draw_edge_mask(image, cv::Vec3b(227, 178, 93), 0);
      }
    }
    else
    {
      p.draw_edge_mask(image, cv::Vec3b(227, 178, 93), 0);
    }
  }
  
  cv::imshow("Image", image);
  cv::waitKey();

  if (has_path_out)
  {
    cv::imwrite(path_out.string(), image);
  }

  return 0;
}