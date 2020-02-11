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

#include <future>
#include <iostream>
#include <list>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "bilateral_filter_future.hpp"
#include "canny_future.hpp"
#include "ez_grid.hpp"
#include "grid.hpp"
#include "morph_grid_future.hpp"
#include "timer.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

static cv::Mat load_image(const po::variables_map& vm, const std::string& key, int imread_flags=1)
{
  fs::path path_image;
  if (vm.count(key))
  {
    path_image = vm[key].as<fs::path>();

    if (!fs::exists(path_image) || !fs::is_regular_file(path_image))
    {
      std::cerr << "Input image does not exist or is no regular file: " << key << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else
  {
    std::cerr << "No input image file specified: " << key << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return cv::imread(path_image.string(), imread_flags);
}

int main(int argc, char* argv[])
{
  /*try
  {*/
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "Show this help message")
      ("image,i", po::value<fs::path>(), "Input image")
      ("filtered,f", po::value<fs::path>(), "Filtered input image")
      ("size,s", po::value<int>(), "Number of pixels per grid cell")
      ("out,o", po::value<fs::path>(), "Output JSON file")
      ("load_full_state", po::value<fs::path>(), "Load full state using json file")
      ("load_partial_state", po::value<fs::path>(), "Load partial state using json file (edges, masks)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
      std::cout << desc << std::endl;
      return 0;
    }

    cv::Mat image = load_image(vm, "image");
    cv::Mat image_filtered = load_image(vm, "filtered", cv::IMREAD_GRAYSCALE);

    fs::path path_partial_state;
    if (vm.count("load_partial_state"))
    {
      path_partial_state = vm["load_partial_state"].as<fs::path>();

      if (!fs::exists(path_partial_state) || !fs::is_regular_file(path_partial_state))
      {
        std::cerr << "Partial state JSON path does not exist or is no regular file." << std::endl
          << desc << std::endl;
        return -1;
      }
    }

    fs::path path_full_state;
    if (vm.count("load_full_state"))
    {
      path_full_state = vm["load_full_state"].as<fs::path>();

      if (!fs::exists(path_full_state) || !fs::is_regular_file(path_full_state))
      {
        std::cerr << "Full state JSON path does not exist or is no regular file." << std::endl
          << desc << std::endl;
        return -1;
      }
    }

    fs::path path_out, path_out_parent;
    if (vm.count("out"))
    {
      path_out = vm["out"].as<fs::path>();
      path_out_parent = path_out.parent_path();
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

    int grid_size;
    if (vm.count("size"))
    {
      grid_size = vm["size"].as<int>();
    }
    else
    {
      std::cerr << "No grid cell size specified." << std::endl;
      return -1;
    }

    EZGrid ez_grid(image, image_filtered, grid_size, path_out);

    if (fs::exists(path_partial_state))
    {
      ez_grid.load_partial_state(path_partial_state);
    }

    if (fs::exists(path_full_state))
    {
      std::cout << "Load full state" << std::endl;
      ez_grid.load_full_state(path_full_state);
    }

    ez_grid.run();
    ez_grid.save();
  /*}
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }*/

  return 0;
}