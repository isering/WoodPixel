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

#include <iomanip>
#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "tree_match.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
  try
  {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "Show this help message")
      ("patches,p", po::value<fs::path>(), "Input folder containing patch images")
      ("json,j", po::value<fs::path>(), "Input JSON file containing texture description")
      ("out,o", po::value<fs::path>(), "Output directory");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
      std::cout << desc << std::endl;
      return 0;
    }

    std::vector<fs::path> files_patches;
    fs::path path_json, path_out;

    if (vm.count("patches"))
    {
      fs::path path_patches = vm["patches"].as<fs::path>();
      if (!fs::exists(path_patches) || !fs::is_directory(path_patches))
      {
        std::cerr << "Path to patches is not a directory." << std::endl
          << desc << std::endl;
        return -1;
      }

      fs::directory_iterator end;
      for (fs::directory_iterator iter(path_patches); iter != end; ++iter)
      {
        if (fs::is_regular_file(*iter))
        {
          if (iter->path().has_extension())
          {
            std::string extension = iter->path().extension().string();
            if (extension == ".jpg" || extension == ".jpeg" || extension == ".png" || extension == ".tiff")
            {
              files_patches.push_back(iter->path());
            }
          }
        }
      }

      if (files_patches.empty())
      {
        std::cerr << "No patches loaded." << std::endl
          << desc << std::endl;
        return -1;
      }
    }
    else
    {
      std::cerr << "No patch path specified." << std::endl
        << desc << std::endl;
      return -1;
    }

    if (vm.count("json"))
    {
      path_json = vm["json"].as<fs::path>();
      if (!fs::exists(path_json) || !fs::is_regular_file(path_json))
      {
        std::cerr << "JSON file does not exist or is no regular file: " << path_json << std::endl;
        return -1;
      }
    }
    else
    {
      std::cerr << "No JSON file specified." << std::endl
        << desc << std::endl;
    }

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
    }
    else
    {
      std::cerr << "No output path specified." << std::endl
        << desc << std::endl;
      return -1;
    }

    TreeMatch matcher = TreeMatch::load(path_json, true);

    std::cout << "Loaded " << files_patches.size() << " patches." << std::endl
      << "Loaded " << matcher.textures().size() << " textures." << std::endl;

    for (const fs::path file : files_patches)
    {
      cv::Mat patch = matcher.fit_single_patch(file);
      cv::imwrite((path_out / file.filename()).string(), patch);
    }
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  return 0;
}