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

#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "tree_match.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

int main(int argc, char *argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "Show this help message")("in,i", po::value<fs::path>(), "Input JSON file")("out,o", po::value<fs::path>(), "Saliency image output filename")("scale,s", po::value<double>(), "Output image scale");

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
      std::cerr << "Input JSON file does not exist or is no regular file." << std::endl
                << desc << std::endl;
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
  bool write_out = false;
  if (vm.count("out"))
  {
    path_out = vm["out"].as<fs::path>();
    write_out = true;

    if (fs::exists(path_out) && !fs::is_regular_file(path_out))
    {
      std::cerr << "Output file already exists but is no regular file." << std::endl
                << desc << std::endl;
      return -1;
    }
  }

  double scale = 1.0;
  if (vm.count("scale"))
  {
    scale = vm["scale"].as<double>();
    if (scale <= 0.0)
    {
      std::cerr << "Output image scale has to be greater than zero." << std::endl
                << desc << std::endl;
      return -1;
    }
  }

  TreeMatch matcher = TreeMatch::load(path_in, false);

  cv::Mat saliency_map, saliency_color_map;
  std::tie(saliency_map, saliency_color_map) = matcher.draw_saliency(0);

  cv::imshow("Saliency map", saliency_map);
  cv::waitKey();

  if (write_out)
  {
    cv::imwrite(path_out.string(), saliency_map);
  }

  return 0;
}