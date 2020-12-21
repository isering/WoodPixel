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

#include "patch.hpp"
#include "texture.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

int main(int argc, char *argv[])
{
  try
  {

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Show this help message")("in,i", po::value<fs::path>(), "Input JSON file")("out,o", po::value<fs::path>(), "Output path")("scale,s", po::value<double>(), "Output image scale");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
      std::cout << desc << std::endl;
      return 0;
    }

    fs::path path_in, base_path;
    if (vm.count("in"))
    {
      path_in = vm["in"].as<fs::path>();
      base_path = path_in.parent_path();

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
        std::cerr << "Unable to create output directory." << std::endl;
        return -1;
      }
    }
    else
    {
      std::cerr << "No output path specified." << std::endl
                << desc << std::endl;
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

    pt::ptree root;
    pt::read_json(path_in.string(), root);

    std::vector<Patch> patches = Serializable::deserialize_vec<Patch>(root, "patches", path_in.parent_path());
    std::vector<Texture> textures = Serializable::deserialize_vec<Texture>(root, "textures_source", path_in.parent_path());

    std::vector<Texture> textures_fullres;
    for (const Texture &t : textures)
    {
      textures_fullres.emplace_back(base_path / t.filename, t.dpi, 1.0);
    }

    cv::Rect bbox = Patch::bounding_box(patches, scale);
    cv::Mat image(bbox.size(), CV_8UC3);
    cv::Mat edge_mask = cv::Mat::zeros(bbox.size(), CV_8UC1);
    for (const Patch &p : patches)
    {
      p.draw(image, scale, textures_fullres[p.source_index], textures[p.source_index].scale);
      p.draw_edge_mask(edge_mask, scale);
    }

    cv::imwrite((path_out / "image.png").string(), image);

    const cv::Vec3b dark_brown(14, 29, 43);
    edge_mask.convertTo(edge_mask, CV_32FC1, 1.0 / 255.0);
    cv::GaussianBlur(edge_mask, edge_mask, cv::Size(5, 5), 0.0);

    for (int y = 0; y < image.rows; ++y)
    {
      cv::Vec3b *ptr_image = reinterpret_cast<cv::Vec3b *>(image.ptr(y));
      const float *ptr_mask = reinterpret_cast<const float *>(edge_mask.ptr(y));
      for (int x = 0; x < image.cols; ++x)
      {
        ptr_image[x] = ptr_mask[x] * dark_brown + (1.0f - ptr_mask[x]) * ptr_image[x];
      }
    }

    cv::imwrite((path_out / "image_boundary.png").string(), image);
  }
  catch (std::exception &e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  return 0;
}