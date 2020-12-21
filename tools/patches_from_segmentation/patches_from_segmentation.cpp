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
#include <random>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <opencv2/opencv.hpp>

#include "bezier_curve.hpp"
#include "blob.hpp"
#include "patch_region.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

static cv::Vec3b hsv_to_bgr(int hue, int saturation, int value)
{
  cv::Mat hsv_mat(1, 1, CV_8UC3, cv::Vec3b(hue, saturation, value));
  cv::Mat rgb_mat;
  cv::cvtColor(hsv_mat, rgb_mat, cv::COLOR_HSV2BGR);
  return cv::Vec3b(rgb_mat.data[0], rgb_mat.data[1], rgb_mat.data[2]);
}

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
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show this help message")
    ("in,i", po::value<fs::path>(), "Input image")
    ("out,o", po::value<fs::path>(), "Output JSON file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 0;
  }

  cv::Mat image = load_image(vm, "in", cv::IMREAD_GRAYSCALE);

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

  double max_val;
  cv::minMaxLoc(image, 0, &max_val);
  
  const int num_level = static_cast<int>(max_val);
  std::cout << "num_level: " << num_level << std::endl;

  std::vector<Blob> blobs = Blob::detect(image, 255);
  std::cout << "num_blobs: " << blobs.size() << std::endl;

  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, 180);

  cv::Mat blob_image(image.size(), CV_8UC3);

  for (const Blob& blob : blobs)
  {
    cv::Vec3b color = hsv_to_bgr(distribution(generator), 128, 240);
    for (const cv::Point& p : blob.points())
    {
      blob_image.at<cv::Vec3b>(p) = color;
    }
  }

  std::vector<PatchRegion> patches;
  for (const Blob& blob : blobs)
  {
    patches.emplace_back(1, blob.points());
  }

  pt::ptree root;
  pt::ptree tree_patches;
  for (const PatchRegion& patch : patches)
  {
    tree_patches.push_back(std::make_pair("", patch.save(path_out_parent, "patches")));
  }
  root.add_child("patches", tree_patches);
  pt::write_json(path_out.string(), root);

  cv::imshow("Image", image *(255.0/max_val));
  cv::imshow("Blob Image", blob_image);
  cv::waitKey();

  return 0;
}