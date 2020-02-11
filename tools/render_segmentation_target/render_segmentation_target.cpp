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
    ("input,i", po::value<fs::path>(), "Input JSON (result.json)")
    ("segmentation,s", po::value<fs::path>(), "Input segmentation for overlay")
    ("out,o", po::value<fs::path>(), "Output image");

  fs::path path_out;
  cv::Mat segmentation;
  std::vector<Patch> patches;
  std::vector<Texture> textures_source;
  Texture texture_target;

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

    if (vm.count("input"))
    {
      fs::path path_input = vm["input"].as<fs::path>();

      if (!fs::exists(path_input) || !fs::is_regular_file(path_input))
      {
        std::cerr << "Input JSON does not exist or is no regular file." << std::endl;
        return -1;
      }

      pt::ptree tree_json;
      pt::read_json(path_input.string(), tree_json);

      patches = Serializable::deserialize_vec<Patch>(tree_json, "patches", path_input.parent_path());
      textures_source = Serializable::deserialize_vec<Texture>(tree_json, "textures_source", path_input.parent_path());
      texture_target.load(path_input.parent_path(), tree_json.get_child("texture_target"));
    }
    else
    {
      std::cerr << "No input image specified." << std::endl
        << desc << std::endl;
      return -1;
    }

    if (vm.count("segmentation"))
    {
      fs::path path_segmentation = vm["segmentation"].as<fs::path>();

      if (!fs::exists(path_segmentation) || !fs::is_regular_file(path_segmentation))
      {
        std::cerr << "Input segmentation does not exist or is no regular file." << std::endl;
        return -1;
      }

      segmentation = cv::imread(path_segmentation.string(), cv::IMREAD_GRAYSCALE);

      if (segmentation.empty())
      {
        std::cerr << "Unable to load segmentation image." << std::endl;
        return -1;
      }
    }

    if (vm.count("out"))
    {
      path_out = vm["out"].as<fs::path>();
      if (fs::exists(path_out) && !fs::is_regular_file(path_out))
      {
        std::cerr << "Output image exists but is no regular file." << std::endl;
        return -1;
      }

      if (!fs::exists(path_out.parent_path()))
      {
        fs::create_directories(path_out.parent_path());
      }
    }
    else
    {
      std::cerr << "No output path specified." << std::endl
        << desc << std::endl;
    }
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  cv::Mat image = cv::Mat::zeros(texture_target.texture.size(), CV_8UC3);
  cv::Mat mask = cv::Mat::zeros(texture_target.texture.size(), CV_8UC1);
  for (const Patch& p : patches)
  {
    p.draw(image, 1.0, textures_source[p.source_index], 1.0);
    
    cv::Mat patch_mask = p.region_target.mask();
    cv::Mat patch_mask_dilated;
    cv::dilate(patch_mask, patch_mask_dilated, cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3)));
    patch_mask = patch_mask_dilated - patch_mask;
    for (int y = 0; y < patch_mask.rows; ++y)
    {
      unsigned char* ptr_patch_mask = patch_mask.ptr(y);
      unsigned char* ptr_mask = mask.ptr(y+p.region_target.anchor().y);

      for (int x = 0; x < patch_mask.cols; ++x)
      {
        if (ptr_patch_mask[x])
        {
          ptr_mask[x+p.region_target.anchor().x] = ptr_patch_mask[x];
        }
      }
    }
  }

  const cv::Vec3b black(0, 0, 0);
  const cv::Vec3b dark_brown(14, 29, 43);
  
  for (int y = 0; y < image.rows; ++y)
  {
    cv::Vec3b* ptr_image = reinterpret_cast<cv::Vec3b*>(image.ptr(y));
    unsigned char* ptr_mask = reinterpret_cast<unsigned char*>(mask.ptr(y));
    for (int x = 0; x < image.cols; ++x)
    {
      if (ptr_image[x] == black)
      {
        ptr_image[x] = dark_brown;
        ptr_mask[x] = 255;
      }
    }
  }

  if (!segmentation.empty())
  {
    std::vector<Blob> blobs = Blob::detect(segmentation, 0);

    for (const Blob& blob : blobs)
    {
      blob.draw_contour(mask);
    }

    cv::ximgproc::thinning(mask, mask, cv::ximgproc::THINNING_GUOHALL);
  }

  mask.convertTo(mask, CV_32FC1, 1.0/255.0);
  cv::GaussianBlur(mask, mask, cv::Size(5, 5), 0.75);

  for (int y = 0; y < image.rows; ++y)
  {
    cv::Vec3b* ptr_image = reinterpret_cast<cv::Vec3b*>(image.ptr(y));
    const float* ptr_mask = reinterpret_cast<const float*>(mask.ptr(y));
    for (int x = 0; x < image.cols; ++x)
    {
      ptr_image[x] = ptr_mask[x] * dark_brown + (1.0f - ptr_mask[x]) * ptr_image[x];
    }
  }

  cv::imwrite(path_out.string(), image);

  return 0;
}