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

#ifndef EZ_GRID_HPP_
#define EZ_GRID_HPP_

#include <boost/filesystem.hpp>
#include <opencv2/opencv.hpp>

#include "bilateral_filter_future.hpp"
#include "canny_future.hpp"
#include "morph_grid_future.hpp"
#include "serializable.hpp"

struct EZGridData
{
  Grid morphed_grid;
  std::vector<PatchRegion> patches;
};

class EZGrid : public Serializable
{
public:
  EZGrid(cv::Mat image, cv::Mat image_filtered, int grid_size, const boost::filesystem::path& path_out);

  void run();
  void load_partial_state(const boost::filesystem::path& path);
  void load_full_state(const boost::filesystem::path& path);
  void save();

  virtual void load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const override;

private:
  std::unique_ptr<BilateralFilterFuture> m_bilateral_filter_future;
  std::unique_ptr<CannyFuture> m_canny_bilateral_future, m_canny_filtered_future;
  std::unique_ptr<MorphGridFuture> m_morph_grid_future;

  std::string m_window_name;
  std::string m_gui_name;
  boost::filesystem::path m_path_out;
};

#endif /* EZ_GRID_HPP_ */
