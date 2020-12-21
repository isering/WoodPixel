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

#ifndef TRLIB_PATCH_GENERATOR_HPP_
#define TRLIB_PATCH_GENERATOR_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

#include "bezier_curve.hpp"
#include "grid.hpp"
#include "patch_region.hpp"
#include "serializable.hpp"
#include "vector.hpp"

std::vector<PatchRegion> generate_patches(int target_index, const Grid &morphed_grid, cv::Mat edge_image);

std::vector<PatchRegion> generate_patches_square(int target_index, const cv::Size &image_size, const cv::Size &patch_size, const cv::Size &filter_kernel_size);

/*
 * Low level stuff.
 */
class CurveData : public Serializable
{
public:
  BezierCurve curve_horiz;
  BezierCurve curve_vert;
  BezierCurve curve_diag;

  virtual void load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const override;
};

DataGrid<Vector<cv::Point2d>> extract_edge_grid(const Grid &grid, cv::Mat edge_mat);
void fit_single_edge(DataGrid<CurveData> &curve_grid, DataGrid<Vector<cv::Point2d>> &edge_grid, cv::Point p_1, cv::Point p_2);
void fill_smooth_curves(DataGrid<CurveData> &curve_grid);
std::vector<PatchRegion> generate_patches_from_curve_grid(int target_index, const DataGrid<CurveData> &curve_grid);

#endif /* TRLIB_PATCH_GENERATOR_HPP_ */