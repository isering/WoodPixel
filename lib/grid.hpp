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

#ifndef TRLIB_GRID_HPP_
#define TRLIB_GRID_HPP_

#include <opencv2/opencv.hpp>

#include "patch.hpp"
#include "patch_region.hpp"
#include "serializable.hpp"

template <typename T>
class DataGrid;

class Grid : public Serializable
{
public:
  Grid() : m_rows(0),
           m_cols(0)
  {
  }

  Grid(int rows, int cols)
  {
    resize(rows, cols);
  }

  Grid(const std::string &morphed_grid_filename)
  {
    load_morphed_grid(morphed_grid_filename);
  }

  virtual void resize(int rows, int cols);

  cv::Point2d &operator()(int row, int col)
  {
    return m_grid_points[row * m_cols + col];
  }

  const cv::Point2d &operator()(int row, int col) const
  {
    return m_grid_points[row * m_cols + col];
  }

  cv::Point2d &operator()(const cv::Point &p)
  {
    return m_grid_points[p.y * m_cols + p.x];
  }

  const cv::Point2d &operator()(const cv::Point &p) const
  {
    return m_grid_points[p.y * m_cols + p.x];
  }

  DataGrid<PatchRegion> generate_patches(const cv::Mat &edge_image);

  bool empty() const
  {
    return (m_rows == 0) || (m_cols == 0);
  }

  int rows() const
  {
    return m_rows;
  }

  int cols() const
  {
    return m_cols;
  }

  void scale(double factor);

  bool load_morphed_grid(const std::string &morphed_grid_filename);

  template <typename T>
  cv::Mat_<cv::Vec<T, 2>> to_mat() const
  {
    cv::Mat_<cv::Vec<T, 2>> mat(m_rows, m_cols);

    int c = 0;
    for (int y = 0; y < m_rows; ++y)
    {
      cv::Vec<T, 2> *ptr = reinterpret_cast<cv::Vec<T, 2> *>(mat.ptr(y));
      for (int x = 0; x < m_cols; ++x)
      {
        ptr[x][0] = static_cast<T>(m_grid_points[c].x);
        ptr[x][1] = static_cast<T>(m_grid_points[c].y);
        ++c;
      }
    }
    return mat;
  }

  template <typename T>
  static Grid from_mat(cv::Mat_<cv::Vec<T, 2>> mat)
  {
    Grid grid(mat.rows, mat.cols);

    int c = 0;
    for (int y = 0; y < mat.rows; ++y)
    {
      const cv::Vec<T, 2> *ptr = reinterpret_cast<const cv::Vec<T, 2> *>(mat.ptr(y));
      for (int x = 0; x < mat.cols; ++x)
      {
        grid.m_grid_points[c++] = cv::Point2d(ptr[x]);
      }
    }

    return grid;
  }

  virtual void load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const override;

private:
  std::vector<cv::Point2d> m_grid_points;
  int m_rows, m_cols;
};

#endif /* TRLIB_GRID_HPP_ */