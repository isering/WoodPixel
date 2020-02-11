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

#ifndef TRLIB_DATA_GRID_HPP_
#define TRLIB_DATA_GRID_HPP_

#include "grid.hpp"
#include "vector.hpp"

template <class T>
class DataGrid : public Grid
{
public:
  DataGrid<T>() :
    Grid()
  {
  }

  DataGrid<T>(int rows, int cols) :
    Grid(rows, cols)
  {
    m_grid_data.resize(rows*cols);
  }

  DataGrid<T>(const Grid& grid) :
    Grid(grid)
  {
    m_grid_data.resize(rows()*cols());
  }

  DataGrid<T>(int rows, int cols, const T& val) :
    Grid(rows, cols)
  {
    m_grid_data.resize(rows*cols, val);
  }

  DataGrid<T>(const Grid& grid, const T& val) :
    Grid(grid)
  {
    m_grid_data.resize(rows()*cols(), val);
  }

  T& data(int row, int col)
  {
    return m_grid_data[row*cols()+col];
  }

  const T& data(int row, int col) const
  {
    return m_grid_data[row*cols()+col];
  }

  T& data(const cv::Point& p)
  {
    return m_grid_data[p.y*cols()+p.x];
  }

  const T& data(const cv::Point& p) const
  {
    return m_grid_data[p.y*cols()+p.x];
  }

  virtual void resize(int rows, int cols)
  {
    Grid::resize(rows, cols);
    m_grid_data.resize(rows*cols);
  }

  virtual boost::property_tree::ptree save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const override
  {
    boost::property_tree::ptree tree;
    tree.add_child("grid", Grid::save(base_path, path));
    serialize(tree, "grid_data", m_grid_data, base_path, path);
    return tree;
  }

  virtual void load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree) override
  {
    Grid::load(base_path, tree.get_child("grid"));
    deserialize(tree, "grid_data", m_grid_data, base_path);
  }

private:
  Vector<T> m_grid_data;
};

#endif /* TRLIB_DATA_GRID_HPP_ */