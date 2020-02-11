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

#include <list>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/tokenizer.hpp>

#include "bezier_curve.hpp"
#include "bezier_ransac.hpp"
#include "grid.hpp"
#include "patch_region.hpp"

namespace fs = boost::filesystem;

bool Grid::load_morphed_grid(const std::string& morphed_grid_filename)
{
  typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;
  std::vector<std::vector<std::string>> separated_vec;

  if (fs::exists(morphed_grid_filename) && fs::is_regular_file(morphed_grid_filename))
  {
    std::string line;

    std::ifstream str(morphed_grid_filename);
    while (std::getline(str, line))
    {
      Tokenizer tokenizer(line);
      separated_vec.emplace_back(tokenizer.begin(), tokenizer.end());
    }
  }

  if (separated_vec.empty())
  {
    return false;
  }

  m_rows = static_cast<int>(separated_vec.size());
  m_cols = static_cast<int>(separated_vec[0].size() / 2);
  m_grid_points.resize(m_rows*m_cols);

  int c = 0;
  for (int y = 0; y < m_rows; ++y)
  {
    for (int x = 0; x < m_cols; ++x)
    {
      m_grid_points[c].x = std::stod(separated_vec[y][x]);
      m_grid_points[c].y = std::stod(separated_vec[y][x+m_cols]);
      ++c;
    }
  }

  return true;
}

void Grid::resize(int rows, int cols)
{
  m_rows = rows;
  m_cols = cols;
  m_grid_points.resize(rows*cols);

  int c = 0;
  for (int y = 0; y < rows; ++y)
  {
    for (int x = 0; x < cols; ++x)
    {
      m_grid_points[c].x = x;
      m_grid_points[c].y = y;
      ++c;
    }
  }
}

void Grid::scale(double factor)
{
  for (cv::Point2d& p : m_grid_points)
  {
    p *= factor;
  }
}

boost::property_tree::ptree Grid::save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const
{
  boost::property_tree::ptree tree;
  serialize(tree, "grid_points", m_grid_points, base_path, path);
  serialize(tree, "rows", m_rows, base_path, path);
  serialize(tree, "cols", m_cols, base_path, path);
  return tree;
}

void Grid::load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree)
{
  deserialize(tree, "grid_points", m_grid_points, base_path);
  deserialize(tree, "rows", m_rows, base_path);
  deserialize(tree, "cols", m_cols, base_path);
}
