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

#ifndef TRLIB_MAT_HPP_
#define TRLIB_MAT_HPP_

#include <vector>

template <typename T>
class mat
{
public:
  mat() : m_height(0),
          m_width(0),
          m_rotate_x(0),
          m_rotate_y(0)
  {
  }

  mat(int height, int width, const T &init = T()) : m_rotate_x(0),
                                                    m_rotate_y(0)
  {
    resize(height, width, init);
  }

  void resize(int height, int width, const T &init = T())
  {
    m_height = height;
    m_width = width;
    m_mat.resize(height * width, init);
  }

  T &operator[](int i)
  {
    const std::pair<int, int> coordinate = index_to_coordinate(i);
    const int index_x = rotate_coordinate_x(coordinate.first, m_rotate_x);
    const int index_y = rotate_coordinate_y(coordinate.second, m_rotate_y);
    i = coordinate_to_index(index_x, index_y);

    return m_mat[i];
  }

  const T &operator[](int i) const
  {
    const std::pair<int, int> coordinate = index_to_coordinate(i);
    const int index_x = rotate_coordinate_x(coordinate.first, m_rotate_x);
    const int index_y = rotate_coordinate_y(coordinate.second, m_rotate_y);
    i = coordinate_to_index(index_x, index_y);

    return m_mat[i];
  }

  T &get(int i, int rotate_rows, int rotate_cols)
  {
    const std::pair<int, int> coordinate = index_to_coordinate(i);
    const int index_x = rotate_coordinate_x(coordinate.first, rotate_cols);
    const int index_y = rotate_coordinate_y(coordinate.second, rotate_rows);
    i = coordinate_to_index(index_x, index_y);

    return m_mat[i];
  }

  const T &get(int i, int rotate_rows, int rotate_cols) const
  {
    const std::pair<int, int> coordinate = index_to_coordinate(i);
    const int index_x = rotate_coordinate_x(coordinate.first, rotate_cols);
    const int index_y = rotate_coordinate_y(coordinate.second, rotate_rows);
    i = coordinate_to_index(index_x, index_y);

    return m_mat[i];
  }

  T &operator()(int row, int col)
  {
    row = rotate_coordinate_y(row, m_rotate_y);
    col = rotate_coordinate_x(col, m_rotate_x);

    return m_mat[row * m_width + col];
  }

  const T &operator()(int row, int col) const
  {
    row = rotate_coordinate_y(row, m_rotate_y);
    col = rotate_coordinate_x(col, m_rotate_x);

    return m_mat[row * m_width + col];
  }

  T &operator()(int row, int col, int rotate_rows, int rotate_cols)
  {
    row = rotate_coordinate_y(row, rotate_rows);
    col = rotate_coordinate_x(col, rotate_cols);

    return m_mat[row * m_width + col];
  }

  const T &operator()(int row, int col, int rotate_rows, int rotate_cols) const
  {
    row = rotate_coordinate_y(row, rotate_rows);
    col = rotate_coordinate_x(col, rotate_cols);

    return m_mat[row * m_width + col];
  }

  mat<T> roi(cv::Rect region) const
  {
    mat<T> mat_roi(m_height, m_width);
    mat_roi.m_rotate_x = m_rotate_x;
    mat_roi.m_rotate_y = m_rotate_y;

    for (int i = 0; i < m_height * m_width; ++i)
    {
      mat_roi.m_mat[i] = m_mat[i](region);
    }

    return mat_roi;
  }

  void rotate_x(int rotate_diff)
  {
    set_rotate_x(m_rotate_x + rotate_diff);
  }

  void rotate_y(int rotate_diff)
  {
    set_rotate_y(m_rotate_y += rotate_diff);
  }

  void set_rotate_x(int rotate_x)
  {
    m_rotate_x = rotate_x;
    while (m_rotate_x < 0)
      m_rotate_x += m_width;
    while (m_rotate_x >= m_width)
      m_rotate_x -= m_width;
  }

  void set_rotate_y(int rotate_y)
  {
    m_rotate_y = rotate_y;
    while (m_rotate_y < 0)
      m_rotate_y += m_height;
    while (m_rotate_y >= m_height)
      m_rotate_y -= m_height;
  }

  int height() const
  {
    return m_height;
  }

  int width() const
  {
    return m_width;
  }

  int size() const
  {
    return m_height * m_width;
  }

  bool empty() const
  {
    return size() == 0;
  }

private:
  std::pair<int, int> index_to_coordinate(int index) const
  {
    const int index_x = index % m_width;
    const int index_y = index / m_width;

    return std::make_pair(index_x, index_y);
  }

  int coordinate_to_index(int index_x, int index_y) const
  {
    return index_y * m_width + index_x;
  }

  int rotate_coordinate_x(int index_x, int rotate_x) const
  {
    index_x += rotate_x;
    while (index_x < 0)
      index_x += m_width;
    while (index_x >= m_width)
      index_x -= m_width;
    return index_x;
  }

  int rotate_coordinate_y(int index_y, int rotate_y) const
  {
    index_y += rotate_y;
    while (index_y < 0)
      index_y += m_height;
    while (index_y >= m_height)
      index_y -= m_height;
    return index_y;
  }

  std::vector<T> m_mat;

  int m_height;
  int m_width;

  int m_rotate_x;
  int m_rotate_y;
};

#endif /* TRLIB_MAT_HPP_ */