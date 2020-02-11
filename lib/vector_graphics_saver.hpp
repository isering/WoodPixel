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

#ifndef TRLIB_VECTOR_GRAPHICS_SAVER_HPP_
#define TRLIB_VECTOR_GRAPHICS_SAVER_HPP_

#include <ostream>
#include <vector>

#include <boost/filesystem.hpp>

#include "bezier_curve.hpp"
#include "texture.hpp"

class SaverPrimitive
{
public:
  virtual void serialize(std::ostream& os) const = 0;
  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const = 0;
};

std::ostream& operator<<(std::ostream& os, const SaverPrimitive& primitive);

class VectorGraphicsSaver
{
public:
  VectorGraphicsSaver(const Texture& texture, const std::string &output_name = "") :
    rows(static_cast<int>(texture.texture.rows)),
    cols(static_cast<int>(texture.texture.cols)),
    dpi(texture.dpi),
    marker(texture.marker),
    output_name(output_name)
  {}
  VectorGraphicsSaver(const std::string &output_name) :
    rows(-1),
    cols(-1),
    dpi(-1),
    output_name(output_name)
  {}

  virtual void add_line(cv::Point2f p_start, cv::Point2f p_end) = 0;
  virtual void add_polygon(const std::vector<cv::Point2f>& points) = 0;
  virtual void add_text(const cv::Point2f& pos, const std::string& text, double rotation_rad) = 0;
  virtual void add_bezier(const std::vector<BezierCurve>& curve) = 0;
  virtual void add_circle(const cv::Point2f& center, double radius) = 0;

  virtual void add_debug_circle(const cv::Point2f& center, double radius)
  {
    add_circle(center, radius);
  }

  virtual void save(const boost::filesystem::path& filename, const cv::Size2d& table_dimensions_mm) const = 0;

  virtual std::string extension() const = 0;
  
  std::string get_output_name() const 
  {
    return output_name;
  }
  const std::vector<std::shared_ptr<SaverPrimitive>> get_primitives() const
  {
    return m_primitives;
  }
  const std::vector<std::shared_ptr<SaverPrimitive>> get_text() const
  {
    return m_text;
  }

protected:
  std::vector<std::shared_ptr<SaverPrimitive>> m_primitives;
  std::vector<std::shared_ptr<SaverPrimitive>> m_debug_primitives;
  std::vector<std::shared_ptr<SaverPrimitive>> m_text;

  int rows, cols;
  double dpi;
  std::string output_name;
  TextureMarker marker;
};

#endif /* TRLIB_VECTOR_GRAPHICS_SAVER_HPP_ */