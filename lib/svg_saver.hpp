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

#ifndef TRLIB_SVG_SAVER_HPP_
#define TRLIB_SVG_SAVER_HPP_

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <opencv2/opencv.hpp>

#include "affine_transformation.hpp"
#include "bezier_curve.hpp"
#include "hershey_font.hpp"
#include "texture_marker.hpp"
#include "vector_graphics_saver.hpp"

const double svg_precision = 10000.0;

class SVGLine : public SaverPrimitive
{
public:
  SVGLine(cv::Point2f p_start, cv::Point2f p_end, const std::string& stroke = "str0", const std::string& fill = "fil0") :
    p_start(p_start),
    p_end(p_end),
	stroke(stroke),
	fill(fill)
  {}

  virtual void serialize(std::ostream& os) const override
  {
    const int x_start = static_cast<int>(svg_precision * p_start.x);
    const int y_start = static_cast<int>(svg_precision * p_start.y);
    const int x_end = static_cast<int>(svg_precision * p_end.x);
    const int y_end = static_cast<int>(svg_precision * p_end.y);
    os << "<line class=\"" << fill << " " << stroke << "\" x1=\"" << x_start << "\" y1=\"" << y_start << "\" x2=\"" << x_end << "\" y2=\"" << y_end << "\" />";
  }

  std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const cv::Point2f p_start_transformed = AffineTransformation::transform(T, p_start);
    const cv::Point2f p_end_transformed = AffineTransformation::transform(T, p_end);
    return std::make_shared<SVGLine>(p_start_transformed, p_end_transformed, stroke, fill);
  }

  cv::Point2f p_start;
  cv::Point2f p_end;
  std::string stroke;
  std::string fill;
};

class SVGPolygon : public SaverPrimitive
{
public:
	SVGPolygon(const std::vector<cv::Point2f>& points, const std::string& stroke = "str0", const std::string& fill = "fil0") :
		points(points),
		stroke(stroke),
		fill(fill)
	{}

  virtual void serialize(std::ostream& os) const override
  {
    os << "<polyline class=\"" << fill << " " << stroke << "\" points=\"";
    for (const cv::Point2f& p : points)
    {
      const int x = static_cast<int>(svg_precision * p.x);
      const int y = static_cast<int>(svg_precision * p.y);
      os << x << "," << y << " ";
    }
    const int x = static_cast<int>(svg_precision * points.front().x);
    const int y = static_cast<int>(svg_precision * points.front().y);
    os << x << "," << y << "\" />";
  }

  std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const std::vector<cv::Point2f> points_transformed = AffineTransformation::transform(T, points);
    return std::make_shared<SVGPolygon>(points_transformed, stroke, fill);
  }
  
  std::vector<cv::Point2f> points;
  std::string stroke;
  std::string fill;
};

class SVGText : public SaverPrimitive
{
public:
  SVGText(const cv::Point2f& pos, const std::string& text, double rotation_rad) :
    pos(pos),
    text(text),
    rotation_deg(180.0 * rotation_rad / boost::math::double_constants::pi)
  {}

  virtual void serialize(std::ostream& os) const override
  {
    HersheyFont font(pos, text, rotation_deg, 4.5, 0.8);
    font.serialize_svg(os, svg_precision);
  }

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const cv::Point2f pos_transformed = AffineTransformation::transform(T, pos);
    const double rotation_deg_transformed = rotation_deg + AffineTransformation::rotation_deg(T);
    return std::make_shared<SVGText>(pos_transformed, text, rotation_deg_transformed / 180.0 * boost::math::double_constants::pi);
  }

  cv::Point2f pos;
  std::string text;
  double rotation_deg;
};

class SVGBezier : public SaverPrimitive
{
public:
  SVGBezier(const std::vector<BezierCurve>& curve, const std::string& stroke = "str0", const std::string& fill = "fil0") :
    curve(curve),
    stroke(stroke),
    fill(fill)
  {}

  virtual void serialize(std::ostream& os) const override
  {
    const int x_front = static_cast<int>(svg_precision * curve[0].control_point(0).x);
    const int y_front = static_cast<int>(svg_precision * curve[0].control_point(0).y);
    os << "<path class=\"" << fill << " " << stroke << "\" d=\"M " << x_front << " " << y_front;
    for (const BezierCurve& c : curve)
    {
      os << " C ";
      for (int i = 1; i < c.degree(); ++i)
      {
        const int x = static_cast<int>(svg_precision * c.control_point(i).x);
        const int y = static_cast<int>(svg_precision * c.control_point(i).y);
        os << x << " " << y << ", ";
      }
      const int x_end = static_cast<int>(svg_precision * c.control_point(c.degree()).x);
      const int y_end = static_cast<int>(svg_precision * c.control_point(c.degree()).y);
      os << x_end << " " << y_end;
    }
    os << "\" />";
  }

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const std::vector<BezierCurve> curve_transformed = BezierCurve::transformed(curve, T);
    return std::make_shared<SVGBezier>(curve_transformed, stroke, fill);
  }

  std::vector<BezierCurve> curve;
  std::string stroke;
  std::string fill;
};

class SVGCircle : public SaverPrimitive
{
public:
  SVGCircle(cv::Point2f center, double radius, const std::string& stroke = "str0", const std::string& fill = "fil0") :
    center(center),
    radius(radius),
    stroke(stroke),
    fill(fill)
  {}

  virtual void serialize(std::ostream& os) const override
  {
    const int x = static_cast<int>(svg_precision * center.x);
    const int y = static_cast<int>(svg_precision * center.y);
    const int r = static_cast<int>(svg_precision * radius);
    os << "<circle class=\"" << fill << " " << stroke << "\" cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << r << "\" />";
  }

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    cv::Point2f center_transformed = AffineTransformation::transform(T, center);
    return std::make_shared<SVGCircle>(center_transformed, radius, stroke, fill);
  }

  cv::Point2f center;
  double radius;
  std::string stroke;
  std::string fill;
};

class SVGSaver : public VectorGraphicsSaver
{
public:
  SVGSaver(const Texture& texture, const std::string &output_name = "") :
    VectorGraphicsSaver(texture, output_name)
  {
  }

  SVGSaver(const std::string &output_name) :
    VectorGraphicsSaver(output_name)
  {
  }

  virtual void add_line(cv::Point2f p_start, cv::Point2f p_end) override
  {
    m_primitives.push_back(std::make_shared<SVGLine>(p_start, p_end));
  }

  virtual void add_polygon(const std::vector<cv::Point2f>& points) override
  {
    m_primitives.push_back(std::make_shared<SVGPolygon>(points, "str0", "fil0"));
  }

  virtual void add_text(const cv::Point2f& pos, const std::string& text, double rotation_rad) override
  {
    m_text.push_back(std::make_shared<SVGText>(pos, text, rotation_rad));
  }

  virtual void add_bezier(const std::vector<BezierCurve>& curve) override
  {
    m_primitives.push_back(std::make_shared<SVGBezier>(curve));
  }

  virtual void add_circle(const cv::Point2f& center, double radius) override
  {
    m_primitives.push_back(std::make_shared<SVGCircle>(center, radius, "str1", "fil0"));
  }

  virtual void add_debug_circle(const cv::Point2f& center, double radius) override
  {
    m_debug_primitives.push_back(std::make_shared<SVGCircle>(center, radius, "str1", "fil0"));
  }

  virtual void add_debug_polygon(const std::vector<cv::Point2f>& points)
  {
    m_primitives.push_back(std::make_shared<SVGPolygon>(points, "str1", "fil0"));
  }

  virtual std::string extension() const
  {
    return ".svg";
  }

  virtual void save(const boost::filesystem::path& filename, const cv::Size2d& table_dimensions_mm) const override;

private:
  void write_header(std::ostream& out, double size_x_mm, double size_y_mm) const;
  void write_footer(std::ostream& out) const;
};
#endif /* TRLIB_SVG_SAVER_HPP_ */