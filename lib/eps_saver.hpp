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

#ifndef TRLIB_EPS_SAVER_HPP_
#define TRLIB_EPS_SAVER_HPP_

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <opencv2/opencv.hpp>

#include "affine_transformation.hpp"
#include "bezier_curve.hpp"
#include "hershey_font.hpp"
#include "texture.hpp"
#include "vector_graphics_saver.hpp"

class EPSLine : public SaverPrimitive
{
public:
	EPSLine(cv::Point2f p_start, cv::Point2f p_end, double pix_to_mm_scale) :
		p_start(p_start),
		p_end(p_end),
		pix_to_mm_scale(pix_to_mm_scale)
	{}

	virtual void serialize(std::ostream& os) const override
	{
		os << "setlinecolor " << p_start.x * pix_to_mm_scale << " " << -p_start.y * pix_to_mm_scale
      << " moveto " << p_end.x * pix_to_mm_scale << " " << -p_end.y * pix_to_mm_scale
      << " lineto stroke";
	}

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const cv::Point2f p_start_transformed = AffineTransformation::transform(T, p_start);
    const cv::Point2f p_end_transformed = AffineTransformation::transform(T, p_end);
    return std::make_shared<EPSLine>(p_start_transformed, p_end_transformed, pix_to_mm_scale);
  }

	cv::Point2f p_start;
	cv::Point2f p_end;
	double pix_to_mm_scale;
};

class EPSPolygon : public SaverPrimitive
{
public:
	EPSPolygon(const std::vector<cv::Point2f>& points, double pix_to_mm_scale) :
		points(points),
		pix_to_mm_scale(pix_to_mm_scale)
	{}

  virtual void serialize(std::ostream& os) const override
	{
		if (points.empty())
		{
			return;
		}

		os << "setpolygoncolor " << points[0].x * pix_to_mm_scale << " " << -points[0].y * pix_to_mm_scale << " moveto ";
		for (size_t i = 1; i < points.size(); ++i)
		{
			os << points[i].x * pix_to_mm_scale << " " << -points[i].y * pix_to_mm_scale << " lineto ";
		}
		os << "closepath stroke";
	}

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const std::vector<cv::Point2f> points_transformed = AffineTransformation::transform(T, points);
    return std::make_shared<EPSPolygon>(points_transformed, pix_to_mm_scale);
  }

	std::vector<cv::Point2f> points;
	double pix_to_mm_scale;
};

class EPSText : public SaverPrimitive
{
public:
	EPSText(const cv::Point2f& pos, const std::string& text, double rotation_rad, double pix_to_mm_scale) :
    pos(pos),
    text(text),
    rotation_deg(180.0*rotation_rad/boost::math::double_constants::pi),
    pix_to_mm_scale(pix_to_mm_scale)
	{}

	virtual void serialize(std::ostream& os) const override
	{
    HersheyFont font(pos * pix_to_mm_scale, text, rotation_deg, 2.5, 0.8);
    os << "settextcolor ";
    font.serialize_ps(os);
	}

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const cv::Point2f pos_transformed = AffineTransformation::transform(T, pos);
    const double rotation_deg_transformed = rotation_deg + AffineTransformation::rotation_deg(T);
    return std::make_shared<EPSText>(pos_transformed, text, rotation_deg_transformed / 180.0 * boost::math::double_constants::pi, pix_to_mm_scale);
  }

  cv::Point2f pos;
  std::string text;
  double rotation_deg;
  double pix_to_mm_scale;
};

class EPSBezier : public SaverPrimitive
{
public:
	EPSBezier(const std::vector<BezierCurve>& curve, double pix_to_mm_scale) :
		curve(curve),
		pix_to_mm_scale(pix_to_mm_scale)
	{}

	virtual void serialize(std::ostream& os) const override
	{
		if (curve.empty())
		{
			return;
		}

    bool moveto = true;

		os << "\nsetbeziercolor % " << curve.size() << "\n";
		for (int i = 0; i < curve.size(); ++i)
		{
      const BezierCurve& c = curve.at(i);
      if (moveto)
      {
        os << c.control_point(0).x * pix_to_mm_scale << " " << -c.control_point(0).y * pix_to_mm_scale << " moveto ";
        moveto = false;
      }
			for (int i = 1; i <= c.degree(); ++i)
			{
				os << c.control_point(i).x * pix_to_mm_scale << " " << -c.control_point(i).y * pix_to_mm_scale << " ";
			}
			os << "curveto ";
		}
		os << "closepath stroke ";
	}

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const std::vector<BezierCurve> curve_transformed = BezierCurve::transformed(curve, T);
    return std::make_shared<EPSBezier>(curve_transformed, pix_to_mm_scale);
  }

private:
	std::vector<BezierCurve> curve;
	double pix_to_mm_scale;
};

class EPSCircle : public SaverPrimitive
{
public:
	EPSCircle(cv::Point2f center, double radius, double pix_to_mm_scale) :
		center(center),
		radius(radius),
		pix_to_mm_scale(pix_to_mm_scale)
	{
	}

	virtual void serialize(std::ostream& os) const override
	{
		os << "setcirclecolor " << center.x * pix_to_mm_scale << " " << -center.y * pix_to_mm_scale << " " << radius * pix_to_mm_scale << " 0 360 arc stroke ";
	}

  virtual std::shared_ptr<SaverPrimitive> transformed(cv::Mat T) const override
  {
    const cv::Point2f center_transformed = AffineTransformation::transform(T, center);
    return std::make_shared<EPSCircle>(center, radius, pix_to_mm_scale);
  }

private:
  cv::Point2f center;
	double radius;
	double pix_to_mm_scale;
};

class EPSSaver : public VectorGraphicsSaver
{
public:
	EPSSaver(const Texture& texture, const std::string &output_name = "") :
    VectorGraphicsSaver(texture, output_name),
		pixel_to_mm_scale(25.4 / texture.dpi)
	{
    for (const cv::Point2d& p : texture.marker.markers_pix)
    {
      add_circle(p, 5.0 / 25.4 * texture.dpi);
    }

    add_polygon({
      cv::Point2f(0.0f, 0.0f),
      cv::Point2f(static_cast<float>(texture.texture.cols), 0.0f),
      cv::Point2f(static_cast<float>(texture.texture.cols), static_cast<float>(texture.texture.rows)),
      cv::Point2f(0.0f, static_cast<float>(texture.texture.rows))});
  }

  EPSSaver(const std::string &output_name) :
     VectorGraphicsSaver(output_name),
		 pixel_to_mm_scale(-1),
     output_name(output_name)
  {}


	virtual void add_line(cv::Point2f p_start, cv::Point2f p_end) override
	{
		m_primitives.push_back(std::make_shared<EPSLine>(p_start, p_end, pixel_to_mm_scale));
	}

  virtual void add_polygon(const std::vector<cv::Point2f>& points) override
	{
		m_primitives.push_back(std::make_shared<EPSPolygon>(points, pixel_to_mm_scale));
	}

  virtual void add_text(const cv::Point2f& pos, const std::string& text, double rotation_rad) override
	{
		m_text.push_back(std::make_shared<EPSText>(pos, text, rotation_rad, pixel_to_mm_scale));
	}

  virtual void add_bezier(const std::vector<BezierCurve>& curve) override
	{
		m_primitives.push_back(std::make_shared<EPSBezier>(curve, pixel_to_mm_scale));
	}

  virtual void add_circle(const cv::Point2f& center, double radius) override
	{
		m_primitives.push_back(std::make_shared<EPSCircle>(center, radius, pixel_to_mm_scale));
	}

  virtual std::string extension() const override
  {
    return ".eps";
  }

  virtual void save(const boost::filesystem::path& filename, const cv::Size2d& table_dimensions_mm) const override;

private:
	double pixel_to_mm_scale;
  std::string output_name;
};

#endif /* TRLIB_EPS_SAVER_HPP_ */