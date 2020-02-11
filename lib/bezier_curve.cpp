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

#include <random>
#include <stack>

#include "sort_pca.hpp"
#include "transformations.hpp"

#include "bezier_curve.hpp"

cv::Point2d BezierCurve::eval(double t) const
{
  std::vector<cv::Point2d> curve = m_control_points;

  for (int i = 1; i <= m_degree; ++i)
  {
    for (int j = 0; j <= m_degree-i; ++j)
    {
      curve[j] = (1.0-t) * curve[j] + t * curve[j+1];
    }
  }

  return curve[0];
}

double BezierCurve::eval_x(double t) const
{
  std::vector<double> curve(m_degree+1);
  for (int i = 0; i <= m_degree; ++i)
  {
    curve[i] = m_control_points[i].x;
  }

  for (int i = 1; i <= m_degree; ++i)
  {
    for (int j = 0; j <= m_degree-i; ++j)
    {
      curve[j] = (1.0-t) * curve[j] + t * curve[j+1];
    }
  }

  return curve[0];
}

double BezierCurve::eval_y(double t) const
{
  std::vector<double> curve(m_degree+1);
  for (int i = 0; i <= m_degree; ++i)
  {
    curve[i] = m_control_points[i].y;
  }

  for (int i = 1; i <= m_degree; ++i)
  {
    for (int j = 0; j <= m_degree-i; ++j)
    {
      curve[j] = (1.0-t) * curve[j] + t * curve[j+1];
    }
  }

  return curve[0];
}

std::vector<BezierCurve> BezierCurve::translated(const std::vector<BezierCurve>& curves, const cv::Point2d& delta)
{
  std::vector<BezierCurve> curves_out(curves.size());
  std::transform(curves.begin(), curves.end(), curves_out.begin(), [delta](const BezierCurve& c) {return c + delta; });
  return curves_out;
}

BezierCurve BezierCurve::transformed(cv::Mat T) const
{
  BezierCurve curve(m_degree);
  for (int i = 0; i <= m_degree; ++i)
  {
    curve.m_control_points[i].x = T.at<double>(0, 0) * m_control_points[i].x + T.at<double>(0, 1) * m_control_points[i].y + T.at<double>(0, 2);
    curve.m_control_points[i].y = T.at<double>(1, 0) * m_control_points[i].x + T.at<double>(1, 1) * m_control_points[i].y + T.at<double>(1, 2);
  }
  return curve;
}

std::vector<BezierCurve> BezierCurve::transformed(const std::vector<BezierCurve>& curves, cv::Mat transformation_matrix)
{
  std::vector<BezierCurve> curves_transformed(curves.size());
  std::transform(curves.begin(), curves.end(), curves_transformed.begin(), [transformation_matrix](const BezierCurve& c) {return c.transformed(transformation_matrix); });
  return curves_transformed;
}

BezierCurve BezierCurve::reversed() const
{
  std::vector<cv::Point2d> control_points_reversed(m_degree+1);
  for (int i = 0; i <= m_degree; ++i)
  {
    control_points_reversed[m_degree-i] = m_control_points[i];
  }
  return BezierCurve(control_points_reversed);
}

BezierCurve BezierCurve::scaled(double factor) const
{
  return (*this) * factor;
}

std::pair<BezierCurve, BezierCurve> BezierCurve::split_cubic(double t) const
{
  if (m_degree != 3)
  {
    throw(std::invalid_argument("BezierCurve::split_cubic called on a non-cubic curve."));
  }

  std::vector<cv::Point2d> d1;
  for (int i = 0; i < 3; ++i)
  {
    d1.push_back(t * m_control_points[i+1] + (1.0-t) * m_control_points[i]);
  }

  std::vector<cv::Point2d> d2;
  for (int i = 0; i < 2; ++i)
  {
    d2.push_back(t * d1[i+1] + (1.0-t) * d1[i]);
  }

  cv::Point2d d3 = t * d2[1] + (1.0-t) * d2[0];

  BezierCurve b1(std::vector<cv::Point2d>({m_control_points[0], d1[0], d2[0], d3}));
  BezierCurve b2(std::vector<cv::Point2d>({d3, d2[1], d1[2], m_control_points[3]}));

  return std::make_pair(b1, b2);
}

bool BezierCurve::is_empty_curve() const
{
  double dist = 0.0;
  for (int i = 1; i <= m_degree; ++i)
  {
    const cv::Point2d diff = m_control_points[0] - m_control_points[i];
    dist += diff.x * diff.x + diff.y * diff.y;
  }
  return dist < 1.0e-2;
}

bool BezierCurve::is_collinear() const
{
  if (m_degree < 2)
  {
    return true;
  }

  for (int i = 2; i <= m_degree; ++i)
  {
    const double triangle_area =
      m_control_points[0].x * (m_control_points[1].y - m_control_points[i].y) +
      m_control_points[1].x * (m_control_points[i].y - m_control_points[0].y) +
      m_control_points[i].x * (m_control_points[0].y - m_control_points[1].y);

    if (std::abs(triangle_area) > 1.0e-6)
    {
      return false;
    }
  }

  return true;
}

BezierCurve BezierCurve::canonical_form_cubic() const
{
  if (m_degree != 3)
  {
    throw(std::invalid_argument("BezierCurve::canonical_form_cubic called for non-cubic curve."));
  }

  if (is_collinear())
  {
    return *this;
  }

  BezierCurve curve_canonical = *this - m_control_points[0];
  curve_canonical = curve_canonical.transformed(transformations::shear_x(-curve_canonical.control_point(1).x / curve_canonical.control_point(1).y));
  curve_canonical = curve_canonical.transformed(transformations::scale(1.0 / curve_canonical.control_point(2).x, 1.0 / curve_canonical.control_point(1).y));
  curve_canonical = curve_canonical.transformed(transformations::shear_y(1.0 - curve_canonical.control_point(2).y / curve_canonical.control_point(2).x));
  
  return curve_canonical;
}

cv::Point2d BezierCurve::canonical_form_cubic_free() const
{
  if (m_degree != 3)
  {
    throw(std::invalid_argument("BezierCurve::canonical_form_cubic_free called for non-cubic curve."));
  }

  if (is_collinear())
  {
    throw(std::invalid_argument("BezierCurve::canonical_form_cubic_free called with collinear curve"));
  }

  const cv::Point2d p2 = m_control_points[1] - m_control_points[0];
  const cv::Point2d p3 = m_control_points[2] - m_control_points[0];
  const cv::Point2d p4 = m_control_points[3] - m_control_points[0];

  const double f32 = p3.y / p2.y;
  const double f42 = p4.y / p2.y;

  const double x = (p4.x - p2.x * f42) / (p3.x - p2.x * f32);
  const double y = f42 + (1.0 - f32) * x;

  return cv::Point2d(x, y);
}

// https://pomax.github.io/bezierinfo/#canonical
bool BezierCurve::is_simple_curve_cubic(double margin) const
{
  if (is_collinear())
  {
    return true;
  }

  const cv::Point2d p = canonical_form_cubic_free();

  // Check if curve has one inflection.
  if (p.y > 1.0 + margin)
  {
    return true;
  }

  // Check if curve is plain.
  if (p.x > 1.0)
  {
    return true;
  }
  else if (p.x > 0.0)
  {
    double loop_curve = 0.5 * (std::sqrt(3.0 * (4.0*p.x - p.x*p.x)) - p.x);
    if (p.y < loop_curve)
    {
      return true;
    }
  }
  else
  {
    double loop_curve = (-p.x*p.x + 3.0*p.x) / 3.0;
    if (p.y < loop_curve)
    {
      return true;
    }
  }

  return false;
}

cv::Mat BezierCurve::compute_separating_curve_mask(cv::Rect region)
{
  struct FloodFillData
  {
    FloodFillData(int p_x, int p_y, unsigned char color) :
      p_x(p_x),
      p_y(p_y),
      color(color)
    {}

    int p_x, p_y;
    unsigned char color;
  };

  cv::Mat mask = cv::Mat::zeros(region.size(), CV_8UC1);
  
  std::vector<BezierCurve> curves({*this});

  const cv::Point2d p_left_1(eval(0.0));
  if (region.contains(p_left_1))
  {
    std::vector<double> dist_vals({std::abs(p_left_1.x - region.tl().x),
      std::abs(p_left_1.y - region.tl().y),
      std::abs(p_left_1.x - region.br().x),
      std::abs(p_left_1.y - region.br().y)});
    size_t index_min = std::distance(dist_vals.begin(), std::min_element(dist_vals.begin(), dist_vals.end()));

    cv::Point2d p_left_2;
    if (index_min == 0)
    {
      p_left_2 = cv::Point2d(region.tl().x, p_left_1.y);
    }
    else if (index_min == 1)
    {
      p_left_2 = cv::Point2d(p_left_1.x, region.tl().y);
    }
    else if (index_min == 2)
    {
      p_left_2 = cv::Point2d(region.br().x, p_left_1.y);
    }
    else
    {
      p_left_2 = cv::Point2d(p_left_1.x, region.br().y);
    }

    const std::vector<cv::Point2d> control_points({
      p_left_2,
      p_left_2 + (p_left_1 - p_left_2) / 3.0,
      p_left_2 + (p_left_1 - p_left_2) * 2.0 / 3.0,
      p_left_1});

    curves.emplace(curves.begin(), control_points);
  }

  const cv::Point2d p_right_1(eval(1.0));
  if (region.contains(p_right_1))
  {
    std::vector<double> dist_vals({std::abs(p_right_1.x - region.tl().x),
                                   std::abs(p_right_1.y - region.tl().y),
                                   std::abs(p_right_1.x - region.br().x),
                                   std::abs(p_right_1.y - region.br().y)});
    size_t index_min = std::distance(dist_vals.begin(), std::min_element(dist_vals.begin(), dist_vals.end()));

    cv::Point2d p_right_2;
    if (index_min == 0)
    {
      p_right_2 = cv::Point2d(region.tl().x, p_right_1.y);
    }
    else if (index_min == 1)
    {
      p_right_2 = cv::Point2d(p_right_1.x, region.tl().y);
    }
    else if (index_min == 2)
    {
      p_right_2 = cv::Point2d(region.br().x, p_right_1.y);
    }
    else
    {
      p_right_2 = cv::Point2d(p_right_1.x, region.br().y);
    }

    const std::vector<cv::Point2d> control_points({
      p_right_1,
      p_right_1 + (p_right_2 - p_right_1) / 3.0,
      p_right_1 + (p_right_2 - p_right_1) * 2.0 / 3.0,
      p_right_2});
    curves.emplace_back(control_points);
  }

  for (BezierCurve& curve : curves)
  {
    curve -= region.tl();
  }
   
  // Compute derivatives.
  std::vector<BezierCurve> curves_deriv;
  for (const BezierCurve& curve : curves)
  {
    curves_deriv.push_back(curve.deriv());
  }

  // Draw curves onto mask.
  for (const BezierCurve& curve : curves)
  {
    curve.draw<unsigned char>(mask, 255);
  }

  std::stack<FloodFillData> fill_stack;
  for (size_t i = 0; i < curves.size(); ++i)
  {
    for (const CurveDrawPoint& d : curves[i].get_draw_points())
    {
      const cv::Vec2d tangent = cv::normalize(cv::Vec2d(curves_deriv[i].eval(d.t)));
      const cv::Vec2d normal(tangent[1], -tangent[0]);

      const int delta_x = normal[0] >= 0.383 ? 1 : normal[0] >= -0.383 ? 0 : -1;
      const int delta_y = normal[1] >= 0.383 ? 1 : normal[1] >= -0.383 ? 0 : -1;

      fill_stack.emplace(d.p.x + delta_x, d.p.y + delta_y, 1);
      fill_stack.emplace(d.p.x - delta_x, d.p.y - delta_y, 2);
    }
  }

  // Perform flood fill.
  while (!fill_stack.empty())
  {
    const FloodFillData d = fill_stack.top();
    fill_stack.pop();

    if (d.p_x >= 0 && d.p_x < mask.cols && d.p_y >= 0 && d.p_y < mask.rows)
    {
      unsigned char& val = mask.at<unsigned char>(d.p_y, d.p_x);
      if (val == 0)
      {
        val = d.color;
        fill_stack.emplace(d.p_x-1, d.p_y, d.color);
        fill_stack.emplace(d.p_x+1, d.p_y, d.color);
        fill_stack.emplace(d.p_x, d.p_y-1, d.color);
        fill_stack.emplace(d.p_x, d.p_y+1, d.color);
      }
    }
  }

  return mask;
}

double BezierCurve::curve_length() const
{
  double dist = 0.0;
  for (int i = 1; i <= m_degree; ++i)
  {
    const cv::Point2d diff = m_control_points[i] - m_control_points[i+1];
    dist += std::sqrt(diff.x * diff.x + diff.y * diff.y);
  }
  return dist;
}

int BezierCurve::degree() const
{
  return m_degree;
}

int BezierCurve::num_control_points() const
{
  return static_cast<int>(m_control_points.size());
}

double BezierCurve::find_x(double x) const
{
  double x_start = eval_x(0.0);
  double x_end = eval_x(1.0);

  double x_min = std::min(x_start, x_end);
  double x_max = std::max(x_start, x_end);

  if (x < x_min || x > x_max)
  {
    return -1.0;
  }

  int dir = x_end > x_start ? 1 : -1;

  double t = 0.5;
  double x_t = eval_x(t);
  double step = 0.25;

  while (std::abs(x_t - x) > 1.0e-6)
  {
    if (x_t < x)
    {
      t += dir * step;
    }
    else
    {
      t -= dir * step;
    }
    x_t = eval_x(t);
    step *= 0.5;
  }

  return t;
}

double BezierCurve::find_y(double y) const
{
  double y_start = eval_y(0.0);
  double y_end = eval_y(1.0);

  double y_min = std::min(y_start, y_end);
  double y_max = std::max(y_start, y_end);

  if (y < y_min || y > y_max)
  {
    return -1.0;
  }

  int dir = y_end > y_start ? 1 : -1;

  double t = 0.5;
  double y_t = eval_y(t);
  double step = 0.25;

  while (std::abs(y_t - y) > 1.0e-6)
  {
    if (y_t < y)
    {
      t += dir * step;
    }
    else
    {
      t -= dir * step;
    }
    y_t = eval_y(t);
    step *= 0.5;
  }

  return t;
}

BezierCurve BezierCurve::average(const BezierCurve& b1, const BezierCurve& b2)
{
  if (b1.m_degree != b2.m_degree)
  {
    throw(std::invalid_argument("Curve degrees differ."));
  }

  std::vector<cv::Point2d> control_points(b1.m_degree+1);
  for (int i = 0; i <= b1.m_degree; ++i)
  {
    control_points[i] = 0.5 * (b1.m_control_points[i] + b2.m_control_points[i]);
  }
  return BezierCurve(control_points);
}

// https://pomax.github.io/bezierinfo/#boundingbox
cv::Rect2d BezierCurve::bounding_box_cubic() const
{
  BezierCurve curve_deriv = deriv(1);
  cv::Rect bbox_deriv = curve_deriv.bounding_box();

  const cv::Vec2d root_a(curve_deriv.m_control_points[0] - 2.0 * curve_deriv.m_control_points[1] + curve_deriv.m_control_points[2]);
  const cv::Vec2d root_b(2.0 * (curve_deriv.m_control_points[1] - curve_deriv.m_control_points[0]));
  const cv::Vec2d root_c(curve_deriv.m_control_points[0]);
  const cv::Vec2d root_sqrt = root_b.mul(root_b) - 4.0 * root_a.mul(root_c);

  std::vector<double> t_x, t_y;
  if (root_sqrt[0] >= 0.0)
  {
    t_x.push_back((-root_b[0] + std::sqrt(root_sqrt[0])) / (2.0 * root_a[0]));
    t_x.push_back((-root_b[0] - std::sqrt(root_sqrt[0])) / (2.0 * root_a[0]));
  }

  if (root_sqrt[1] >= 0.0)
  {
    t_y.push_back((-root_b[1] + std::sqrt(root_sqrt[1])) / (2.0 * root_a[1]));
    t_y.push_back((-root_b[1] - std::sqrt(root_sqrt[1])) / (2.0 * root_a[1]));
  }
  
  double x1 = std::min(m_control_points[0].x, m_control_points[3].x);
  double x2 = std::max(m_control_points[0].x, m_control_points[3].x);
  double y1 = std::min(m_control_points[0].y, m_control_points[3].y);
  double y2 = std::max(m_control_points[0].y, m_control_points[3].y);

  for (double t : t_x)
  {
    if (t > 0.0 && t < 1.0)
    {
      const cv::Point2d p = eval(t);
      x1 = std::min(x1, p.x);
      x2 = std::max(x2, p.x);
      y1 = std::min(y1, p.y);
      y2 = std::max(y2, p.y);
    }
  }

  for (double t : t_y)
  {
    if (t > 0.0 && t < 1.0)
    {
      const cv::Point2d p = eval(t);
      x1 = std::min(x1, p.x);
      x2 = std::max(x2, p.x);
      y1 = std::min(y1, p.y);
      y2 = std::max(y2, p.y);
    }
  }

  return cv::Rect2d(x1, y1, x2 - x1, y2 - y1);
}

cv::Rect BezierCurve::bounding_box() const
{
  if (m_draw_points.empty())
  {
    get_draw_points();
  }

  std::vector<cv::Point> points(m_draw_points.size());
  for (size_t i = 0; i < m_draw_points.size(); ++i)
  {
    points[i] = m_draw_points[i].p;
  }

  return cv::boundingRect(points);
}

std::vector<BezierCurve> BezierCurve::fit(const std::vector<double>& points, int degree, int subpatch_size, double weight_reg, double weight_slope)
{
  if (points.empty())
  {
    throw(std::invalid_argument("Empty point vector."));
  }

  std::vector<BezierCurve> curves;

  const int num_segments = (static_cast<int>(points.size()) - 1) / subpatch_size;
  cv::Mat A = cv::Mat::zeros((subpatch_size + degree + 2) * num_segments - 1, degree * num_segments + 1, CV_64FC1);
  std::vector<double> t_vals = linspace(0.0, 1.0, subpatch_size + 1, true);
  cv::Mat basis = bernstein_basis(degree, t_vals);

  for (int i = 0; i < num_segments; ++i)
  {

    basis.copyTo(A(cv::Rect(degree*i, (subpatch_size + 1)*i, basis.cols, basis.rows)));
    if (i > 0)
    {
      A.at<double>(A.rows - i, degree*i - 1) = -1.0 * weight_reg;
      A.at<double>(A.rows - i, degree*i) = 2.0 * weight_reg;
      A.at<double>(A.rows - i, degree*i + 1) = -1.0 * weight_reg;
    }
  }

  for (int i = 0; i < degree * num_segments; ++i)
  {
    A.at<double>(A.rows-num_segments-i, i+1) = -1.0 * weight_slope;
    A.at<double>(A.rows-num_segments-i, i) = 1.0 * weight_slope;
  }

  cv::Mat Ainv = pseudo_inverse(A);

  cv::Mat rhs((subpatch_size + degree + 2) * num_segments - 1, 1, CV_64FC1);

  int c = 0;
  for (int i = 0; i < num_segments; ++i)
  {
    for (auto iter = points.begin() + i*subpatch_size; iter != points.begin() + (i + 1)*subpatch_size + 1; ++iter)
    {
      rhs.at<double>(c++, 0) = *iter;
    }
  }
  for (int i = 0; i < (degree + 1) * num_segments - 1; ++i)
  {
    rhs.at<double>(c++, 0) = 0.0;
  }

  cv::Mat coeffs = Ainv * rhs;

  for (int i = 0; i < num_segments; ++i)
  {
    const std::vector<double> x_vals = linspace(
      static_cast<double>(i * subpatch_size),
      static_cast<double>((i+1) * subpatch_size),
      degree+1, true);

    const std::vector<double> y_vals = std::vector<double>(
      coeffs.begin<double>() + i*degree,
      coeffs.begin<double>() + (i+1) * degree + 1);

    curves.emplace_back(x_vals, y_vals);
  }

  return curves;
}

void BezierCurve::draw(cv::Mat image, int subpatch_size, double factor) const
{
  subpatch_size = static_cast<int>(subpatch_size * factor);

  for (int i = 0; i < subpatch_size; ++i)
  {
    const double t = static_cast<double>(2 * i + 1) / static_cast<double>(2 * subpatch_size);
    const int x = static_cast<int>(factor * eval_x(t));
    const int y = static_cast<int>(factor * eval_y(t));
    if (0 <= y && y < image.rows && 0 <= x && x < image.cols)
    {
      image.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 255, 255);
    }
  }
}

static bool is_8_connected(cv::Point p1, cv::Point p2)
{
  const cv::Point d = p1 - p2;
  return std::abs(d.x) <= 1 && std::abs(d.y) <= 1;
}

BezierCurve BezierCurve::deriv(int order) const
{
  BezierCurve curve = *this;

  for (int i = 0; i < order; ++i)
  {
    if (curve.degree() <= 0)
    {
      return BezierCurve();
    }

    std::vector<cv::Point2d> control_points_deriv(curve.degree());
    for (int i = 0; i < curve.degree(); ++i)
    {
      control_points_deriv[i] = curve.degree() * (curve.control_point(i+1) - curve.control_point(i));
    }
    
    curve = BezierCurve(control_points_deriv);
  }

  return curve;
}

cv::Vec2d BezierCurve::tangent(double t) const
{
  return cv::normalize(cv::Vec2d(this->deriv().eval(t)));
}

cv::Vec2d BezierCurve::normal(double t) const
{
  cv::Vec2d tangent_vector = tangent(t);
  return cv::Vec2d(tangent_vector[1], -tangent_vector[0]);
}

static void intersect_cubic(std::vector<std::pair<double, double>>& t_vals, const BezierCurve& b1, const BezierCurve& b2, double t11, double t12, double t21, double t22, double tolerance = 1.0e-3)
{
  if (b1.degree() != 3 || b2.degree() != 3)
  {
    throw(std::invalid_argument("Bezier::intersect_cubic called with non-cubic curve."));
  }

  cv::Rect2d bbox_1 = b1.bounding_box_cubic();
  cv::Rect2d bbox_2 = b2.bounding_box_cubic();

  if (bbox_1.tl().x > bbox_2.br().x ||
      bbox_1.br().x < bbox_2.tl().x ||
      bbox_1.tl().y > bbox_2.br().y ||
      bbox_1.br().y < bbox_2.tl().y)
  {
    return;
  }

  const double t1 = 0.5 * (t11 + t12);
  const double t2 = 0.5 * (t21 + t22);

  cv::Rect2d bbox = bbox_1 | bbox_2;

  if (bbox.width < tolerance && bbox.height < tolerance)
  {
    t_vals.emplace_back(t1, t2);
    return;
  }

  BezierCurve b11, b12, b21, b22;
  std::tie(b11, b12) = b1.split_cubic(0.5);
  std::tie(b21, b22) = b2.split_cubic(0.5);

  intersect_cubic(t_vals, b11, b21, t11, t1, t21, t2, tolerance);
  intersect_cubic(t_vals, b12, b21, t1, t12, t21, t2, tolerance);
  intersect_cubic(t_vals, b11, b22, t11, t1, t2, t22, tolerance);
  intersect_cubic(t_vals, b12, b22, t1, t12, t2, t22, tolerance);
}

struct ClusterData
{
  ClusterData(const std::pair<double, double>& t, const cv::Point2d& p) :
    t_vals(1, t),
    centroid(p),
    num_points(1)
  {
  }

  void insert(const std::pair<double, double>& t, const cv::Point2d& p)
  {
    t_vals.push_back(t);
    centroid = (num_points * centroid + p) / (num_points + 1);
    ++num_points;
  }

  std::vector<std::pair<double, double>> t_vals;
  cv::Point2d centroid;
  int num_points;
};

static bool intersect_cubic(BezierCurve b1, BezierCurve b2, double& t_result_1, double& t_result_2)
{
  std::vector<std::pair<double, double>> t_vals;
  intersect_cubic(t_vals, b1, b2, 0.0, 1.0, 0.0, 1.0, 1.0e-3);

  if (!t_vals.empty())
  {
    std::vector<ClusterData> clusters;

    for (const std::pair<double, double>& t : t_vals)
    {
      bool inserted = false;
      cv::Point2d p = b1.eval(t.first);
      for (ClusterData& cluster : clusters)
      {
        if (cv::norm(cluster.centroid - p) < 1.0e-2)
        {
          cluster.insert(t, p);
          inserted = true;
          break;
        }
      }
      if (!inserted)
      {
        clusters.emplace_back(t, p);
      }
    }

    if (clusters.size() > 1)
    {
      throw(std::runtime_error("intersect_cubic: Multiple intersection points not implemented"));
    }

    double t1 = 0.0;
    double t2 = 0.0;
    for (const std::pair<double, double>& t : clusters[0].t_vals)
    {
      t1 += t.first;
      t2 += t.second;
    }

    t_result_1 = t1 / clusters[0].num_points;
    t_result_2 = t2 / clusters[0].num_points;

    return true;
  }

  return false;
}

std::pair<double, double> BezierCurve::intersect_curves_cubic(const BezierCurve& b1, const BezierCurve& b2)
{
  if (b1.degree() != 3 || b2.degree() != 3)
  {
    throw(std::invalid_argument("Bezier::intersect_cubic called with non-cubic curve."));
  }

  double t1, t2;
  if (intersect_cubic(b1, b2, t1, t2))
  {
    return std::make_pair(t1, t2);
  }
  else
  {
    return std::make_pair(-1.0, -1.0);
  }
}

std::pair<CurveIntersection, CurveIntersection> BezierCurve::intersect_curves_cubic(const std::vector<BezierCurve>& curve_1, bool from_left_1, const std::vector<BezierCurve>& curve_2, bool from_left_2)
{
  const int size_1 = static_cast<int>(curve_1.size());
  const int size_2 = static_cast<int>(curve_2.size());

  double t1 = -1.0;
  double t2 = -1.0;

  if (from_left_1 && from_left_2)
  {
    for (int i = 0; i < size_1; ++i)
    {
      for (int j = 0; j < size_2; ++j)
      {
        if (intersect_cubic(curve_1[i], curve_2[j], t1, t2))
        {
          return std::pair<CurveIntersection, CurveIntersection>({ i, t1 }, { j, t2 });
        }
      }
    }
  }
  else if (!from_left_1 && from_left_2)
  {
    for (int i = size_1-1; i >= 0; --i)
    {
      for (int j = 0; j < size_2; ++j)
      {
        if (intersect_cubic(curve_1[i], curve_2[j], t1, t2))
        {
          return std::pair<CurveIntersection, CurveIntersection>({ i, t1 }, { j, t2 });
        }
      }
    }
  }
  else if (from_left_1 && !from_left_2)
  {
    for (int i = 0; i < size_1; ++i)
    {
      for (int j = size_2-1; j >= 0; --j)
      {
        if (intersect_cubic(curve_1[i], curve_2[j], t1, t2))
        {
          return std::pair<CurveIntersection, CurveIntersection>({ i, t1 }, { j, t2 });
        }
      }
    }
  }
  else
  {
    for (int i = size_1-1; i >= 0; --i)
    {
      for (int j = size_2-1; j >= 0; --j)
      {
        if (intersect_cubic(curve_1[i], curve_2[j], t1, t2))
        {
          return std::pair<CurveIntersection, CurveIntersection>({ i, t1 }, { j, t2 });
        }
      }
    }
  }

  return std::pair<CurveIntersection, CurveIntersection>({ -1, t1 }, { -1, t2 });
}

void BezierCurve::trim_curves_cubic(std::vector<BezierCurve>& curve_1, bool trim_left_1, std::vector<BezierCurve>& curve_2, bool trim_left_2)
{
  std::pair<CurveIntersection, CurveIntersection> crossing = BezierCurve::intersect_curves_cubic(curve_1, !trim_left_1, curve_2, !trim_left_2);

  if (crossing.first.index < 0 || crossing.second.index < 0)
  {
    throw(std::runtime_error("Unable to find cut."));
  }

  if (trim_left_1)
  {
    curve_1.erase(curve_1.begin(), curve_1.begin() + crossing.first.index);
    curve_1.front() = curve_1.front().split_cubic(crossing.first.t).second;
  }
  else
  {
    curve_1.erase(curve_1.begin() + crossing.first.index + 1, curve_1.end());
    curve_1.back() = curve_1.back().split_cubic(crossing.first.t).first;
  }

  if (trim_left_2)
  {
    curve_2.erase(curve_2.begin(), curve_2.begin() + crossing.second.index);
    curve_2.front() = curve_2.front().split_cubic(crossing.second.t).second;
  }
  else
  {
    curve_2.erase(curve_2.begin() + crossing.second.index + 1, curve_2.end());
    curve_2.back() = curve_2.back().split_cubic(crossing.second.t).first;
  }

  // Finally, just to be absolutely on the safe side, set first point of second to last point of first.
  //curve_2.front().x_vals[0] = curve_1.back().x_vals[3];
  //curve_2.front().y_vals[0] = curve_1.back().y_vals[3];
}

/*
 * Returns a straight segment that extends this curve at t=0.
 */
BezierCurve BezierCurve::extend_curve_left() const
{
  if (m_degree < 1)
  {
    return BezierCurve();
  }

  const cv::Point2d delta = cv::normalize(cv::Vec2d(m_control_points[1] - m_control_points[0]));
  const std::vector<cv::Point2d> curve({
    m_control_points[0] - 48.0 * delta,
    m_control_points[0] - 32.0 * delta,
    m_control_points[0] - 16.0 * delta,
    m_control_points[0]});

  return BezierCurve(curve);
}

/*
 * Returns a straight segment that extends this curve at t=1.
 */
BezierCurve BezierCurve::extend_curve_right() const
{
  if (m_degree < 1)
  {
    return BezierCurve();
  }

  const cv::Point2d delta = cv::normalize(cv::Vec2d(m_control_points[m_degree] - m_control_points[m_degree-1]));
  const std::vector<cv::Point2d> curve({
    m_control_points[m_degree],
    m_control_points[m_degree] + 16.0 * delta,
    m_control_points[m_degree] + 32.0 * delta,
    m_control_points[m_degree] + 48.0 * delta});

  return BezierCurve(curve);
}

double BezierCurve::max_curvature_norm_approx() const
{
  if (m_draw_points.empty())
  {
    generate_draw_points();
  }
  
  BezierCurve dd = deriv(2);
  double max_curvature = 0.0;
  
  for (const CurveDrawPoint& p : m_draw_points)
  {
    max_curvature = std::max(max_curvature, cv::norm(cv::Vec2d(dd.eval(p.t)), cv::NORM_L2SQR));
  }

  return std::sqrt(max_curvature);
}

double BezierCurve::max_gradient_norm_approx() const
{
  if (m_draw_points.empty())
  {
    generate_draw_points();
  }

  BezierCurve dd = deriv(1);
  double max_gradient = 0.0;

  for (const CurveDrawPoint& p : m_draw_points)
  {
    max_gradient = std::max(max_gradient, cv::norm(cv::Vec2d(dd.eval(p.t)), cv::NORM_L2SQR));
  }

  return std::sqrt(max_gradient);
}

double BezierCurve::min_dist_approx(const cv::Point2d& p) const
{
  double min_dist = std::numeric_limits<double>::max();

  if (m_draw_points.empty())
  {
    generate_draw_points();
  }

  for (const CurveDrawPoint& d : m_draw_points)
  {
    double dist = cv::norm(cv::Vec2d(d.p.x, d.p.y) - cv::Vec2d(p), cv::NORM_L2SQR);
    if (dist < min_dist)
    {
      min_dist = dist;
    }
  }

  return std::sqrt(min_dist);
}

/*
* ComputeLeftTangent, ComputeRightTangent, ComputeCenterTangent :
* Approximate unit tangents at endpoints and "center" of digitized curve.
*/
static cv::Vec2d compute_left_tangent(const std::vector<cv::Point2d>& points, int end)
{
  return cv::normalize(cv::Vec2d(points[end+1] - points[end]));
}

static cv::Vec2d compute_right_tangent(const std::vector<cv::Point2d>& points, int end)
{
  return cv::normalize(cv::Vec2d(points[end-1] - points[end]));
}

/*
*  ChordLengthParameterize :
*	Assign parameter values to digitized points
*	using relative distances between points.
*/
static std::vector<double> chord_length_parameterize(const std::vector<cv::Point2d> points, int first, int last)
{
  std::vector<double> u(last - first + 1);

  u[0] = 0.0;
  for (int i = first+1; i <= last; ++i)
  {
    u[i-first] = u[i-first-1] + cv::norm(cv::Vec2d(points[i] - points[i-1]), cv::NORM_L2);
  }

  for (size_t i = 1; i < u.size(); ++i)
  {
    u[i] /= u.back();
  }

  return(u);
}

/*
*  B0, B1, B2, B3 :
*	 Bezier multipliers
*/
static double B0(double u)
{
  const double tmp = 1.0 - u;
  return tmp * tmp * tmp;
}

static double B1(double u)
{
  const double tmp = 1.0 - u;
  return 3.0 * u * (tmp * tmp);
}

static double B2(double u)
{
  double tmp = 1.0 - u;
  return 3.0 * u * u * tmp;
}

static double B3(double u)
{
  return u * u * u;
}

/*
*  GenerateBezier :
*  Use least-squares method to find Bezier control points for region.
*
*/
static BezierCurve generate_bezier(const std::vector<cv::Point2d>& points, int first, int last, const std::vector<double> u_prime, const cv::Vec2d& t_hat_1, const cv::Vec2d& t_hat_2)
{
  const int num_points = last - first + 1;

  std::vector<cv::Point2d> curve(4);
  std::vector<double[2][2]> A(num_points);
  double C[2][2];
  double X[2];
  
  for (int i = 0; i < num_points; ++i)
  {
    const cv::Vec2d v1 = t_hat_1 * B1(u_prime[i]);
    const cv::Vec2d v2 = t_hat_2 * B2(u_prime[i]);
    A[i][0][0] = v1[0];
    A[i][0][1] = v1[1];
    A[i][1][0] = v2[0];
    A[i][1][1] = v2[1];
  }

  C[0][0] = 0.0;
  C[0][1] = 0.0;
  C[1][0] = 0.0;
  C[1][1] = 0.0;
  X[0] = 0.0;
  X[1] = 0.0;

  for (int i = 0; i < num_points; ++i)
  {
    C[0][0] += A[i][0][0] * A[i][0][0] + A[i][0][1] * A[i][0][1];
    C[0][1] += A[i][0][0] * A[i][1][0] + A[i][0][1] * A[i][1][1];
    C[1][0] = C[0][1];
    C[1][1] += A[i][1][0] * A[i][1][0] + A[i][1][1] * A[i][1][1];

    const cv::Vec2d tmp = points[first+i]
      - points[first] * B0(u_prime[i])
      - points[first] * B1(u_prime[i])
      - points[last] * B2(u_prime[i])
      - points[last] * B3(u_prime[i]);

    X[0] += A[i][0][0] * tmp[0] + A[i][0][1] * tmp[1];
    X[1] += A[i][1][0] * tmp[0] + A[i][1][1] * tmp[1];
  }

  /* Compute the determinants of C and X	*/
  double det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
  double det_C0_X  = C[0][0] * X[1]    - C[0][1] * X[0];
  double det_X_C1  = X[0]    * C[1][1] - X[1]    * C[0][1];

  /* Finally, derive alpha values	*/
  if (det_C0_C1 == 0.0)
  {
    det_C0_C1 = (C[0][0] * C[1][1]) * 10e-12;
  }

  double alpha_l = det_X_C1 / det_C0_C1;
  double alpha_r = det_C0_X / det_C0_C1;

  /*  If alpha negative, use the Wu/Barsky heuristic (see text) */
  /* (if alpha is 0, you get coincident control points that lead to
  * divide by zero in any subsequent NewtonRaphsonRootFind() call. */
  if (alpha_l < 1.0e-6 || alpha_r < 1.0e-6)
  {
    double	dist = cv::norm(cv::Vec2d(points[last] - points[first]), cv::NORM_L2) / 3.0;
    curve[0] = points[first];
    curve[1] = points[first] + cv::Point2d(dist * t_hat_1);
    curve[2] = points[last] + cv::Point2d(dist * t_hat_2);
    curve[3] = points[last];
    return curve;
  }

  /*  First and last control points of the Bezier curve are */
  /*  positioned exactly at the first and last data points */
  /*  Control points 1 and 2 are positioned an alpha distance out */
  /*  on the tangent vectors, left and right, respectively */
  curve[0] = points[first];
  curve[1] = points[first] + cv::Point2d(alpha_l * t_hat_1);
  curve[2] = points[last] + cv::Point2d(alpha_r * t_hat_2);
  curve[3] = points[last];

  return curve;
}

/*
*  NewtonRaphsonRootFind :
*	Use Newton-Raphson iteration to find better root.
*/
static double newton_raphson_root_find(const BezierCurve& curve, const cv::Point2d& P, double u)
{
  /* Compute Q(u)	*/
  const cv::Point2d Q_u = curve.eval(u);

  /* Generate control vertices for Q'	*/
  BezierCurve Q1(2);
  for (int i = 0; i <= 2; i++)
  {
    Q1.set_control_point(i, (curve.control_point(i+1) - curve.control_point(i)) * 3.0);
  }

  /* Generate control vertices for Q'' */
  BezierCurve Q2(1);
  for (int i = 0; i <= 1; i++)
  {
    Q2.set_control_point(i, (Q1.control_point(i+1) - Q1.control_point(i)) * 2.0);
  }

  /* Compute Q'(u) and Q''(u)	*/
  const cv::Point2d Q1_u = Q1.eval(u);
  const cv::Point2d Q2_u = Q2.eval(u);

  /* Compute f(u)/f'(u) */
  const double numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
  const double denominator = (Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
    (Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y);

  /* u = u - f(u)/f'(u) */
  return u - numerator / denominator;
}

/*
*  Reparameterize:
*	Given set of points and their parameterization, try to find
*   a better parameterization.
*
*/
static std::vector<double> reparameterize(const std::vector<cv::Point2d>& points, int first, int last, const std::vector<double>& u, const BezierCurve& curve)
{
  const int num_points = last - first + 1;
  std::vector<double> u_prime(num_points);

  for (int i = first; i <= last; ++i)
  {
    u_prime[i-first] = newton_raphson_root_find(curve, points[i], u[i-first]);
  }

  return u_prime;
}

/*
*  ComputeMaxError :
*	Find the maximum squared distance of digitized points
*	to fitted curve.
*/
static double compute_max_error(const std::vector<cv::Point2d>& points, int first, int last, const BezierCurve& curve, const std::vector<double>& u)
{
  double max_dist = 0.0;

  for (int i = first + 1; i < last; ++i)
  {
    const cv::Point2d P = curve.eval(u[i-first]);
    const double dist = cv::norm(cv::Vec2d(P - points[i]), cv::NORM_L2SQR);
    if (dist > max_dist)
    {
      max_dist = dist;
    }
  }

  return max_dist;
}

BezierCurve fit_cubic(const std::vector<cv::Point2d>& points, int first, int last, const cv::Vec2d& t_hat_1, const cv::Vec2d& t_hat_2)
{
  BezierCurve curve(3);
  const int num_points = last - first + 1;

  if (num_points == 2)
  {
    double dist = cv::norm(cv::Vec2d(points[last] - points[first]), cv::NORM_L2) / 3.0;
    curve.set_control_point(0, points[first]);
    curve.set_control_point(1, points[first] + cv::Point2d(dist * t_hat_1));
    curve.set_control_point(2, points[last] + cv::Point2d(dist * t_hat_2));
    curve.set_control_point(3, points[last]);
    return curve;
  }

  std::vector<double> u = chord_length_parameterize(points, first, last);
  curve = generate_bezier(points, first, last, u, t_hat_1, t_hat_2);
  double error = compute_max_error(points, first, last, curve, u);

  const int max_iterations = 8;
  for (int i = 0; i < max_iterations; ++i)
  {
    const std::vector<double> u_prime = reparameterize(points, first, last, u, curve);
    const BezierCurve curve_prime = generate_bezier(points, first, last, u_prime, t_hat_1, t_hat_2);
    const double error_prime = compute_max_error(points, first, last, curve_prime, u_prime);

    if (error_prime >= error)
    {
      break;
    }

    u = u_prime;
    curve = curve_prime;
    error = error_prime;
  }

  return curve;
}

/*
An Algorithm for Automatically Fitting Digitized Curves
by Philip J. Schneider
from "Graphics Gems", Academic Press, 1990
*/
BezierCurve BezierCurve::fit_cubic(std::vector<cv::Point2d> points)
{
  const int num_points = static_cast<int>(points.size());

  sort_pca(points);

  cv::Vec2d t_hat_1 = compute_left_tangent(points, 0);
  cv::Vec2d t_hat_2 = compute_right_tangent(points, num_points-1);
  
  return ::fit_cubic(points, 0, num_points-1, t_hat_1, t_hat_2);
}

BezierCurve BezierCurve::fit_cubic(const cv::Point2d& p1, const cv::Point2d& p2, const std::vector<cv::Point2d>& points)
{
  std::vector<cv::Point2d> points_copy;
  points_copy.push_back(p1);
  points_copy.insert(points_copy.end(), points.begin(), points.end());
  points_copy.push_back(p2);

  const int num_points = static_cast<int>(points_copy.size());

  cv::Vec2d t_hat_1 = compute_left_tangent(points_copy, 0);
  cv::Vec2d t_hat_2 = compute_right_tangent(points_copy, num_points-1);

  return ::fit_cubic(points_copy, 0, num_points-1, t_hat_1, t_hat_2);
}

void BezierCurve::generate_draw_points() const
{
  m_draw_points.clear();

  if (m_degree < 1)
  {
    return;
  }

  m_draw_points.push_back({0.0, eval(0.0)});
  m_draw_points.push_back({1.0, eval(1.0)});

  // Generate points to draw using bisection.
  bool finished;
  do
  {
    finished = true;
    for (size_t i = 0; i < m_draw_points.size() - 1; ++i)
    {
      if (!is_8_connected(m_draw_points[i].p, m_draw_points[i+1].p))
      {
        const double t_mid = 0.5 * (m_draw_points[i].t + m_draw_points[i+1].t);
        m_draw_points.insert(m_draw_points.begin()+i+1, {t_mid, eval(t_mid)});
        finished = false;
        break;
      }
    }
  } while (!finished);

  // Clear duplicates.
  do
  {
    finished = true;
    for (size_t i = 0; i < m_draw_points.size() - 1; ++i)
    {
      if (m_draw_points[i].p == m_draw_points[i+1].p)
      {
        m_draw_points.erase(m_draw_points.begin()+i+1);
        finished = false;
        break;
      }
    }
  } while (!finished);
}

double BezierCurve::bernstein_polynomial(int degree, int index, double t)
{
  return math::n_choose_k(degree, index) * std::pow(1.0 - t, degree - index) * std::pow(t, index);
}

cv::Mat BezierCurve::bernstein_basis(int degree, const std::vector<double>& t_vals)
{
  cv::Mat basis(static_cast<int>(t_vals.size()), degree + 1, CV_64FC1);
  for (int y = 0; y < basis.rows; ++y)
  {
    double* ptr = reinterpret_cast<double*>(basis.ptr(y));
    for (int x = 0; x < basis.cols; ++x)
    {
      ptr[x] = bernstein_polynomial(degree, x, t_vals[y]);
    }
  }

  return basis;
}

cv::Mat BezierCurve::pseudo_inverse(cv::Mat A)
{
  cv::Mat At;
  cv::transpose(A, At);
  cv::Mat AtAinv;
  cv::invert(At*A, AtAinv);
  return AtAinv * At;
}

void BezierCurve::draw_debug() const
{
  if (m_draw_points.empty())
  {
    generate_draw_points();
  }

  std::vector<cv::Point> points;
  for (const CurveDrawPoint& p : m_draw_points)
  {
    points.push_back(p.p);
  }

  if (!points.empty())
  {
    cv::Rect bbox = cv::boundingRect(points);
    cv::Mat image_debug = cv::Mat::zeros(bbox.size(), CV_8UC1);
    for (const cv::Point& p : points)
    {
      image_debug.at<unsigned char>(p-bbox.tl()) = 255;
    }
    cv::resize(image_debug, image_debug, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);
    cv::imshow("Image Debug", image_debug);
    cv::waitKey();
  }

}

BezierCurve::BezierCurve() :
  m_degree(-1)
{}

BezierCurve::BezierCurve(int degree) :
  m_degree(degree)
{
  m_control_points.resize(degree+1);
}

BezierCurve::BezierCurve(const std::vector<double>& control_points_x, const std::vector<double>& control_points_y)
{
  m_degree = static_cast<int>(control_points_x.size()) - 1;
  m_control_points.resize(control_points_x.size());
  for (int i = 0; i <= m_degree; ++i)
  {
    m_control_points[i] = cv::Point2d(control_points_x[i], control_points_y[i]);
  }
}

BezierCurve::BezierCurve(const std::vector<cv::Point2d>& control_points) :
  m_degree(static_cast<int>(control_points.size())-1),
  m_control_points(control_points)
{}

BezierCurve::BezierCurve(const cv::Point2d& p_start, const cv::Point2d& p_end, int degree) :
  m_degree(degree)
{
  const cv::Point2d& diff = p_end - p_start;
  for (int i = 0; i <= degree; ++i)
  {
    const double t = static_cast<double>(i) / static_cast<double>(degree);
    m_control_points.push_back(p_start + t * diff);
  }
}

std::ostream& operator<<(std::ostream& os, const BezierCurve& curve)
{
  os << "[";
  for (int i = 0; i <= curve.m_degree; ++i)
  {
    os << " (" << curve.m_control_points[i].x << " " << curve.m_control_points[i].y << ")";
  }
  os << " ]";
  return os;
}

boost::property_tree::ptree BezierCurve::save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const
{
  boost::property_tree::ptree tree;
  serialize(tree, "control_points", m_control_points, base_path, path);
  serialize(tree, "degree", m_degree, base_path, path);
  return tree;
}

void BezierCurve::load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree)
{
  deserialize(tree, "control_points", m_control_points, base_path);
  deserialize(tree, "degree", m_degree, base_path);
}
