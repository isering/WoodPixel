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

#ifndef TRLIB_BEZIER_CURVE_
#define TRLIB_BEZIER_CURVE_

#include <numeric>
#include <random>
#include <vector>

#include <opencv2/opencv.hpp>

#include "linspace.hpp"
#include "math.hpp"
#include "serializable.hpp"

struct CurveIntersection
{
  int index;
  double t;
};

struct CurveDrawPoint
{
  double t;
  cv::Point p;
};

class BezierCurve : public Serializable
{
public:
  /*
   * Constructors
   */
  BezierCurve();
  explicit BezierCurve(int degree);
  BezierCurve(const std::vector<double> &control_points_x, const std::vector<double> &control_points_y);
  BezierCurve(const std::vector<cv::Point2d> &control_points);
  BezierCurve(const cv::Point2d &p_start, const cv::Point2d &p_end, int degree);

  /*
   * Evaluate curve
   */
  cv::Point2d eval(double t) const;
  double eval_x(double t) const;
  double eval_y(double t) const;
  cv::Vec2d normal(double t) const;
  cv::Vec2d tangent(double t) const;

  /*
   * Derivatives
   */
  BezierCurve deriv(int order = 1) const;

  /*
   * Get control points
   */
  const cv::Point2d &control_point(int index) const
  {
    return m_control_points[index];
  }

  cv::Point2d &control_point(int index)
  {
    return m_control_points[index];
  }

  /*
   * Manipulate control points
   */
  void set_control_point(int index, const cv::Point2d &p)
  {
    m_control_points[index] = p;
    clear_draw_points();
  }

  /*
   * Curve fitting
   */
  static std::vector<BezierCurve> fit(const std::vector<double> &points, int degree, int subpatch_size, double weight_reg, double weight_slope);
  static BezierCurve fit_cubic(std::vector<cv::Point2d> points);
  static BezierCurve fit_cubic(const cv::Point2d &p1, const cv::Point2d &p2, const std::vector<cv::Point2d> &points);

  /*
   * Curve intersection
   */

  /*
   * Curve information
   */
  cv::Rect bounding_box() const;
  cv::Rect2d bounding_box_cubic() const;
  double curve_length() const;
  int degree() const;
  bool is_collinear() const;
  bool is_empty_curve() const;
  double max_curvature_norm_approx() const;
  double max_gradient_norm_approx() const;
  double min_dist_approx(const cv::Point2d &p) const;
  int num_control_points() const;

  /*
   * Curve manpulations
   */
  template <typename T>
  void operator+=(const cv::Point_<T> &p);
  template <typename T>
  void operator-=(const cv::Point_<T> &p);
  template <typename T>
  void operator*=(const T &scalar);
  template <typename T>
  void operator/=(const T &scalar);
  BezierCurve reversed() const;
  BezierCurve scaled(double factor) const;
  BezierCurve transformed(cv::Mat T) const;

  /*
   * Batch manipulation
   */
  static std::vector<BezierCurve> transformed(const std::vector<BezierCurve> &curves, cv::Mat transformation_matrix);
  static std::vector<BezierCurve> translated(const std::vector<BezierCurve> &curves, const cv::Point2d &delta);

  /*
   * Drawing functions
   */
  void draw(cv::Mat image, int subpatch_size, double factor) const;
  template <typename T>
  void draw(cv::Mat image, const T &color, const cv::Point &offset = cv::Point(0, 0)) const;
  void draw_debug() const;

  std::vector<CurveDrawPoint> get_draw_points() const
  {
    if (m_draw_points.empty())
    {
      generate_draw_points();
    }
    return m_draw_points;
  }

  /*
   * Extent curve linearly
   */
  BezierCurve extend_curve_left() const;
  BezierCurve extend_curve_right() const;

  /*
   * Find points on curve
   */
  double find_x(double x) const;
  double find_y(double y) const;

  void swap_x_y()
  {
    for (cv::Point2d &p : m_control_points)
    {
      std::swap(p.x, p.y);
    }
    clear_draw_points();
  }

  /*
   * Operations for cubic Bezier curves
   */
  static std::pair<double, double> intersect_curves_cubic(const BezierCurve &b1, const BezierCurve &b2);
  static std::pair<CurveIntersection, CurveIntersection> intersect_curves_cubic(const std::vector<BezierCurve> &curve_1, bool from_left_1, const std::vector<BezierCurve> &curve_2, bool from_left_2);
  std::pair<BezierCurve, BezierCurve> split_cubic(double t) const;
  static void trim_curves_cubic(std::vector<BezierCurve> &curve_1, bool trim_front_1, std::vector<BezierCurve> &curve_2, bool trim_front_2);

  BezierCurve canonical_form_cubic() const;
  cv::Point2d canonical_form_cubic_free() const;
  bool is_simple_curve_cubic(double margin = 0.0) const;

  static BezierCurve average(const BezierCurve &b1, const BezierCurve &b2);

  friend std::ostream &operator<<(std::ostream &os, const BezierCurve &curve);

  cv::Mat compute_separating_curve_mask(cv::Rect region);

  virtual void load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const override;

private:
  std::vector<cv::Point2d> m_control_points;
  int m_degree;

  mutable std::vector<CurveDrawPoint> m_draw_points;

  static BezierCurve fit_cubic(const std::vector<cv::Point2d> &points, int first, int last, const cv::Vec2d &t_hat_1, const cv::Vec2d &t_hat_2);

  void generate_draw_points() const;
  void clear_draw_points() const
  {
    m_draw_points.clear();
  }

  static double bernstein_polynomial(int degree, int index, double t);
  static cv::Mat bernstein_basis(int degree, const std::vector<double> &t_vals);
  static cv::Mat pseudo_inverse(cv::Mat A);
};

template <typename T>
BezierCurve operator+(const BezierCurve &curve, const cv::Point_<T> &p);
template <typename T>
BezierCurve operator+(const cv::Point_<T> &p, const BezierCurve &curve);
template <typename T>
BezierCurve operator+(const BezierCurve &curve, const T &p);
template <typename T>
BezierCurve operator+(const T &p, const BezierCurve &curve);
template <typename T>
BezierCurve operator-(const BezierCurve &curve, const cv::Point_<T> &p);
template <typename T>
BezierCurve operator-(const cv::Point_<T> &p, const BezierCurve &curve);
template <typename T>
BezierCurve operator-(const BezierCurve &curve, const T &p);
template <typename T>
BezierCurve operator-(const T &p, const BezierCurve &curve);
template <typename T>
BezierCurve operator*(const BezierCurve &curve, const T &scalar);
template <typename T>
BezierCurve operator*(const T &scalar, const BezierCurve &curve);
template <typename T>
BezierCurve operator/(const BezierCurve &curve, const T &scalar);

#include "bezier_curve_impl.hpp"

#endif /* TRLIB_BEZIER_CURVE_ */
