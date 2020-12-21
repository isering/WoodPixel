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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

#include "flood_fill.hpp"
#include "patch_region.hpp"

const std::string output_characters = "0123456789ABCDEFGHIJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

void PatchRegion::init_curves(cv::Rect bounding_box)
{
  std::vector<cv::Rect> bounding_box_vec;

  /*
  if (m_curves_top.empty(), m_curves_left.empty() && m_curves_diag.empty())
  {
    throw(std::invalid_argument("Top, left, and diag curves empty in PatchRegion::init_curves."));
  }
  */

  if (!m_curves_top.empty())
  {
    m_tl = m_curves_top[0].eval(0.0);
  }
  else if (!m_curves_left.empty())
  {
    m_tl = m_curves_left[0].eval(0.0);
  }
  else if (!m_curves_diag.empty())
  {
    m_tl = m_curves_diag[0].eval(0.0);
  }
  else
  {
    m_tl.x = -1;
    m_tl.y = -1;
  }

  update_bounding_box(bounding_box);
}

void PatchRegion::update_bounding_box(cv::Rect bounding_box)
{
  std::vector<std::reference_wrapper<BezierCurve>> curves;
  curves.insert(curves.end(), m_curves_top.begin(), m_curves_top.end());
  curves.insert(curves.end(), m_curves_bot.begin(), m_curves_bot.end());
  curves.insert(curves.end(), m_curves_left.begin(), m_curves_left.end());
  curves.insert(curves.end(), m_curves_right.begin(), m_curves_right.end());
  curves.insert(curves.end(), m_curves_diag.begin(), m_curves_diag.end());

  if (curves.empty())
  {
    std::vector<cv::Point> draw_points;
    cv::findNonZero(m_mask, draw_points);
    return;
  }

  std::vector<cv::Point> draw_points;
  for (auto iter = curves.begin(); iter != curves.end(); ++iter)
  {
    std::vector<CurveDrawPoint> curve_draw_points = iter->get().get_draw_points();
    for (const CurveDrawPoint &p : curve_draw_points)
    {
      draw_points.push_back(p.p);
    }
  }
  m_bounding_box = bounding_box == cv::Rect() ? cv::boundingRect(draw_points) : bounding_box;

  m_mask = cv::Mat::zeros(m_bounding_box.size(), CV_8UC1);

  for (const BezierCurve &c : curves)
  {
    (c - m_bounding_box.tl()).draw<unsigned char>(m_mask, 255);
  }
  flood_fill_patch<unsigned char>(m_mask, 255);

  m_rectangular = (cv::countNonZero(255 - m_mask) == 0 ? true : false);
}

void PatchRegion::init_rect()
{
  m_mask = cv::Mat(m_bounding_box.size(), CV_8UC1, 255);
  m_rectangular = true;
}

void PatchRegion::set_sub_regions(const BezierCurve &separating_curve)
{
  if (separating_curve.degree() != 3)
  {
    throw(std::invalid_argument("PatchRegion::set_sub_regions: Only cubic Bezier curves supported."));
  }

  const std::vector<BezierCurve> no_curves;
  const std::vector<BezierCurve> curves_diag(1, separating_curve);

  const double d1 = cv::norm(m_curves_top[0].control_point(0) - separating_curve.control_point(0));
  const double d2 = cv::norm(m_curves_bot[0].control_point(0) - separating_curve.control_point(0));
  const bool origin_tl = d1 < d2 ? true : false;

  m_sub_regions.clear();
  if (origin_tl)
  {
    m_sub_regions.emplace_back(m_target_index, m_coordinate, m_curves_top, no_curves, no_curves, m_curves_right, curves_diag);
    m_sub_regions.emplace_back(m_target_index, m_coordinate, no_curves, m_curves_bot, m_curves_left, no_curves, curves_diag);
  }
  else
  {
    m_sub_regions.emplace_back(m_target_index, m_coordinate, m_curves_top, no_curves, m_curves_left, no_curves, curves_diag);
    m_sub_regions.emplace_back(m_target_index, m_coordinate, no_curves, m_curves_bot, no_curves, m_curves_right, curves_diag);
  }
}

cv::Mat PatchRegion::draw(double scale)
{
  cv::Mat image = cv::Mat::zeros(static_cast<int>(m_mask.rows * scale), static_cast<int>(m_mask.cols * scale), CV_8UC3);

  for (const BezierCurve &c : m_curves_top)
  {
    c.scaled(scale).draw<cv::Vec3b>(image, cv::Vec3b(255, 255, 0), scale * anchor());
  }
  for (const BezierCurve &c : m_curves_bot)
  {
    c.scaled(scale).draw<cv::Vec3b>(image, cv::Vec3b(255, 255, 0), scale * anchor());
  }
  for (const BezierCurve &c : m_curves_left)
  {
    c.scaled(scale).draw<cv::Vec3b>(image, cv::Vec3b(255, 255, 0), scale * anchor());
  }
  for (const BezierCurve &c : m_curves_right)
  {
    c.scaled(scale).draw<cv::Vec3b>(image, cv::Vec3b(255, 255, 0), scale * anchor());
  }
  for (const BezierCurve &c : m_curves_diag)
  {
    c.scaled(scale).draw<cv::Vec3b>(image, cv::Vec3b(255, 255, 0), scale * anchor());
  }

  return image;
}

bool PatchRegion::diag_origin_tl() const
{
  if (m_curves_diag.empty() || m_curves_diag[0].is_empty_curve())
  {
    return false;
  }

  cv::Point2d p_tl, p_bl;

  bool has_tl = true;
  if (!m_curves_top.empty() && !m_curves_top[0].is_empty_curve())
  {
    p_tl = m_curves_top[0].control_point(0);
  }
  else if (!m_curves_left.empty() && !m_curves_left[0].is_empty_curve())
  {
    p_tl = m_curves_left[0].control_point(0);
  }
  else
  {
    has_tl = false;
  }

  bool has_bl = true;
  if (!m_curves_bot.empty() && !m_curves_bot[0].is_empty_curve())
  {
    p_bl = m_curves_bot[0].control_point(0);
  }
  else if (!m_curves_left.empty() && !m_curves_left.back().is_empty_curve())
  {
    p_bl = m_curves_left.back().control_point(m_curves_left.back().num_control_points() - 1);
  }
  else
  {
    has_bl = false;
  }

  if (!has_tl && !has_bl)
  {
    throw(std::invalid_argument("PatchRegion::diag_origin_tl encountered invalid patch."));
  }
  else if (!has_tl)
  {
    return false;
  }
  else if (!has_bl)
  {
    return true;
  }
  else
  {
    const double d1 = cv::norm(p_tl - m_curves_diag[0].control_point(0));
    const double d2 = cv::norm(p_bl - m_curves_diag[0].control_point(0));

    return d1 < d2;
  }
}

bool PatchRegion::diag_origin_bl() const
{
  if (m_curves_diag.empty() || m_curves_diag[0].is_empty_curve())
  {
    return false;
  }

  cv::Point2d p_tl, p_bl;

  bool has_tl = true;
  if (!m_curves_top.empty() && !m_curves_top[0].is_empty_curve())
  {
    p_tl = m_curves_top[0].control_point(0);
  }
  else if (!m_curves_left.empty() && !m_curves_left[0].is_empty_curve())
  {
    p_tl = m_curves_left[0].control_point(0);
  }
  else
  {
    has_tl = false;
  }

  bool has_bl = true;
  if (!m_curves_bot.empty() && !m_curves_bot[0].is_empty_curve())
  {
    p_bl = m_curves_bot[0].control_point(0);
  }
  else if (!m_curves_left.empty() && !m_curves_left.back().is_empty_curve())
  {
    p_bl = m_curves_left.back().control_point(m_curves_left.back().num_control_points() - 1);
  }
  else
  {
    has_bl = false;
  }

  if (!has_tl && !has_bl)
  {
    throw(std::invalid_argument("PatchRegion::diag_origin_bl encountered invalid patch."));
  }
  else if (!has_tl)
  {
    return true;
  }
  else if (!has_bl)
  {
    return false;
  }
  else
  {
    const double d1 = cv::norm(p_tl - m_curves_diag[0].control_point(0));
    const double d2 = cv::norm(p_bl - m_curves_diag[0].control_point(0));

    return d2 < d1;
  }
}

static double get_relative_angle(const cv::Point2d &p, const cv::Point2d &cp_1, const cv::Point2d &cp_2)
{
  const double pi = boost::math::constants::pi<double>();

  const cv::Point2d dcp_1 = cp_1 - p;
  const cv::Point2d dcp_2 = cp_2 - p;

  double a_1 = std::atan2(dcp_1.y, dcp_1.x);
  double a_2 = std::atan2(dcp_2.y, dcp_2.x);

  double rel_a = a_2 - a_1;

  while (rel_a < -pi)
  {
    rel_a += 2.0 * pi;
  }

  while (rel_a > pi)
  {
    rel_a -= 2.0 * pi;
  }

  return rel_a;
}

void PatchRegion::fix_fabricability(PatchRegion &patch_tr, PatchRegion &patch_br, PatchRegion &patch_bl, PatchRegion &patch_tl)
{
  std::vector<std::reference_wrapper<cv::Point2d>> control_points;

  control_points.emplace_back(patch_tr.m_curves_bot[0].control_point(1));

  if (patch_br.has_sub_regions() && patch_br.sub_regions()[0].diag_origin_tl())
  {
    control_points.emplace_back(patch_br.sub_regions()[0].m_curves_diag[0].control_point(1));
  }

  control_points.emplace_back(patch_br.m_curves_left[0].control_point(1));

  if (patch_bl.has_sub_regions() && patch_bl.sub_regions()[0].diag_origin_bl())
  {
    control_points.emplace_back(patch_bl.sub_regions()[0].m_curves_diag.back().control_point(2));
  }

  control_points.emplace_back(patch_bl.m_curves_top.back().control_point(2));

  if (patch_tl.has_sub_regions() && patch_tl.sub_regions()[0].diag_origin_tl())
  {
    control_points.emplace_back(patch_tl.sub_regions()[0].m_curves_diag.back().control_point(2));
  }

  control_points.emplace_back(patch_tl.m_curves_right.back().control_point(2));

  if (patch_tr.has_sub_regions() && patch_tr.sub_regions()[0].diag_origin_bl())
  {
    control_points.emplace_back(patch_tr.sub_regions()[0].m_curves_diag[0].control_point(1));
  }

  const cv::Point2d p = patch_tr.m_curves_bot[0].control_point(0);

  auto cp_iter = control_points.begin();
  while (cp_iter != control_points.end())
  {
    auto cp_iter_2 = cp_iter + 1;
    if (cp_iter_2 == control_points.end())
    {
      cp_iter_2 = control_points.begin();
    }

    const double rel_angle = get_relative_angle(p, cp_iter->get(), cp_iter_2->get());

    if (rel_angle < 0.0 && rel_angle > -0.5 * boost::math::constants::pi<double>())
    {
      std::swap(cp_iter->get(), cp_iter_2->get());
      cp_iter = control_points.begin();
    }
    else
    {
      ++cp_iter;
    }
  }

  patch_br.m_curves_top[0].control_point(1) = patch_tr.m_curves_bot[0].control_point(1);
  for (PatchRegion &sub_patch : patch_br.sub_regions())
  {
    if (!sub_patch.m_curves_top.empty())
    {
      sub_patch.m_curves_top[0].control_point(1) = patch_tr.m_curves_bot[0].control_point(1);
    }
    if (!sub_patch.m_curves_left.empty())
    {
      sub_patch.m_curves_left[0].control_point(1) = patch_br.m_curves_left[0].control_point(1);
    }
  }

  patch_bl.m_curves_right[0].control_point(1) = patch_br.m_curves_left[0].control_point(1);
  for (PatchRegion &sub_patch : patch_bl.sub_regions())
  {
    if (!sub_patch.m_curves_right.empty())
    {
      sub_patch.m_curves_right[0].control_point(1) = patch_br.m_curves_left[0].control_point(1);
    }
    if (!sub_patch.m_curves_top.empty())
    {
      sub_patch.m_curves_top.back().control_point(2) = patch_bl.m_curves_top.back().control_point(2);
    }
  }

  patch_tl.m_curves_bot.back().control_point(2) = patch_bl.m_curves_top.back().control_point(2);
  for (PatchRegion &sub_patch : patch_tl.sub_regions())
  {
    if (!sub_patch.m_curves_bot.empty())
    {
      sub_patch.m_curves_bot.back().control_point(2) = patch_bl.m_curves_top.back().control_point(2);
    }
    if (!sub_patch.m_curves_right.empty())
    {
      sub_patch.m_curves_right.back().control_point(2) = patch_tl.m_curves_right.back().control_point(2);
    }
  }

  patch_tr.m_curves_left.back().control_point(2) = patch_tl.m_curves_right.back().control_point(2);
  for (PatchRegion &sub_patch : patch_tr.sub_regions())
  {
    if (!sub_patch.m_curves_left.empty())
    {
      sub_patch.m_curves_left.back().control_point(2) = patch_tl.m_curves_right.back().control_point(2);
    }
    if (!sub_patch.m_curves_bot.empty())
    {
      sub_patch.m_curves_bot[0].control_point(1) = patch_tr.m_curves_bot[0].control_point(1);
    }
  }

  for (size_t i = 1; i < patch_br.sub_regions().size(); ++i)
  {
    if (patch_br.sub_regions()[0].diag_origin_tl())
    {
      patch_br.sub_regions()[i].m_curves_diag[0].control_point(1) = patch_br.sub_regions()[0].m_curves_diag[0].control_point(1);
    }
  }

  for (size_t i = 1; i < patch_bl.sub_regions().size(); ++i)
  {
    if (patch_bl.sub_regions()[0].diag_origin_bl())
    {
      patch_bl.sub_regions()[i].m_curves_diag.back().control_point(2) = patch_bl.sub_regions()[0].m_curves_diag.back().control_point(2);
    }
  }

  for (size_t i = 1; i < patch_tl.sub_regions().size(); ++i)
  {
    if (patch_tl.sub_regions()[0].diag_origin_tl())
    {
      patch_tl.sub_regions()[i].m_curves_diag.back().control_point(2) = patch_tl.sub_regions()[0].m_curves_diag.back().control_point(2);
    }
  }

  for (size_t i = 1; i < patch_tr.sub_regions().size(); ++i)
  {
    if (patch_tr.sub_regions()[0].diag_origin_bl())
    {
      patch_tr.sub_regions()[i].m_curves_diag[0].control_point(1) = patch_tr.sub_regions()[0].m_curves_diag[0].control_point(1);
    }
  }
}

std::vector<cv::Point> PatchRegion::get_draw_points() const
{
  std::vector<cv::Point> points;
  for (const BezierCurve &curve : m_curves_top)
  {
    const std::vector<CurveDrawPoint> curve_points = curve.get_draw_points();
    for (const CurveDrawPoint &c : curve_points)
    {
      points.push_back(c.p);
    }
  }
  for (const BezierCurve &curve : m_curves_bot)
  {
    const std::vector<CurveDrawPoint> curve_points = curve.get_draw_points();
    for (const CurveDrawPoint &c : curve_points)
    {
      points.push_back(c.p);
    }
  }
  for (const BezierCurve &curve : m_curves_left)
  {
    const std::vector<CurveDrawPoint> curve_points = curve.get_draw_points();
    for (const CurveDrawPoint &c : curve_points)
    {
      points.push_back(c.p);
    }
  }
  for (const BezierCurve &curve : m_curves_right)
  {
    const std::vector<CurveDrawPoint> curve_points = curve.get_draw_points();
    for (const CurveDrawPoint &c : curve_points)
    {
      points.push_back(c.p);
    }
  }
  for (const BezierCurve &curve : m_curves_diag)
  {
    const std::vector<CurveDrawPoint> curve_points = curve.get_draw_points();
    for (const CurveDrawPoint &c : curve_points)
    {
      points.push_back(c.p);
    }
  }

  if (points.empty())
  {
    cv::findNonZero(m_mask, points);
  }

  return points;
}

PatchRegion PatchRegion::scaled(double scale) const
{
  PatchRegion p(*this);
  p.scale(scale);
  return p;
}

void PatchRegion::scale(double scale)
{
  for (PatchRegion &p : m_sub_regions)
  {
    p.scale(scale);
  }

  std::transform(m_curves_top.begin(), m_curves_top.end(), m_curves_top.begin(), std::bind(&BezierCurve::scaled, std::placeholders::_1, scale));
  std::transform(m_curves_bot.begin(), m_curves_bot.end(), m_curves_bot.begin(), std::bind(&BezierCurve::scaled, std::placeholders::_1, scale));
  std::transform(m_curves_left.begin(), m_curves_left.end(), m_curves_left.begin(), std::bind(&BezierCurve::scaled, std::placeholders::_1, scale));
  std::transform(m_curves_right.begin(), m_curves_right.end(), m_curves_right.begin(), std::bind(&BezierCurve::scaled, std::placeholders::_1, scale));
  std::transform(m_curves_diag.begin(), m_curves_diag.end(), m_curves_diag.begin(), std::bind(&BezierCurve::scaled, std::placeholders::_1, scale));

  if (m_curves_top.empty() && m_curves_bot.empty() && m_curves_left.empty() && m_curves_right.empty() && m_curves_diag.empty())
  {
    cv::resize(m_mask, m_mask, cv::Size(), scale, scale, cv::INTER_NEAREST);
    m_bounding_box.x = static_cast<int>(m_bounding_box.x * scale);
    m_bounding_box.y = static_cast<int>(m_bounding_box.y * scale);
    m_bounding_box.width = m_mask.cols;
    m_bounding_box.height = m_mask.rows;
  }
  else
  {
    update_bounding_box();
  }
}

cv::Mat PatchRegion::mask(cv::Rect bounding_box) const
{
  const cv::Rect bbox_union = bounding_box & m_bounding_box;

  cv::Mat mask_out = mask_out.zeros(bounding_box.size(), CV_8UC1);
  if (bbox_union.area() > 0)
  {
    cv::Rect rect_source(bbox_union.tl() - m_bounding_box.tl(), bbox_union.size());
    cv::Rect rect_target(bbox_union.tl() - bounding_box.tl(), bbox_union.size());
    mask_out(rect_target).setTo(m_mask(rect_source));
  }

  return mask_out;
}

/*
 * Check PatchRegion for validity, i.e. whether its position and bounding box are valid and 
 * whether it contains a valid configuration of boundary patches.
*/
bool PatchRegion::valid() const
{
  if (m_tl.x < 0.0 || m_tl.y < 0.0 || m_bounding_box.width <= 0 || m_bounding_box.height <= 0)
  {
    return false;
  }

  if (!m_curves_top.empty() && !m_curves_bot.empty() && !m_curves_left.empty() && !m_curves_right.empty() && m_curves_diag.empty())
  {
    return true;
  }

  if (!m_curves_diag.empty() && diag_origin_tl() && !m_curves_right.empty() && !m_curves_top.empty())
  {
    return true;
  }

  if (!m_curves_diag.empty() && diag_origin_tl() && !m_curves_left.empty() && !m_curves_bot.empty())
  {
    return true;
  }

  if (!m_curves_diag.empty() && diag_origin_bl() && !m_curves_right.empty() && !m_curves_bot.empty())
  {
    return true;
  }

  if (!m_curves_diag.empty() && diag_origin_bl() && !m_curves_left.empty() && !m_curves_top.empty())
  {
    return true;
  }

  return false;
}

PatchRegion::PatchRegion(int target_index, const std::vector<cv::Point> &points) : m_target_index(target_index),
                                                                                   m_coordinate(0, 0)
{
  m_bounding_box = cv::boundingRect(points);
  m_tl = m_bounding_box.tl();
  m_mask = cv::Mat::zeros(m_bounding_box.size(), CV_8UC1);

  for (const cv::Point &p : points)
  {
    m_mask.at<unsigned char>(p - m_bounding_box.tl()) = 255;
  }

  m_rectangular = cv::countNonZero(m_mask) == m_mask.size().area() ? true : false;
}

PatchRegion::PatchRegion(int target_index, const cv::Point &coordinate, const BezierCurve &curve_top, const BezierCurve &curve_bot, const BezierCurve &curve_left, const BezierCurve &curve_right, const BezierCurve &curve_diag, cv::Rect bounding_box) : m_target_index(target_index),
                                                                                                                                                                                                                                                           m_coordinate(coordinate)
{
  if (!curve_top.is_empty_curve())
  {
    m_curves_top.push_back(curve_top);
  }

  if (!curve_bot.is_empty_curve())
  {
    m_curves_bot.push_back(curve_bot);
  }

  if (!curve_left.is_empty_curve())
  {
    m_curves_left.push_back(curve_left);
  }

  if (!curve_right.is_empty_curve())
  {
    m_curves_right.push_back(curve_right);
  }

  if (!curve_diag.is_empty_curve())
  {
    m_curves_diag.push_back(curve_diag);
  }

  init_curves(bounding_box);
}

PatchRegion::PatchRegion(int target_index, const cv::Point &coordinate, const std::vector<BezierCurve> &curves_top, const std::vector<BezierCurve> &curves_bot, const std::vector<BezierCurve> &curves_left, const std::vector<BezierCurve> &curves_right, const std::vector<BezierCurve> &curves_diag, cv::Rect bounding_box) : m_target_index(target_index),
                                                                                                                                                                                                                                                                                                                                 m_coordinate(coordinate),
                                                                                                                                                                                                                                                                                                                                 m_curves_top(curves_top),
                                                                                                                                                                                                                                                                                                                                 m_curves_bot(curves_bot),
                                                                                                                                                                                                                                                                                                                                 m_curves_left(curves_left),
                                                                                                                                                                                                                                                                                                                                 m_curves_right(curves_right),
                                                                                                                                                                                                                                                                                                                                 m_curves_diag(curves_diag)
{
  init_curves(bounding_box);
}

PatchRegion::PatchRegion(int target_index, const cv::Point &coordinate, const cv::Rect &bounding_box) : m_target_index(target_index),
                                                                                                        m_coordinate(coordinate),
                                                                                                        m_bounding_box(bounding_box)
{
  init_rect();
}

std::vector<BezierCurve> PatchRegion::get_bezier_curves() const
{
  std::vector<BezierCurve> curves;
  curves.insert(curves.end(), m_curves_top.begin(), m_curves_top.end());
  curves.insert(curves.end(), m_curves_right.begin(), m_curves_right.end());
  curves.insert(curves.end(), m_curves_bot.begin(), m_curves_bot.end());
  curves.insert(curves.end(), m_curves_left.begin(), m_curves_left.end());
  curves.insert(curves.end(), m_curves_diag.begin(), m_curves_diag.end());
  return curves;
}

std::string PatchRegion::get_id_string() const
{
  return (boost::format("%02d%02d%02d") % m_target_index % m_coordinate.x % m_coordinate.y).str();
}

std::string PatchRegion::get_short_id_string() const
{
  std::string id_string;
  id_string += output_characters[m_target_index % output_characters.size()];
  id_string += ' ';
  id_string += output_characters[m_coordinate.x % output_characters.size()];
  id_string += output_characters[m_coordinate.y % output_characters.size()];
  return id_string;
}

boost::property_tree::ptree PatchRegion::save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const
{
  boost::property_tree::ptree tree;
  serialize(tree, "curves_top", m_curves_top, base_path, path);
  serialize(tree, "curves_bot", m_curves_bot, base_path, path);
  serialize(tree, "curves_left", m_curves_left, base_path, path);
  serialize(tree, "curves_right", m_curves_right, base_path, path);
  serialize(tree, "curves_diag", m_curves_diag, base_path, path);
  serialize(tree, "bounding_box", m_bounding_box, base_path, path);
  serialize_image(tree, "mask", m_mask, base_path, path);
  serialize(tree, "sub_regions", m_sub_regions, base_path, path);
  serialize(tree, "rectangular", m_rectangular, base_path, path);
  serialize(tree, "target_index", m_target_index, base_path, path);
  serialize(tree, "coordinate", m_coordinate, base_path, path);
  return tree;
}

void PatchRegion::load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree)
{
  deserialize(tree, "curves_top", m_curves_top, base_path);
  deserialize(tree, "curves_bot", m_curves_bot, base_path);
  deserialize(tree, "curves_left", m_curves_left, base_path);
  deserialize(tree, "curves_right", m_curves_right, base_path);
  deserialize(tree, "curves_diag", m_curves_diag, base_path);
  deserialize(tree, "bounding_box", m_bounding_box, base_path);
  deserialize_image(tree, "mask", m_mask, base_path);
  deserialize(tree, "sub_regions", m_sub_regions, base_path);
  deserialize(tree, "rectangular", m_rectangular, base_path);
  deserialize(tree, "target_index", m_target_index, base_path);
  deserialize(tree, "coordinate", m_coordinate, base_path);
}
