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

#ifndef TRLIB_PATCH_REGION_HPP_
#define TRLIB_PATCH_REGION_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

#include "bezier_curve.hpp"
#include "serializable.hpp"

class PatchRegion : public Serializable
{
public:
  PatchRegion() = default;
  PatchRegion(const PatchRegion&) = default;
  PatchRegion(int target_index, const std::vector<cv::Point>& points);
  PatchRegion(int target_index, const cv::Point& coordinate, const BezierCurve& curve_top, const BezierCurve& curve_bot, const BezierCurve& curve_left, const BezierCurve& curve_right, const BezierCurve& curve_diag=BezierCurve(), cv::Rect bounding_box=cv::Rect());
  PatchRegion(int target_index, const cv::Point& coordinate, const std::vector<BezierCurve>& curves_top, const std::vector<BezierCurve>& curves_bot, const std::vector<BezierCurve>& curves_left, const std::vector<BezierCurve>& curves_right, const std::vector<BezierCurve>& curves_diag=std::vector<BezierCurve>(), cv::Rect bounding_box=cv::Rect());
  PatchRegion(int target_index, const cv::Point& coordinate, const cv::Rect& bounding_box);

  void set_sub_regions(const BezierCurve& separating_curve);
  bool diag_origin_tl() const;
  bool diag_origin_bl() const;

  cv::Rect bounding_box() const
  {
    return m_bounding_box;
  }

  cv::Point anchor() const
  {
    return m_bounding_box.tl();
  }

  bool is_rectangular() const
  {
    return m_rectangular;
  }

  bool has_diag_curve() const
  {
    return !m_curves_diag.empty();
  }

  bool has_sub_regions() const
  {
    return !m_sub_regions.empty();
  }

  static void fix_fabricability(PatchRegion& patch_tr, PatchRegion& patch_br, PatchRegion& patch_bl, PatchRegion& patch_tl);

  cv::Mat mask() const
  {
    return m_mask;
  }

  cv::Mat mask(cv::Rect bounding_box) const;

  cv::Size size() const
  {
    return m_bounding_box.size();
  }

  std::vector<PatchRegion>& sub_regions()
  {
    return m_sub_regions;
  }

  const std::vector<PatchRegion>& sub_regions() const
  {
    return m_sub_regions;
  }

  cv::Point bbox_tl() const
  {
    return m_bounding_box.tl();
  }

  cv::Point2f curve_tl() const
  {
    return m_tl;
  }

  void remove_sub_regions()
  {
    m_sub_regions.clear();
  }

  int target_index() const
  {
    return m_target_index;
  }

  void set_target_index(int target_index)
  {
    m_target_index = target_index;
    for (PatchRegion& sub_region : m_sub_regions)
    {
      sub_region.m_target_index = target_index;
    }
  }

  cv::Point coordinate() const
  {
    return m_coordinate;
  }

  const std::vector<BezierCurve>& curves_top() const
  {
    return m_curves_top;
  }

  std::vector<BezierCurve> curves_top_reversed() const
  {
    std::vector<BezierCurve> curves = m_curves_top;
    std::reverse(curves.begin(), curves.end());
    for (BezierCurve& c : curves)
    {
      c = c.reversed();
    }
    return curves;
  }

    const std::vector<BezierCurve>& curves_bot() const
  {
    return m_curves_bot;
  }

  std::vector<BezierCurve> curves_bot_reversed() const
  {
    std::vector<BezierCurve> curves = m_curves_bot;
    std::reverse(curves.begin(), curves.end());
    for (BezierCurve& c : curves)
    {
      c = c.reversed();
    }
    return curves;
  }

    const std::vector<BezierCurve>& curves_left() const
  {
    return m_curves_left;
  }

  std::vector<BezierCurve> curves_left_reversed() const
  {
    std::vector<BezierCurve> curves = m_curves_left;
    std::reverse(curves.begin(), curves.end());
    for (BezierCurve& c : curves)
    {
      c = c.reversed();
    }
    return curves;
  }

    const std::vector<BezierCurve>& curves_right() const
  {
    return m_curves_right;
  }

  std::vector<BezierCurve> curves_right_reversed() const
  {
    std::vector<BezierCurve> curves = m_curves_right;
    std::reverse(curves.begin(), curves.end());
    for (BezierCurve& c : curves)
    {
      c = c.reversed();
    }
    return curves;
  }

    const std::vector<BezierCurve>& curves_diag() const
  {
    return m_curves_diag;
  }

  std::vector<BezierCurve> curves_diag_reversed() const
  {
    std::vector<BezierCurve> curves = m_curves_diag;
    std::reverse(curves.begin(), curves.end());
    for (BezierCurve& c : curves)
    {
      c = c.reversed();
    }
    return curves;
  }

  bool valid() const;

  PatchRegion scaled(double scale) const;
  void scale(double scale);

  template <typename T> void draw_edge_mask(cv::Mat image, const T& color, int dilate, cv::Point offset=cv::Point(0, 0)) const;
  template<typename T> void draw(cv::Mat image, const T& color, double scale=1.0) const;
  template <typename T, typename Iter> static void draw(cv::Mat image, const T& color, double scale, Iter begin, Iter end);

  std::vector<cv::Point> get_draw_points() const;

  std::vector<BezierCurve> get_bezier_curves() const;

  std::string get_id_string() const;
  std::string get_short_id_string() const;

  virtual void load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const override;

private:
  void init_curves(cv::Rect bounding_box);
  void init_rect();
  cv::Mat draw(double scale);

  void update_bounding_box(cv::Rect bounding_box=cv::Rect());
  
  int m_target_index;
  cv::Point m_coordinate;
  
  std::vector<BezierCurve> m_curves_top, m_curves_bot, m_curves_left, m_curves_right, m_curves_diag;
  cv::Rect m_bounding_box;
  cv::Mat m_mask;

  cv::Point2d m_tl;

  std::vector<PatchRegion> m_sub_regions;

  bool m_rectangular;
};

#include "patch_region_impl.hpp"

#endif /* TRLIB_PATCH_REGION_HPP_ */