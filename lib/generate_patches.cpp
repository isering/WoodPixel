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

#include <list>

#include "bezier_ransac.hpp"
#include "data_grid.hpp"
#include "generate_patches.hpp"
#include "patch_region.hpp"

struct dist_to_curve_less
{
  dist_to_curve_less(const BezierCurve& curve, double max_dist) :
    curve(curve),
    max_dist(max_dist)
  {
  }

  bool operator()(const cv::Point2d& p)
  {
    return curve.min_dist_approx(p) < max_dist;
  }

  BezierCurve curve;
  double max_dist;
};

struct dist_to_line_greater
{
  dist_to_line_greater(const cv::Point2d& p1, const cv::Point2d& p2, double max_dist) :
    p1(p1),
    p2(p2),
    max_dist(max_dist)
  {
    fac_1 = p2.y - p1.y;
    fac_2 = p2.x - p1.x;
    sum_1 = p2.x * p1.y - p2.y * p1.x;
    denom = cv::norm(p2 - p1);
  }

  bool operator()(const cv::Point2d& p)
  {
    const double dist = std::abs(fac_1 * p.x - fac_2 * p.y + sum_1) / denom;
    return dist > max_dist;
  }

  cv::Point2d p1, p2;
  double max_dist;
  double fac_1, fac_2, sum_1, denom;
};

struct projected_line_miss
{
  projected_line_miss(const cv::Point2d& p1, const cv::Point2d& p2) :
    p1(p1),
    p2(p2)
  {
    s = p2 - p1;
    denom = s.dot(s);
  }

  bool operator()(const cv::Point2d& p)
  {
    const double t = s.dot(p - p1) / denom;
    return t < 0.0 || t > 1.0;
  }

  cv::Point2d p1, p2;
  cv::Point2d s;
  double denom;
};

static void append_vector(std::vector<cv::Point2d>& vec1, const std::vector<cv::Point2d>& vec2)
{
  vec1.insert(vec1.end(), vec2.begin(), vec2.end());
}

static void erase_if_close_to_curve(std::vector<cv::Point2d>& vec, const BezierCurve& curve, double max_dist)
{
  vec.erase(
    std::remove_if(vec.begin(), vec.end(), dist_to_curve_less(curve, max_dist)),
    vec.end());
}

static BezierCurve fit_curve(DataGrid<Vector<cv::Point2d>>& grid, cv::Point c1, cv::Point c2)
{
  cv::Point d = c2 - c1;
  if (c1 == c2 || (std::abs(d.x) > 1 || std::abs(d.y) > 1))
  {
    throw(std::invalid_argument("Grid::FitCurve invalid points supplied."));
  }

  if ((c1.x > c2.x) || ((c1.x == c2.x) && (c1.y > c2.y)))
  {
    std::swap(c1, c2);
  }

  // Generate vector of points.
  Vector<cv::Point2d> points_to_fit;

  if (d.x == 1 && d.y == 0)
  {
    append_vector(points_to_fit, grid.data(c1.y-1, c1.x));
    append_vector(points_to_fit, grid.data(c1.y, c1.x));
  }
  else if (d.x == 1 && d.y == 1)
  {
    append_vector(points_to_fit, grid.data(c1.y, c1.x));
  }
  else if (d.x == 0 && d.y == 1)
  {
    append_vector(points_to_fit, grid.data(c1.y, c1.x-1));
    append_vector(points_to_fit, grid.data(c1.y, c1.x));
  }
  else if (d.x == 1 && d.y == -1)
  {
    append_vector(points_to_fit, grid.data(c1.y-1, c1.x));
  }

  /*
  * Remove points that are either too far away from the connecting line
  * or that when projected onto the line don't fall onto
  * the segment connecting c1 and c2.
  */
  const int num_points_total = static_cast<int>(points_to_fit.size());
  points_to_fit.erase(
    std::remove_if(
    points_to_fit.begin(), points_to_fit.end(),
    dist_to_line_greater(grid(c1), grid(c2), 0.2*cv::norm(grid(c1)-grid(c2)))),
    points_to_fit.end());
  points_to_fit.erase(
    std::remove_if(points_to_fit.begin(), points_to_fit.end(), projected_line_miss(grid(c1), grid(c2))), points_to_fit.end());
  const int num_points_local = static_cast<int>(points_to_fit.size());

  if (points_to_fit.empty())
  {
    return BezierCurve();
  }

  /*
  std::vector<cv::Point> draw_points(points_to_fit.begin(), points_to_fit.end());
  draw_points.push_back(grid(c1));
  draw_points.push_back(grid(c2));
  cv::Rect bbox = cv::boundingRect(draw_points);
  cv::Mat image_debug = cv::Mat::zeros(bbox.size(), CV_8UC1);
  for (size_t i = 0; i < draw_points.size()-2; ++i)
  {
    image_debug.at<unsigned char>(draw_points[i]-bbox.tl()) = 128;
  }
  image_debug.at<unsigned char>(draw_points[draw_points.size()-2]-bbox.tl()) = 255;
  image_debug.at<unsigned char>(draw_points[draw_points.size()-1]-bbox.tl()) = 255;
  cv::resize(image_debug, image_debug, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);
  cv::imshow("DEBUG", image_debug);
  cv::waitKey();
  */

  double inlier_ratio;
  //BezierCurve curve = fit_bezier_cubic_ransac(grid(c1), grid(c2), points_to_fit, 4, 50, &inlier_ratio);
  BezierCurve curve = fit_bezier_cubic_ransac(grid(c1), grid(c2), points_to_fit, 6, 100, &inlier_ratio);

  //std::cout << inlier_ratio << " " << (inlier_ratio / num_points_total) * num_points_local << std::endl;

  //if ((inlier_ratio / num_points_total) * num_points_local < 0.5)
  if (inlier_ratio < 0.7)
  {
    return BezierCurve();
  }

  if (d.x == 1 && d.y == 0)
  {
    erase_if_close_to_curve(grid.data(c1.y-1, c1.x), curve, 1.0);
    erase_if_close_to_curve(grid.data(c1.y, c1.x), curve, 1.0);
  }
  else if (d.x == 1 && d.y == 1)
  {
    erase_if_close_to_curve(grid.data(c1.y, c1.x), curve, 1.0);
  }
  else if (d.x == 0 && d.y == 1)
  {
    erase_if_close_to_curve(grid.data(c1.y, c1.x-1), curve, 1.0);
    erase_if_close_to_curve(grid.data(c1.y, c1.x), curve, 1.0);
  }
  else if (d.x == 1 && d.y == -1)
  {
    erase_if_close_to_curve(grid.data(c1.y-1, c1.x), curve, 1.0);
  }

  return curve;
}

static std::vector<BezierCurve> generate_smooth_curves(std::vector<cv::Point2d> knots, BezierCurve curve_left, BezierCurve curve_right)
{
  if (knots.size() < 2)
  {
    return std::vector<BezierCurve>();
  }

  const int n = static_cast<int>(knots.size())-1;

  /*
  * Assemble system of equations.
  */
  std::vector<double> a(n), b(n), c(n);
  std::vector<cv::Point2d> d(n);

  //if (curve_left.is_empty_curve())
  if (true)
  {
    a[0] = 0.0;
    b[0] = 2.0;
    c[0] = 1.0;
    d[0] = knots[0] + 2.0 * knots[1];
  }
  else
  {
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = 0.0;
    d[0] = 2.0 * knots[0] - curve_left.control_point(2);
  }

  for (int i = 1; i < n-1; ++i)
  {
    a[i] = 1.0;
    b[i] = 4.0;
    c[i] = 1.0;
    d[i] = 4.0 * knots[i] + 2.0 * knots[i+1];
  }

  //if (curve_right.is_empty_curve())
  if (true)
  {
    a[n-1] = 2.0;
    b[n-1] = 7.0;
    c[n-1] = 0.0;
    d[n-1] = 8.0 * knots[n-1] + knots[n];
  }
  else
  {
    a[n-1] = 1.0;
    b[n-1] = 4.0;
    c[n-1] = 0.0;
    d[n-1] = 4.0 * knots[n-1] + 2.0 * knots[n] - curve_right.control_point(1);
  }

  // Solve tridiagonal system of equations using the Thomas method.
  c[0] /= b[0];
  d[0] /= b[0];
  for (int i = 1; i < n; ++i)
  {
    c[i] /= b[i] - c[i-1] * a[i];
    d[i] = (d[i] - d[i-1] * a[i]) / (b[i] - c[i-1] * a[i]);
  }

  for (int i = n-2; i >= 0; --i)
  {
    d[i] -= c[i] * d[i+1];
  }

  std::vector<BezierCurve> curves(n, BezierCurve(3));
  for (int i = 0; i < n-1; ++i)
  {
    curves[i].set_control_point(0, knots[i]);
    curves[i].set_control_point(1, d[i]);
    curves[i].set_control_point(2, 2.0 * knots[i+1] - d[i+1]);
    curves[i].set_control_point(3, knots[i+1]);
  }

  //if (curve_right.is_empty_curve())
  if (true)
  {
    curves[n-1].set_control_point(0, knots[n-1]);
    curves[n-1].set_control_point(1, d[n-1]);
    curves[n-1].set_control_point(2, 0.5 * (d[n-1] + knots[n]));
    curves[n-1].set_control_point(3, knots[n]);
  }
  else
  {
    curves[n-1].set_control_point(0, knots[n-1]);
    curves[n-1].set_control_point(1, d[n-1]);
    curves[n-1].set_control_point(2, 2.0 * knots[n] - curve_right.control_point(1));
    curves[n-1].set_control_point(3, knots[n]);
  }

  return curves;
}

std::vector<PatchRegion> generate_patches_square(int target_index, const cv::Size& image_size, const cv::Size& patch_size, const cv::Size& filter_kernel_size)
{
  const cv::Size subpatch_size = patch_size / 4;

  const int rows = (image_size.height - filter_kernel_size.height - subpatch_size.height) / (patch_size.height - subpatch_size.height);
  const int cols = (image_size.width - filter_kernel_size.width - subpatch_size.width) / (patch_size.width - subpatch_size.width);

  std::vector<PatchRegion> reconstruction_regions;
  for (int y = 0; y < rows; ++y)
  {
    for (int x = 0; x < cols; ++x)
    {
      reconstruction_regions.emplace_back(target_index, cv::Point(x, y), cv::Rect(
        filter_kernel_size.width / 2 + x * (patch_size.width - subpatch_size.width),
        filter_kernel_size.height / 2 + y * (patch_size.height - subpatch_size.height),
        patch_size.width,
        patch_size.height));
    }
  }

  return reconstruction_regions;
}

DataGrid<Vector<cv::Point2d>> extract_edge_grid(const Grid& grid, cv::Mat edge_mat)
{
  // Extract edge points.
  std::list<cv::Point2f> edge_points;
  for (int y = 0; y < edge_mat.rows; ++y)
  {
    const unsigned char* ptr = edge_mat.ptr(y);
    for (int x = 0; x < edge_mat.cols; ++x)
    {
      if (ptr[x])
      {
        edge_points.emplace_back(static_cast<float>(x), static_cast<float>(y));
      }
    }
  }

  // Sort edge points into grid.
  DataGrid<Vector<cv::Point2d>> edge_grid(grid);
  for (int y = 0; y < edge_grid.rows()-1; ++y)
  {
    for (int x = 0; x < edge_grid.cols()-1; ++x)
    {
      std::vector<cv::Point2f> cell_contour({
        edge_grid(y, x),
        edge_grid(y, x+1),
        edge_grid(y+1, x+1),
        edge_grid(y+1, x)});

      std::list<cv::Point2f>::const_iterator iter = edge_points.begin();
      while (iter != edge_points.end())
      {
        if (cv::pointPolygonTest(cell_contour, *iter, false) >= 0.0)
        {
          edge_grid.data(y, x).push_back(*iter);
          edge_points.erase(iter++);
        }
        else
        {
          ++iter;
        }
      }
    }
  }

  return edge_grid;
}

void fit_single_edge(DataGrid<CurveData>& curve_grid, DataGrid<Vector<cv::Point2d>>& edge_grid, cv::Point p_1, cv::Point p_2)
{
  if (p_1.x > p_2.x || (p_1.x == p_2.x && p_1.y > p_2.y))
  {
    std::swap(p_1, p_2);
  }

  const cv::Point d = p_2 - p_1;
  if (d.x < 0 || d.y < -1 || d.x > 1 || d.y > 1 || (d.x == 0 && d.y == 0) || (d.x == 0 && d.y < 0))
  {
    throw(std::invalid_argument("fit_single_edge: Invalid grid points specified."));
  }

  BezierCurve c = fit_curve(edge_grid, p_1, p_2);
  if (!c.is_empty_curve())
  {
    if (d.x == 1 && d.y == 1)
    {
      curve_grid.data(p_1).curve_diag = c;
    }
    else if (d.x == 1 && d.y == -1)
    {
      curve_grid.data(p_1.y-1, p_1.x).curve_diag = c;
    }
    else if (d.x == 1)
    {
      curve_grid.data(p_1).curve_horiz = c;
    }
    else if (d.y == 1)
    {
      curve_grid.data(p_1).curve_vert = c;
    }
  }
}

void fill_smooth_curves(DataGrid<CurveData>& curve_grid)
{
  for (int y = 0; y < curve_grid.rows(); ++y)
  {
    int x = 0;

    while (x < curve_grid.cols() - 1)
    {
      if (!curve_grid.data(y, x).curve_horiz.is_empty_curve())
      {
        ++x;
        continue;
      }

      int x_left = x;

      const BezierCurve curve_left = (x == 0 ? BezierCurve() : curve_grid.data(y, x-1).curve_horiz);

      std::vector<cv::Point2d> knots;
      do
      {
        knots.push_back(curve_grid(y, x));
        ++x;
      } while (x < curve_grid.cols()-1 && curve_grid.data(y, x).curve_horiz.is_empty_curve());
      knots.push_back(curve_grid(y, x));

      int x_right = x;

      const BezierCurve curve_right = (x == curve_grid.cols()-1 ? BezierCurve() : curve_grid.data(y, x).curve_horiz);

      std::vector<BezierCurve> curves = generate_smooth_curves(knots, curve_left, curve_right);

      for (int i = 0; i < static_cast<int>(curves.size()); ++i)
      {
        curve_grid.data(y, static_cast<int>(x-curves.size()+i)).curve_horiz = curves[i];
      }
    }
  }

  for (int x = 0; x < curve_grid.cols(); ++x)
  {
    int y = 0;

    while (y < curve_grid.rows()-1)
    {
      if (!curve_grid.data(y, x).curve_vert.is_empty_curve())
      {
        ++y;
        continue;
      }

      const BezierCurve curve_left = (y == 0 ? BezierCurve() : curve_grid.data(y-1, x).curve_vert);

      std::vector<cv::Point2d> knots;

      do
      {
        knots.push_back(curve_grid(y, x));
        ++y;
      } while (y < curve_grid.rows()-1 && curve_grid.data(y, x).curve_vert.is_empty_curve());
      knots.push_back(curve_grid(y, x));

      const BezierCurve curve_right = (y == curve_grid.rows()-1 ? BezierCurve() : curve_grid.data(y, x).curve_vert);

      std::vector<BezierCurve> curves = generate_smooth_curves(knots, curve_left, curve_right);

      for (int i = 0; i < static_cast<int>(curves.size()); ++i)
      {
        curve_grid.data(static_cast<int>(y-curves.size()+i), x).curve_vert = curves[i];
      }
    }
  }
}

std::vector<PatchRegion> generate_patches_from_curve_grid(int target_index, const DataGrid<CurveData>& curve_grid)
{
  DataGrid<PatchRegion> reconstruction_regions_grid(curve_grid.rows()-1, curve_grid.cols()-1);
  for (int y = 0; y < curve_grid.rows()-1; ++y)
  {
    for (int x = 0; x < curve_grid.cols()-1; ++x)
    {
      PatchRegion& region = reconstruction_regions_grid.data(y, x);

      reconstruction_regions_grid(y, x) = curve_grid(y, x);
      region = PatchRegion(target_index, cv::Point(x, y),
        curve_grid.data(y, x).curve_horiz, curve_grid.data(y+1, x).curve_horiz,
        curve_grid.data(y, x).curve_vert, curve_grid.data(y, x+1).curve_vert);

      if (!curve_grid.data(y, x).curve_diag.is_empty_curve())
      {
        region.set_sub_regions(curve_grid.data(y, x).curve_diag);
      }
    }
  }

  /*
  // FIXME: Consider boundary cases as well if necessary.
  for (int y = 1; y < reconstruction_regions_grid.rows(); ++y)
  {
  for (int x = 1; x < reconstruction_regions_grid.cols(); ++x)
  {
  PatchRegion::fix_fabricability(
  reconstruction_regions_grid.data(y-1, x),
  reconstruction_regions_grid.data(y, x),
  reconstruction_regions_grid.data(y, x-1),
  reconstruction_regions_grid.data(y-1, x-1));
  }
  }
  */

  std::vector<PatchRegion> reconstruction_regions;
  for (int y = 0; y < reconstruction_regions_grid.rows(); ++y)
  {
    for (int x = 0; x < reconstruction_regions_grid.cols(); ++x)
    {
      if (reconstruction_regions_grid.data(y, x).valid())
      {
        reconstruction_regions.push_back(reconstruction_regions_grid.data(y, x));
      }
    }
  }

  return reconstruction_regions;
}

std::vector<PatchRegion> generate_patches(int target_index, const Grid& morphed_grid, cv::Mat edge_image)
{
  if (morphed_grid.empty() || !edge_image.data)
  {
    return std::vector<PatchRegion>();
  }

  const int grid_rows = morphed_grid.rows();
  const int grid_cols = morphed_grid.cols();

  DataGrid<Vector<cv::Point2d>> edge_grid = extract_edge_grid(morphed_grid, edge_image);

  DataGrid<CurveData> curve_grid(morphed_grid);
  for (int y = 1; y < grid_rows; ++y)
  {
    for (int x = 0; x < grid_cols-1; ++x)
    {
      fit_single_edge(curve_grid, edge_grid, cv::Point(x, y), cv::Point(x+1, y));
    }
  }

  for (int y = 0; y < grid_rows-1; ++y)
  {
    for (int x = 1; x < grid_cols; ++x)
    {
      fit_single_edge(curve_grid, edge_grid, cv::Point(x, y), cv::Point(x, y+1));
    }
  }

  for (int y = 0; y < grid_rows-1; ++y)
  {
    for (int x = 0; x < grid_cols-1; ++x)
    {
      fit_single_edge(curve_grid, edge_grid, cv::Point(x, y+1), cv::Point(x+1, y));
    }
  }

  for (int y = 0; y < grid_rows-1; ++y)
  {
    for (int x = 0; x < grid_cols-1; ++x)
    {
      fit_single_edge(curve_grid, edge_grid, cv::Point(x, y), cv::Point(x+1, y+1));
    }
  }

  fill_smooth_curves(curve_grid);
  std::vector<PatchRegion> reconstruction_regions = generate_patches_from_curve_grid(target_index, curve_grid);

  /*
  cv::Mat image_debug = edge_image.clone();
  image_debug *= 0.25;
  cv::resize(image_debug, image_debug, cv::Size(), 4.0, 4.0, cv::INTER_NEAREST);
  cv::merge(std::vector<cv::Mat>(3, image_debug), image_debug);
  for (const PatchRegion& region : reconstruction_regions)
  {
    region.draw(image_debug, cv::Vec3b(255, 0, 255), 4.0);
  }
  cv::imshow("image_debug_grid", image_debug);
  cv::waitKey(0);
  std::exit(EXIT_SUCCESS);
  */

  return reconstruction_regions;
}

boost::property_tree::ptree CurveData::save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const
{
  boost::property_tree::ptree tree;
  serialize(tree, "curve_horiz", curve_horiz, base_path, path);
  serialize(tree, "curve_vert", curve_vert, base_path, path);
  serialize(tree, "curve_diag", curve_diag, base_path, path);
  return tree;
}

void CurveData::load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree)
{
  deserialize(tree, "curve_horiz", curve_horiz, base_path);
  deserialize(tree, "curve_vert", curve_vert, base_path);
  deserialize(tree, "curve_diag", curve_diag, base_path);
}
