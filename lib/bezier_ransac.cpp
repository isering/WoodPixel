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

#include "bezier_ransac.hpp"

static std::random_device random_device;

/*
* Samples a vector of random points.
*/
template <typename RNG>
static std::vector<cv::Point2d> get_random_points(const std::vector<cv::Point2d>& points, int num_points, RNG& generator)
{
  if (points.size() == num_points)
  {
    return points;
  }

  if (points.size() < num_points)
  {
    return std::vector<cv::Point2d>();
  }

  std::uniform_int_distribution<> dist(0, static_cast<int>(points.size()-1));
  std::vector<int> random_index_vec;
  while (random_index_vec.size() < num_points)
  {
    const int random_index = dist(generator);
    if (std::find(random_index_vec.begin(), random_index_vec.end(), random_index) == random_index_vec.end())
    {
      random_index_vec.push_back(random_index);
    }
  }

  std::vector<cv::Point2d> points_random(num_points);
  for (int i = 0; i < num_points; ++i)
  {
    points_random[i] = points[random_index_vec[i]];
  }

  return points_random;
}

template <typename RNG>
static std::vector<cv::Point2d> get_random_points_stable(const std::vector<cv::Point2d>& points, int num_points, RNG& generator)
{
  if (points.size() == num_points)
  {
    return points;
  }

  if (points.size() < num_points)
  {
    return std::vector<cv::Point2d>();
  }

  std::uniform_int_distribution<> dist(0, static_cast<int>(points.size()-1));
  std::vector<int> random_index_vec;
  while (random_index_vec.size() < num_points)
  {
    const int random_index = dist(generator);
    if (std::find(random_index_vec.begin(), random_index_vec.end(), random_index) == random_index_vec.end())
    {
      random_index_vec.push_back(random_index);
    }
  }

  std::sort(random_index_vec.begin(), random_index_vec.end());

  std::vector<cv::Point2d> points_random(num_points);
  for (int i = 0; i < num_points; ++i)
  {
    points_random[i] = points[random_index_vec[i]];
  }

  return points_random;
}

static int compute_num_inliers(const BezierCurve& curve, const std::vector<cv::Point2d>& points, double inlier_dist)
{
  //const double inlier_dist = std::sqrt(2);
//  const double inlier_dist = 1.0;
  int num_inliers = 0;

  for (const cv::Point2d& p : points)
  {
    if (curve.min_dist_approx(p) < inlier_dist)
    {
      ++num_inliers;
    }
  }

  return num_inliers;
}

static double compute_error(const BezierCurve& curve, const std::vector<cv::Point2d>& points)
{
  double error = 0.0;
  for (const cv::Point2d& p : points)
  {
    error += curve.min_dist_approx(p);
  }
  return error;
}

BezierCurve fit_bezier_cubic_ransac(const std::vector<cv::Point2d>& points, cv::Rect region, int num_estimate, int max_iterations)
{
  if (points.size() < num_estimate)
  {
    return BezierCurve();
  }

  BezierCurve curve_best;
  int num_inliers_max = 0;
  double error_min = std::numeric_limits<double>::max();

  std::default_random_engine engine(random_device());
  
  for (int i = 0; i < max_iterations; ++i)
  {
    const std::vector<cv::Point2d> points_random = get_random_points(points, num_estimate, random_device);

    BezierCurve curve = BezierCurve::fit_cubic(points_random);

    if (!curve.is_simple_curve_cubic(0.0))
    {
      continue;
    }

    cv::Mat mask = curve.compute_separating_curve_mask(region);

    /*
    cv::Mat mask_draw = mask.clone();
    mask_draw.setTo(64, mask_draw == 1);
    mask_draw.setTo(128, mask_draw == 2);
    cv::resize(mask_draw, mask_draw, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);
    cv::imshow("mask_draw", mask_draw);
    cv::waitKey();
    std::exit(EXIT_SUCCESS);
    */

    int num_pixels_left = cv::countNonZero(mask == 1);
    int num_pixels_right = cv::countNonZero(mask == 2);
    int num_pixels_curve = cv::countNonZero(mask == 255);

    const double ratio_left = static_cast<double>(num_pixels_left) / static_cast<double>(num_pixels_left + num_pixels_right + num_pixels_curve);
    const double ratio_right = static_cast<double>(num_pixels_left) / static_cast<double>(num_pixels_right + num_pixels_right + num_pixels_curve);
    if (ratio_left < 0.3 || ratio_right < 0.3)
    {
      continue;
    }

    const int num_inliers = compute_num_inliers(curve, points, 1.0);
    //const double error = compute_error(curve, points);

    if (num_inliers < num_inliers_max)
    {
      num_inliers_max = num_inliers;
      curve_best = curve;
    }
  }

  return curve_best;
}

struct proj_length_sort
{
  proj_length_sort(const cv::Point2d& c1, const cv::Point2d& c2) :
    c(c1),
    s(c2 - c1)
  {
    denom = s.dot(s);
  }

  bool operator()(const cv::Point2d& p1, const cv::Point2d& p2)
  {
    return proj_length(p1) < proj_length(p2);
  }

  double proj_length(const cv::Point2d& p)
  {
    return s.dot(p - c) / denom;
  }

  cv::Point2d c, s;
  double denom;
};

BezierCurve fit_bezier_cubic_ransac(const cv::Point2d& p1, const cv::Point2d& p2, std::vector<cv::Point2d> points, int num_estimate, int max_iterations, double* inlier_ratio)
{
  if (points.size() < num_estimate)
  {
    return BezierCurve();
  }

  std::sort(points.begin(), points.end(), proj_length_sort(p1, p2));

  /*
  std::copy(points.begin(), points.end(), std::ostream_iterator<cv::Point2d>(std::cout, "\n"));
  std::cout << p1 << " " << p2 << std::endl;
  */

  BezierCurve curve_best;
  int num_inliers_max = 0;

  std::default_random_engine engine(random_device());

  for (int i = 0; i < max_iterations; ++i)
  {
    const std::vector<cv::Point2d> points_random = get_random_points_stable(points, num_estimate, random_device);
    /*
    std::copy(points_random.begin(), points_random.end(), std::ostream_iterator<cv::Point2d>(std::cout, "\n"));
    std::cout << std::endl;
    */

    BezierCurve curve = BezierCurve::fit_cubic(p1, p2, points_random);
    //curve.draw_debug();

    if (!curve.is_simple_curve_cubic(0.0))
    {
      continue;
    }

    //std::cout << "is_simple" << std::endl;

    /*
    cv::Mat mask = curve.compute_separating_curve_mask(region);
    mask = mask(cv::Rect(mask.cols/8, mask.rows/8, mask.cols-mask.cols/4, mask.rows-mask.rows/4));

    cv::Mat mask_draw = mask.clone();
    mask_draw.setTo(64, mask_draw == 1);
    mask_draw.setTo(128, mask_draw == 2);
    cv::resize(mask_draw, mask_draw, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);
    cv::imshow("mask_draw", mask_draw);
    cv::waitKey();
    std::exit(EXIT_SUCCESS);
    
    int num_pixels_left = cv::countNonZero(mask == 1);
    int num_pixels_right = cv::countNonZero(mask == 2);
    int num_pixels_curve = cv::countNonZero(mask == 255);

    const double ratio_left = static_cast<double>(num_pixels_left) / static_cast<double>(num_pixels_left + num_pixels_right + num_pixels_curve);
    const double ratio_right = static_cast<double>(num_pixels_left) / static_cast<double>(num_pixels_right + num_pixels_right + num_pixels_curve);
    if (ratio_left < 0.3 || ratio_right < 0.3)
    {
      continue;
    }
    */

    /*
    const double len_1 = cv::norm(curve.control_point(0)-curve.control_point(1));
    const double len_2 = cv::norm(curve.control_point(2)-curve.control_point(3));
    std::cout << len_1 << " "  << len_2 << std::endl;
    
    if (len_1 >= 15.0 || len_2 >= 15.0)
    {
      continue;
    }
    */

    
    
    //const int num_inliers = compute_num_inliers(curve, points, std::sqrt(2.0));
    const int num_inliers = compute_num_inliers(curve, points, std::sqrt(2.0));
    if (num_inliers > num_inliers_max)
    {
      num_inliers_max = num_inliers;
      curve_best = curve;
    }
  }

  if (inlier_ratio != 0)
  {
    *inlier_ratio = static_cast<double>(num_inliers_max) / static_cast<double>(points.size());
  }

  std::vector<cv::Point2f> points_float(points.begin(), points.end());
  points_float.push_back(p1);
  points_float.push_back(p2);
  cv::Rect bbox = cv::boundingRect(points_float);
  bbox.width += 1;
  bbox.height += 1;

  /*
  cv::Mat im = cv::Mat::zeros(bbox.size(), CV_8UC3);
  for (const cv::Point2d& p : points)
  {
    const cv::Point p_int(p);
    if (p_int.x >= bbox.x && p_int.x < bbox.x+bbox.width && p_int.y >= bbox.y && p_int.y < bbox.y+bbox.height)
    {
      cv::Vec3b& c = im.at<cv::Vec3b>(p_int-bbox.tl());
       c = cv::Vec3b(255, 255, 255);
    }
  }

  std::cout << "num_points: " << points.size() << std::endl;
  std::cout << "inlier ratio: " << *inlier_ratio << std::endl;
  std::cout << "max_gradient: " << curve_best.max_gradient_norm_approx() << std::endl;
  std::cout << "max_curvature: " << curve_best.max_curvature_norm_approx() << std::endl;;

  im.at<cv::Vec3b>(cv::Point(p1)-bbox.tl()) = cv::Vec3b(255, 0, 255);
  im.at<cv::Vec3b>(cv::Point(p2)-bbox.tl()) = cv::Vec3b(255, 0, 255);

  cv::resize(im, im, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);

  curve_best.scaled(8.0).draw(im, cv::Vec3b(255, 255, 0), 8*bbox.tl());

  cv::imshow("RANSAC", im);
  cv::waitKey();
  */
  return curve_best;
}