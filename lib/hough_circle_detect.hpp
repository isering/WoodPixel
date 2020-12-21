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

#ifndef TRLIB_HOUGH_CIRCLE_DETECT_HPP_
#define TRLIB_HOUGH_CIRCLE_DETECT_HPP_

#include <future>
#include <vector>

#include <opencv2/opencv.hpp>

class HoughCircleDetect
{
public:
  HoughCircleDetect(double scale) : m_scale(scale)
  {
  }

  std::vector<cv::Vec3f> compute(const cv::Mat &image, double param_1, double param_2, double dpi, double circle_size_mm);
  std::future<std::vector<cv::Vec3f>> compute_background(const cv::Mat &image, double param_1, double param_2, double dpi, double circle_size_mm);
  std::vector<cv::Vec3f> compute_interactive(const std::vector<cv::Vec3f> &segmentation_data_in, const cv::Mat &image, const cv::Mat &image_background, double dpi, double circle_size_mm, std::string window_name);

  static void draw_to_file(const std::string &fname, const cv::Mat &im, const std::vector<cv::Vec3f> &seg);

private:
  void draw();

  static void on_mouse(int event, int x, int y, int, void *ptr);

  std::future<std::vector<cv::Vec3f>> m_future;

  std::vector<cv::Vec3f> m_data, m_data_in;

  cv::Mat m_image_background;
  std::string m_window_name;

  int m_param_1, m_param_2;
  int m_param_old_1, m_param_old_2;

  double m_scale;

  static double raw_to_param_1(int val);
  static double raw_to_param_2(int val);
};

#endif /* TRLIB_HOUGH_CIRCLE_DETECT_HPP_ */