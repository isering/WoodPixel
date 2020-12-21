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

#ifndef TRLIB_HERSHEY_FONT_HPP_
#define TRLIB_HERSHEY_FONT_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

class HersheyFont
{
public:
  HersheyFont(cv::Point2f pos, const std::string &str, double rotation_deg, double scale, double stretch_x = 1.0) : pos(pos),
                                                                                                                    str(str),
                                                                                                                    rotation_deg(rotation_deg),
                                                                                                                    scale(scale),
                                                                                                                    stretch_x(stretch_x)
  {
  }

  void serialize_ps(std::ostream &os) const;
  void serialize_svg(std::ostream &os, double svg_precision) const;

private:
  static int char_width(const char c);
  static cv::Point2i get_vertex(int char_index, int vertex_index);
  static std::vector<std::pair<cv::Point2i, bool>> get_char_path(const char c);

  cv::Point2f pos;
  std::string str;
  double rotation_deg;
  double scale;
  double stretch_x;

  static const int num_chars;
  static const int first_char;
  static const double base_scale;
  static short simplex[95][112];
};

#endif /* TRLIB_HERSHEY_FONT_HPP_ */