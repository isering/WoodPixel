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

#ifndef TRLIB_BEZIER_RANSAC_HPP_
#define TRLIB_BEZIER_RANSAC_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

#include "bezier_curve.hpp"

BezierCurve fit_bezier_cubic_ransac(const std::vector<cv::Point2d> &points, cv::Rect region, int num_estimate, int max_iterations);
BezierCurve fit_bezier_cubic_ransac(const cv::Point2d &p1, const cv::Point2d &p2, std::vector<cv::Point2d> points, int num_estimate, int max_iterations, double *inlier_ratio = 0);

#endif /* TRLIB_BEZIER_RANSAC_HPP_ */