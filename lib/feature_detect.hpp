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

#ifndef TRLIB_FEATURE_DETECT_HPP_
#define TRLIB_FEATURE_DETECT_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> feature_detect(cv::Mat texture, int max_corners=1500, double quality_level=0.01, double min_distance=25.0);

void draw_features(const std::string& window_name, cv::Mat texture, const std::vector<cv::Point2f>& features);

cv::Mat features_to_binary_image(const std::vector<cv::Point2f>& features, int height, int width);

#endif /* TRLIB_FEATURE_DETECT_HPP_ */