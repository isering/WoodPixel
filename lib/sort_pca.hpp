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

#ifndef TRLIB_SORT_PCA_HPP_
#define TRLIB_SORT_PCA_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

template <typename T>
void sort_pca(std::vector<cv::Point_<T>> &points)
{
  struct sort_pca_data
  {
    double comp_1;
    double comp_2;
    int index;
  };

  struct sort_pca_comp
  {
    bool operator()(const sort_pca_data &lhs, const sort_pca_data &rhs)
    {
      return (lhs.comp_1 < rhs.comp_1) || (lhs.comp_1 == rhs.comp_1 && lhs.comp_2 < rhs.comp_2);
    }
  };

  const int num_points = static_cast<int>(points.size());

  // Copy PCA data to matrix.
  cv::Mat pca_data(num_points, 2, CV_64FC1);
  for (int i = 0; i < num_points; ++i)
  {
    double *ptr = reinterpret_cast<double *>(pca_data.ptr(i));
    ptr[0] = points[i].x;
    ptr[1] = points[i].y;
  }

  // Apply PCA.
  cv::PCA pca(pca_data, cv::noArray(), cv::PCA::DATA_AS_ROW, 2);
  cv::Mat projected = pca.project(pca_data);

  // Sort points according to PCA.
  std::vector<sort_pca_data> sort_data(points.size());
  for (int i = 0; i < num_points; ++i)
  {
    const double *ptr = reinterpret_cast<const double *>(projected.ptr(i));
    sort_data[i].comp_1 = ptr[0];
    sort_data[i].comp_2 = ptr[1];
    sort_data[i].index = i;
  }
  std::sort(sort_data.begin(), sort_data.end(), sort_pca_comp());

  // Sort original points.
  std::vector<cv::Point_<T>> points_sorted(points.size());
  for (int i = 0; i < num_points; ++i)
  {
    points_sorted[i] = points[sort_data[i].index];
  }
  points.swap(points_sorted);
}

#endif /* TRLIB_SORT_PCA_HPP_ */