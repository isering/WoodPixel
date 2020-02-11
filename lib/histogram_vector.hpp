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

#ifndef TRLIB_HISTOGRAM_VECTOR_HPP_
#define TRLIB_HISTOGRAM_VECTOR_HPP_

#include "feature_vector.hpp"

class HistogramVector : public FeatureVector
{
public:
  HistogramVector() = default;
  HistogramVector(const std::vector<cv::Mat>& histogram_vector, const std::vector<cv::Mat>& circular_rotation) :
    FeatureVector(histogram_vector),
    m_circular_rotation(circular_rotation)
  {
  }

protected:
  std::vector<cv::Mat> m_circular_rotation;
};

#endif /* TRLIB_HISTOGRAM_VECTOR_HPP_ */