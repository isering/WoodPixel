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

#ifndef TRLIB_OPENCV_EXTRA_HPP_
#define TRLIB_OPENCV_EXTRA_HPP_

#include <opencv2/opencv.hpp>

namespace cvx
{

  void addWithMask(cv::InputArray _src1, cv::InputArray _src2, cv::OutputArray _dst, cv::InputArray _mask)
  {
    //CV_Assert(_src1.sameSize(_src2) && _src1.type() == _src2.type());
    CV_Assert(_src1.sameSize(_mask));
    cv::Mat mask = _mask.getMat();
    CV_Assert(mask.type() == CV_8U);
    cv::Mat dst(_src1.size(), _src1.type());
    cv::Mat src1 = _src1.getMat(), src2 = _src2.getMat();

    cv::Mat mask_neg;
    cv::add(src1, src2, dst);
    cv::bitwise_not(mask, mask_neg);
    src1.copyTo(dst, mask_neg);

    dst.copyTo(_dst);
  }

  void multiplywithMask(cv::InputArray _src1, cv::InputArray _src2, cv::OutputArray _dst, cv::InputArray _mask)
  {
    //CV_Assert(_src1.sameSize(_src2) && _src1.type() == _src2.type());
    CV_Assert(_src1.sameSize(_mask));
    cv::Mat mask = _mask.getMat();
    CV_Assert(mask.type() == CV_8U);
    cv::Mat dst(_src1.size(), _src1.type());
    cv::Mat src1 = _src1.getMat(), src2 = _src2.getMat();

    cv::Mat mask_neg;
    cv::multiply(src1, src2, dst);
    cv::bitwise_not(mask, mask_neg);
    src1.copyTo(dst, mask_neg);

    dst.copyTo(_dst);
  }

} // namespace cvx

#endif /* TRLIB_OPENCV_EXTRA_HPP_ */