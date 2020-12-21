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

#ifndef TRLIB_TEXTURE_MARKER_HPP_
#define TRLIB_TEXTURE_MARKER_HPP_

#include <vector>

#include <opencv2/opencv.hpp>

#include "serializable.hpp"

class TextureMarker : public Serializable
{
public:
  TextureMarker scaled(double scale) const
  {
    TextureMarker marker_return = *this;
    for (cv::Point2d &p : marker_return.markers_pix)
    {
      p *= scale;
    }
    return marker_return;
  }

  std::vector<cv::Point2d> markers_pix;

  virtual void load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const override;
};

#endif /* TRLIB_TEXTURE_MARKER_HPP_ */
