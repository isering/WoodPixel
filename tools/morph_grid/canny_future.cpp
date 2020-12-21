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

#include "canny_future.hpp"

boost::property_tree::ptree CannyFuture::save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const
{
  boost::property_tree::ptree tree;
  serialize_image(tree, "image", m_image, base_path, path);
  serialize_image(tree, "edge_image", m_edge_image, base_path, path);
  serialize(tree, "canny_threshold_1", m_canny_threshold_1, base_path, path);
  serialize(tree, "canny_threshold_2", m_canny_threshold_2, base_path, path);
  return tree;
}

void CannyFuture::load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree)
{
  deserialize_image(tree, "image", m_image, base_path);
  deserialize_image(tree, "edge_image", m_edge_image, base_path);
  deserialize(tree, "canny_threshold_1", m_canny_threshold_1, base_path);
  deserialize(tree, "canny_threshold_2", m_canny_threshold_2, base_path);
}
