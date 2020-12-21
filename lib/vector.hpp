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

#ifndef TRLIB_VECTOR_HPP_
#define TRLIB_VECTOR_HPP_

#include <vector>

#include "serializable.hpp"

template <typename T>
class Vector : public std::vector<T>, public Serializable
{
public:
  virtual boost::property_tree::ptree save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const override
  {
    boost::property_tree::ptree tree;
    serialize(tree, "data", this->begin(), this->end(), base_path, path);
    return tree;
  }

  virtual void load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree) override
  {
    std::vector<T> data;
    deserialize(tree, "data", data, base_path);
    this->swap(data);
  }
};

#endif /* TRLIB_VECTOR_HPP_ */
