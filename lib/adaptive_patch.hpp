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

#ifndef TRLIB_ADAPTIVE_PATCH_HPP_
#define TRLIB_ADAPTIVE_PATCH_HPP_

#include <vector>

#include "patch.hpp"

struct AdaptivePatch
{
  AdaptivePatch(const std::vector<Patch> &patches, int level) : patches(patches),
                                                                level(level)
  {
  }

  double cost() const
  {
    double cost_val = 0.0;
    for (const Patch &p : patches)
    {
      cost_val += p.cost();
    }
    return cost_val;
  }

  int num_pixels() const
  {
    int n = 0;
    for (const Patch &p : patches)
    {
      n += p.size().area();
    }
    return n;
  }

  int target_index() const
  {
    return patches.empty() ? -1 : patches[0].region_target.target_index();
  }

  cv::Point coordinate() const
  {
    return patches.empty() ? cv::Point(-1, -1) : patches[0].region_target.coordinate();
  }

  std::vector<Patch> patches;
  int level;
};

#endif /* TRLIB_ADAPTIVE_PATCH_HPP_ */