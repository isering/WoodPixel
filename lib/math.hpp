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

#ifndef TRLIB_MATH_HPP_
#define TRLIB_MATH_HPP_

namespace math
{

  inline int factorial(int n)
  {
    int result = 1;
    for (int i = 1; i <= n; ++i)
    {
      result *= i;
    }
    return result;
  }

  inline int n_choose_k(int n, int k)
  {
    return factorial(n) / (factorial(k) * factorial(n - k));
  }

}; // namespace math

#endif /* TRLIB_MATH_HPP_ */