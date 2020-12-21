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

#ifndef TRLIB_GH_INVARIANT_HPP_
#define TRLIB_GH_INVARIANT_HPP_

#include <cmath>
#include <map>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <opencv2/opencv.hpp>

class GHBasis
{
public:
  GHBasis(int degree, double size_mm) : m_degree(degree),
                                        m_size_mm(size_mm)
  {
    /*
     * Construct (unnormalized) Hermite polynomials
     */
    hermite_polynomials.push_back(boost::math::tools::polynomial<double>{{1.0}});
    hermite_polynomials.push_back(boost::math::tools::polynomial<double>{{0.0, 2.0}});

    for (int i = 2; i <= degree; ++i)
    {
      hermite_polynomials.push_back(hermite_polynomials[1] * hermite_polynomials[i - 1] - 2.0 * (i - 1.0) * hermite_polynomials[i - 2]);
    }
  }

  double hermite_polynomial(int p, double x)
  {
    return boost::math::tools::evaluate_polynomial(hermite_polynomials[p].data().data(), x, p + 1);
  }

  double hermite_polynomial_normalized(int p, double x, double sigma)
  {
    return hermite_polynomial(p, x / sigma) * std::exp(-0.5 * x * x / (sigma * sigma));
  }

  /*
  double moment(int p, int q, cv::Mat im, int x, int y)
  {

  }
  */

private:
  std::map<std::pair<int, int>, boost::math::tools::polynomial<double>> polynomials;
  std::vector<boost::math::tools::polynomial<double>> hermite_polynomials;

  int m_degree;
  double m_size_mm;
};

#endif /* TRLIB_GH_INVARIANT_HPP_ */