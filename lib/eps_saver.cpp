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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/filesystem.hpp>

#include "eps_saver.hpp"

namespace fs = boost::filesystem;

void EPSSaver::save(const boost::filesystem::path& filename, const cv::Size2d& table_dimensions_mm) const
{
  const double box_x = table_dimensions_mm.width * 72.0 / 25.4;
  const double box_y = table_dimensions_mm.height * 72.0 / 25.4;

  std::stringstream out;
  out << "%!PS-Adobe-3.1 EPSF-3.0\n"
    << "%%BoundingBox: 0 0 " << " " << static_cast<int>(std::ceil(box_x)) <<" "<<  static_cast<int>(std::floor(box_y)) << "\n"
    << "%%PageBoundingBox: 0 0 " << " " << static_cast<int>(std::ceil(box_x)) <<" "<< static_cast<int>(std::floor(box_y)) << "\n"
    << "%%HiResBoundingBox: 0 0 " << " " << box_x <<" "<<  box_y << "\n"
    << "%%CropBox: 0 0 "<< " " << box_x <<" "<<  box_y  << "\n"
    << "0 " << box_y << " translate\n\n"
    << "72 25.4 div dup scale\n"
    << "/Helvetica findfont\n"
    << "3 scalefont setfont %3mm font size?\n\n"
    << "/setstr0 {0.0762 setlinewidth 1 0 0 setrgbcolor} def\n"
    << "/setstr1 {0.0762 setlinewidth 0 0 1 setrgbcolor} def\n"
    << "/fill1 {0 0 0 setrgbcolor} def\n"
    << "/fill2 {0 1 0 setrgbcolor} def\n\n"
    << "/setbeziercolor {1 0 0 setrgbcolor} def\n"
    << "/setcirclecolor {0 0 0 setrgbcolor} def\n"
    << "/setpolygoncolor {0 1 0 setrgbcolor} def\n"
    << "/setlinecolor {0 1 0 setrgbcolor} def\n\n"    
    << "/settextcolor {0 0 1 setrgbcolor} def\n\n"    
    << "/centertext { % x y string\n"
    << "/string exch def\n"
    << "/y exch def\n"
    << "/x exch def\n"
    << "0 0 moveto\n"
    << "string dup stringwidth pop\n"
    << "2 div\n"
    << "gsave 0.7 1 scale\n"
    << "x exch sub\n"
    << "y moveto\n"
    << "false charpath stroke grestore\n"
    << "} def \n"
    << "setstr0 \n";

  out << " " << table_dimensions_mm.width << " 0 translate -1 1 scale \n";
  for (std::shared_ptr<SaverPrimitive> p : m_primitives)
  {
    out << *p << "\n";
  }

  out << " " << table_dimensions_mm.width << " 0 translate -1 1 scale \n";
  for (std::shared_ptr<SaverPrimitive> p : m_text)
  {
    EPSText text = *reinterpret_cast<const EPSText*>(p.get());
    text.pos.x = cols - text.pos.x;
    out << *p << "\n";
  }

  if (!fs::exists(filename.parent_path()))
  {
    fs::create_directories(filename.parent_path());
  }

  std::ofstream f_out(filename.string(), std::ios::trunc);
  f_out << out.rdbuf();
}