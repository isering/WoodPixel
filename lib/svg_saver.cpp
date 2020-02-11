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

#include "svg_saver.hpp"

namespace fs = boost::filesystem;

void SVGSaver::save(const boost::filesystem::path& filename, const cv::Size2d& table_dimensions_mm) const
{
  std::stringstream out;

  write_header(out, table_dimensions_mm.width, table_dimensions_mm.height);

  out << "  <g id=\"Layer_x0020_1\">\n"
    << "    <metadata id=\"CorelCorpID_0Corel-Layer\"/>\n"
    << "    <rect class=\"str1\" x=\"0\" y=\"0\" width=\"" << svg_precision * table_dimensions_mm.width << "\" height=\"" << svg_precision * table_dimensions_mm.height << "\" />" << std::endl;
    for (std::shared_ptr<SaverPrimitive> p : m_debug_primitives)
    {
      out << "    " << *p << "\n";
    }
  out << "  </g>\n";


  out << "  <g id=\"Layer_x0020_2\">\n"
    << "    <metadata id=\"CorelCorpID_1Corel-Layer\"/>\n";

  for (std::shared_ptr<SaverPrimitive> p : m_primitives)
  {
    out << "    " << *p << "\n";
  }

  out << "  </g>\n"
    << "  <g id=\"Layer_x0020_3\">\n"
    << "    <metadata id=\"CorelCorpID_2Corel-Layer\"/>\n";

  for (std::shared_ptr<SaverPrimitive> p : m_text)
    {
      SVGText text = *reinterpret_cast<const SVGText*>(p.get());
      out << "  " << text << "\n";
    }

  out << "  </g>\n";
  
  write_footer(out);

  if (!fs::exists(filename.parent_path()))
  {
    fs::create_directories(filename.parent_path());
  }

  std::ofstream f_out(filename.string(), std::ios::trunc);
  f_out << out.rdbuf();
}

void SVGSaver::write_header(std::ostream& out, double size_x_mm, double size_y_mm) const
{
  const int size_x_scaled = static_cast<int>(svg_precision * size_x_mm);
  const int size_y_scaled = static_cast<int>(svg_precision * size_y_mm);
  const double stroke_width_hairline = 1;
  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
    << "<svg xmlns=\"http://www.w3.org/2000/svg\"\n"
    << "     xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
    << "     version=\"1.1\" baseProfile=\"full\"\n"
    << "     width=\"" << size_x_mm << "mm\" height=\"" << size_y_mm << "mm\"\n"
    << "     viewBox=\"0 0 " << size_x_scaled << " " << size_y_scaled << "\">\n"
    << "  <defs>\n"
    << "    <style type=\"text/css\">\n"
    << "      <![CDATA[\n"
    << "        .str0 {stroke:red; stroke-width:" << stroke_width_hairline << "}\n"
    << "        .str1 {stroke:green; stroke-width:" << stroke_width_hairline << "}\n"
    << "        .str2 {stroke:black; stroke-width:" << stroke_width_hairline << "}\n"
    << "        .fil0 {fill:none}\n"
    << "        .fil1 {fill:black}\n"
    << "        .fil2 {fill:green}\n"
    << "        .fnt0 {text-anchor:'middle'; dominant-baseline:'central'; font-weight:normal; font-size:3pt; font-family:'Arial'; text-decoration:underline;}\n"
    << "      ]]>\n"
    << "    </style>\n"
    << "  </defs>\n";
}

void SVGSaver::write_footer(std::ostream& out) const
{
  out << "</svg>\n";
}
