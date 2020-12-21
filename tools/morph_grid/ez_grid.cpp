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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "ez_grid.hpp"
#include "timer.hpp"

namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

EZGrid::EZGrid(cv::Mat image, cv::Mat image_filtered, int grid_size, const boost::filesystem::path &path_out) : m_window_name("ezGrid 0.1"),
                                                                                                                m_gui_name("Parameters"),
                                                                                                                m_path_out(path_out)
{
  cv::namedWindow(m_window_name, cv::WINDOW_AUTOSIZE);
  cv::namedWindow(m_gui_name, cv::WINDOW_AUTOSIZE);

  m_bilateral_filter_future = std::make_unique<BilateralFilterFuture>(m_gui_name, image);
  m_canny_bilateral_future = std::make_unique<CannyFuture>(m_gui_name, "cannyBilateral");
  m_canny_filtered_future = std::make_unique<CannyFuture>(m_gui_name, "cannyFiltered");
  m_morph_grid_future = std::make_unique<MorphGridFuture>(m_window_name, image, grid_size);

  m_canny_filtered_future->set_image(image_filtered);
  m_morph_grid_future->set_image_filtered(image_filtered);
}

void EZGrid::run()
{
  int key = -1;
  Timer<boost::milli> t;
  int remaining = 1;

  do
  {
    /*
    * Keyboard callback.
    */
    key = cv::waitKey(remaining);
    t.restart();
    m_morph_grid_future->key_press(key);

    /*
    *  Bilateral filter loop.
    */
    m_bilateral_filter_future->loop();
    if (m_bilateral_filter_future->ready())
    {
      cv::Mat image_bilateral = m_bilateral_filter_future->filtered_image();
      m_canny_bilateral_future->set_image(image_bilateral);
      m_morph_grid_future->set_image_bilateral(image_bilateral);
    }

    /*
    * Canny edge detection loop.
    */
    m_canny_bilateral_future->loop();
    if (m_canny_bilateral_future->ready())
    {
      cv::Mat edge_image = m_canny_bilateral_future->edge_image();
      m_morph_grid_future->set_edge_bilateral(edge_image);
    }

    m_canny_filtered_future->loop();
    if (m_canny_filtered_future->ready())
    {
      cv::Mat edge_image = m_canny_filtered_future->edge_image();
      m_morph_grid_future->set_edge_filtered(edge_image);
    }

    /*
    * Morph grid loop.
    */
    m_morph_grid_future->loop();

    /*
     * Save state if necessary.
     */
    if (m_morph_grid_future->save_state())
    {
      std::cout << "Saving state. This can take a while..." << std::endl;
      save();
      std::cout << "Done." << std::endl;
    }

    /*
    * Draw output image.
    */
    cv::Mat image_out = m_morph_grid_future->draw();
    cv::imshow(m_window_name, image_out);

    remaining = std::max(static_cast<int>(16.0f - t.duration().count()), 1);
  } while (key != ' ');
}

void EZGrid::load_partial_state(const boost::filesystem::path &path)
{
  const boost::filesystem::path base_path = path.parent_path();

  pt::ptree tree;
  pt::read_json(path.string(), tree);
  pt::ptree tree_grid = tree.get_child("ez_grid");

  deserialize(tree_grid, "bilateral_filter", *m_bilateral_filter_future, base_path);
  deserialize(tree_grid, "canny_bilateral", *m_canny_bilateral_future, base_path);
  deserialize(tree_grid, "canny_filtered", *m_canny_filtered_future, base_path);
  m_morph_grid_future->load_input_partial(base_path, tree_grid.get_child("morph_grid"));
  deserialize(tree_grid, "window_name", m_window_name, base_path);
}

void EZGrid::load_full_state(const boost::filesystem::path &path)
{
  pt::ptree tree;
  pt::read_json(path.string(), tree);
  load(path.parent_path(), tree.get_child("ez_grid"));
}

void EZGrid::save()
{
  const std::vector<PatchRegion> patches = m_morph_grid_future->patches();
  const DataGrid<unsigned char> grid = m_morph_grid_future->grid();

  pt::ptree root;
  pt::ptree tree_patches;
  for (const PatchRegion &p : patches)
  {
    tree_patches.push_back(std::make_pair("", p.save(m_path_out.parent_path(), "patches")));
  }
  root.add_child("patches", tree_patches);
  root.add_child("grid", grid.save(m_path_out.parent_path(), "grid"));
  root.add_child("ez_grid", save(m_path_out.parent_path(), "ez_grid"));
  pt::write_json(m_path_out.string(), root);
}

boost::property_tree::ptree EZGrid::save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const
{
  boost::property_tree::ptree tree;
  serialize(tree, "bilateral_filter", *m_bilateral_filter_future, base_path, path);
  serialize(tree, "canny_bilateral", *m_canny_bilateral_future, base_path, path);
  serialize(tree, "canny_filtered", *m_canny_filtered_future, base_path, path);
  serialize(tree, "morph_grid", *m_morph_grid_future, base_path, path);
  serialize(tree, "window_name", m_window_name, base_path, path);
  return tree;
}

void EZGrid::load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree)
{
  deserialize(tree, "bilateral_filter", *m_bilateral_filter_future, base_path);
  deserialize(tree, "canny_bilateral", *m_canny_bilateral_future, base_path);
  deserialize(tree, "canny_filtered", *m_canny_filtered_future, base_path);
  deserialize(tree, "morph_grid", *m_morph_grid_future, base_path);
  deserialize(tree, "window_name", m_window_name, base_path);
}
