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

#ifndef MORPH_GRID_FUTURE_HPP_
#define MORPH_GRID_FUTURE_HPP_

#include <deque>
#include <future>

#include "data_grid.hpp"
#include "edge.hpp"
#include "future_base.hpp"
#include "generate_patches.hpp"
#include "grid.hpp"
#include "patch_region.hpp"
#include "serializable.hpp"

class MorphGridFuture : public FutureBase, public Serializable
{
public:
  MorphGridFuture(const std::string &window_name, cv::Mat image, int grid_size) : FutureBase(window_name),
                                                                                  m_image(image),
                                                                                  m_output_image(OutputImage::BILATERAL),
                                                                                  m_output_edge(OutputEdge::EDGE_DRAW),
                                                                                  m_mode(Mode::EDIT_GRID),
                                                                                  m_mouse_pos(-1, -1),
                                                                                  m_mouse_tool_size(15.0),
                                                                                  m_mouse_l_pressed(false),
                                                                                  m_mouse_shift_pressed(false),
                                                                                  m_grid_size(grid_size),
                                                                                  m_is_grid_point_selected(false),
                                                                                  m_edge_select_last(-1, -1),
                                                                                  m_is_grid_updated(false),
                                                                                  m_update_patches(false),
                                                                                  m_draw_grid(true),
                                                                                  m_visibility_mode(false),
                                                                                  m_save_state(false)
  {
    if (m_image.channels() == 1)
    {
      cv::cvtColor(m_image, m_image, cv::COLOR_GRAY2BGR);
    }

    m_grid = DataGrid<unsigned char>(image.rows / grid_size + 1, image.cols / grid_size + 1, 1);
    m_grid.scale(grid_size);
    m_edge_bilateral = cv::Mat::zeros(image.size(), CV_8UC1);
    m_edge_filtered = cv::Mat::zeros(image.size(), CV_8UC1);
    m_edge_drawn = cv::Mat::zeros(image.size(), CV_8UC1);
    m_mask_bilateral = cv::Mat::zeros(image.size(), CV_8UC1);
    m_mask_filtered = cv::Mat::zeros(image.size(), CV_8UC1);
    m_grid_density = cv::Mat::zeros(image.size(), CV_32FC1);
    m_distance = cv::Mat::zeros(image.size(), CV_32FC1);
    m_distance_grad_x = cv::Mat::zeros(image.size(), CV_32FC1);
    m_distance_grad_y = cv::Mat::zeros(image.size(), CV_32FC1);
  }

  enum class OutputImage
  {
    NONE,
    IMAGE,
    BILATERAL,
    FILTERED,
    MASK_FILTERED,
    MASK_BILATERAL,
    DISTANCE,
    DISTANCE_GRAD,
    DENSITY
  };

  enum class OutputEdge
  {
    NONE,
    EDGE,
    EDGE_DRAW,
    FILTERED,
    BILATERAL,
    DRAWN,
    PATCHES,
    DRRAW_GRID
  };

  enum class Mode
  {
    NONE,
    EDIT_GRID,
    EDGE_SELECT,
    EDIT_FILTERED_EDGES,
    EDIT_BILATERAL_EDGES,
    EDIT_DRAW_EDGES,
    EDIT_MAKE_DENSER,
    EDIT_MAKE_NORMAL,
    EDIT_MAKE_WIDER
  };

  struct MorphGridData
  {
    DataGrid<unsigned char> grid;
    cv::Mat im_distance;
    cv::Mat im_label;
    cv::Mat gradient_distance_x, gradient_distance_y;
  };

  struct PatchData
  {
    std::vector<PatchRegion> patches;
    std::vector<BezierCurve> curves_fitted;
    DataGrid<CurveData> curve_grid;
  };

  virtual void loop() override;

  cv::Mat edge_image() const;

  cv::Mat draw() const;
  void draw_edge_image(cv::Mat image) const;

  DataGrid<unsigned char> relax_grid(DataGrid<unsigned char> grid, int num_iterations, cv::Mat im_dist, cv::Mat gradient_dist_x, cv::Mat gradient_dist_y, cv::Mat density_image, int grid_size) const;

  MorphGridData update_morphed_grid(DataGrid<unsigned char> grid, cv::Mat edge_image, cv::Mat density_image, int grid_size);
  PatchData update_patch_data(DataGrid<unsigned char> grid, cv::Mat edge_image, std::vector<Edge> selected_edges);

  void handle_key_press();

  void set_image_bilateral(cv::Mat image_bilateral)
  {
    m_image_bilateral = image_bilateral.clone();
  }

  void set_edge_bilateral(cv::Mat edge_bilateral)
  {
    m_edge_bilateral = edge_bilateral.clone();
    m_update = true;
  }

  void set_image_filtered(cv::Mat image_filtered)
  {
    m_image_filtered = image_filtered.clone();
  }

  void set_edge_filtered(cv::Mat edge_filtered)
  {
    m_edge_filtered = edge_filtered.clone();
    m_update = true;
  }

  void select_grid_point(int p_x, int p_y);
  void deselect_grid_point(int p_x, int p_y, bool shift_pressed);

  void handle_edge_select(int p_x, int p_y, bool button_l_pressed, bool delete_modifier);

  const DataGrid<unsigned char> &grid() const
  {
    return m_grid;
  }

  const std::vector<PatchRegion> &patches() const
  {
    return m_patches;
  }

  bool save_state()
  {
    bool save_state_last = m_save_state;
    m_save_state = false;
    return save_state_last;
  }

  virtual void load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree) override;
  virtual void load_input_partial(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree);
  virtual boost::property_tree::ptree save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const override;

protected:
  virtual void on_mouse(int event, int x, int y, int flags) override;

private:
  OutputImage m_output_image;
  OutputEdge m_output_edge;
  Mode m_mode;

  cv::Mat m_image;

  cv::Mat m_image_bilateral;
  cv::Mat m_image_filtered;

  cv::Mat m_mask_bilateral;
  cv::Mat m_mask_filtered;

  cv::Mat m_edge_bilateral;
  cv::Mat m_edge_filtered;
  cv::Mat m_edge_drawn;
  cv::Mat m_grid_density;

  cv::Mat m_distance, m_label;
  cv::Mat m_distance_grad_x, m_distance_grad_y;

  std::string m_message;

  std::vector<PatchRegion> m_patches;
  std::vector<BezierCurve> m_curves_fitted;
  DataGrid<CurveData> m_curve_grid;

  std::future<MorphGridData> m_morph_grid_data_future;
  std::future<PatchData> m_patch_data_future;

  /*
   * Mouse stuff.
   */
  cv::Point m_mouse_pos, m_mouse_pos_old;
  double m_mouse_tool_size;
  bool m_mouse_l_pressed, m_mouse_shift_pressed;

  DataGrid<unsigned char> m_grid;
  int m_grid_size;

  bool m_is_grid_point_selected;
  cv::Point m_selected_grid_point;

  std::vector<Edge> m_selected_edges;
  cv::Point m_edge_select_last;

  cv::Point get_nearest_grid_point(int p_x, int p_y);
  int get_nearest_selected_edge(int p_x, int p_y);

  void draw_new_edge(cv::Mat edge_image, cv::Point p_1, cv::Point p_2, double size_draw, double size_delete, int event);
  void draw_soft_blob(cv::Mat mask, cv::Point p_1, double target_value, double size, double amount, int event);

  bool m_is_grid_updated;
  bool m_update_patches;
  bool m_draw_grid;
  bool m_visibility_mode;
  bool m_save_state;
};

#endif /* MORPH_GRID_FUTURE_HPP_ */