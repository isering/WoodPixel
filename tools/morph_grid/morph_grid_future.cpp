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

#include <boost/format.hpp>

#include "generate_patches.hpp"
#include "line.hpp"
#include "morph_grid_future.hpp"

void MorphGridFuture::loop()
{
  handle_key_press();
  
  // DEBUG
  //m_update = true;

  if (m_update && !m_morph_grid_data_future.valid())
  {
    m_morph_grid_data_future = std::async(&MorphGridFuture::update_morphed_grid, this, m_grid, edge_image(), m_grid_density, m_grid_size);
    m_update = false;
  }

  if (m_morph_grid_data_future.valid() && m_morph_grid_data_future._Is_ready())
  {
    if (m_is_grid_updated)
    {
      m_morph_grid_data_future.get();
      m_is_grid_updated = false;
      m_update = true;
    }
    else
    {
      MorphGridData data = m_morph_grid_data_future.get();
      m_grid = data.grid;
      m_distance = data.im_distance;
      m_distance_grad_x = data.gradient_distance_x;
      m_distance_grad_y = data.gradient_distance_y;
      m_label = data.im_label;
      m_update_patches = true;
    }
  }

  if (m_update_patches && !m_patch_data_future.valid())
  {
    m_patch_data_future = std::async(&MorphGridFuture::update_patch_data, this, m_grid, edge_image(), m_selected_edges);
    m_update_patches = false;
  }

  if (m_patch_data_future.valid() && m_patch_data_future._Is_ready())
  {
    PatchData data = m_patch_data_future.get();
    m_patches = data.patches;
    m_curves_fitted = data.curves_fitted;
    m_curve_grid = data.curve_grid;

    /*
     * Check patch validity.
     */
    for (const PatchRegion& patch : m_patches)
    {
      if (patch.has_sub_regions())
      {
        for (const PatchRegion& sub_region : patch.sub_regions())
        {
          if (!sub_region.valid())
          {
            throw(std::runtime_error(
              (boost::format("MorphGridFuture: Encountered invalid subregion from patch %d %d %d")
              % sub_region.target_index() % sub_region.coordinate().x % sub_region.coordinate().y).str()));
          }
        }
      }
      else
      {
        if (!patch.valid())
        {
          throw(std::runtime_error(
            (boost::format("TreeMatch::load: Tried to load invalid patch %d %d %d")
            % patch.target_index() % patch.coordinate().x % patch.coordinate().y).str()));
        }
      }
    }
  }
}

cv::Mat MorphGridFuture::edge_image() const
{
  cv::Mat image = cv::Mat::zeros(m_image.size(), CV_8UC1);
  if (!m_edge_bilateral.empty() && !m_mask_bilateral.empty())
  {
    image.setTo(255, m_edge_bilateral & m_mask_bilateral);
  }
  if (!m_edge_filtered.empty() && !m_mask_filtered.empty())
  {
    image.setTo(255, m_edge_filtered & m_mask_filtered);
  }
  image.setTo(255, m_edge_drawn);
  return image;
}

void MorphGridFuture::draw_edge_image(cv::Mat image) const
{
  cv::Mat edge_mat;
  if (!m_edge_bilateral.empty())
  {
    if (m_visibility_mode)
    {
      cv::dilate(m_edge_bilateral, edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
      image.setTo(cv::Vec3b(22, 97, 184), edge_mat);
    }
    else
    {
      image.setTo(cv::Vec3b(22, 97, 184), m_edge_bilateral);
    }
  }
  if (!m_edge_filtered.empty())
  {
    if (m_visibility_mode)
    {
      cv::dilate(m_edge_filtered, edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
      image.setTo(cv::Vec3b(12, 204, 127), edge_mat);
    }
    else
    {
      image.setTo(cv::Vec3b(12, 204, 127), m_edge_filtered);
    }
  }

  if (m_visibility_mode)
  {
    cv::dilate(m_edge_drawn, edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
    image.setTo(cv::Vec3b(124, 153, 6), edge_mat);

    cv::dilate(edge_image(), edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
    image.setTo(cv::Vec3b(255, 255, 255), edge_mat);
  }
  else
  {
    image.setTo(cv::Vec3b(124, 153, 6), m_edge_drawn);
    image.setTo(cv::Vec3b(255, 255, 255), edge_image());
  }
}

cv::Mat MorphGridFuture::draw() const
{
  double max_val;

  cv::Mat image_in = m_image;
  switch (m_output_image)
  {
  case OutputImage::IMAGE:
    image_in = m_image;
    break;
  case OutputImage::BILATERAL:
    if (!m_image_bilateral.empty())
    {
      image_in = m_image_bilateral;
    }
    break;
  case OutputImage::FILTERED:
    if (!m_image_filtered.empty())
    {
      image_in = m_image_filtered;
    }
    break;
  case OutputImage::MASK_BILATERAL:
    image_in = m_mask_bilateral;
    break;
  case OutputImage::MASK_FILTERED:
    image_in = m_mask_filtered;
    break;
  case OutputImage::DENSITY:
    m_grid_density.convertTo(image_in, CV_8UC1, 127.0, 128.0);
    break;
  case OutputImage::DISTANCE:
    cv::minMaxLoc(m_distance, 0, &max_val);
    m_distance.convertTo(image_in, CV_8UC1, 255.0/max_val);
    cv::applyColorMap(image_in, image_in, cv::COLORMAP_PARULA);
    image_in *= 0.5;
    break;
  case OutputImage::DISTANCE_GRAD:
    cv::magnitude(m_distance_grad_x, m_distance_grad_y, image_in);
    cv::minMaxLoc(image_in, 0, &max_val);
    image_in.convertTo(image_in, CV_8UC1, 255.0/max_val);
    cv::applyColorMap(image_in, image_in, cv::COLORMAP_PARULA);
    image_in *= 0.5;
    break;
  }

  cv::Mat image_out = image_in.clone();
  if (image_out.channels() == 1)
  {
    cv::cvtColor(image_out, image_out, cv::COLOR_GRAY2BGR);
  }

  cv::Mat edge_mat;
  switch (m_output_edge)
  {
  case OutputEdge::EDGE:
    
    if (m_visibility_mode)
    {
      cv::dilate(edge_image(), edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
      image_out.setTo(cv::Vec3b(255, 255, 255), edge_mat);
    }
    else
    {
      image_out.setTo(cv::Vec3b(255, 255, 255), edge_image());
    }
    
    break;
  case OutputEdge::EDGE_DRAW:
    draw_edge_image(image_out);
    break;
  case OutputEdge::BILATERAL:
    if (!m_edge_bilateral.empty())
    {
      if (m_visibility_mode)
      {
        cv::dilate(m_edge_bilateral, edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
        image_out.setTo(cv::Vec3b(255, 255, 255), edge_mat);
      }
      else
      {
        image_out.setTo(cv::Vec3b(255, 255, 255), m_edge_bilateral);
      }
    }
    break;
  case OutputEdge::FILTERED:
    if (!m_edge_filtered.empty())
    {
      if (m_visibility_mode)
      {
        cv::dilate(m_edge_filtered, edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
        image_out.setTo(cv::Vec3b(255, 255, 255), edge_mat);
      }
      else
      {
        image_out.setTo(cv::Vec3b(255, 255, 255), m_edge_filtered);
      }
    }
    break;
  case OutputEdge::DRAWN:
    if (m_visibility_mode)
    {
      cv::dilate(m_edge_drawn, edge_mat, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
      image_out.setTo(cv::Vec3b(255, 255, 255), edge_mat);
    }
    else
    {
      image_out.setTo(cv::Vec3b(255, 255, 255), m_edge_drawn);
    }
    break;
  }

  if (m_draw_grid)
  {
    Grid grid_draw = m_grid;
    if (m_mode == Mode::EDIT_GRID && m_is_grid_point_selected)
    {
      grid_draw(m_selected_grid_point) = m_mouse_pos;
    }

    if (m_output_edge != OutputEdge::PATCHES)
    {
      for (int y = 0; y < m_grid.rows(); ++y)
      {
        for (int x = 1; x < m_grid.cols(); ++x)
        {
          cv::line(image_out, cv::Point(grid_draw(y, x-1)), cv::Point(grid_draw(y, x)), m_visibility_mode ? cv::Scalar(199, 155, 34) : cv::Scalar(200, 200, 200), m_visibility_mode ? 2 : 1);
        }
      }

      for (int y = 1; y < m_grid.rows(); ++y)
      {
        for (int x = 0; x < m_grid.cols(); ++x)
        {
          cv::line(image_out, cv::Point(grid_draw(y-1, x)), cv::Point(grid_draw(y, x)), m_visibility_mode ? cv::Scalar(199, 155, 34) : cv::Scalar(200, 200, 200), m_visibility_mode ? 2 : 1);
        }
      }
    }

    for (int y = 0; y < m_grid.rows(); ++y)
    {
      for (int x = 0; x < m_grid.cols(); ++x)
      {
        cv::circle(image_out, grid_draw(y, x), m_visibility_mode ? 5 : 3, cv::Scalar(196, 15, 219), -1);
      }
    }

    for (const Edge& edge : m_selected_edges)
    {
      cv::line(image_out, cv::Point(grid_draw(edge.p1)), cv::Point(grid_draw(edge.p2)), cv::Scalar(104, 138, 54), m_visibility_mode ? 3 : 2);
    }
  }

  if (m_output_edge == OutputEdge::PATCHES)
  {
    cv::Mat mask = cv::Mat::zeros(image_out.size(), CV_8UC1);
    for (const PatchRegion& p : m_patches)
    {
      p.draw<unsigned char>(mask, 255);
    }

    if (m_visibility_mode)
    {
      cv::dilate(mask, mask, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
    }

    image_out.setTo(cv::Vec3b(255, 255, 255), mask != 0);

    mask = 0;
    for (const BezierCurve& c : m_curves_fitted)
    {
      c.draw<unsigned char>(mask, 255);
    }

    if (m_visibility_mode)
    {
      cv::dilate(mask, mask, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));
    }

    image_out.setTo(cv::Vec3b(50, 8, 156), mask);
  }

  if (m_mode == Mode::EDIT_BILATERAL_EDGES || m_mode == Mode::EDIT_FILTERED_EDGES || m_mode == Mode::EDIT_DRAW_EDGES)
  {
    if (m_mouse_pos.x >= 0 && m_mouse_pos.y >= 0)
    {
      cv::circle(image_out, m_mouse_pos, static_cast<int>(m_mouse_tool_size), cv::Scalar(128, 128, 128));
    }
  }
  if (m_mode == Mode::EDIT_MAKE_DENSER || m_mode == Mode::EDIT_MAKE_NORMAL || m_mode == Mode::EDIT_MAKE_WIDER)
  {
    if (m_mouse_pos.x >= 0 && m_mouse_pos.y >= 0)
    {
      cv::circle(image_out, m_mouse_pos, static_cast<int>(4.0 * m_mouse_tool_size), cv::Scalar(128, 128, 128));
    }
  }

  const double font_scale = 0.8;
  const int font = cv::FONT_HERSHEY_TRIPLEX;
  const int thickness = 1;
  cv::copyMakeBorder(image_out, image_out, 0, 32, 0, 0, cv::BORDER_CONSTANT, cv::Scalar(255, 255, 255));
  cv::putText(image_out, m_message, cv::Point(0, image_out.rows-10), font, font_scale, cv::Scalar(0, 0, 0), thickness);

  return image_out;
}

MorphGridFuture::MorphGridData MorphGridFuture::update_morphed_grid(DataGrid<unsigned char> grid, cv::Mat edge_image, cv::Mat density_image, int grid_size)
{
  MorphGridData data;

  const double attraction_radius = 0.5 * grid_size;

  cv::Mat dist;
  cv::distanceTransform(255-edge_image, dist, data.im_label, cv::DIST_L2, cv::DIST_MASK_PRECISE);  
  dist = cv::min(dist, attraction_radius-dist);
  dist = cv::max(dist, 0.0f);
  data.im_distance = dist;

  cv::Sobel(data.im_distance, data.gradient_distance_x, CV_32FC1, 1, 0);
  cv::Sobel(data.im_distance, data.gradient_distance_y, CV_32FC1, 0, 1);

  data.grid = relax_grid(grid, 300, data.im_distance, data.gradient_distance_x, data.gradient_distance_y, density_image, edge_image.rows / grid.rows());

  for (int y = 0; y < data.grid.rows(); ++y)
  {
    for (int x = 0; x < data.grid.cols(); ++x)
    {
      cv::Point2d& p = data.grid(y, x);
      if (p.x < 0.0)
      {
        p.x = 0.0;
      }
      if (p.x > edge_image.cols-1)
      {
        p.x = edge_image.cols-1;
      }
      if (p.y < 0.0)
      {
        p.y = 0.0;
      }
      if (p.y > edge_image.rows-1)
      {
        p.y = edge_image.rows-1;
      }
    }
  }

  return data;
}

MorphGridFuture::PatchData MorphGridFuture::update_patch_data(DataGrid<unsigned char> grid, cv::Mat edge_image, std::vector<Edge> selected_edges)
{
  PatchData patch_data;
  DataGrid<Vector<cv::Point2d>> edge_grid = extract_edge_grid(grid, edge_image);
  DataGrid<CurveData> curve_grid(grid);
  for (const Edge& p : selected_edges)
  {
    fit_single_edge(curve_grid, edge_grid, p.p1, p.p2);
  }

  std::vector<BezierCurve> curves;
  for (int y = 0; y < curve_grid.rows(); ++y)
  {
    for (int x = 0; x < curve_grid.cols(); ++x)
    {
      if (!curve_grid.data(y, x).curve_diag.is_empty_curve())
      {
        curves.push_back(curve_grid.data(y, x).curve_diag);
      }
      if (!curve_grid.data(y, x).curve_horiz.is_empty_curve())
      {
        curves.push_back(curve_grid.data(y, x).curve_horiz);
      }
      if (!curve_grid.data(y, x).curve_vert.is_empty_curve())
      {
        curves.push_back(curve_grid.data(y, x).curve_vert);
      }
    }
  }
  patch_data.curves_fitted = curves;

  fill_smooth_curves(curve_grid);
  patch_data.curve_grid = curve_grid;

  patch_data.patches = generate_patches_from_curve_grid(0, curve_grid);

  return patch_data;
}

void MorphGridFuture::handle_key_press()
{
  if (key_update())
  {
    if (key() == 0x17)
    {
      m_save_state =true;
      m_message = "Save state.";
      return;
    }

    switch (std::tolower(key()))
    {
    case '1':
      m_output_image = OutputImage::IMAGE;
      m_message = "Switch to image.";
      break;
    case '2':
      m_output_image = OutputImage::BILATERAL;
      m_message = "Switch to bilateral filter image.";
      break;
    case '3':
      m_output_image = OutputImage::FILTERED;
      m_message = "Switch to rolling guidance filtered image.";
      break;
    case '4':
      m_output_image = OutputImage::MASK_FILTERED;
      m_message = "Switch to rolling guidance filter mask.";
      break;
    case '5':
      m_output_image = OutputImage::MASK_BILATERAL;
      m_message = "Switch to bilateral filter mask.";
      break;
    case '6':
      m_output_image = OutputImage::DISTANCE;
      m_message = "Switch to distance map.";
      break;
    case '7':
      m_output_image = OutputImage::DISTANCE_GRAD;
      m_message = "Switch to distance gradient.";
      break;
    case '8':
      m_output_image = OutputImage::DENSITY;
      m_message = "Switch to grid density.";
      break;
    case 'q':
      m_output_edge = OutputEdge::NONE;
      m_message = "Don't show edge pixels.";
      break;
    case 'w':
      m_output_edge = OutputEdge::EDGE;
      m_message = "Show computed edge pixels.";
      break;
    case 'e':
      m_output_edge = OutputEdge::EDGE_DRAW;
      m_message = "Show edge visualization.";
      break;
    case 'r':
      m_output_edge = OutputEdge::FILTERED;
      m_message = "Show rolling guidance filtered edge image.";
      break;
    case 't':
      m_output_edge = OutputEdge::BILATERAL;
      m_message = "Show bilateral filter edge image";
      break;
    case 'y':
      m_output_edge = OutputEdge::DRAWN;
      m_message = "Show drawn edge image.";
      break;
    case 'u':
      m_output_edge = OutputEdge::PATCHES;
      m_message = "Show computed patches.";
      break;
    case 'i':
      m_draw_grid = !m_draw_grid;
      m_message = "Toggle draw grid.";
      break;
    case 'o':
      m_visibility_mode = !m_visibility_mode;
      m_message = std::string("Toggle high visibility mode ") + (m_visibility_mode ? "on" : "off");
      break;
    case 'a':
      m_mode = Mode::NONE;
      m_message = "Disable editing.";
      break;
    case 's':
      m_mode = Mode::EDIT_GRID;
      m_message = "Enter grid edit mode.";
      break;
    case 'd':
      m_mode = Mode::EDGE_SELECT;
      m_message = "Enter edge select mode.";
      break;
    case 'f':
      m_mode = Mode::EDIT_FILTERED_EDGES;
      m_message = "Enter filtered edge edit mode.";
      break;
    case 'g':
      m_mode = Mode::EDIT_BILATERAL_EDGES;
      m_message = "Enter bilateral edge edit mode.";
      break;
    case 'h':
      m_mode = Mode::EDIT_DRAW_EDGES;
      m_message = "Enter manual edge draw mode.";
      break;
    case 'z':
      m_mode = Mode::EDIT_MAKE_DENSER;
      m_message = "Locally increase (Shift: decrease) grid density.";
      break;
    case 'x':
      m_mode = Mode::EDIT_MAKE_NORMAL;
      m_message = "Locally restore original grid cell size.";
      break;
    case 'c':
      m_mode = Mode::EDIT_MAKE_WIDER;
      m_message = "Locally decrease (Shift: increase) grid density.";
      break;
    }
  }
}

cv::Point MorphGridFuture::get_nearest_grid_point(int p_x, int p_y)
{
  const cv::Point2d p(p_x, p_y);
  double min_dist = std::numeric_limits<double>::max();
  cv::Point selected_grid_point;
  for (int y = 1; y < m_grid.rows()-1; ++y)
  {
    for (int x = 1; x < m_grid.cols()-1; ++x)
    {
      double dist = cv::norm(p - m_grid(y, x));
      if (dist < min_dist)
      {
        min_dist = dist;
        selected_grid_point.x = x;
        selected_grid_point.y = y;
      }
    }
  }

  if (min_dist < m_grid_size / 4.0)
  {
    return selected_grid_point;
  }
  else
  {
    return cv::Point(-1, -1);
  }
}

int MorphGridFuture::get_nearest_selected_edge(int p_x, int p_y)
{
  const cv::Point2d p(p_x, p_y);

  int min_index = -1;
  double min_dist = std::numeric_limits<double>::max();

  for (int i = 0; i < static_cast<int>(m_selected_edges.size()); ++i)
  {
    const cv::Point2d& l_1 = m_grid(m_selected_edges[i].p1);
    const cv::Point2d& l_2 = m_grid(m_selected_edges[i].p2);
    const double t = project_point_line(p, l_1, l_2);

    if (t >= 0.0 && t <= 1.0)
    {
      const double dist = cv::norm(l_1 + t * (l_2 - l_1) - p);
      if (dist < min_dist)
      {
        min_index = i;
        min_dist = dist;
      }
    }
  }

  if (min_dist < m_grid_size / 4.0)
  {
    return min_index;
  }
  else
  {
    return -1;
  }
}

void MorphGridFuture::select_grid_point(int p_x, int p_y)
{
  cv::Point p = get_nearest_grid_point(p_x, p_y);
  if (p.x != -1 && p.y != -1)
  {
    m_is_grid_point_selected = true;
    m_selected_grid_point = p;
  }
}

void MorphGridFuture::deselect_grid_point(int p_x, int p_y, bool shift_pressed)
{
  if (m_mode == Mode::EDIT_GRID && m_is_grid_point_selected)
  {
    p_x = std::max(p_x, 1);
    p_x = std::min(p_x, m_image.cols-2);
    p_y = std::max(p_y, 1);
    p_y = std::min(p_y, m_image.rows-2);
    m_grid(m_selected_grid_point) = cv::Point2d(p_x, p_y);
    m_grid.data(m_selected_grid_point) = (shift_pressed ? 0 : 1);
    m_update = true;
  }
  m_is_grid_updated = true;
  m_is_grid_point_selected = false;
}

void MorphGridFuture::handle_edge_select(int p_x, int p_y, bool button_l_pressed, bool delete_modifier)
{
  if (m_mode != Mode::EDGE_SELECT)
  {
    return;
  }

  if (delete_modifier)
  {
    if (button_l_pressed)
    {
      const int delete_index = get_nearest_selected_edge(p_x, p_y);
      if (delete_index >= 0)
      {
        m_selected_edges.erase(m_selected_edges.begin() + delete_index);
        m_update = true;
      }
    }
  }
  else
  {
    cv::Point p = get_nearest_grid_point(p_x, p_y);
    if (p == m_edge_select_last || p.x == -1 || p.y == -1)
    {
      return;
    }

    if (button_l_pressed)
    {
      if (m_edge_select_last.x != -1 && m_edge_select_last.y != -1 &&
        std::abs(m_edge_select_last.x - p.x) <= 1 &&
        std::abs(m_edge_select_last.y - p.y) <= 1)
      {
        m_selected_edges.emplace_back(m_edge_select_last, p);
        m_update = true;
      }
      m_edge_select_last = p;
    }
    else
    {
      m_edge_select_last.x = -1;
      m_edge_select_last.y = -1;
    }
  }
}

cv::Mat get_gaussian_blob(int rows, int cols, double x0, double y0, double sigma_x, double sigma_y=-1.0)
{
  if (sigma_y < 0.0)
  {
    sigma_y = sigma_x;
  }
  
  const double x_spread = 1.0 / (2.0 * sigma_x * sigma_x);
  const double y_spread = 1.0 / (2.0 * sigma_y * sigma_y);

  std::vector<double> gauss_x, gauss_y;

  gauss_x.reserve(cols);
  for (int i = 0; i < cols; ++i)
  {
    double x = i - x0;
    gauss_x.push_back(std::exp(-x*x * x_spread));
  }

  gauss_y.reserve(rows);
  for (int i = 0; i < rows; ++i)
  {
    double y = i - y0;
    gauss_y.push_back(std::exp(-y*y * y_spread));
  }

  cv::Mat kernel = cv::Mat::zeros(rows, cols, CV_32FC1);
  for (int y = 0; y < rows; ++y)
    for (int x = 0; x < cols; ++x)
    {
      kernel.at<float>(y, x) = static_cast<float>(gauss_x[x] * gauss_y[y]);
    }

  return kernel;
}

void MorphGridFuture::draw_new_edge(cv::Mat mask, cv::Point p_1, cv::Point p_2, double size_draw, double size_delete, int event)
{
  if (m_mouse_l_pressed)
  {
    if (m_mouse_shift_pressed)
    {
      cv::line(mask, p_1, p_2, cv::Scalar(0), static_cast<int>(size_delete));
    }
    else
    {
      cv::line(mask, p_1, p_2, cv::Scalar(255), static_cast<int>(size_draw));
    }
  }
  if (event == cv::EVENT_LBUTTONUP)
  {
    m_update = true;
  }
}

void MorphGridFuture::draw_soft_blob(cv::Mat mask, cv::Point p_1, double target_value, double size, double amount, int event)
{
  cv::Mat blob = amount * get_gaussian_blob(mask.rows, mask.cols, p_1.x, p_1.y, size, size);
  cv::Mat target_value_mat = target_value * cv::Mat(mask.size(), mask.type(), cv::Scalar(target_value));

  if (m_mouse_l_pressed)
  {
    if (m_mouse_shift_pressed) 
    {
      mask = (-1.0 * target_value_mat.mul(blob)) + mask.mul(1.0 - blob);
    }
    else
    {
      mask = target_value_mat.mul(blob) + mask.mul(1.0 - blob);
    }
    m_update = true;
  }

  if (event == cv::EVENT_LBUTTONUP)
  {
    m_update = true;
  }
}

void MorphGridFuture::on_mouse(int event, int x, int y, int flags)
{
  switch (event)
  {
  case cv::EVENT_MOUSEWHEEL:
    m_mouse_tool_size += static_cast<double>(cv::getMouseWheelDelta(flags)) / 120.0;
    m_mouse_tool_size = std::max(1.0, m_mouse_tool_size);
    m_mouse_tool_size = std::min(60.0, m_mouse_tool_size);
    break;
  case cv::EVENT_LBUTTONDOWN:
    m_mouse_l_pressed = true;
    m_mouse_shift_pressed = ((flags & cv::EVENT_FLAG_SHIFTKEY) != 0);
    select_grid_point(x, y);
    break;
  case cv::EVENT_LBUTTONUP:
    m_mouse_l_pressed = false;
    deselect_grid_point(x, y, ((flags & cv::EVENT_FLAG_SHIFTKEY) != 0));
    break;
  }

  handle_edge_select(x, y, m_mouse_l_pressed, m_mouse_shift_pressed);

  if (event != cv::EVENT_MOUSEWHEEL)
  {
    m_mouse_pos_old = m_mouse_pos;
    m_mouse_pos.x = x;
    m_mouse_pos.y = y;

    if (m_mode == Mode::EDIT_BILATERAL_EDGES)
    {
      draw_new_edge(m_mask_bilateral, m_mouse_pos_old, m_mouse_pos, 2.0*m_mouse_tool_size, 2.0*m_mouse_tool_size, event);
    }
    else if (m_mode == Mode::EDIT_FILTERED_EDGES)
    {
      draw_new_edge(m_mask_filtered, m_mouse_pos_old, m_mouse_pos, 2.0*m_mouse_tool_size, 2.0*m_mouse_tool_size, event);
    }
    else if (m_mode == Mode::EDIT_DRAW_EDGES)
    {
      draw_new_edge(m_edge_drawn, m_mouse_pos_old, m_mouse_pos, 1.0, 2.0*m_mouse_tool_size, event);
    }
    else if (m_mode == Mode::EDIT_MAKE_DENSER || m_mode == Mode::EDIT_MAKE_WIDER)
    {
      draw_soft_blob(m_grid_density, m_mouse_pos, m_mode == Mode::EDIT_MAKE_DENSER ? 1.0 : -1.0, 3.0 * m_mouse_tool_size, 0.05, event);
    }
    else if (m_mode == Mode::EDIT_MAKE_NORMAL)
    {
      draw_soft_blob(m_grid_density, m_mouse_pos, 0.0, 3.0 * m_mouse_tool_size, 0.2, event);
    }
  }
}

template <typename T, typename TPoint>
T interpolate(const cv::Mat_<T> image, cv::Point_<TPoint> p)
{
  cv::Point p_int(p);
  p -= cv::Point_<TPoint>(p_int);
  
  if (p_int.x >= 0 && p_int.x < image.cols-1 && p_int.y >= 0 && p_int.y < image.rows-1)
  {
    T f_1 = static_cast<T>((1.0 - p.x) * image.at<T>(p_int.y, p_int.x) + p.x * image.at<T>(p_int.y, p_int.x+1));
    T f_2 = static_cast<T>((1.0 - p.x) * image.at<T>(p_int.y+1, p_int.x) + p.x * image.at<T>(p_int.y+1, p_int.x+1));
    return static_cast<T>((1.0 - p.y) * f_1 + p.y * f_2);
  }

  if (p_int.x < 0)
  {
    p_int.x = 0;
  }

  if (p_int.x > image.cols-1)
  {
    p_int.x = image.cols-1;
  }

  if (p_int.y < 0)
  {
    p_int.y = 0;
  }

  if (p_int.y > image.rows-1)
  {
    p_int.y = image.rows-1;
  }

  return image.at<T>(p_int);
}

DataGrid<unsigned char> MorphGridFuture::relax_grid(DataGrid<unsigned char> grid, int num_iterations, cv::Mat im_dist, cv::Mat gradient_dist_x, cv::Mat gradient_dist_y, cv::Mat density, int grid_size) const
{
  const float w_1 = 100.0f;
  const float w_2 = 1.5f;
  const float d_t = 0.02f;
  const float damping = 0.2f;

  cv::Mat X = grid.to_mat<float>();
  cv::Mat V = cv::Mat::zeros(grid.rows(),  grid.cols(), CV_32FC2);

  for (int i = 0; i < num_iterations; ++i)
  {
    cv::Mat F = cv::Mat::zeros(grid.rows(), grid.cols(), CV_32FC2);

    for (int y = 1; y < grid.rows()-1; ++y)
    {
      const cv::Vec2f* X_ptr = reinterpret_cast<const cv::Vec2f*>(X.ptr(y));
      const cv::Vec2f* V_ptr = reinterpret_cast<const cv::Vec2f*>(V.ptr(y));
      cv::Vec2f* F_ptr = reinterpret_cast<cv::Vec2f*>(F.ptr(y));
      
      for (int x = 1; x < grid.cols()-1; ++x)
      {
        if (grid.data(y, x))
        {
          cv::Point2f p(X_ptr[x]);

          const float dist = interpolate<float>(im_dist, p);
          const float d_x = interpolate<float>(gradient_dist_x, p);
          const float d_y = interpolate<float>(gradient_dist_y, p);

          F_ptr[x][0] -= w_1 * d_x;
          F_ptr[x][1] -= w_1 * d_y;

          for (int del_y = -1; del_y <= 1; ++del_y)
          {
            for (int del_x = -1; del_x <= 1; ++del_x)
            {
              if (x != 0 || y != 0)
              {
                const cv::Vec2f& p_other = X.at<cv::Vec2f>(y+del_y, x+del_x);
                float spring_constant = std::powf(1.5f, interpolate<float>(density, cv::Point2f(p_other[0], p_other[1])));
                spring_constant *= del_x * del_x + del_y * del_y;
                F_ptr[x] += spring_constant * w_2 * (p_other - cv::Vec2f(p));
              }
            }
          }

          F_ptr[x] -= damping * V_ptr[x];
        }
      }
    }

    V += d_t * F;
    X += d_t * V;

    for (int y = 1; y < X.rows-1; ++y)
    {
      cv::Vec2f* X_ptr = reinterpret_cast<cv::Vec2f*>(X.ptr(y));
      for (int x = 1; x < X.cols-1; ++x)
      {
        X_ptr[x][0] = std::max(X_ptr[x][0], 1.0f);
        X_ptr[x][0] = std::min(X_ptr[x][0], im_dist.cols-2.0f);
        X_ptr[x][1] = std::max(X_ptr[x][1], 1.0f);
        X_ptr[x][1] = std::min(X_ptr[x][1], im_dist.rows-2.0f);
      }
    }
  }

  for (int y = 1; y < X.rows-1; ++y)
  {
    const cv::Vec2f* X_ptr = reinterpret_cast<const cv::Vec2f*>(X.ptr(y));
    for (int x = 1; x < X.cols-1; ++x)
    {
      grid(y, x) = cv::Point2d(X_ptr[x]);
    }
  }

  return grid;
}

boost::property_tree::ptree MorphGridFuture::save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const
{
  boost::property_tree::ptree tree;
  serialize_enum(tree, "output_image", m_output_image, base_path, path);
  serialize_enum(tree, "output_edge", m_output_edge, base_path, path);
  serialize_enum(tree, "mode", m_mode, base_path, path);
  serialize_image(tree, "image", m_image, base_path, path);
  serialize_image(tree, "image_bilateral", m_image_bilateral, base_path, path);
  serialize_image(tree, "image_filtered", m_image_filtered, base_path, path);
  serialize_image(tree, "mask_bilateral", m_mask_bilateral, base_path, path);
  serialize_image(tree, "mask_filtered", m_mask_filtered, base_path, path);
  serialize_image(tree, "edge_bilateral", m_edge_bilateral, base_path, path);
  serialize_image(tree, "edge_filtered", m_edge_filtered, base_path, path);
  serialize_image(tree, "edge_drawn", m_edge_drawn, base_path, path);
  serialize_image(tree, "distance", m_distance, base_path, path);
  serialize_image(tree, "density", m_grid_density, base_path, path);
  serialize_image(tree, "label", m_label, base_path, path);
  serialize_image(tree, "distance_grad_x", m_distance_grad_x, base_path, path);
  serialize_image(tree, "distance_grad_y", m_distance_grad_y, base_path, path);
  serialize(tree, "message", m_message, base_path, path);
  serialize(tree, "patches", m_patches, base_path, path);
  serialize(tree, "curves_fitted", m_curves_fitted, base_path, path);
  serialize(tree, "curve_grid", m_curve_grid, base_path, path);
  serialize(tree, "mouse_tool_size", m_mouse_tool_size, base_path, path);
  serialize(tree, "grid", m_grid, base_path, path);
  serialize(tree, "grid_size", m_grid_size, base_path, path);
  serialize(tree, "selected_edges", m_selected_edges, base_path, path);
  serialize(tree, "edge_select_last", m_edge_select_last, base_path, path);
  serialize(tree, "visibility_mode", m_visibility_mode, base_path, path);
  return tree;
}

void MorphGridFuture::load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree)
{
  deserialize_enum(tree, "output_image", m_output_image, base_path);
  deserialize_enum(tree, "output_edge", m_output_edge, base_path);
  deserialize_enum(tree, "mode", m_mode, base_path);
  deserialize_image(tree, "image", m_image, base_path);
  deserialize_image(tree, "image_bilateral", m_image_bilateral, base_path);
  deserialize_image(tree, "image_filtered", m_image_filtered, base_path);
  deserialize_image(tree, "mask_bilateral", m_mask_bilateral, base_path);
  deserialize_image(tree, "mask_filtered", m_mask_filtered, base_path);
  deserialize_image(tree, "edge_bilateral", m_edge_bilateral, base_path);
  deserialize_image(tree, "edge_filtered", m_edge_filtered, base_path);
  deserialize_image(tree, "edge_drawn", m_edge_drawn, base_path);
  deserialize_image(tree, "distance", m_distance, base_path);
  deserialize_image(tree, "density", m_grid_density, base_path);
  deserialize_image(tree, "label", m_label, base_path);
  deserialize_image(tree, "distance_grad_x", m_distance_grad_x, base_path);
  deserialize_image(tree, "distance_grad_y", m_distance_grad_y, base_path);
  deserialize(tree, "message", m_message, base_path);
  deserialize(tree, "patches", m_patches, base_path);
  deserialize(tree, "curves_fitted", m_curves_fitted, base_path);
  deserialize(tree, "curve_grid", m_curve_grid, base_path);
  deserialize(tree, "mouse_tool_size", m_mouse_tool_size, base_path);
  deserialize(tree, "grid", m_grid, base_path);
  deserialize(tree, "grid_size", m_grid_size, base_path);
  deserialize(tree, "selected_edges", m_selected_edges, base_path);
  deserialize(tree, "edge_select_last", m_edge_select_last, base_path);
  deserialize(tree, "visibility_mode", m_visibility_mode, base_path);
}

void MorphGridFuture::load_input_partial(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree)
{
  deserialize_image(tree, "mask_bilateral", m_mask_bilateral, base_path);
  deserialize_image(tree, "mask_filtered", m_mask_filtered, base_path);
  deserialize_image(tree, "edge_drawn", m_edge_drawn, base_path);
}
