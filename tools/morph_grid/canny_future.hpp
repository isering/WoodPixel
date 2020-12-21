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

#ifndef CANNY_FUTURE_HPP_
#define CANNY_FUTURE_HPP_

#include <future>

#include "future_base.hpp"
#include "serializable.hpp"

class CannyFuture : public FutureBase, public Serializable
{
public:
  CannyFuture(const std::string &window_name, const std::string &class_name = "canny") : FutureBase(window_name),
                                                                                         m_canny_threshold_1(100),
                                                                                         m_canny_threshold_2(200)
  {
    add_trackbar(class_name + "Thresh1", &m_canny_threshold_1, 500);
    add_trackbar(class_name + "Thresh2", &m_canny_threshold_2, 500);
  }

  cv::Mat canny(cv::Mat image, int threshold_1, int threshold_2)
  {
    cv::Mat edge_image;
    cv::Canny(image, edge_image, threshold_1, threshold_2, 3);
    return edge_image;
  }

  void set_image(cv::Mat image)
  {
    m_image = image.clone();
    m_update = true;
  }

  cv::Mat edge_image() const
  {
    return m_edge_image.clone();
  }

  virtual void loop() override
  {
    if (m_update && !m_edge_image_future.valid() && !m_image.empty())
    {
      m_edge_image_future = std::async(&CannyFuture::canny, this, m_image, m_canny_threshold_1, m_canny_threshold_2);
      m_update = false;
    }

    if (m_edge_image_future.valid() && m_edge_image_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
      m_edge_image = m_edge_image_future.get();
      m_ready = true;
    }
  }

  virtual void load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const override;

private:
  cv::Mat m_image;
  cv::Mat m_edge_image;
  std::future<cv::Mat> m_edge_image_future;

  int m_canny_threshold_1;
  int m_canny_threshold_2;
};

#endif /* CANNY_FUTURE_HPP_ */