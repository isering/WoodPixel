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

#ifndef BILATERAL_FILTER_FUTURE_HPP_
#define BILATERAL_FILTER_FUTURE_HPP_

#include <chrono>
#include <future>

#include "future_base.hpp"
#include "serializable.hpp"

class BilateralFilterFuture : public FutureBase, public Serializable
{
public:
  BilateralFilterFuture(const std::string &window_name, cv::Mat image) : FutureBase(window_name),
                                                                         m_image(image),
                                                                         m_d(7),
                                                                         m_sigma_color(50),
                                                                         m_sigma_space(50)
  {
    if (m_image.channels() == 3)
    {
      cv::cvtColor(m_image, m_image, cv::COLOR_BGR2GRAY);
    }

    add_trackbar("d", &m_d, 25);
    add_trackbar("sigmaColor", &m_sigma_color, 200);
    add_trackbar("sigmaSpace", &m_sigma_space, 200);
  }

  cv::Mat filter(int d, int sigma_color, int sigma_space)
  {
    cv::Mat image_filtered;
    cv::bilateralFilter(m_image, image_filtered, d, sigma_color, sigma_space);
    return image_filtered;
  }

  cv::Mat filtered_image() const
  {
    return m_filtered_image.clone();
  }

  virtual void loop() override
  {
    if (m_update && !m_filtered_image_future.valid())
    {
      m_filtered_image_future = std::async(&BilateralFilterFuture::filter, this, m_d, m_sigma_color, m_sigma_space);
      m_update = false;
    }

    if (m_filtered_image_future.valid() && m_filtered_image_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
      m_filtered_image = m_filtered_image_future.get();
      m_ready = true;
    }
  }

  virtual void load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree) override;
  virtual boost::property_tree::ptree save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const override;

private:
  cv::Mat m_image;
  int m_d, m_sigma_color, m_sigma_space;

  cv::Mat m_filtered_image;
  std::future<cv::Mat> m_filtered_image_future;
};

#endif /* BILATERAL_FILTER_FUTURE_HPP_ */