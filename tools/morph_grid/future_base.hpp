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

#ifndef FUTURE_BASE_HPP_
#define FUTURE_BASE_HPP_

#include <opencv2/opencv.hpp>

class FutureBase
{
public:
  FutureBase(const std::string &window_name) : m_window_name(window_name),
                                               m_update(true),
                                               m_ready(false),
                                               m_key(-1),
                                               m_key_update(false)
  {
    cv::setMouseCallback(m_window_name, &FutureBase::mouse_callback, this);
  }

  virtual void loop() = 0;

  void add_trackbar(const std::string &trackbar_name, int *value, int count)
  {
    cv::createTrackbar(trackbar_name, m_window_name, value, count, FutureBase::on_trackbar_change, this);
  }

  bool ready()
  {
    const bool rval = m_ready;
    m_ready = false;
    return rval;
  }

  void key_press(int key)
  {
    if (key != -1)
    {
      m_key = key;
      m_key_update = true;
    }
  }

protected:
  virtual void on_mouse(int event, int x, int y, int flags) {}

  static void on_trackbar_change(int, void *ptr)
  {
    reinterpret_cast<FutureBase *>(ptr)->m_update = true;
  }

  bool key_update()
  {
    bool rval = m_key_update;
    m_key_update = false;
    return rval;
  }

  int key() const
  {
    return m_key;
  }

  std::string m_window_name;
  bool m_update;
  bool m_ready;

private:
  static void mouse_callback(int event, int x, int y, int flags, void *ptr)
  {
    reinterpret_cast<FutureBase *>(ptr)->on_mouse(event, x, y, flags);
  }

  bool m_key_update;
  int m_key;
};

#endif /* FUTURE_BASE_HPP_ */