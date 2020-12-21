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

#include <list>

#include <opencv2/highgui/highgui_c.h>

#include "hough_circle_detect.hpp"

std::vector<cv::Vec3f> HoughCircleDetect::compute(const cv::Mat &image, double param_1, double param_2, double dpi, double circle_size_mm)
{
	std::vector<cv::Vec3f> circle_vec;

	double circle_radius_pixel = 0.5 * circle_size_mm * dpi / 25.4;
	int circle_radius_min = static_cast<int>(0.9 * circle_radius_pixel);
	int circle_radius_max = static_cast<int>(1.1 * circle_radius_pixel);

	cv::HoughCircles(image, circle_vec, cv::HOUGH_GRADIENT, 2.0, 25.0, param_1, param_2, circle_radius_min, circle_radius_max);

	return circle_vec;
}

std::future<std::vector<cv::Vec3f>> HoughCircleDetect::compute_background(const cv::Mat &image, double param_1, double param_2, double dpi, double circle_size_mm)
{
	return std::async(std::launch::async, &HoughCircleDetect::compute, this, image, param_1, param_2, dpi, circle_size_mm);
}

std::vector<cv::Vec3f> HoughCircleDetect::compute_interactive(const std::vector<cv::Vec3f> &segmentation_data_in, const cv::Mat &image, const cv::Mat &image_background, double dpi, double circle_size_mm, std::string window_name)
{
	cv::Mat image_background_scaled;
	cv::resize(image_background, image_background_scaled, cv::Size(), m_scale, m_scale);

	if (!cvGetWindowHandle(window_name.c_str()))
	{
		m_param_1 = 100;
		m_param_2 = 100;
		cv::namedWindow(window_name);
		cv::imshow(window_name, image_background_scaled);
		cv::createTrackbar("Canny", window_name, &m_param_1, 200);
		cv::createTrackbar("Accumulator", window_name, &m_param_2, 200);
		cv::setMouseCallback(window_name, HoughCircleDetect::on_mouse, this);
	}

	m_param_old_1 = -1;
	m_param_old_2 = -1;

	m_data.clear();
	m_data_in = segmentation_data_in;

	m_image_background = image_background;
	m_window_name = window_name;

	draw();

	do
	{
		if (segmentation_data_in.empty())
		{
			if (m_future.valid() && m_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
			{
				m_data = m_future.get();
				draw();
			}
			else if (!m_future.valid() && (m_param_1 != m_param_old_1 || m_param_2 != m_param_old_2))
			{

				m_param_old_1 = m_param_1;
				m_param_old_2 = m_param_2;

				m_future = compute_background(image, raw_to_param_1(m_param_1), raw_to_param_2(m_param_2), dpi, circle_size_mm);
			}
		}

		if (cv::waitKey(1) > 0 || !cvGetWindowHandle(window_name.c_str()))
			break;
	} while (true);

	if (cvGetWindowHandle(window_name.c_str()))
		cv::destroyWindow(window_name);

	return m_data;
}

void HoughCircleDetect::draw()
{
	const cv::Vec3b green(0, 255, 0);
	const cv::Vec3b red(0, 0, 255);
	const cv::Vec3b cyan(255, 255, 0);

	cv::Mat im_out = m_image_background.clone();

	for (const cv::Vec3f &s : m_data_in)
	{
		cv::Point center(static_cast<int>(m_scale * s[0]), static_cast<int>(m_scale * s[1]));
		int radius = static_cast<int>(m_scale * s[2]);
		cv::circle(im_out, center, radius, cyan, -1);
	}

	for (const cv::Vec3f &s : m_data)
	{
		cv::Point center(static_cast<int>(m_scale * s[0]), static_cast<int>(m_scale * s[1]));
		int radius = static_cast<int>(m_scale * s[2]);
		cv::circle(im_out, center, radius, green, -1);
	}

	cv::imshow(m_window_name, im_out);
	cv::waitKey(1);
}

void HoughCircleDetect::draw_to_file(const std::string &fname, const cv::Mat &im, const std::vector<cv::Vec3f> &seg)
{
	const cv::Vec3b green(0, 255, 0);
	const cv::Vec3b red(0, 0, 255);
	const cv::Vec3b cyan(255, 255, 0);

	cv::Mat im_out = im.clone();

	for (const cv::Vec3f &s : seg)
	{
		cv::Point center(static_cast<int>(s[0]), static_cast<int>(s[1]));
		int radius = static_cast<int>(s[2]);
		cv::circle(im_out, center, radius, green, -1);
	}

	cv::imwrite(fname, im_out);
}

static inline int dist_sqr(int x1, int y1, int x2, int y2)
{
	const int x_dist = x1 - x2;
	const int y_dist = y1 - y2;
	return x_dist * x_dist + y_dist * y_dist;
}

static inline int line_dist_sqr(int x, int y, const cv::Point &l1, const cv::Point &l2)
{
	const cv::Vec2i M = l2 - l1;
	const cv::Point P(x, y);

	double t = M.dot(P - l2) / static_cast<double>(M.dot(M));
	t = std::max(t, 0.0);
	t = std::min(t, 1.0);

	const cv::Point P2 = l1 + t * l2;

	const int x_dist = P.x - P2.x;
	const int y_dist = P.y - P2.y;

	return x_dist * x_dist + y_dist * y_dist;
}

void HoughCircleDetect::on_mouse(int event, int x, int y, int, void *ptr)
{
	HoughCircleDetect *p = static_cast<HoughCircleDetect *>(ptr);

	if (event == cv::EVENT_RBUTTONUP)
	{
		p->m_data.erase(
			std::remove_if(p->m_data.begin(), p->m_data.end(), [x, y, p](const cv::Vec3f &s) { return dist_sqr(static_cast<int>(x / p->m_scale), static_cast<int>(y / p->m_scale), static_cast<int>(s[0]), static_cast<int>(s[1])) <= s[2] * s[2]; }),
			p->m_data.end());

		p->draw();
	}
	else if (event == cv::EVENT_LBUTTONUP)
	{
		/*
		p->m_data.emplace_back(static_cast<float>(x) / m_scale, static_cast<float>(y) / m_scale, 25.0f);
		p->draw();
	*/
	}
}

double HoughCircleDetect::raw_to_param_1(int val)
{
	return 0.5 * static_cast<double>(val) + 1.0;
}

double HoughCircleDetect::raw_to_param_2(int val)
{
	return static_cast<double>(val) + 1.0;
}