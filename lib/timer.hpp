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

#ifndef TRLIB_TIMER_HPP_
#define TRLIB_TIMER_HPP_

#include <iostream>

#include <boost/chrono.hpp>

template < class ratio=boost::ratio<1LL, 1LL> >
class Timer
{
public:
	Timer() :
		m_time_start(boost::chrono::process_real_cpu_clock::now()),
		m_time_stop(m_time_start),
		m_running(true) {}

	~Timer() {}

	void restart()
	{
		m_time_start = boost::chrono::process_real_cpu_clock::now();
		m_running = true;
	}

	void stop()
	{
		if (m_running) {
			m_time_stop = boost::chrono::process_real_cpu_clock::now();
			m_running = false;
		} else {
			m_time_stop = m_time_start;
		}
	}

	boost::chrono::duration<float, ratio> duration()
	{
		return (m_running ? boost::chrono::process_real_cpu_clock::now() : m_time_stop) - m_time_start;
	}

	friend std::ostream& operator<<(std::ostream& str, const Timer& rhs)
	{
		boost::chrono::duration<float, ratio> duration = (rhs.m_running ? boost::chrono::process_real_cpu_clock::now() : rhs.m_time_stop) - rhs.m_time_start;
		return str << boost::chrono::duration_short <<  duration;
	}

private:
	boost::chrono::process_real_cpu_clock::time_point m_time_start;
	boost::chrono::process_real_cpu_clock::time_point m_time_stop;
	bool m_running;
};

#endif /* TRLIB_TIMER_HPP_ */