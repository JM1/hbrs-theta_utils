/* Copyright (c) 2016-2018 Jakob Meng, <jakobmeng@web.de>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#ifndef HBRS_THETA_UTILS_FWD_FN_VISUALIZE_HPP
#define HBRS_THETA_UTILS_FWD_FN_VISUALIZE_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/fn/vtk.hpp>
#include <string>
#include <vector>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

struct visualize_options {
	std::string path;
	std::string grid_prefix;
	std::string input_prefix;
	std::string output_prefix;
	vtk_file_format output_format;
	bool overwrite;
	bool simple_numbering;
	std::vector<std::string> includes;
	std::vector<std::string> excludes;
};

void
visualize(
	visualize_options v_opts,
	bool verbose
);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_FN_VISUALIZE_HPP
