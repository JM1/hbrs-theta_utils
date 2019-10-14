/* Copyright (c) 2016-2019 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_DT_COMMAND_OPTION_IMPL_HPP
#define HBRS_THETA_UTILS_DT_COMMAND_OPTION_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/detail/vtk.hpp>
#include <vector>
#include <string>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

struct generic_options {
	bool verbose = false;
	bool debug = false;
};

struct theta_input_options {
	std::string path;
	std::string grid_prefix;
	std::string pval_prefix;
};

struct theta_output_options {
	std::string path;
	std::string prefix;
	bool overwrite;
};

struct visualize_options {
	std::vector<std::string> includes;
	std::vector<std::string> excludes;
	vtk_file_format format;
	bool simple_numbering;
};

struct pca_options {
	std::vector<std::string> pc_nr_seqs;
	pca_backend backend;
	bool center;
	bool normalize;
};

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_COMMAND_OPTION_IMPL_HPP
