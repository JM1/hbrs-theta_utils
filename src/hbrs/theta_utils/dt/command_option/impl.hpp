/* Copyright (c) 2016-2020 Jakob Meng, <jakobmeng@web.de>
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

struct HBRS_THETA_UTILS_API generic_options {
	std::size_t verbosity = 0;
	bool debug = false;
};

struct HBRS_THETA_UTILS_API theta_input_options {
	std::string path;
	std::string grid_prefix;
	std::string pval_prefix;
};

struct HBRS_THETA_UTILS_API theta_output_options {
	std::string path;
	std::string prefix;
	bool overwrite;
};

struct HBRS_THETA_UTILS_API visualize_options {
	std::vector<std::string> includes;
	std::vector<std::string> excludes;
	vtk_file_format format;
	bool simple_numbering;
};

struct HBRS_THETA_UTILS_API pca_options {
	std::vector<std::string> pc_nr_seqs;
	pca_backend backend;
	bool center;
	bool normalize;
	bool keep_centered;
};

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_COMMAND_OPTION_IMPL_HPP
