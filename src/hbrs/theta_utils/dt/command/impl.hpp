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

#ifndef HBRS_THETA_UTILS_DT_COMMAND_IMPL_HPP
#define HBRS_THETA_UTILS_DT_COMMAND_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/dt/command_option.hpp>
#include <string>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

struct version_cmd {
	generic_options g_opts;
	std::string version;
};

struct help_cmd {
	generic_options g_opts;
	std::string help;
};

struct visualize_cmd {
	generic_options g_opts;
	theta_input_options i_opts;
	theta_output_options o_opts;
	visualize_options v_opts;
};

struct pca_cmd {
	generic_options g_opts;
	theta_input_options i_opts;
	theta_output_options o_opts;
	pca_options pca_opts;
};

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_COMMAND_IMPL_HPP
