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

#ifndef HBRS_THETA_UTILS_FN_EXECUTE_FWD_HPP
#define HBRS_THETA_UTILS_FN_EXECUTE_FWD_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/dt/command/fwd.hpp>

//TODO: Make execute() a generic function using HBRS_MPL_DEC_F1 macros!

HBRS_THETA_UTILS_NAMESPACE_BEGIN

HBRS_THETA_UTILS_API
void
execute(version_cmd);

HBRS_THETA_UTILS_API
void
execute(help_cmd);

HBRS_THETA_UTILS_API
void
execute(visualize_cmd cmd);

HBRS_THETA_UTILS_API
void
execute(pca_cmd cmd);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FN_EXECUTE_FWD_HPP
