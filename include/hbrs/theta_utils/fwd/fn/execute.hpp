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

#ifndef HBRS_THETA_UTILS_FWD_FN_EXECUTE_HPP
#define HBRS_THETA_UTILS_FWD_FN_EXECUTE_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/dt/command.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

void
execute(version_cmd);

void
execute(help_cmd);

void
execute(visualize_cmd cmd);

void
execute(pca_cmd cmd);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_FN_EXECUTE_HPP
