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

#include "../impl.hpp"

#include <hbrs/theta_utils/dt/command.hpp>
#include <iostream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

HBRS_THETA_UTILS_API
void
execute(help_cmd cmd) {
	std::cout << cmd.help << std::endl;
}

HBRS_THETA_UTILS_NAMESPACE_END
