/* Copyright (c) 2016 Jakob Meng, <jakobmeng@web.de>
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

#include <hbrs/theta_utils/dt/nc_dimension.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

nc_dimension::nc_dimension(
		std::string name,
		std::size_t length
	) : name_{name}, length_{length}
	{}

HBRS_THETA_UTILS_DEFINE_ATTR(name, std::string, nc_dimension)
HBRS_THETA_UTILS_DEFINE_ATTR(length, std::size_t, nc_dimension)

HBRS_THETA_UTILS_NAMESPACE_END