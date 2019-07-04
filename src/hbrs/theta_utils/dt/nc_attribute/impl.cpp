/* Copyright (c) 2019 Jakob Meng, <jakobmeng@web.de>
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

#include <hbrs/theta_utils/dt/nc_attribute.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

nc_attribute::nc_attribute(std::string name, array value) 
: name_{name}, value_{value} {}

HBRS_THETA_UTILS_DEFINE_ATTR(name, std::string, nc_attribute)
HBRS_THETA_UTILS_DEFINE_ATTR(value, nc_attribute::array, nc_attribute)

HBRS_THETA_UTILS_NAMESPACE_END