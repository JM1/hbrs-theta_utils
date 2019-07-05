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

#include "impl.hpp"

HBRS_THETA_UTILS_NAMESPACE_BEGIN

nc_variable::nc_variable(std::string name, std::vector<nc_dimension> dims, array data) 
: name_{name}, dimensions_{dims}, data_{data} {}

HBRS_THETA_UTILS_DEFINE_ATTR(name, std::string, nc_variable)
HBRS_THETA_UTILS_DEFINE_ATTR(dimensions, std::vector<nc_dimension>, nc_variable)
HBRS_THETA_UTILS_DEFINE_ATTR(data, nc_variable::array, nc_variable)

HBRS_THETA_UTILS_NAMESPACE_END
