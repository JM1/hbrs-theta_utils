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
namespace detail {

HBRS_THETA_UTILS_API
bool
iff(bool lhs, bool rhs) {
	/* p | q | p != q | !(p != q)
	 * --+---+--------+----------
	 * F | F |   F    |     T
	 * T | F |   T    |     F
	 * F | T |   T    |     F
	 * T | T |   F    |     T
	 */
	return !(lhs != rhs);
}

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END
