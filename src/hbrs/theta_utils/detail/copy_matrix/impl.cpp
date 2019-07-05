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

mpl::matrix_size<std::size_t, std::size_t>
size(std::vector<theta_field> const& series) {
	BOOST_ASSERT(series.size() > 0);
	std::size_t xv_sz = series[0].x_velocity().size();
	std::size_t yv_sz = series[0].y_velocity().size();
	std::size_t zv_sz = series[0].z_velocity().size();
	BOOST_ASSERT(xv_sz == yv_sz);
	BOOST_ASSERT(xv_sz == zv_sz);
	
	std::size_t m = xv_sz + yv_sz + zv_sz;
	std::size_t n = series.size();
	
	BOOST_ASSERT(m > 0);
	BOOST_ASSERT(n > 0);
	
	return {m,n};
}

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END
