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

#include <hbrs/theta_utils/dt/exception.hpp>
#include <hbrs/mpl/detail/mpi.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpi = hbrs::mpl::detail::mpi;
namespace detail {

mpl::matrix_size<std::size_t, std::size_t>
local_size(std::vector<theta_field> const& series) {
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

mpl::matrix_size<std::size_t, std::size_t>
global_size(std::vector<theta_field> const& series) {
	auto lcl_sz = local_size(series);
	
	if (mpi::comm_size() == 1) {
		return lcl_sz;
	}
	
	std::size_t lcl_m = lcl_sz.m();
	std::size_t lcl_n = lcl_sz.n();
	std::size_t gbl_m;
	mpi::allreduce(&lcl_m, &gbl_m, 1, MPI_SUM, MPI_COMM_WORLD);
	
	std::size_t gbl_min_n, gbl_max_n;
	mpi::allreduce(&lcl_n, &gbl_min_n, 1, MPI_MIN, MPI_COMM_WORLD);
	mpi::allreduce(&lcl_n, &gbl_max_n, 1, MPI_MAX, MPI_COMM_WORLD);
	
	if(gbl_min_n != gbl_max_n) {
		BOOST_THROW_EXCEPTION((mpl::incompatible_matrix_exception{} << mpl::errinfo_matrix_size{lcl_sz}));
	}
	
	return { gbl_m, gbl_min_n };
}

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
