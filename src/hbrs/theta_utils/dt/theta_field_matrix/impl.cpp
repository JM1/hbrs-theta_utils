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

#include <hbrs/mpl/dt/exception.hpp>
#include <hbrs/mpl/detail/log.hpp>
#include <boost/assert.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

theta_field_matrix::theta_field_matrix(std::vector<theta_field> data) : data_{data} {}

theta_field_matrix::theta_field_matrix(mpl::matrix_size<std::size_t, std::size_t> sz) {
	using namespace hbrs::mpl;
	
	if (sz.m()%3 != 0) {
		HBRS_MPL_LOG_TRIVIAL(debug) << "theta_field_matrix.size() == {" << sz.m() << "," << sz.n() << "}";
		BOOST_THROW_EXCEPTION((incompatible_matrix_exception{} << errinfo_matrix_size{sz}));
	}
	
	data_ = std::vector<theta_field>(
		sz.n(),
		theta_field{
			{} /*density*/,
			std::vector<double>(sz.m()/3) /*x_velocity*/,
			std::vector<double>(sz.m()/3) /*y_velocity*/,
			std::vector<double>(sz.m()/3) /*z_velocity*/,
			{} /*pressure*/,
			{} /*residual*/,
			{} /* global_id */,
			{} /* ndomains */
		}
	);
}

HBRS_THETA_UTILS_DEFINE_ATTR(data, std::vector<theta_field>, theta_field_matrix)

mpl::matrix_size<std::size_t, std::size_t>
theta_field_matrix::size() const {
	BOOST_ASSERT(data_.size() > 0);
	std::size_t xv_sz = data_[0].x_velocity().size();
	std::size_t yv_sz = data_[0].y_velocity().size();
	std::size_t zv_sz = data_[0].z_velocity().size();
	BOOST_ASSERT(xv_sz == yv_sz);
	BOOST_ASSERT(xv_sz == zv_sz);
	
	std::size_t m = xv_sz + yv_sz + zv_sz;
	std::size_t n = data_.size();
	
	BOOST_ASSERT(m > 0);
	BOOST_ASSERT(n > 0);
	
	return {m,n};
}

HBRS_THETA_UTILS_NAMESPACE_END
