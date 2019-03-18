/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_DETAIL_MAKE_THETA_FIELDS_HPP
#define HBRS_THETA_UTILS_DETAIL_MAKE_THETA_FIELDS_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>

#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/detail/copy_matrix.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>

#include <boost/assert.hpp>
#include <vector>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace detail {

std::vector<theta_field>
make_theta_fields(mpl::rtsam<double, mpl::storage_order::row_major> const& from) {
	auto sz = from.size();
	BOOST_ASSERT(sz.m()%3 == 0);
	
	auto to = std::vector<theta_field>(
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
	
	return detail::copy_matrix(from, to);
};

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_MAKE_THETA_FIELDS_HPP
