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

#ifndef HBRS_THETA_UTILS_DETAIL_GATHER_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_GATHER_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/theta_utils/dt/theta_field_matrix.hpp>
#include <hbrs/theta_utils/detail/scatter.hpp>

#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace detail {

template<typename Algorithm, typename LocalSize>
struct gather_control {
	Algorithm algorithm;
	LocalSize local_size;
};

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
theta_field_matrix
gather(
	mpl::el_dist_matrix<
		double, El::VC, El::STAR, El::ELEMENT
	> const& from,
	gather_control<
		theta_field_distribution_1,
		mpl::matrix_size<std::size_t, std::size_t>
	> const& ctrl
);
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_GATHER_IMPL_HPP
