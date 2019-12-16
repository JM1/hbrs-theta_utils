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

#ifndef HBRS_THETA_UTILS_DETAIL_SCATTER_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_SCATTER_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/dt/theta_field_matrix.hpp>

#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL
#include <hbrs/mpl/dt/matrix_size.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace detail {

struct theta_field_distribution_1{};

template<typename Algorithm>
struct scatter_control{
	Algorithm algorithm;
};

mpl::matrix_size<std::size_t, std::size_t>
distributed_size(
	theta_field_matrix const& series,
	theta_field_distribution_1
);

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
hbrs::mpl::el_dist_matrix<double, El::VC, El::STAR, El::ELEMENT>
scatter(
	theta_field_matrix const& from,
	scatter_control<theta_field_distribution_1>
);
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_SCATTER_IMPL_HPP
