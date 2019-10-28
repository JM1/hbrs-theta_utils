/* Copyright (c) 2018-2019 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_DETAIL_MAKE_THETA_FIELDS_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_MAKE_THETA_FIELDS_IMPL_HPP

#include "fwd.hpp"

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/ctsam.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/exception.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/detail/copy_matrix.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <type_traits>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace detail {

template<
	typename From,
	typename std::enable_if_t<
		(
			(
				std::is_same_v< hana::tag_of_t<From>, mpl::rtsam_tag > ||
				std::is_same_v< hana::tag_of_t<From>, mpl::ctsam_tag > ||
				std::is_same_v< hana::tag_of_t<From>, mpl::sm_tag >
			) &&
			std::is_same_v<
				/* Ring == double? */
				std::decay_t<
					decltype( std::declval<From>().at(mpl::matrix_index<std::size_t, std::size_t>{0u,0u}))
				>, double
			>
		)
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
		|| (
			(
				std::is_same_v< hana::tag_of_t<From>, mpl::el_matrix_tag >
			) &&
			std::is_same_v<
				/* Ring == double? */
				std::decay_t<
					decltype( std::declval<From>().at(mpl::matrix_index<El::Int, El::Int>{0,0}))
				>, double
			>
		)
		#endif // !HBRS_MPL_ENABLE_ELEMENTAL
	>* = nullptr
>
std::vector<theta_field>
make_theta_fields(
	From && from
) {
	using namespace hbrs::mpl;
	auto sz_ = (*mpl::size)(from);
	auto m_ = boost::numeric_cast<std::size_t> ((*m)(sz_));
	auto n_ = boost::numeric_cast<std::size_t> ((*n)(sz_));
	
	if (m_%3 != 0) {
		BOOST_THROW_EXCEPTION((incompatible_matrix_exception{} << errinfo_matrix_size{{m_, n_}}));
	}
	
	std::vector<theta_field> to(
		n_,
		theta_field{
			{} /*density*/,
			std::vector<double>(m_/3) /*x_velocity*/,
			std::vector<double>(m_/3) /*y_velocity*/,
			std::vector<double>(m_/3) /*z_velocity*/,
			{} /*pressure*/,
			{} /*residual*/,
			{} /* global_id */,
			{} /* ndomains */
		}
	);
	
	return copy_matrix(HBRS_MPL_FWD(from), to);
}

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_MAKE_THETA_FIELDS_IMPL_HPP
