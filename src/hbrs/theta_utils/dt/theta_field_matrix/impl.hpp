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

#ifndef HBRS_THETA_UTILS_DT_THETA_FIELD_MATRIX_IMPL_HPP
#define HBRS_THETA_UTILS_DT_THETA_FIELD_MATRIX_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/core/preprocessor.hpp>

#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

#ifdef HBRS_MPL_ENABLE_MATLAB
    #include <hbrs/mpl/dt/ml_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_MATLAB

#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/ctsam.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>

#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/detail/matrix.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/hana/core.hpp>
#include <vector>
#include <type_traits>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace mpl = hbrs::mpl;

struct theta_field_matrix {
public:
	theta_field_matrix(std::vector<theta_field> data = {});
	theta_field_matrix(mpl::matrix_size<std::size_t, std::size_t> sz);
	
	theta_field_matrix(theta_field_matrix const&) = default;
	theta_field_matrix(theta_field_matrix &&) = default;
	
	theta_field_matrix&
	operator=(theta_field_matrix const&) = default;
	theta_field_matrix&
	operator=(theta_field_matrix &&) = default;
	
	mpl::matrix_size<std::size_t, std::size_t>
	size() const;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(data, std::vector<theta_field>)
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::theta_field_matrix > {
	using type = hbrs::theta_utils::theta_field_matrix_tag;
};

template <>
struct make_impl<hbrs::theta_utils::theta_field_matrix_tag> {
	template<
		typename From,
		typename std::enable_if_t<
			(
				(
					std::is_same_v< hana::tag_of_t<From>, hbrs::mpl::rtsam_tag > ||
					std::is_same_v< hana::tag_of_t<From>, hbrs::mpl::ctsam_tag > ||
					std::is_same_v< hana::tag_of_t<From>, hbrs::mpl::sm_tag >
				) &&
				std::is_same_v<
					/* Ring == double? */
					std::decay_t<
						decltype( std::declval<From>().at(hbrs::mpl::matrix_index<std::size_t, std::size_t>{0u,0u}))
					>, double
				>
			)
			#ifdef HBRS_MPL_ENABLE_MATLAB
			|| (
				(
					std::is_same_v< hana::tag_of_t<From>, hbrs::mpl::ml_matrix_tag >
				) &&
				std::is_same_v<
					/* Ring == double? */
					std::decay_t<
						decltype( std::declval<From>().at(hbrs::mpl::matrix_index<int, int>{0,0}))
					>, double
				>
			)
			#endif // !HBRS_MPL_ENABLE_MATLAB
			#ifdef HBRS_MPL_ENABLE_ELEMENTAL
			|| (
				(
					std::is_same_v< hana::tag_of_t<From>, hbrs::mpl::el_matrix_tag >
				) &&
				std::is_same_v<
					/* Ring == double? */
					std::decay_t<
						decltype( std::declval<From>().at(hbrs::mpl::matrix_index<El::Int, El::Int>{0,0}))
					>, double
				>
			)
			#endif // !HBRS_MPL_ENABLE_ELEMENTAL
		>* = nullptr
	>
	static hbrs::theta_utils::theta_field_matrix
	apply(From && from) {
		namespace mpl = hbrs::mpl;
		
		auto sz_ = (*mpl::size)(from);
		auto m_ = boost::numeric_cast<std::size_t> ((*mpl::m)(sz_));
		auto n_ = boost::numeric_cast<std::size_t> ((*mpl::n)(sz_));
		
		hbrs::theta_utils::theta_field_matrix to{{m_, n_}};
		
		return hbrs::theta_utils::detail::copy_matrix(HBRS_MPL_FWD(from), to);
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_THETA_FIELD_MATRIX_IMPL_HPP
