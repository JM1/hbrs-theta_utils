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

#ifndef HBRS_THETA_UTILS_DETAIL_COPY_MATRIX_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_COPY_MATRIX_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>

#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>

#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/less_equal.hpp>
#include <hbrs/mpl/fn/greater_equal.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>

#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_MATLAB
    #include <hbrs/mpl/dt/ml_matrix.hpp>
#endif
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif

#include <boost/assert.hpp>
#include <vector>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;

namespace detail {

mpl::matrix_size<std::size_t, std::size_t>
size(std::vector<theta_field> const& series);

template<
	typename To,
	typename std::enable_if_t<
		#ifdef HBRS_MPL_ENABLE_MATLAB
			std::is_same< hana::tag_of_t<To>, mpl::ml_matrix_tag >::value ||
		#endif
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
			std::is_same< hana::tag_of_t<To>, mpl::el_matrix_tag >::value ||
		#endif
		std::is_same< hana::tag_of_t<To>, mpl::rtsam_tag >::value
	>* = nullptr
>
decltype(auto)
copy_matrix(std::vector<theta_field> const& from, To && to) {
	auto from_sz = size(from);
	auto from_m = (*mpl::m)(from_sz);
	auto from_n = (*mpl::n)(from_sz);
	
	auto to_sz = (*mpl::size)(to);
	auto to_m = (*mpl::m)(to_sz);
	auto to_n = (*mpl::n)(to_sz);
	
	BOOST_ASSERT((*mpl::less_equal)(from_m, to_m));
	BOOST_ASSERT((*mpl::equal)(from_n, to_n));
	
	for (std::size_t j = 0; j < from_n; ++j) {
		BOOST_ASSERT(from[j].density().empty());
		BOOST_ASSERT(from[j].pressure().empty());
		BOOST_ASSERT(from[j].residual().empty());
	
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i, j)) = from[j].x_velocity().begin()[i];
		}
		
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i+from_m/3, j)) = from[j].y_velocity().begin()[i];
		}
		
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i+from_m/3*2, j)) = from[j].z_velocity().begin()[i];
		}
	}
	
	return HBRS_MPL_FWD(to);
}

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	inline decltype(auto)
	copy_matrix(
		std::vector<theta_field> const& from,
		hbrs::mpl::el_dist_matrix<double, El::VC, El::STAR, El::ELEMENT> && to
	) {
		El::Zero(to.data()); // if to-matrix is larger than from-field then unused rows will just be zero.
		to.data().Matrix() = copy_matrix(from, hbrs::mpl::el_matrix<double>(to.data().Matrix())).data();
		return HBRS_MPL_FWD(to);
	}
#endif

template<
	typename From,
	typename std::enable_if_t<
		#ifdef HBRS_MPL_ENABLE_MATLAB
			std::is_same< hana::tag_of_t<From>, mpl::ml_matrix_tag >::value ||
		#endif
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
			std::is_same< hana::tag_of_t<From>, mpl::el_matrix_tag >::value ||
		#endif
		std::is_same< hana::tag_of_t<From>, mpl::rtsam_tag >::value
	>* = nullptr
>
decltype(auto)
copy_matrix(From const& from, std::vector<theta_field> & to) {
	auto from_sz = (*mpl::size)(from);
	auto from_m = (*mpl::m)(from_sz);
	auto from_n = (*mpl::n)(from_sz);
	
	auto to_sz = size(to);
	auto to_m = (*mpl::m)(to_sz);
	auto to_n = (*mpl::n)(to_sz);
	
	BOOST_ASSERT((*mpl::greater_equal)(from_m, to_m));
	BOOST_ASSERT((*mpl::equal)(from_n, to_n));
	
	for (std::size_t j = 0; j < to_n; ++j) {
		to[j].density().clear();
		to[j].pressure().clear();
		to[j].residual().clear();
		BOOST_ASSERT(to[j].x_velocity().size() == to_m/3);
		BOOST_ASSERT(to[j].y_velocity().size() == to_m/3);
		BOOST_ASSERT(to[j].z_velocity().size() == to_m/3);
	
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to[j].x_velocity().begin()[i] = (*mpl::at)(from, mpl::make_matrix_index(i, j));
		}
		
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to[j].y_velocity().begin()[i] = (*mpl::at)(from, mpl::make_matrix_index(i+to_m/3, j));
		}
		
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to[j].z_velocity().begin()[i] = (*mpl::at)(from, mpl::make_matrix_index(i+to_m/3*2, j));
		}
	}
	
	return HBRS_MPL_FWD(to);
}

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	inline decltype(auto)
	copy_matrix(
		hbrs::mpl::el_dist_matrix<double, El::VC, El::STAR, El::ELEMENT> const& from,
		std::vector<theta_field> & to
	) {
		return copy_matrix(hbrs::mpl::el_matrix<double const>(from.data().LockedMatrix()), to);
	}
#endif

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_COPY_MATRIX_IMPL_HPP
