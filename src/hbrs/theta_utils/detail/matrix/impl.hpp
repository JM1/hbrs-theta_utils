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

#ifndef HBRS_THETA_UTILS_DETAIL_MATRIX_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_MATRIX_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>

#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

#ifdef HBRS_MPL_ENABLE_MATLAB
    #include <hbrs/mpl/dt/ml_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_MATLAB

#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/ctsam.hpp>
#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/theta_utils/dt/theta_field_matrix/fwd.hpp>

#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/less_equal.hpp>
#include <hbrs/mpl/fn/greater_equal.hpp>

#include <boost/hana/core/to.hpp>
#include <boost/hana/ext/std/vector.hpp>
#include <boost/assert.hpp>
#include <algorithm>

namespace boost { namespace hana {

template<typename From>
struct to_impl<ext::std::vector_tag, From, when<
	#ifdef HBRS_MPL_ENABLE_MATLAB
		std::is_same_v< tag_of_t<From>, hbrs::mpl::ml_column_vector_tag > ||
	#endif // !HBRS_MPL_ENABLE_MATLAB
	#ifdef HBRS_MPL_ENABLE_ELEMENTAL
		std::is_same_v< tag_of_t<From>, hbrs::mpl::el_column_vector_tag > ||
	#endif // !HBRS_MPL_ENABLE_ELEMENTAL
	false
>> {
	template <typename From_>
	static constexpr auto
	apply(From_ const& from) {
		using namespace hbrs::mpl;
		typedef decltype(hbrs::mpl::at(from, 0)) Ring;
		typedef std::decay_t<Ring> _Ring_;
		
		std::size_t from_sz = boost::numeric_cast<std::size_t>((*hbrs::mpl::size)(from));
		std::vector<_Ring_> to;
		to.resize(from_sz);
		BOOST_ASSERT(to.size() == from_sz);
		
		for (std::size_t i = 0; i < from_sz; ++i) {
			to.at(i) = (*hbrs::mpl::at)(from, i);
		}
		
		return to;
	}
};


#ifdef HBRS_MPL_ENABLE_MATLAB
template<>
struct to_impl<hbrs::mpl::ml_matrix_tag, hbrs::theta_utils::theta_field_matrix_tag> {
	HBRS_THETA_UTILS_API
	static hbrs::mpl::ml_matrix<double>
	apply(hbrs::theta_utils::theta_field_matrix const& rhs);
};
#endif // !HBRS_MPL_ENABLE_MATLAB

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
template<>
struct to_impl<ext::std::vector_tag, hbrs::mpl::el_dist_column_vector_tag> {
	template<typename Ring, El::Dist Columnwise, El::Dist Rowwise, El::DistWrap Wrapping>
	static constexpr auto
	apply(hbrs::mpl::el_dist_column_vector<Ring, Columnwise, Rowwise, Wrapping> const& from) {
		typedef std::decay_t<Ring> _Ring_;
		
		El::DistMatrixReadProxy<Ring, _Ring_, El::STAR, El::STAR, Wrapping> from_pxy = {from.data()};
		
		auto from_local = from_pxy.GetLocked();
		std::size_t from_sz = boost::numeric_cast<std::size_t>(from_local.Height());
		BOOST_ASSERT(from_local.Width() == 1);
		BOOST_ASSERT(from_local.Height() == from.data().Height());
		
		std::vector<_Ring_> to;
		to.resize(from_sz);
		BOOST_ASSERT(to.size() == from_sz);
		
		for (std::size_t i = 0; i < from_sz; ++i) {
			to.at(i) = from_local.Get(boost::numeric_cast<El::Int>(i), 0);
		}
		
		return to;
	}
};

template<>
struct to_impl<hbrs::mpl::el_matrix_tag, hbrs::theta_utils::theta_field_matrix_tag> {
	HBRS_THETA_UTILS_API
	static hbrs::mpl::el_matrix<double>
	apply(hbrs::theta_utils::theta_field_matrix const& rhs);
};
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

/* namespace hana */ } /* namespace boost */ }

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;

namespace detail {

template<
	typename From,
	typename To,
	typename std::enable_if_t<
		(
			#ifdef HBRS_MPL_ENABLE_MATLAB
				std::is_same_v< hana::tag_of_t<To>, mpl::ml_matrix_tag > ||
			#endif
			#ifdef HBRS_MPL_ENABLE_ELEMENTAL
				std::is_same_v< hana::tag_of_t<To>, mpl::el_matrix_tag > ||
			#endif
			std::is_same_v< hana::tag_of_t<To>, mpl::rtsam_tag > ||
			std::is_same_v< hana::tag_of_t<To>, mpl::ctsam_tag > ||
			std::is_same_v< hana::tag_of_t<To>, mpl::sm_tag >
		) &&
		std::is_same_v<From, theta_field_matrix>
	>* = nullptr
>
decltype(auto)
copy_matrix(From const& from, To && to) {
	auto from_sz = from.size();
	auto from_m = (*mpl::m)(from_sz);
	auto from_n = (*mpl::n)(from_sz);
	
	auto to_sz = (*mpl::size)(to);
	auto to_m = (*mpl::m)(to_sz);
	auto to_n = (*mpl::n)(to_sz);
	
	BOOST_ASSERT((*mpl::less_equal)(from_m, to_m));
	BOOST_ASSERT((*mpl::equal)(from_n, to_n));
	
	for (std::size_t j = 0; j < from_n; ++j) {
		BOOST_ASSERT(from.data().at(j).density().empty());
		BOOST_ASSERT(from.data().at(j).pressure().empty());
		BOOST_ASSERT(from.data().at(j).residual().empty());
	
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i, j)) = from.data().at(j).x_velocity().at(i);
		}
		
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i+from_m/3, j)) = from.data().at(j).y_velocity().at(i);
		}
		
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i+from_m/3*2, j)) = from.data().at(j).z_velocity().at(i);
		}
	}
	
	return HBRS_MPL_FWD(to);
}

template<
	typename To,
	typename From,
	typename std::enable_if_t<
		(
			#ifdef HBRS_MPL_ENABLE_MATLAB
				std::is_same< hana::tag_of_t<From>, mpl::ml_matrix_tag >::value ||
			#endif
			#ifdef HBRS_MPL_ENABLE_ELEMENTAL
				std::is_same< hana::tag_of_t<From>, mpl::el_matrix_tag >::value ||
			#endif
			std::is_same< hana::tag_of_t<From>, mpl::rtsam_tag >::value ||
			std::is_same< hana::tag_of_t<From>, mpl::ctsam_tag >::value ||
			std::is_same< hana::tag_of_t<From>, mpl::sm_tag >::value
		) &&
		std::is_same_v<To, theta_field_matrix>
	>* = nullptr
>
decltype(auto)
copy_matrix(From const& from, To & to) {
	auto from_sz = (*mpl::size)(from);
	auto from_m = (*mpl::m)(from_sz);
	auto from_n = (*mpl::n)(from_sz);
	
	auto to_sz = to.size();
	auto to_m = (*mpl::m)(to_sz);
	auto to_n = (*mpl::n)(to_sz);
	
	BOOST_ASSERT((*mpl::greater_equal)(from_m, to_m));
	BOOST_ASSERT((*mpl::equal)(from_n, to_n));
	
	for (std::size_t j = 0; j < to_n; ++j) {
		to.data().at(j).density().clear();
		to.data().at(j).pressure().clear();
		to.data().at(j).residual().clear();
		BOOST_ASSERT(to.data().at(j).x_velocity().size() * 3 == to_m);
		BOOST_ASSERT(to.data().at(j).y_velocity().size() * 3 == to_m);
		BOOST_ASSERT(to.data().at(j).z_velocity().size() * 3 == to_m);
	
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to.data().at(j).x_velocity().at(i) = (*mpl::at)(from, mpl::make_matrix_index(i, j));
		}
		
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to.data().at(j).y_velocity().at(i) = (*mpl::at)(from, mpl::make_matrix_index(i+to_m/3, j));
		}
		
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to.data().at(j).z_velocity().at(i) = (*mpl::at)(from, mpl::make_matrix_index(i+to_m/3*2, j));
		}
	}
	
	return to;
}


/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_MATRIX_IMPL_HPP
