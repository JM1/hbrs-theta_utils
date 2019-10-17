/* Copyright (c) 2016-2019 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_DETAIL_MAKE_MATRIX_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_MAKE_MATRIX_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/detail/copy_matrix.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/detail/mpi.hpp>

#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_MATLAB
    #include <hbrs/mpl/dt/ml_matrix.hpp>
#endif
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif

#include <boost/numeric/conversion/cast.hpp>
#include <boost/hana/core/tag_of.hpp>

#include <type_traits>
#include <vector>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace mpl = hbrs::mpl;
namespace mpi = hbrs::mpl::detail::mpi;

namespace detail {

#if defined(HBRS_MPL_ENABLE_MATLAB) || defined(HBRS_MPL_ENABLE_ELEMENTAL)
	template<
		typename From,
		typename std::enable_if_t<
			#ifdef HBRS_MPL_ENABLE_MATLAB
				std::is_same< hana::tag_of_t<From>, mpl::ml_column_vector_tag >::value ||
			#endif
			#ifdef HBRS_MPL_ENABLE_ELEMENTAL
				std::is_same< hana::tag_of_t<From>, mpl::el_column_vector_tag >::value ||
			#endif
			false
		>* = nullptr
	>
	auto
	to_vector(From const& from) {
		std::size_t from_sz = boost::numeric_cast<std::size_t>((*mpl::size)(from));
		std::vector<double> to;
		to.resize(from_sz);
		BOOST_ASSERT(to.size() == from_sz);
		
		for (std::size_t i = 0; i < from_sz; ++i) {
			to.at(i) = boost::numeric_cast<double>((*mpl::at)(from, i));
		}
		
		return to;
	}
#endif

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	template<typename Ring, El::Dist Columnwise, El::Dist Rowwise, El::DistWrap Wrapping>
	auto
	to_vector(mpl::el_dist_column_vector<Ring, Columnwise, Rowwise, Wrapping> const& from) {
		typedef std::decay_t<Ring> _Ring_;
		
		El::DistMatrixReadProxy<Ring, _Ring_, El::STAR, El::STAR, Wrapping> from_pxy = {from.data()};
		
		auto from_local = from_pxy.GetLocked();
		std::size_t from_sz = boost::numeric_cast<std::size_t>(from_local.Height());
		BOOST_ASSERT(from_local.Width() == 1);
		BOOST_ASSERT(from_local.Height() == from.data().Height());
		
		std::vector<double> to;
		to.resize(from_sz);
		BOOST_ASSERT(to.size() == from_sz);
		
		for (std::size_t i = 0; i < from_sz; ++i) {
			to.at(i) = boost::numeric_cast<double>(from_local.Get(i, 0));
		}
		
		return to;
	}
#endif

#ifdef HBRS_MPL_ENABLE_MATLAB
	auto
	make_matrix(hana::basic_type<mpl::ml_matrix_tag>, mpl::matrix_size<int, int> sz) {
		return mpl::make_ml_matrix(hana::type_c<double>, sz);
	}
#endif

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	auto
	make_matrix(hana::basic_type<mpl::el_matrix_tag>, mpl::matrix_size<El::Int, El::Int> sz) {
		return mpl::make_el_matrix(hana::type_c<double>, sz);
	}

	auto
	make_matrix(hana::basic_type<mpl::el_dist_matrix_tag>, mpl::matrix_size<El::Int, El::Int> sz) {
		static El::Grid const grid{El::mpi::COMM_WORLD};
		
		auto min_m = El::mpi::AllReduce(sz.m(), El::mpi::MIN, grid.Comm());
		auto max_m = El::mpi::AllReduce(sz.m(), El::mpi::MAX, grid.Comm());
		BOOST_ASSERT(min_m <= max_m || min_m > max_m);
		auto full_m = max_m * mpi::size(); //TODO: Is this correct?

	// 	auto full_m = El::mpi::AllReduce(sz.m(), El::mpi::SUM, grid.Comm());
	// 	BOOST_ASSERT(min_m * mpi::size() == full_m);	
		
		auto min_n = El::mpi::AllReduce(sz.n(), El::mpi::MIN, grid.Comm());
		auto max_n = El::mpi::AllReduce(sz.n(), El::mpi::MAX, grid.Comm());
		BOOST_ASSERT(min_n == max_n);
		
		El::DistMatrix<double, El::VC, El::STAR, El::ELEMENT> a{grid};
		a.Resize(full_m, sz.n());
		
		BOOST_ASSERT(a.Matrix().Width() == sz.n());
		return mpl::make_el_dist_matrix(a);
	}
#endif

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_MAKE_MATRIX_IMPL_HPP
