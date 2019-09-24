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

#define BOOST_TEST_MODULE make_matrix_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>

#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/mpl/dt/ctsav.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/detail/test.hpp>

#include <hbrs/theta_utils/detail/copy_matrix.hpp>
#include <hbrs/theta_utils/detail/make_theta_fields.hpp>
#include <hbrs/theta_utils/detail/make_matrix.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>

#include <boost/hana/tuple.hpp>
#include <boost/hana/for_each.hpp>
#include <vector>
#include <boost/test/results_collector.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(make_matrix_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(make_matrix_, * utf::tolerance(0.000000001)) {
	using namespace hbrs::theta_utils;
	
	static constexpr auto datasets = hana::make_tuple(
		mpl::make_sm(
			mpl::make_ctsav(mpl::detail::mat_a), 
			mpl::make_matrix_size(hana::size_c<mpl::detail::mat_a_m>, hana::size_c<mpl::detail::mat_a_n>),
			mpl::row_major_c
		),
		mpl::make_sm(
			mpl::make_ctsav(mpl::detail::mat_c),
			mpl::make_matrix_size(hana::size_c<mpl::detail::mat_c_m>, hana::size_c<mpl::detail::mat_c_n>),
			mpl::row_major_c
		),
		mpl::make_sm(
			mpl::make_ctsav(mpl::detail::mat_e),
			mpl::make_matrix_size(hana::size_c<mpl::detail::mat_e_m>, hana::size_c<mpl::detail::mat_e_n>),
			mpl::row_major_c
		)
	);
	
	#if defined(HBRS_MPL_ENABLE_MATLAB) || defined(HBRS_MPL_ENABLE_ELEMENTAL)
	std::size_t dataset_nr = 0;
	hana::for_each(datasets, [&dataset_nr](auto const& dataset) {
		BOOST_TEST_MESSAGE("dataset_nr=" << dataset_nr);
		
		hana::for_each(
			hana::make_tuple(
				#ifdef HBRS_MPL_ENABLE_MATLAB
					hana::type_c<mpl::ml_matrix_tag>
				#endif
				#if defined(HBRS_MPL_ENABLE_MATLAB) && defined(HBRS_MPL_ENABLE_ELEMENTAL)
					,
				#endif
				#ifdef HBRS_MPL_ENABLE_ELEMENTAL
					hana::type_c<mpl::el_matrix_tag>,
					hana::type_c<mpl::el_dist_matrix_tag>
				#endif
			),
			[&dataset](auto tag) {
				std::vector<theta_field> series = hbrs::theta_utils::detail::make_theta_fields(
					mpl::rtsam<double, mpl::storage_order::row_major>{(*mpl::size)(dataset)}
				);
				
				#ifdef HBRS_MPL_ENABLE_MATLAB
					if constexpr(decltype(tag){} == hana::type_c<mpl::ml_matrix_tag>) {
						BOOST_TEST_MESSAGE("matlab_lapack");
						hbrs::theta_utils::detail::copy_matrix(mpl::make_ml_matrix(dataset), series);
					} else 
				#endif
				#ifdef HBRS_MPL_ENABLE_ELEMENTAL
					if constexpr(decltype(tag){} == hana::type_c<mpl::el_matrix_tag>) {
						BOOST_TEST_MESSAGE("elemental_openmp");
						hbrs::theta_utils::detail::copy_matrix(mpl::make_el_matrix(dataset), series);
					} else if constexpr(decltype(tag){} == hana::type_c<mpl::el_dist_matrix_tag>) {
						BOOST_TEST_MESSAGE("elemental_mpi");
						hbrs::theta_utils::detail::copy_matrix(
							mpl::make_el_dist_matrix(
								El::Grid{El::mpi::COMM_WORLD},
								mpl::make_el_matrix(dataset)
							), 
							series);
					} else 
				#endif
				{
					BOOST_ASSERT(false);
				}
				
				auto sz_ = hbrs::theta_utils::detail::size(series);
				auto r = hbrs::theta_utils::detail::copy_matrix(
					series,
					hbrs::theta_utils::detail::make_matrix(tag, sz_)
				);
				
				HBRS_MPL_TEST_MMEQ(dataset, r, false);
			}
		);
		++dataset_nr;
	});
	#endif
}

BOOST_AUTO_TEST_SUITE_END()
