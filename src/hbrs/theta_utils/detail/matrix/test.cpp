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

#define BOOST_TEST_MODULE detail_matrix_test
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

#ifdef HBRS_MPL_ENABLE_MATLAB
    #include <hbrs/mpl/dt/ml_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_MATLAB
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

#include <hbrs/mpl/detail/test.hpp>
#include <hbrs/mpl/detail/log.hpp>
#include <hbrs/mpl/detail/mpi.hpp>

#include <hbrs/theta_utils/detail/matrix.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/dt/theta_field_matrix.hpp>

#include <boost/hana/tuple.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/less_equal.hpp>
#include <boost/test/results_collector.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <random>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(detail_matrix_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(matrix_, * utf::tolerance(0.000000001)) {
	using namespace hbrs::theta_utils;
	using namespace hbrs::theta_utils::detail;
	using hbrs::mpl::detail::loggable;
	namespace mpi = hbrs::mpl::detail::mpi;
	
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
		),
		mpl::make_sm(
			mpl::make_ctsav(mpl::detail::mat_l),
			mpl::make_matrix_size(hana::size_c<mpl::detail::mat_l_m>, hana::size_c<mpl::detail::mat_l_n>),
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
				#endif // !HBRS_MPL_ENABLE_MATLAB
				#if defined(HBRS_MPL_ENABLE_MATLAB) && defined(HBRS_MPL_ENABLE_ELEMENTAL)
					,
				#endif //!( defined(HBRS_MPL_ENABLE_MATLAB) && defined(HBRS_MPL_ENABLE_ELEMENTAL) )
				#ifdef HBRS_MPL_ENABLE_ELEMENTAL
					hana::type_c<mpl::el_matrix_tag>
				#endif // !HBRS_MPL_ENABLE_ELEMENTAL
			),
			[&dataset](auto tag) {
				BOOST_ASSERT(3u <= (*mpl::m)(mpl::size(dataset)));
				
				#ifdef HBRS_MPL_ENABLE_MATLAB
					if constexpr(decltype(tag){} == hana::type_c<mpl::ml_matrix_tag>) {
						if (mpi::comm_size() > 1 && mpi::comm_rank() != MPI_ROOT) { return; }
						BOOST_TEST_MESSAGE("matlab_lapack");
						theta_field_matrix series = make_theta_field_matrix(
							mpl::rtsam<double, mpl::storage_order::row_major>{(*mpl::size)(dataset)}
						);
						detail::copy_matrix(mpl::make_ml_matrix(dataset), series);
						
						auto r = detail::copy_matrix(
							series,
							mpl::make_ml_matrix(
								hana::type_c<double>,
								mpl::matrix_size<int, int>{ series.size() }
							)
						);
						
						HBRS_MPL_LOG_TRIVIAL(trace) << "r:" << loggable{r};
						
						HBRS_MPL_TEST_MMEQ(dataset, r, false);
					} else 
				#endif // !HBRS_MPL_ENABLE_MATLAB
				#ifdef HBRS_MPL_ENABLE_ELEMENTAL
					if constexpr(decltype(tag){} == hana::type_c<mpl::el_matrix_tag>) {
						if (mpi::comm_size() > 1 && mpi::comm_rank() != MPI_ROOT) { return; }
						BOOST_TEST_MESSAGE("elemental_openmp");
						theta_field_matrix series = make_theta_field_matrix(
							mpl::rtsam<double, mpl::storage_order::row_major>{(*mpl::size)(dataset)}
						);
						detail::copy_matrix(mpl::make_el_matrix(dataset), series);
						
						auto r = detail::copy_matrix(
							series,
							mpl::make_el_matrix(
								hana::type_c<double>,
								mpl::matrix_size<El::Int, El::Int>{ series.size() }
							)
						);
						
						HBRS_MPL_LOG_TRIVIAL(trace) << "r:" << loggable{r};
						
						HBRS_MPL_TEST_MMEQ(dataset, r, false);
					} else 
				#endif // !HBRS_MPL_ENABLE_ELEMENTAL
				{
					BOOST_ASSERT(false);
				}
			}
		);
		++dataset_nr;
	});
	#endif //!( defined(HBRS_MPL_ENABLE_MATLAB) || defined(HBRS_MPL_ENABLE_ELEMENTAL) )
}

BOOST_AUTO_TEST_SUITE_END()
