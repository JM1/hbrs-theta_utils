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

#define BOOST_TEST_MODULE detail_gather_test
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
#include <hbrs/mpl/detail/log.hpp>
#include <hbrs/mpl/detail/mpi.hpp>

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

#include <hbrs/theta_utils/detail/gather.hpp>
#include <hbrs/theta_utils/detail/scatter.hpp>
#include <hbrs/theta_utils/detail/matrix.hpp>
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

BOOST_AUTO_TEST_SUITE(detail_gather_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(gather_, * utf::tolerance(0.000000001)) {
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
	
	#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	std::size_t dataset_nr = 0;
	hana::for_each(datasets, [&dataset_nr](auto const& dataset) {
		BOOST_TEST_MESSAGE("dataset_nr=" << dataset_nr);
		
		hana::for_each(
			hana::make_tuple(
				#ifdef HBRS_MPL_ENABLE_ELEMENTAL
					hana::type_c<mpl::el_dist_matrix_tag>
				#endif // !HBRS_MPL_ENABLE_ELEMENTAL
			),
			[&dataset](auto tag) {
				BOOST_ASSERT(3u <= (*mpl::m)(mpl::size(dataset)));
				
				#ifdef HBRS_MPL_ENABLE_ELEMENTAL
					if constexpr(decltype(tag){} == hana::type_c<mpl::el_dist_matrix_tag>) {
						BOOST_TEST_MESSAGE("elemental_mpi");
						{
							/* equally-sized submatrices */
							theta_field_matrix series = make_theta_field_matrix(
								mpl::rtsam<double, mpl::storage_order::row_major>{(*mpl::size)(dataset)}
							);
							detail::copy_matrix(mpl::make_el_matrix(dataset), series);
							
							static El::Grid const grid{El::mpi::COMM_WORLD};
							auto dist_vc_star =
								detail::scatter(series, detail::scatter_control<detail::theta_field_distribution_1>{{}});
							
							static_assert(std::is_same_v<
								decltype(dist_vc_star),
								mpl::el_dist_matrix<double, El::VC, El::STAR>
							>, "");
							
							HBRS_MPL_LOG_TRIVIAL(trace) << "dist_vc_star:" << loggable{dist_vc_star};
							
							/* compare complete matrix */
							mpl::el_matrix<double> single = detail::copy_matrix(
								series,
								mpl::make_el_matrix(
									hana::type_c<double>,
									mpl::matrix_size<El::Int, El::Int>{ series.size() }
								)
							);
							
							mpl::el_matrix<double> concat = mpl::make_el_matrix(
								hana::type_c<double>,
								mpl::matrix_size<El::Int, El::Int>{
									single.data().Height() * mpi::comm_size(),
									single.data().Width()
								}
							);
							for(El::Int i = 0; i < mpi::comm_size(); ++i) {
								auto view = concat.data()(
									El::IR(i*single.data().Height(), (i+1)*single.data().Height()),
									El::ALL
								);
								El::Copy(single.data(), view);
							}
							
							HBRS_MPL_LOG_TRIVIAL(trace) << "concat:" << loggable{concat};
							
							HBRS_MPL_TEST_MMEQ(concat, dist_vc_star, false);
							
							/* compare local matrix */
							series = detail::gather(
								dist_vc_star,
								detail::gather_control<
									detail::theta_field_distribution_1,
									mpl::matrix_size<std::size_t, std::size_t>
								>{{}, series.size()}
							);
							mpl::el_matrix<double> local = detail::copy_matrix(
								series,
								mpl::make_el_matrix(
									hana::type_c<double>,
									mpl::matrix_size<El::Int, El::Int>{ series.size() }
								)
							);
							
							HBRS_MPL_LOG_TRIVIAL(trace) << "local:" << loggable{local};
							
							HBRS_MPL_TEST_MMEQ(dataset, local, false);
						}
						{
							/* differently-sized submatrices */
							mpl::el_matrix<double> untouched = mpl::make_el_matrix(dataset);
							
							HBRS_MPL_LOG_TRIVIAL(trace) << "untouched at mpi_rank " << mpi::comm_rank() << ":" << loggable{untouched};
							
							std::random_device rd;
							std::mt19937 mt(rd());
							
							std::uniform_int_distribution<El::Int> m_dist(3., untouched.size().m());
							El::Int rand_m = m_dist(mt);
							rand_m += 3-(rand_m%3);
							if (rand_m > untouched.size().m()) {
								rand_m = untouched.size().m();
							}
							
							std::uniform_int_distribution<El::Int> offset_dist(0., (untouched.size().m() - rand_m));
							El::Int rand_offset = offset_dist(mt);
							
							mpl::el_matrix<double> truncated = {rand_m, untouched.size().n()};
							truncated.data() = untouched.data()(
								El::IR(rand_offset,rand_offset+rand_m),
								El::ALL
							);
							
							HBRS_MPL_LOG_TRIVIAL(trace) << "truncated at mpi_rank " << mpi::comm_rank() << ":" << loggable{truncated};
							theta_field_matrix series = make_theta_field_matrix(
								mpl::rtsam<double, mpl::storage_order::row_major>{(*mpl::size)(truncated)}
							);
							detail::copy_matrix(truncated, series);
							
							static El::Grid const grid{El::mpi::COMM_WORLD};
							mpl::el_dist_matrix<double, El::VC, El::STAR> dist_vc_star =
								detail::scatter(series, detail::scatter_control<detail::theta_field_distribution_1>{{}});
							
							HBRS_MPL_LOG_TRIVIAL(trace) << "dist_vc_star:" << loggable{dist_vc_star};
							
							/* compare local matrix */
							series = detail::gather(
								dist_vc_star,
								detail::gather_control<
									detail::theta_field_distribution_1,
									mpl::matrix_size<std::size_t, std::size_t>
								>{{}, series.size()}
							);
							mpl::el_matrix<double> local = detail::copy_matrix(
								series,
								mpl::make_el_matrix(
									hana::type_c<double>,
									mpl::matrix_size<El::Int, El::Int>{ series.size() }
								)
							);
							
							HBRS_MPL_LOG_TRIVIAL(trace) << "local:" << loggable{local};
							
							HBRS_MPL_TEST_MMEQ(truncated, local, false);
						}
					} else 
				#endif // !HBRS_MPL_ENABLE_ELEMENTAL
				{
					BOOST_ASSERT(false);
				}
			}
		);
		++dataset_nr;
	});
	#endif // !HBRS_MPL_ENABLE_ELEMENTAL
}

BOOST_AUTO_TEST_SUITE_END()
