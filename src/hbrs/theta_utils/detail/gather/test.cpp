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

#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/mpl/dt/ctsav.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

#include <hbrs/mpl/detail/test.hpp>
#include <hbrs/mpl/detail/log.hpp>
#include <hbrs/mpl/detail/mpi.hpp>

#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>

#include <hbrs/theta_utils/detail/gather.hpp>
#include <hbrs/theta_utils/detail/scatter.hpp>
#include <hbrs/theta_utils/detail/matrix.hpp>
#include <hbrs/theta_utils/detail/test.hpp>

#include <hbrs/theta_utils/dt/theta_field_matrix.hpp>


#include <boost/hana/filter.hpp>
#include <boost/hana/zip.hpp>
#include <boost/hana/first.hpp>
#include <boost/hana/second.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana/at.hpp>
#include <boost/hana/plus.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/transform.hpp>
#include <boost/hana/cartesian_product.hpp>
#include <boost/hana/drop_back.hpp>
#include <boost/hana/drop_front.hpp>
#include <boost/hana/front.hpp>
#include <boost/hana/back.hpp>
#include <boost/hana/unpack.hpp>
#include <boost/hana/greater_equal.hpp>
#include <boost/hana/range.hpp>
#include <boost/hana/length.hpp>
#include <boost/hana/mult.hpp>
#include <boost/hana/mod.hpp>
#include <boost/hana/functional/id.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/hana/less_equal.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

#define _TOL 0.000000001

BOOST_AUTO_TEST_SUITE(detail_gather_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(gather_,
	* utf::precondition(hbrs::theta_utils::detail::mpi_world_size_condition{{0,5}})
	* utf::tolerance(_TOL)
) {
	using namespace hbrs::mpl;
	namespace detail = hbrs::theta_utils::detail;
	namespace mpi = hbrs::mpl::detail::mpi;
	namespace hana = boost::hana;
	using namespace hbrs::theta_utils;
	using hbrs::mpl::detail::loggable;
	
	static constexpr auto datasets = hana::make_tuple(
		make_sm(
			make_ctsav(mpl::detail::mat_a), 
			make_matrix_size(hana::size_c<mpl::detail::mat_a_m>, hana::size_c<mpl::detail::mat_a_n>),
			row_major_c
		),
		make_sm(
			make_ctsav(mpl::detail::mat_c),
			make_matrix_size(hana::size_c<mpl::detail::mat_c_m>, hana::size_c<mpl::detail::mat_c_n>),
			row_major_c
		),
		make_sm(
			make_ctsav(mpl::detail::mat_e),
			make_matrix_size(hana::size_c<mpl::detail::mat_e_m>, hana::size_c<mpl::detail::mat_e_n>),
			row_major_c
		),
		make_sm(
			make_ctsav(mpl::detail::mat_l),
			make_matrix_size(hana::size_c<mpl::detail::mat_l_m>, hana::size_c<mpl::detail::mat_l_n>),
			row_major_c
		)
	);
	
	static constexpr auto dimensions = hana::make_tuple(
		hana::to_tuple(hana::make_range(hana::size_c<0>, hana::length(datasets))),
		hana::make_tuple(
			#ifdef HBRS_MPL_ENABLE_ELEMENTAL
				hana::type_c<el_dist_matrix_tag>
			#endif // !HBRS_MPL_ENABLE_ELEMENTAL
		),
		hana::make_tuple(
			hana::type_c<detail::theta_field_distribution_1>, hana::type_c<detail::theta_field_distribution_2>
		)
	);
	
	hana::for_each(hana::cartesian_product(dimensions), [](auto const& cfg) {
		auto const& dataset_nr = hana::at_c<0>(cfg);
		auto const& dist_tag = hana::at_c<1>(cfg);
		auto const& dist_alg = hana::at_c<2>(cfg);
		BOOST_TEST_MESSAGE("dataset_nr=" << dataset_nr);
		
		auto const& dataset = hana::at(datasets, dataset_nr);
		
		using dist_tag_t = typename std::decay_t<decltype(dist_tag)>::type;
		using dist_alg_t = typename std::decay_t<decltype(dist_alg)>::type;
		
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
		if constexpr(std::is_same_v<dist_tag_t, el_dist_matrix_tag>) {
			BOOST_TEST_MESSAGE("elemental_mpi");
			if constexpr(
				std::is_same_v<dist_alg_t, detail::theta_field_distribution_1>
			) {
				/* equally-sized submatrices */
				theta_field_matrix series = make_theta_field_matrix(
					rtsam<double, storage_order::row_major>{(*size)(dataset)}
				);
				detail::copy_matrix(make_el_matrix(dataset), series);
				
				static El::Grid grid{}; // grid is static because reference to grid is required by
				auto dist_vc_star =
					detail::scatter(series, detail::scatter_control<dist_alg_t>{{}});
				
				static_assert(std::is_same_v<
					decltype(dist_vc_star),
					el_dist_matrix<double, El::VC, El::STAR>
				>, "");
				
				HBRS_MPL_LOG_TRIVIAL(trace) << "dist_vc_star:" << loggable{dist_vc_star};
				
				/* compare complete matrix */
				el_matrix<double> single = detail::copy_matrix(
					series,
					make_el_matrix(
						hana::type_c<double>,
						matrix_size<El::Int, El::Int>{ series.size() }
					)
				);
				
				el_matrix<double> concat = make_el_matrix(
					hana::type_c<double>,
					matrix_size<El::Int, El::Int>{
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
						dist_alg_t,
						matrix_size<std::size_t, std::size_t>
					>{{}, series.size()}
				);
				el_matrix<double> local = detail::copy_matrix(
					series,
					make_el_matrix(
						hana::type_c<double>,
						matrix_size<El::Int, El::Int>{ series.size() }
					)
				);
				
				HBRS_MPL_LOG_TRIVIAL(trace) << "local:" << loggable{local};
				
				HBRS_MPL_TEST_MMEQ(dataset, local, false);
			}
			
			if constexpr(
					std::is_same_v<dist_alg_t, detail::theta_field_distribution_1> ||
					std::is_same_v<dist_alg_t, detail::theta_field_distribution_2>
				) {
				/* differently-sized submatrices */
				el_matrix<double> untouched = make_el_matrix(dataset);
				
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
				
				el_matrix<double> truncated = {rand_m, untouched.size().n()};
				truncated.data() = untouched.data()(
					El::IR(rand_offset,rand_offset+rand_m),
					El::ALL
				);
				
				HBRS_MPL_LOG_TRIVIAL(trace) << "truncated at mpi_rank " << mpi::comm_rank() << ":" << loggable{truncated};
				theta_field_matrix series = make_theta_field_matrix(
					rtsam<double, storage_order::row_major>{(*size)(truncated)}
				);
				detail::copy_matrix(truncated, series);
				
				static El::Grid grid{}; // grid is static because reference to grid is required by
				
				el_dist_matrix<double, El::VC, El::STAR> dist_vc_star =
					detail::scatter(series, detail::scatter_control<dist_alg_t>{{}});
				
				HBRS_MPL_LOG_TRIVIAL(trace) << "dist_vc_star:" << loggable{dist_vc_star};
				
				/* compare local matrix */
				series = detail::gather(
					dist_vc_star,
					detail::gather_control<
						dist_alg_t,
						matrix_size<std::size_t, std::size_t>
					>{{}, series.size()}
				);
				el_matrix<double> local = detail::copy_matrix(
					series,
					make_el_matrix(
						hana::type_c<double>,
						matrix_size<El::Int, El::Int>{ series.size() }
					)
				);
				HBRS_MPL_LOG_TRIVIAL(trace) << "local:" << loggable{local};
				HBRS_MPL_TEST_MMEQ(truncated, local, false);
			}
		}
		#endif // !HBRS_MPL_ENABLE_ELEMENTAL
	});
}

BOOST_AUTO_TEST_SUITE_END()
