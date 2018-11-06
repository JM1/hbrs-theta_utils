/* Copyright (c) 2016-2018 Jakob Meng, <jakobmeng@web.de>
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

#define BOOST_TEST_MODULE pca_filter_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/expand.hpp>
#include <hbrs/mpl/fn/columns.hpp>
#include <hbrs/mpl/fn/mean.hpp>

#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>

#include <hbrs/theta_utils/detail/copy_matrix.hpp>
#include <hbrs/theta_utils/fn/pca_filter.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/for_each.hpp>
#include <vector>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/test/results_collector.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

#define _TOL 0.000000001
#define _REQUIRE_TEST_PASSED BOOST_REQUIRE(boost::unit_test::results_collector.results(boost::unit_test::framework::current_test_case().p_id).passed());

BOOST_AUTO_TEST_SUITE(pca_filter_test)

BOOST_AUTO_TEST_CASE(pca_filter_, * utf::tolerance(_TOL)) {
	using namespace hbrs::theta_utils;
	
	auto make_theta_field = [](mpl::rtsam<std::size_t, mpl::storage_order::row_major> from) {
		auto sz = from.size();
		
		auto to = std::vector<theta_field>(
			sz.n(), 
			theta_field{
				{} /*density*/,
				std::vector<double>(sz.m()/3) /*x_velocity*/,
				{} /*x_velocity_old*/,
				std::vector<double>(sz.m()/3) /*y_velocity*/,
				{} /*y_velocity_old*/,
				std::vector<double>(sz.m()/3) /*z_velocity*/,
				{} /*z_velocity_old*/,
				{} /*pressure*/,
				{} /*pressure_old*/,
				{} /*residual*/
			}
		);
		
		BOOST_ASSERT(sz.m()%3 == 0);
		return detail::copy_matrix(from, to);
	};
	
	std::vector< std::vector<theta_field> > datasets = {
		make_theta_field({
			{
				1,2,3,
				4,5,6,
				7,8,9
			},
			mpl::matrix_size<std::size_t, std::size_t>{3u,3u}
		}),
		make_theta_field({
			{
				1,2,3,4,5,6,
				7,8,9,10,11,12,
				13,14,15,16,17,18
			},
			mpl::matrix_size<std::size_t, std::size_t>{3u,6u}
		}),
		make_theta_field({
			{
				1,2,3,
				4,5,6,
				7,8,9,
				10,11,12,
				13,14,15,
				16,17,18
			},
			mpl::matrix_size<std::size_t, std::size_t>{6u,3u}
		})
	};
	
	std::size_t dataset_nr = 0;
	for(auto && dataset : datasets) {
		BOOST_TEST_MESSAGE("dataset_nr=" << dataset_nr);
		auto ds_sz = detail::size(dataset);
		
		for(auto && keep : { true, false }) {
			BOOST_TEST_MESSAGE("keep=" << (keep ? "true" : "false"));
			hana::for_each(
				hana::make_tuple(matlab_lapack_backend_c, elemental_openmp_backend_c /*, elemental_mpi_backend_c*/), //TODO: Use once implemented!
				[&dataset, &ds_sz, &keep](auto && backend) {
					switch(backend.value) {
						case pca_backend::matlab_lapack:
							BOOST_TEST_MESSAGE("pca_backend=matlab_lapack");
							break;
						case pca_backend::elemental_openmp:
							BOOST_TEST_MESSAGE("pca_backend=elemental_openmp");
							break;
						case pca_backend::elemental_mpi:
							BOOST_TEST_MESSAGE("pca_backend=elemental_mpi");
							break;
						default:
							BOOST_ASSERT(false);
					}
					
					auto ds_m = (*mpl::m)(ds_sz);
					auto ds_n = (*mpl::n)(ds_sz);
					
					auto r = pca_filter(
						dataset,
						std::vector<bool>(ds_m-1<ds_n ? ds_m-1 : std::min(ds_m, ds_n), keep),
						backend
					);
					auto & data = (*mpl::at)(r, mpl::pca_filter_data{});
// 					auto & latent = (*mpl::at)(r, mpl::pca_filter_latent{});
					
					auto data_sz = detail::size(data);
					auto data_m = (*mpl::m)(data_sz);
					auto data_n = (*mpl::n)(data_sz);
					BOOST_TEST((*mpl::equal)(ds_sz, data_sz));
// 					_REQUIRE_TEST_PASSED
					
					auto data2_m = boost::numeric_cast<int>(data_m);
					auto data2_n = boost::numeric_cast<int>(data_n);
					auto data2 = detail::copy_matrix(data, El::Matrix<double>{data2_m,data2_n});
					auto data2_sz = (*mpl::size)(data2);
					BOOST_TEST((*mpl::equal)(ds_sz, data2_sz));
// 					_REQUIRE_TEST_PASSED
					
					El::Matrix<double> reference = detail::copy_matrix(dataset, El::Matrix<double>{data2_m,data2_n});
					if (!keep) {
						reference = (*mpl::expand)(mpl::mean(mpl::columns(reference)), mpl::size(reference));
					}
					
					for(std::size_t i = 0; i < data_m; ++i) {
						for(std::size_t j = 0; j < data_n; ++j) {
							auto x = (*mpl::at)(reference, mpl::make_matrix_index(i, j));
							auto y = (*mpl::at)(data2, mpl::make_matrix_index(i, j));
							BOOST_TEST(x == y, "@[" << i << "][" << j << "] := Left: " << x << " Right: " << y);
						}
					}
// 					_REQUIRE_TEST_PASSED
				}
			);
		}
		
		++dataset_nr;
	}
}

BOOST_AUTO_TEST_SUITE_END()
