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
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/expand.hpp>
#include <hbrs/mpl/fn/columns.hpp>
#include <hbrs/mpl/fn/mean.hpp>

#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>

#include <hbrs/mpl/detail/test.hpp>
#include <hbrs/theta_utils/detail/copy_matrix.hpp>
#include <hbrs/theta_utils/detail/make_theta_fields.hpp>
#include <hbrs/theta_utils/fn/pca_filter.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/for_each.hpp>
#include <vector>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/test/results_collector.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(pca_filter_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(pca_filter_, * utf::tolerance(0.000000001)) {
	using namespace hbrs::theta_utils;
	
	std::vector< mpl::rtsam<double, mpl::storage_order::row_major> > datasets = {
		{
			{
				1,2,3,
				4,5,6,
				7,8,9
			},
			mpl::matrix_size<std::size_t, std::size_t>{3u,3u}
		},
		{
			{
				1,2,3,4,5,6,
				7,8,9,10,11,12,
				13,14,15,16,17,18
			},
			mpl::matrix_size<std::size_t, std::size_t>{3u,6u}
		},
		{
			{
				1,2,3,
				4,5,6,
				7,8,9,
				10,11,12,
				13,14,15,
				16,17,18
			},
			mpl::matrix_size<std::size_t, std::size_t>{6u,3u}
		}
	};
	
	std::size_t dataset_nr = 0;
	for(auto const& dataset : datasets) {
		BOOST_TEST_MESSAGE("dataset_nr=" << dataset_nr);
		
		for(auto && keep : { true, false }) {
			BOOST_TEST_MESSAGE("keep=" << (keep ? "true" : "false"));
			hana::for_each(
				hana::make_tuple(matlab_lapack_backend_c, elemental_openmp_backend_c, elemental_mpi_backend_c),
				[&dataset, &keep](auto && backend) {
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
					
					auto ds_sz = (*mpl::size)(dataset);
					auto ds_m = (*mpl::m)(ds_sz);
					auto ds_n = (*mpl::n)(ds_sz);
					
					hbrs::theta_utils::detail::int_ranges<std::size_t> filter;
					if (keep) {
						filter.push_back({ 0, ds_m-1<ds_n ? ds_m-1 : std::min(ds_m, ds_n) });
					}
					
					auto r = pca_filter(hbrs::theta_utils::detail::make_theta_fields(dataset), filter, backend);
					
					auto & data = (*mpl::at)(r, mpl::pca_filter_data{});
// 					auto & latent = (*mpl::at)(r, mpl::pca_filter_latent{});
					
					auto data_sz = hbrs::theta_utils::detail::size(data);
					auto data_m = (*mpl::m)(data_sz);
					auto data_n = (*mpl::n)(data_sz);
					auto data2_m = boost::numeric_cast<El::Int>(data_m);
					auto data2_n = boost::numeric_cast<El::Int>(data_n);
					auto data2 = hbrs::theta_utils::detail::copy_matrix(data, El::Matrix<double>{data2_m,data2_n});
					
					if (keep) {
						HBRS_MPL_TEST_MMEQ(dataset, data2, false);
					} else {
						auto mean_ = (*mpl::expand)(mpl::mean(mpl::columns(dataset)), mpl::size(dataset));
						HBRS_MPL_TEST_MMEQ(mean_, data2, false);
					}
				}
			);
		}
		
		++dataset_nr;
	}
}

BOOST_AUTO_TEST_SUITE_END()
