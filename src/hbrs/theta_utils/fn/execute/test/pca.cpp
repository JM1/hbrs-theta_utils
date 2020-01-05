/* Copyright (c) 2018-2020 Jakob Meng, <jakobmeng@web.de>
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

#define BOOST_TEST_MODULE fn_execute_pca_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/detail/test.hpp>
#include <hbrs/mpl/detail/gather.hpp>
#include <hbrs/mpl/detail/not_supported.hpp>

#include <hbrs/theta_utils/dt/command.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/fn/execute.hpp>
#include <hbrs/theta_utils/detail/test.hpp>
#include <hbrs/theta_utils/detail/matrix.hpp>
#include <hbrs/mpl/fn/zip.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/plus.hpp>
#include <hbrs/mpl/fn/expand.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/mean.hpp>
#include <hbrs/mpl/fn/columns.hpp>
#include <hbrs/mpl/fn/transpose.hpp>

#include <hbrs/mpl/dt/sm.hpp>
#include <hbrs/mpl/dt/ctsav.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/dt/ctsam.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

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

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

#define _TOL 0.000000001

BOOST_AUTO_TEST_SUITE(fn_execute_pca_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(write_read,
	* utf::precondition(hbrs::theta_utils::detail::mpi_world_size_condition{{0,3}})
	* utf::tolerance(_TOL)
) {
	using namespace hbrs::mpl;
	namespace detail = hbrs::mpl::detail;
	namespace mpi = hbrs::mpl::detail::mpi;
	namespace hana = boost::hana;
	using namespace hbrs::theta_utils;
	
	static constexpr auto datasets = hana::make_tuple(
		make_sm(
			make_ctsav(detail::mat_k), make_matrix_size(hana::size_c<detail::mat_k_m>, hana::size_c<detail::mat_k_n>), row_major_c
		),
		make_sm(
			make_ctsav(detail::mat_l), make_matrix_size(hana::size_c<detail::mat_l_m>, hana::size_c<detail::mat_l_n>), row_major_c
		),
		make_sm(
			make_ctsav(detail::mat_m), make_matrix_size(hana::size_c<detail::mat_m_m>, hana::size_c<detail::mat_m_n>), row_major_c
		),
		make_sm(
			make_ctsav(detail::mat_p), make_matrix_size(hana::size_c<detail::mat_p_m>, hana::size_c<detail::mat_p_n>), row_major_c
		)
	);
	
	static constexpr auto dimensions = hana::make_tuple(
		hana::to_tuple(hana::make_range(hana::size_c<0>, hana::length(datasets))),
		hana::make_tuple(theta_field_path::naming_scheme::theta, theta_field_path::naming_scheme::tau_unsteady),
		hana::make_tuple(hana::true_c, hana::false_c) /* center */,
		hana::make_tuple(hana::true_c, hana::false_c) /* normalize */,
		hana::make_tuple(hana::true_c, hana::false_c) /* keep_centered */
	);
	
	static constexpr auto factories = hana::drop_back(hana::make_tuple(
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
		hana::make_pair(
			/* scatter */[](auto && dataset) {
				typedef decltype(dataset.at(matrix_index<std::size_t, std::size_t>{0u,0u})) Ring;
				typedef std::decay_t<Ring> _Ring_;
				
				static El::Grid grid {}; // grid is static because reference to grid is required by El::DistMatrix<...>
				
				el_dist_matrix<_Ring_, El::STAR, El::STAR> dist_star_star = make_el_dist_matrix(grid, make_el_matrix(dataset));
				el_dist_matrix<_Ring_, El::VC,   El::STAR> dist_vc_star   = {dist_star_star.data()};
				
				return make_el_matrix(dist_vc_star.data().Matrix());
			},
			/* gather */ [](auto && dataset, auto && dist_size) {
				typedef decltype(dataset.at(matrix_index<std::size_t, std::size_t>{0u,0u})) Ring;
				typedef std::decay_t<Ring> _Ring_;
				
				static El::Grid grid {}; // grid is static because reference to grid is required by El::DistMatrix<...>
				
				El::DistMatrix<_Ring_, El::VC, El::STAR> dist_vc_star{
					boost::numeric_cast<El::Int>(dist_size.m().value),
					boost::numeric_cast<El::Int>(dist_size.n().value),
					grid
				};
				dist_vc_star.Matrix() = make_el_matrix(dataset).data();
				
				El::DistMatrix<_Ring_, El::STAR, El::STAR> dist_star_star{dist_vc_star};
				return make_el_matrix(dist_star_star.Matrix());
			}
		),
		#endif // !HBRS_MPL_ENABLE_ELEMENTAL
		"SEQUENCE_TERMINATOR___REMOVED_BY_DROP_BACK"
	));
	
	hana::for_each(hana::cartesian_product(dimensions), [](auto const& cfg) {
		auto const& dataset_nr = hana::at_c<0>(cfg);
		auto const& scheme = hana::at_c<1>(cfg);
		auto const& center = hana::at_c<2>(cfg);
		auto const& normalize = hana::at_c<3>(cfg);
		auto const& keep_centered = hana::at_c<4>(cfg);
		
		BOOST_TEST_MESSAGE("dataset_nr=" << dataset_nr);
		BOOST_TEST_MESSAGE(
			"scheme=" <<
				((scheme == theta_field_path::naming_scheme::theta) ? "theta" : "tau_unsteady")
		);
		BOOST_TEST_MESSAGE("center=" << (center ? "true" : "false"));
		BOOST_TEST_MESSAGE("normalize=" << (normalize ? "true" : "false"));
		BOOST_TEST_MESSAGE("keep_centered=" << (keep_centered ? "true" : "false"));
		
		auto const& dataset = hana::at(datasets, dataset_nr);
		
		auto sz_ = (*size)(dataset);
		auto m_ = (*m)(sz_);
		auto n_ = (*n)(sz_);
		
		if (m_%(3u * boost::numeric_cast<std::size_t>(mpi::comm_size())) != 0u) {
			BOOST_TEST_MESSAGE(
				"Skipping dataset_nr=" << dataset_nr 
				<< " because incompatible for executing with " << mpi::comm_size() << " MPI processes"
			);
			return;
		}
		
		auto testcases = hana::transform(
			factories,
			hana::id
		);
		
		auto supported_indices = hana::transform(
			hana::filter(
				hana::zip(
					hana::transform(testcases, detail::is_supported), // tuple of hana::true_ and hana::false_
					hana::to_tuple(hana::make_range(hana::size_c<0>,hana::length(testcases))) // indices of testcases
				),
				hana::front
			),
			hana::back
		);
		
		auto results = hana::transform(
			supported_indices,
			[&](auto i) {
				BOOST_TEST_MESSAGE("Running impl nr " << i);
				
				auto testcase = hana::at(testcases, i);
				
				auto scatter = hana::first(testcase);
				auto gather = hana::second(testcase);
				
				auto local_dataset = scatter(dataset);
				
				hbrs::theta_utils::detail::io_fixture fxo{"pca_output"};
				
				{
					/* When fxi goes out of scope, then this input directory is removed.
					 * This prevents buggy code from accidentally reading the input directory.
					 */
					hbrs::theta_utils::detail::io_fixture fxi{"pca_input"};
					BOOST_TEST_MESSAGE("PCA input directory: " << fxi.wd().path().string());
					BOOST_TEST_MESSAGE("PCA output directory: " << fxo.wd().path().string());
					
					auto domain_num = mpi::comm_size() > 1
						? boost::optional<int>{mpi::comm_rank()}
						: boost::optional<int>{boost::none};
					
					theta_field_matrix local_series = hbrs::theta_utils::make_theta_field_matrix(local_dataset);
					for(theta_field & field : local_series.data()) {
						if (mpi::comm_size() > 1) {
							field.global_id() = std::vector<int>(local_series.size().m()/3, 0);
						}
						field.ndomains() = mpi::comm_size();
					}
					
					std::vector<theta_field_path> pca_input_paths =
						hbrs::theta_utils::detail::make_theta_field_paths(
							fxi.wd().path(), fxi.prefix(), local_series, scheme
						);
					
					for(theta_field_path & path : pca_input_paths) {
						path.domain_num() = domain_num;
					}
					
					write_theta_fields(
						mpl::detail::zip_impl_std_tuple_vector{}(local_series.data(), pca_input_paths),
						false
					);
					
					pca_cmd cmd;
					cmd.i_opts.path = fxi.wd().path().string();
					cmd.i_opts.pval_prefix = fxi.prefix();
					cmd.o_opts.path = fxo.wd().path().string();
					cmd.o_opts.prefix = fxo.prefix();
					cmd.o_opts.overwrite = false;
					cmd.pca_opts.pc_nr_seqs = {/* all */};
					cmd.pca_opts.backend = pca_backend::elemental_mpi;
					cmd.pca_opts.center = center;
					cmd.pca_opts.normalize = normalize;
					cmd.pca_opts.keep_centered = keep_centered;
					execute(cmd);
				}
				
				auto all_paths = find_theta_fields(fxo.wd().path(), fxo.prefix() + "_all");
				auto paths = filter_theta_fields_by_domain_num(
					all_paths,
					mpi::comm_size() > 1
						? boost::optional<int>{mpi::comm_rank()}
						: boost::optional<int>{boost::none}
				);
				
				BOOST_TEST((*equal)(all_paths.size(), n_)); // because every MPI process creates its own temporary directory
				BOOST_TEST((*equal)(paths.size(), n_));
				
				theta_field_matrix local_series = theta_field_matrix{ read_theta_fields(paths) };
				rtsam<double, storage_order::row_major> local_series_as_matrix =
					hbrs::theta_utils::detail::copy_matrix(
						local_series, 
						rtsam<double, storage_order::row_major>{
							local_series.size()
						}
					);
				
				return gather(local_series_as_matrix, sz_);
			}
		);
		
		auto results_indices = hana::to_tuple(hana::make_range(
			hana::size_c<0>,
			hana::length(results)
		));
		
		BOOST_TEST_MESSAGE("All runs finished.");
		
		BOOST_TEST_PASSPOINT();
		
		hana::for_each(
			results_indices,
			[&](auto i) {
				auto impl_idx = hana::at(supported_indices, i);
				auto result = hana::at(results, i);
				
				BOOST_TEST_MESSAGE("Comparing original data and reconstructed data computed by impl nr " << impl_idx);
				
				if (center && keep_centered) {
					// if matrix was centered and mean was not readded after pca,
					// then we have to add mean now to be able to do a comparison
					auto testcase = hana::at(testcases, impl_idx);
					auto scatter = hana::first(testcase);
					auto gather = hana::second(testcase);
					auto dataset_ = gather(scatter(dataset), sz_);
					// matrix will be transposed before pca_cmd, thus
					// mean must be computed from transposed matrix.
					auto dataset_t = transpose(dataset_);
					auto transpose_t = transpose(result);
					auto result_w_mean_t = (*plus)(transpose_t, expand(mean(columns(dataset_t)), size(dataset_t)));
					auto result_w_mean = transpose(result_w_mean_t);
					HBRS_MPL_TEST_MMEQ(dataset, result_w_mean, false);
				} else {
					HBRS_MPL_TEST_MMEQ(dataset, result, false);
				}
			}
		);
		
		BOOST_TEST_MESSAGE("Comparing original and reconstructed datasets done.");
		
		BOOST_TEST_PASSPOINT();
	});
	
		
				
				
	
	#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	
	#else // HBRS_MPL_ENABLE_ELEMENTAL
	
	#endif // !HBRS_MPL_ENABLE_ELEMENTAL
	
}



BOOST_AUTO_TEST_SUITE_END()
