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

#include "../impl.hpp"

#include <hbrs/theta_utils/dt/command.hpp>
#include <hbrs/theta_utils/dt/command_option.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/detail/int_ranges.hpp>
#include <hbrs/theta_utils/detail/make_matrix.hpp>
#include <hbrs/mpl/dt/pca_filter_result.hpp>
#include <hbrs/mpl/dt/pca_control.hpp>

#include <hbrs/mpl/fn/pca_filter.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/zip.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>

#include <boost/numeric/conversion/cast.hpp>
#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/format.hpp>
#include <boost/system/error_code.hpp>
#include <boost/assert.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/combine.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <sstream>
#include <limits>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace mpi = hbrs::mpl::detail::mpi;
namespace io = boost::iostreams;

namespace {

fs::path
make_stats_output_path(
	fs::path const& output_folder,
	std::string const& output_prefix,
	std::string const& output_tag,
	boost::optional<int> const& domain_num
) {
	return 
		output_folder /
			(output_prefix + '_' + output_tag + ".stat" + 
				(domain_num ? std::string{".domain_"} + boost::lexical_cast<std::string>(*domain_num) : "") 
			);
}

typedef boost::integer_range<std::size_t> number_range;
typedef std::vector<number_range> number_ranges;
typedef std::vector<number_ranges> number_ranges_sequence;

number_ranges_sequence
parse_nr_sequences(std::vector<std::string> const& rngs_seq) {
	static constexpr auto min = 0;
	static constexpr auto max = std::numeric_limits<std::size_t>::max();
	
	if (rngs_seq.empty()) {
		return { number_ranges{ number_range{min, max} } };
	}

	number_ranges_sequence nr_rngs_seq;
	nr_rngs_seq.reserve(rngs_seq.size());
	
	for(auto rngs : rngs_seq) {
		boost::replace_all(rngs, "first", boost::lexical_cast<std::string>(min));
		boost::replace_all(rngs, "last", boost::lexical_cast<std::string>(max));
		number_ranges nr_rngs;
		bool valid = detail::parse_int_ranges(rngs.begin(), rngs.end(), nr_rngs);
		
		if (valid) {
			for(auto & nr_rng : nr_rngs) {
				if (nr_rng.empty() || (*nr_rng.begin() < min)) {
					valid = false;
					break;
				}
			}
		}
		
		if (!valid) {
			BOOST_THROW_EXCEPTION(invalid_number_range_spec_exception{} << errinfo_number_range_spec{rngs});
		}
		
		nr_rngs_seq.push_back(nr_rngs);
	}
	
	return nr_rngs_seq;
}

void
write_stats(
	std::vector<double> const& latent,
	fs::path const& file_path
) {
	io::stream_buffer<io::file_sink> buf(file_path.string(), BOOST_IOS::binary);
	std::ostream out(&buf);
	
	out << "Nr of principal components: " << latent.size() << std::endl;
	
	double sum = std::accumulate(latent.begin(), latent.end(), 0., std::plus<>{});
	double cumsum{0.};
	std::size_t i{0};
	out << "PC#     |    Variance | Variance/Sum (%) |   Cumsum/Sum (%)" << std::endl;
	for(auto && variance : latent) {
		cumsum += variance;
		out 
			<< boost::format("PC %05i:  %10.05f - %15.4f%% - %15.4f%%")
				% (i+1) % variance % (variance/sum*100.) % (cumsum/sum*100.)
			<< std::endl;
			
		++i;
	}
}

// template<typename T1, typename T2>
// auto
// unzip2(std::vector<std::tuple<T1, T2>> && zipped) {
// 	std::vector<T1> unzipped1;
// 	unzipped1.reserve(zipped.size());
// 	
// 	std::vector<T2> unzipped2;
// 	unzipped2.reserve(zipped.size());
// 	
// 	for(auto && [ zipped1, zipped2 ] : zipped) {
// 		unzipped1.push_back(std::move(HBRS_MPL_FWD(zipped1)));
// 		unzipped2.push_back(std::move(HBRS_MPL_FWD(zipped2)));
// 	}
// 	
// 	return std::make_tuple(unzipped1, unzipped2);
// }

#if defined(HBRS_MPL_ENABLE_MATLAB) || defined(HBRS_MPL_ENABLE_ELEMENTAL)
	template<typename ToTag>
	auto
	convert_to(std::vector<theta_field> const& series, hana::basic_type<ToTag>) {
		return detail::copy_matrix(
			series,
			detail::make_matrix(hana::type_c<ToTag>, detail::size(series))
		);
	}

	template<typename From, typename To>
	mpl::pca_filter_result<
		std::vector<theta_field> /* data */,
		std::vector<double> /* latent*/
	>
	convert_from(From from, To && to) {
		return {
			detail::copy_matrix(from.data(), to),
			detail::to_vector(from.latent())
		};
	}
	
	template<typename Matrix>
	decltype(auto)
	reduce(
		Matrix && data,
		std::function<bool(std::size_t)> keep,
		mpl::pca_control<bool,bool,bool> ctrl
	) {
		// NOTE: in our data matrix, rows correspond to variables and columns correspond to observations,
		//       but in pca it is vice versa.
		decltype(auto) r = mpl::pca_filter(mpl::transpose(HBRS_MPL_FWD(data)), keep, ctrl);
		
		return mpl::make_pca_filter_result(
			mpl::transpose(HBRS_MPL_FWD(r).data()),
			HBRS_MPL_FWD(r).latent()
		);
	}
	
	template<typename ToTag>
	decltype(auto)
	reduce(
		std::vector<theta_field> series,
		hana::basic_type<ToTag>,
		std::function<bool(std::size_t)> keep,
		mpl::pca_control<bool,bool,bool> ctrl
	) {
		return convert_from(
			reduce(
				convert_to(series, hana::type_c<ToTag>), 
				keep, 
				ctrl
			),
			series
		);
	}
#endif

void
decompose_with_pca(
	std::vector<theta_field_path> paths,
	detail::int_ranges<std::size_t> const& includes,
	pca_backend const& backend,
	bool center,
	bool normalize,
	fs::path const& output_folder,
	std::string const& output_prefix,
	std::string const& output_tag,
	bool const& overwrite
) {
	BOOST_ASSERT(paths.size() > 0);
	
	std::vector<theta_field> series = read_theta_fields(paths, {".*_velocity"});
	
	boost::optional<int> domain_num = paths[0].domain_num();
	BOOST_ASSERT(series.at(0).ndomains() == mpi::size()); //TODO: Turn assertion into exception?
	// TODO: Warn user that his theta_field was computed with n domains but his current number of mpi processes is different!
	
	// test for existing files before starting decomposition
	std::vector<theta_field_path> output_paths = paths;
	for(auto & path : output_paths) {
		BOOST_ASSERT(path.domain_num() == domain_num);
		
		// transform input paths to output paths
		path.folder() = output_folder;
		path.prefix() = output_prefix + '_' + output_tag;
		
		if (fs::exists(path.full_path()) && !overwrite) {
			BOOST_THROW_EXCEPTION((
				fs::filesystem_error{
					(boost::format("output file %s already exists in folder %s") % path.filename().string() % path.full_path().parent_path().string()).str(),
					make_error_code(boost::system::errc::file_exists)
				}
			));
		}
	}
	
	auto stats_path = make_stats_output_path(output_folder, output_prefix, output_tag, domain_num);
	if (fs::exists(stats_path) && !overwrite) {
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("stats file %s already exists in folder %s") % stats_path.string() % stats_path.parent_path().string()).str(),
				make_error_code(boost::system::errc::file_exists)
			}
		));
	}
	
	mpl::pca_filter_result<
		std::vector<theta_field>,
		std::vector<double>
	> reduced;
	
	BOOST_ASSERT(mpi::size() > 1 ? backend == pca_backend::elemental_mpi : true);
	
	std::function<bool(std::size_t)> keep = [&includes](std::size_t i) {
		return detail::in_int_ranges(includes, i); 
	};
	
	mpl::pca_control<bool,bool,bool> ctrl {true /* economy */, center, normalize};
	
	switch (backend) {
		case pca_backend::matlab_lapack:
			#ifdef HBRS_MPL_ENABLE_MATLAB
				reduced = reduce(std::move(series), hana::type_c<mpl::ml_matrix_tag>, keep, ctrl);
			#else
				BOOST_THROW_EXCEPTION(invalid_backend_exception{} << errinfo_pca_backend{matlab_lapack_backend_c});
			#endif
			break;
		case pca_backend::elemental_openmp:
			#ifdef HBRS_MPL_ENABLE_ELEMENTAL
				reduced = reduce(std::move(series), hana::type_c<mpl::el_matrix_tag>, keep, ctrl);
			#else
				BOOST_THROW_EXCEPTION(invalid_backend_exception{} << errinfo_pca_backend{elemental_openmp_backend_c});
			#endif
			break;
		case pca_backend::elemental_mpi:
			#ifdef HBRS_MPL_ENABLE_ELEMENTAL
				reduced = reduce(std::move(series), hana::type_c<mpl::el_dist_matrix_tag>, keep, ctrl);
			#else
				BOOST_THROW_EXCEPTION(invalid_backend_exception{} << errinfo_pca_backend{elemental_mpi_backend_c});
			#endif
			break;
		default:
			BOOST_ASSERT_MSG(false, "unknown pca backend");
	};
	
	#if !defined(NDEBUG)
	{
		auto latent_sz = (*mpl::size)(reduced.latent());
		auto data_sz = detail::size(reduced.data());
		auto data_n = (*mpl::n)(data_sz);
		//NOTE: pca was applied to transposed data matrix
		auto DOF = data_n - (ctrl.center() ? 1 : 0);
		BOOST_ASSERT(latent_sz == DOF);
	}
	#endif
	
	{
		// we need global_id field if distributed, e.g. for visualization
		std::vector<theta_field> global_ids = read_theta_fields(paths, {"global_id"});
		BOOST_ASSERT(reduced.data().size() == global_ids.size());
		BOOST_ASSERT(mpi::size() > 1
			? global_ids[0].global_id().size() > 0
			: global_ids[0].global_id().size() == 0
		);
		
		for(std::size_t i = 0; i < reduced.data().size(); ++i) {
			auto & src = global_ids[i].global_id();
			auto & tgt = reduced.data()[i].global_id();
			
			BOOST_ASSERT(mpi::size() > 1
				? src.size() == reduced.data()[i].x_velocity().size() /* distributed */
				: src.size() == 0
			);
			
			tgt = std::move(src);
		}
	}
	write_theta_fields(
		mpl::detail::zip_impl_std_tuple_vector{}(std::move(reduced.data()), std::move(output_paths)),
		overwrite
	);
	write_stats(std::move(reduced.latent()), stats_path);
}

/* unnamed namespace */ }

void
execute(pca_cmd cmd) {
	BOOST_ASSERT(mpi::initialized());
	
	auto paths = filter_theta_fields_by_domain_num(
		find_theta_fields(cmd.i_opts.path, cmd.i_opts.pval_prefix),
		mpi::size() > 1
			? boost::optional<int>{mpi::rank()}
			: boost::optional<int>{boost::none}
	);
	
	if (paths.empty()) { 
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("no theta field files (*.pval.*) with prefix %s found in folder %s") % cmd.i_opts.pval_prefix % cmd.i_opts.path).str(),
				boost::system::errc::make_error_code(boost::system::errc::no_such_file_or_directory)
			}
		));
	}
	
	{
		auto const& first = paths[0];
		for(auto const& path : paths) {
			if (path.naming_scheme() != first.naming_scheme()) {
				BOOST_THROW_EXCEPTION((
					ambiguous_naming_scheme_exception{}
					<< errinfo_ambiguous_field_paths{{first.full_path(), path.full_path()}}
				));
			}
		}
	}
	
	auto includes_seqs = parse_nr_sequences(cmd.pca_opts.pc_nr_seqs);
	auto tags = cmd.pca_opts.pc_nr_seqs.empty() == false ? cmd.pca_opts.pc_nr_seqs : std::vector<std::string>{"all"};
	BOOST_ASSERT(includes_seqs.size() == tags.size());
	
	for(auto && [ includes, tag ] : mpl::detail::zip_impl_std_tuple_vector{}(std::move(includes_seqs), std::move(tags))) {
		decompose_with_pca(
			paths,
			HBRS_MPL_FWD(includes),
			cmd.pca_opts.backend,
			cmd.pca_opts.center,
			cmd.pca_opts.normalize,
			{ cmd.o_opts.path },
			cmd.o_opts.prefix,
			HBRS_MPL_FWD(tag),
			cmd.o_opts.overwrite
		);
	}
}

HBRS_THETA_UTILS_NAMESPACE_END
