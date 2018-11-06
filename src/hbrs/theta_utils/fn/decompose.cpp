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

#include <hbrs/theta_utils/fn/decompose.hpp>

#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/dt/theta_grid.hpp>
#include <hbrs/theta_utils/dt/number_sequence.hpp>
#include <hbrs/theta_utils/fn/pca_filter.hpp>
#include <hbrs/theta_utils/fn/vtk.hpp>
#include <boost/filesystem.hpp>

#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/format.hpp>
#include <boost/system/error_code.hpp>
#include <boost/assert.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/combine.hpp>
#include <boost/iostreams/device/file.hpp>
#include <sstream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace io = boost::iostreams;

static
fs::path
make_vtk_output_path(
	theta_field_path input_path, 
	std::string const& output_prefix, 
	std::string const& tag, 
	bool simple_numbering,
	vtk_file_format frmt
) {
	input_path.prefix() = output_prefix;
	fs::path folder = input_path.folder();
	
	fs::path filename;
	if (simple_numbering) {
		std::stringstream fn;
		fn  << input_path.prefix() << ".pval";
		
		if (input_path.domain_num()) {
			fn << ".domain_" << *input_path.domain_num();
		}
		
		fn << "." << tag;
		fn << ".vtk." << input_path.step();
		
		filename = { fn.str() };
	} else {
		filename = { input_path.filename().string() + '.' + tag + '.' + vtk_file_extension(frmt) };
	}
	
	return folder / filename;
}

static
fs::path
make_pvd_output_path(fs::path const& output_folder, std::string const& output_prefix, std::string const& output_tag) {
	return output_folder / (output_prefix + '.' + output_tag + ".pvd");
}

static fs::path
make_stats_output_path(
	fs::path const& output_path,
	std::string const& output_prefix,
	std::string const& output_tag
) {
	return output_path / (output_prefix + '.' + output_tag + ".stat");
}

static std::vector<std::string>
make_pc_nr_tags(std::vector<std::string> const& nr_seqs, std::size_t max_pc_nr) {
	if (nr_seqs.empty()) {
		return { "all" };
	}

	std::vector<std::string> tags;
	tags.reserve(nr_seqs.size());
	
	for(auto tag : nr_seqs) {
		boost::replace_all(tag, "first", "1");
		boost::replace_all(tag, "last", boost::lexical_cast<std::string>(max_pc_nr));
		
		tags.push_back(tag);
	}
	
	return tags;
}

static std::vector<std::vector<bool>>
parse_pc_nr_sequences(std::vector<std::string> const& nr_seqs, std::size_t max_pc_nr) {
	if (nr_seqs.empty()) {
		return { std::vector<bool>(max_pc_nr, true) };
	}

	std::vector<std::vector<bool>> includes_seq;
	includes_seq.reserve(nr_seqs.size());
	
	for(auto nr_seq : nr_seqs) {
		boost::replace_all(nr_seq, "first", "1");
		boost::replace_all(nr_seq, "last", boost::lexical_cast<std::string>(max_pc_nr));
		std::vector<unsigned long> range;
		bool valid = parse_number_sequence(nr_seq.begin(), nr_seq.end(), range);
		
		if (!valid || range.empty() || (range.front() == 0) || (range.back() > max_pc_nr)) {
			BOOST_THROW_EXCEPTION(invalid_number_range_spec_exception{} << errinfo_number_range_spec{nr_seq});
		}
		
		std::vector<bool> includes(max_pc_nr, false);
		BOOST_ASSERT(includes.size() == max_pc_nr);
		
		for (auto && nr : range) {
			includes[nr-1] = true;
		}
		
		includes_seq.push_back(includes);
	}
	
	return includes_seq;
}

static void
write_stats(
	std::vector<double> const& latent,
	fs::path const& output_filename
) {
	io::stream_buffer<io::file_sink> buf(output_filename.string(), BOOST_IOS::binary);
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

static void
decompose_with_pca(
	theta_grid const& grid,
	std::vector<std::tuple<theta_field, theta_field_path>> const& input_series,
	std::vector<bool> const& includes,
	pca_backend const& backend,
	fs::path const& output_path,
	std::string const& output_prefix,
	std::string const& output_tag,
	vtk_file_format const& output_format,
	bool const& simple_numbering,
	bool const& overwrite
) {
	BOOST_ASSERT(includes.size() == input_series.size());
	
	std::vector<theta_field> output_fields;
	output_fields.reserve(input_series.size());
	std::vector<fs::path> output_field_paths;
	output_field_paths.reserve(input_series.size());
	
	for(auto && fnp : input_series) {
		auto && field = std::get<0>(fnp);
		auto && input_path = std::get<1>(fnp);
		
		output_fields.push_back(std::move(field));
		
		fs::path output_filename = make_vtk_output_path(input_path, output_prefix, output_tag, simple_numbering, output_format);
		
		if (fs::exists(output_filename) && !overwrite) {
			BOOST_THROW_EXCEPTION((
				fs::filesystem_error{
					(boost::format("output file %s already exists in folder %s") % output_filename.string() % output_path).str(),
					make_error_code(boost::system::errc::file_exists)
				}
			));
		}
		
		output_field_paths.push_back(std::move(output_filename));
	}
	
	auto pvd_path = make_pvd_output_path(output_path, output_prefix, output_tag);
	if (fs::exists(pvd_path) && !overwrite) {
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("output file %s already exists in folder %s") % pvd_path.string() % output_path.string()).str(),
				make_error_code(boost::system::errc::file_exists)
			}
		));
	}
	
	auto stats_path = make_stats_output_path(output_path, output_prefix, output_tag);
	if (fs::exists(stats_path) && !overwrite) {
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("stats file %s already exists in folder %s") % stats_path.string() % output_path).str(),
				make_error_code(boost::system::errc::file_exists)
			}
		));
	}
	
	mpl::pca_filter_result<
		std::vector<theta_field>,
		std::vector<double>
	> reduced;
	
	switch (backend) {
		case pca_backend::matlab_lapack:
			reduced = pca_filter(std::move(output_fields), includes, matlab_lapack_backend_c);
			break;
		case pca_backend::elemental_openmp:
			reduced = pca_filter(std::move(output_fields), includes, elemental_openmp_backend_c);
			break;
		case pca_backend::elemental_mpi:
			reduced = pca_filter(std::move(output_fields), includes, elemental_mpi_backend_c);
			break;
		default:
			BOOST_ASSERT_MSG(false, "unknown pca backend");
	};
	
	output_fields = std::move(reduced.data());
	auto && latent = reduced.latent();
	
	BOOST_ASSERT((*mpl::size)(latent) == output_fields.size());
	
	std::vector<std::tuple<theta_field, fs::path>> output_series;
	output_series.reserve(input_series.size());
	
	for(std::size_t i = 0; i < input_series.size(); ++i) {
		output_series.push_back({ std::move(output_fields[i]), output_field_paths[i] });
	}
	
	write_vtk(grid, output_series, output_format);
	if (output_format == vtk_file_format::xml_binary) {
		write_pvd(output_field_paths, pvd_path.string());
	}
	write_stats(latent, stats_path);
}

void
decompose_with_pca(
	visualize_options v_opts,
	pca_decomposition_options pca_opts,
	bool verbose
) {
	auto grid = read_theta_grid(v_opts.path, v_opts.grid_prefix);
	auto input_series = read_theta_fields(v_opts.path, v_opts.input_prefix, v_opts.includes, v_opts.excludes);
	
	auto includes_seqs = parse_pc_nr_sequences(pca_opts.pc_nr_seqs, input_series.size());
	auto tags = make_pc_nr_tags(pca_opts.pc_nr_seqs, input_series.size());
	
	BOOST_ASSERT(includes_seqs.size() == tags.size());
	
	for(auto && pc_nr_spec : boost::combine(includes_seqs, tags)) {
		auto && includes = boost::get<0>(pc_nr_spec);
		auto && tag = boost::get<1>(pc_nr_spec);
		
		BOOST_ASSERT(includes.size() == input_series.size());
		decompose_with_pca(grid, input_series, includes, pca_opts.backend, { v_opts.path },
			v_opts.output_prefix, tag, v_opts.output_format, v_opts.simple_numbering, v_opts.overwrite);
	}
}



HBRS_THETA_UTILS_NAMESPACE_END
