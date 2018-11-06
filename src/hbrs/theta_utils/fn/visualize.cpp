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

#include <hbrs/theta_utils/fn/visualize.hpp>

#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/dt/theta_grid.hpp>
#include <hbrs/theta_utils/fn/vtk.hpp>
#include <boost/filesystem.hpp>

#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/format.hpp>
#include <boost/system/error_code.hpp>

#include <sstream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace fs = boost::filesystem;

static
fs::path
make_vtk_output_path(
	theta_field_path input_path,
	std::string const& output_prefix,
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
		
		fn << ".vtk." << input_path.step();
		
		filename = { fn.str() };
	} else {
		filename = { input_path.filename().string() + '.' + vtk_file_extension(frmt) };
	}
	
	return folder / filename;
}

static
fs::path
make_pvd_output_path(fs::path const& output_folder, std::string const& output_prefix) {
	return output_folder / (output_prefix + ".pvd");
}

void
visualize(
	visualize_options v_opts,
	bool verbose
) {
	auto grid = read_theta_grid(v_opts.path, v_opts.grid_prefix);
	auto input_series = read_theta_fields(v_opts.path, v_opts.input_prefix, v_opts.includes, v_opts.excludes);
	
	std::vector<std::tuple<theta_field, fs::path>> output_series;
	output_series.reserve(input_series.size());
	
	for(auto && fnp : input_series) {
		auto && field = std::get<0>(fnp);
		auto && input_path = std::get<1>(fnp);
		fs::path output_path = make_vtk_output_path(input_path, v_opts.output_prefix, v_opts.simple_numbering, v_opts.output_format);
		if (fs::exists(output_path) && !v_opts.overwrite) {
			BOOST_THROW_EXCEPTION((
				fs::filesystem_error{
					(boost::format("output file %s already exists in folder %s") % output_path.string() % v_opts.path).str(),
					make_error_code(boost::system::errc::file_exists)
				}
			));
		}
		
		output_series.push_back({std::move(field), std::move(output_path)});
	}
	
	write_vtk(grid, output_series, v_opts.output_format);
	if (v_opts.output_format == vtk_file_format::xml_binary) {
		auto pvd_file = make_pvd_output_path({ v_opts.path }, v_opts.output_prefix);
		if (fs::exists(pvd_file) && !v_opts.overwrite) {
			BOOST_THROW_EXCEPTION((
				fs::filesystem_error{
					(boost::format("output file %s already exists in folder %s") % pvd_file.string() % v_opts.path).str(),
					make_error_code(boost::system::errc::file_exists)
				}
			));
		}
		std::vector<fs::path> output_paths;
		output_paths.reserve(output_series.size());
		for(auto && fnp : output_series) {
			output_paths.push_back(std::get<1>(fnp));
		}
		
		write_pvd(output_paths, pvd_file.string());
	}
}

HBRS_THETA_UTILS_NAMESPACE_END