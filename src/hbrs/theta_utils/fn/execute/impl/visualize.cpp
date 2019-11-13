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
#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/dt/theta_grid.hpp>
#include <hbrs/theta_utils/detail/vtk.hpp>
#include <boost/filesystem.hpp>

#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/format.hpp>
#include <boost/system/error_code.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace fs = boost::filesystem;
namespace mpi = hbrs::mpl::detail::mpi;

void
execute(visualize_cmd cmd) {
	BOOST_ASSERT(mpi::initialized());
	
	boost::optional<theta_grid_path> grid_path = find_theta_grid(cmd.i_opts.path, cmd.i_opts.grid_prefix);
	if (!grid_path) {
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("No *.grid file with prefix %s found in folder %s") % cmd.i_opts.grid_prefix % cmd.i_opts.path).str(),
				boost::system::errc::make_error_code(boost::system::errc::no_such_file_or_directory)
			}
		));
	}
	
	std::vector<theta_field_path> field_paths = filter_theta_fields_by_domain_num(
		find_theta_fields(cmd.i_opts.path, cmd.i_opts.pval_prefix),
		mpi::comm_size() > 1
			? boost::optional<int>{mpi::comm_rank()}
			: boost::optional<int>{boost::none}
	);
	
	if (field_paths.empty() && mpi::comm_size() == 1) {
		// if simulation was split among several processes but only one process was choosen for visalization
		field_paths = filter_theta_fields_by_domain_num(
			find_theta_fields(cmd.i_opts.path, cmd.i_opts.pval_prefix),
			boost::optional<int>{mpi::comm_rank()}
		);
	}
	
	if (field_paths.empty()) {
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("No *.pval.* file with prefix %s found in folder %s") % cmd.i_opts.pval_prefix % cmd.i_opts.path).str(),
				boost::system::errc::make_error_code(boost::system::errc::no_such_file_or_directory)
			}
		));
	}
	
	{
		auto const& first = field_paths[0];
		for(auto const& path : field_paths) {
			if (path.naming_scheme() != first.naming_scheme()) {
				BOOST_THROW_EXCEPTION((
					ambiguous_naming_scheme_exception{}
					<< errinfo_ambiguous_field_paths{{first.full_path(), path.full_path()}}
				));
			}
		}
	}
	
	convert_to_vtk(
		*grid_path,
		field_paths,
		cmd.o_opts.path,
		cmd.o_opts.prefix,
		cmd.v_opts.includes,
		cmd.v_opts.excludes,
		cmd.v_opts.simple_numbering,
		cmd.v_opts.format,
		cmd.o_opts.overwrite);
}

HBRS_THETA_UTILS_NAMESPACE_END
