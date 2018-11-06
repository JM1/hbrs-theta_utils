/* Copyright (c) 2016 Jakob Meng, <jakobmeng@web.de>
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

#include <hbrs/theta_utils/dt/nc_cntr.hpp>
#include <hbrs/theta_utils/dt/nc_exception.hpp>
#include <boost/throw_exception.hpp>

#include <netcdf.h>
#include <algorithm>
#include <numeric>
#include <regex>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

nc_cntr::nc_cntr(
	std::vector<nc_dimension> dims,
	std::vector<nc_variable> vars
) : dimensions_{dims}, variables_{vars} {}

HBRS_THETA_UTILS_DEFINE_ATTR(dimensions, std::vector<nc_dimension>, nc_cntr)
HBRS_THETA_UTILS_DEFINE_ATTR(variables, std::vector<nc_variable>, nc_cntr)

#define __find_if(__cntr)                                                                                              \
	auto const it = std::find_if(                                                                                      \
		__cntr.begin(),                                                                                                \
		__cntr.end(),                                                                                                  \
		[&name](auto const& x) { return x.name() == name; }                                                            \
	);                                                                                                                 \
	if (it == __cntr.end()) { return {}; }                                                                             \
	else { return *it; }

boost::optional<nc_dimension &> nc_cntr::dimension(std::string const& name) { __find_if(dimensions_) }
boost::optional<nc_dimension const&> nc_cntr::dimension(std::string const& name) const { __find_if(dimensions_) }

boost::optional<nc_variable &> nc_cntr::variable(std::string const& name) { __find_if(variables_) }
boost::optional<nc_variable const&> nc_cntr::variable(std::string const& name) const { __find_if(variables_) }

nc_cntr read_nc_cntr(
	std::string const& path,
	std::vector<std::string> const& includes /*regex filter*/,
	std::vector<std::string> const& excludes /*regex filter*/
) {
	int ncid, ndims, nvars, status;
	
	std::vector<nc_dimension> dimensions;
	std::vector<nc_variable> variables;
	
	struct var {
		int id;
		std::string name;
		nc_type type;
		std::vector<int> dimids;
		int natts;
	};
	
	std::vector<var> vars;
	std::vector<bool> use_vars;
	
	auto const handle_error = [&ncid, &status, &path]() {
		if(status != NC_NOERR) {
			nc_close(ncid);
			BOOST_THROW_EXCEPTION(
				nc_exception{} \
				<< errinfo_nc_status(status) \
				<< boost::errinfo_file_name(path)
			);
		}
	};
	
	status = nc_open(path.data(), NC_NOWRITE, &ncid);
	if(status != NC_NOERR) { 
		BOOST_THROW_EXCEPTION(
			nc_exception{} \
			<< errinfo_nc_status(status) \
			<< boost::errinfo_file_name(path)
		);
	}
	
	status = nc_inq_ndims(ncid, &ndims);
	handle_error();
	dimensions.reserve(ndims);
	{
		std::array<char, NC_MAX_NAME+1> name;
		size_t length;
		for(int dimid = 0; dimid < ndims; ++dimid) {
			status = nc_inq_dim(ncid, dimid, name.data(), &length);
			handle_error();
			dimensions.push_back({name.data(), length});
		}
	}
	
	status = nc_inq_nvars(ncid, &nvars);
	handle_error();
	vars.reserve(nvars);
	{
		std::array<char, NC_MAX_NAME+1> name;
		nc_type type;
		int ndims;
		int dimids[NC_MAX_VAR_DIMS];
		int natts;
		
		for(int i = 0; i < nvars; ++i) {
			status = nc_inq_var(ncid, i, name.data(), &type, &ndims, dimids, &natts);
			handle_error();
			vars.push_back({i, name.data(), type, {dimids, dimids+ndims}, natts});
		}
	}
	
	if (includes.empty()) {
		use_vars.resize(nvars, true);
	} else {
		use_vars.resize(nvars, false);
		for(auto const& include : includes) {
			std::regex regex{include};
			for(int i = 0; i < nvars; ++i) {
				if (use_vars[i] == false) {
					use_vars[i] = std::regex_search(vars[i].name, regex);
				}
			}
		}
	}
	
	if (excludes.size()) {
		for(auto const& exclude : excludes) {
			std::regex regex{exclude};
			for(int i = 0; i < nvars; ++i) {
				if (use_vars[i] == true) {
					use_vars[i] = !std::regex_search(vars[i].name, regex);
				}
			}
		}
	}
	
	auto const total = [&dimensions](std::vector<int> const& dimids) {
		return std::accumulate(dimids.begin(), dimids.end(), 1, 
			[&dimensions](std::size_t const& length, int const& dim) {
				return length * dimensions[dim].length();
			}
		);
	};
	
	auto const get_dims = [&dimensions](std::vector<int> const& dimids) {
		return std::accumulate(dimids.begin(), dimids.end(), std::vector<nc_dimension>{}, 
			[&dimensions](std::vector<nc_dimension> dims, int const& dim) {
				dims.push_back(dimensions[dim]);
				return dims;
			}
		);
	};
	
	variables.reserve(vars.size());
#define __nc_type_case(__nc_type, __type)                                                                              \
	if (var.type == __nc_type) {                                                                                       \
		std::vector<__type> data(total(var.dimids));                                                                   \
		status = nc_get_var(ncid, var.id, data.data());                                                                \
		variables.push_back({var.name, get_dims(var.dimids), {data}});                                                 \
		handle_error();                                                                                                \
	} else

	for(int i = 0; i < nvars; ++i) {
		if (!use_vars[i]) { continue; }
		
		auto const& var = vars[i];
		
		__nc_type_case(NC_DOUBLE, double)
		__nc_type_case(NC_FLOAT, float)
		__nc_type_case(NC_INT, int)
		__nc_type_case(NC_LONG, long)
		__nc_type_case(NC_SHORT, short)
		__nc_type_case(NC_CHAR, char)
		__nc_type_case(NC_UINT, unsigned int)
		__nc_type_case(NC_USHORT, unsigned short)
		{
			status = NC_EBADTYPID;
			handle_error();
		}
	}
	
	status = nc_close(ncid);
	if (status != NC_NOERR) { 
		BOOST_THROW_EXCEPTION(
			nc_exception{} \
			<< errinfo_nc_status(status) \
			<< boost::errinfo_file_name(std::string{path})
		);
	}
	
	return nc_cntr{dimensions, variables};
}

HBRS_THETA_UTILS_NAMESPACE_END