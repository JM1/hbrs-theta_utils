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
	std::vector<nc_variable> vars,
	std::vector<nc_attribute> attrs
) : dimensions_{dims}, variables_{vars}, attributes_{attrs} {
	// each dimension a variable uses must be in dims
	for(nc_variable const& var : vars) {
		for(nc_dimension const& dim : var.dimensions()) {
			auto equals_dim = [&dim](auto const& other){ return dim == other; };
			if (std::none_of(dims.begin(), dims.end(), equals_dim)) {
				BOOST_THROW_EXCEPTION(
					nc_exception{}
					<< errinfo_nc_status(NC_EBADDIM)
				);
			}
		}
	}
}

HBRS_THETA_UTILS_DEFINE_ATTR(dimensions, std::vector<nc_dimension>, nc_cntr)
HBRS_THETA_UTILS_DEFINE_ATTR(variables, std::vector<nc_variable>, nc_cntr)
HBRS_THETA_UTILS_DEFINE_ATTR(attributes, std::vector<nc_attribute>, nc_cntr)

#define __find_if(__cntr)                                                                                              \
	auto const it = std::find_if(                                                                                      \
		__cntr.begin(),                                                                                                \
		__cntr.end(),                                                                                                  \
		[&name](auto const& x) { return x.name() == name; }                                                            \
	);                                                                                                                 \
	if (it == __cntr.end()) { return {}; }                                                                             \
	else { return *it; }

boost::optional<nc_dimension &>
nc_cntr::dimension(std::string const& name) {
	__find_if(dimensions_)
}

boost::optional<nc_dimension const&>
nc_cntr::dimension(std::string const& name) const {
	__find_if(dimensions_) 
}

boost::optional<nc_variable &>
nc_cntr::variable(std::string const& name) {
	__find_if(variables_)
}

boost::optional<nc_variable const&>
nc_cntr::variable(std::string const& name) const {
	__find_if(variables_)
}

boost::optional<nc_attribute &>
nc_cntr::attribute(std::string const& name) {
	__find_if(attributes_)
}

boost::optional<nc_attribute const&>
nc_cntr::attribute(std::string const& name) const {
	__find_if(attributes_)
}

#undef __find_if

namespace {

void
throw_if_error(int ncid, int status, std::string path, bool close_fd) {
	if(status != NC_NOERR) {
		if (close_fd) {
			nc_close(ncid);
		}
		BOOST_THROW_EXCEPTION(
			nc_exception{} \
			<< errinfo_nc_status(status) \
			<< boost::errinfo_file_name(path)
		);
	}
}

/* unnamed namespace */ }

nc_cntr
read_nc_cntr(
	std::string const& path,
	std::vector<std::string> const& includes /*regex filter*/,
	std::vector<std::string> const& excludes /*regex filter*/
) {
	int ncid, ndims, nvars, ngatts, status;
	
	std::vector<nc_dimension> dimensions;
	std::vector<nc_variable> variables;
	std::vector<nc_attribute> attributes;
	
	struct var {
		int id;
		std::string name;
		nc_type type;
		std::vector<int> dimids;
		int natts;
	};
	
	std::vector<var> vars;
	std::vector<bool> use_vars;
	
	status = nc_open(path.data(), NC_NOWRITE, &ncid);
	throw_if_error(ncid, status, path, false);
	
	status = nc_inq_ndims(ncid, &ndims);
	throw_if_error(ncid, status, path, true);
	dimensions.reserve(ndims);
	{
		std::array<char, NC_MAX_NAME+1> name;
		size_t length;
		for(int dimid = 0; dimid < ndims; ++dimid) {
			status = nc_inq_dim(ncid, dimid, name.data(), &length);
			throw_if_error(ncid, status, path, true);
			dimensions.push_back({name.data(), length});
		}
	}
	
	status = nc_inq_nvars(ncid, &nvars);
	throw_if_error(ncid, status, path, true);
	vars.reserve(nvars);
	{
		std::array<char, NC_MAX_NAME+1> name;
		nc_type type;
		int ndims;
		int dimids[NC_MAX_VAR_DIMS];
		int natts;
		
		for(int i = 0; i < nvars; ++i) {
			status = nc_inq_var(ncid, i, name.data(), &type, &ndims, dimids, &natts);
			throw_if_error(ncid, status, path, true);
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
		throw_if_error(ncid, status, path, true);                                                                      \
		variables.push_back({var.name, get_dims(var.dimids), {data}});                                                 \
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
			throw_if_error(ncid, status, path, true);
		}
	}
	
#undef __nc_type_case

#define __nc_type_case(__nc_type, __type)                                                                              \
	if (xtype == __nc_type) {                                                                                          \
		std::vector<__type> value(len);                                                                                \
		status = nc_get_att (ncid, NC_GLOBAL, name.data(), value.data());                                              \
		throw_if_error(ncid, status, path, true);                                                                      \
		attributes.push_back({name.data(), value});                                                                    \
	} else

	status = nc_inq_natts(ncid, &ngatts);
	throw_if_error(ncid, status, path, true);
	attributes.reserve(ngatts);
	{
		std::array<char, NC_MAX_NAME+1> name;
		size_t length;
		for(int gattrid = 0; gattrid < ngatts; ++gattrid) {
			name.fill(0);
			status = nc_inq_attname(ncid, NC_GLOBAL, gattrid, name.data());
			throw_if_error(ncid, status, path, true);
			
			nc_type xtype;
			size_t len;
			status = nc_inq_att(ncid, NC_GLOBAL, name.data(), &xtype, &len);
			throw_if_error(ncid, status, path, true);
			
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
				throw_if_error(ncid, status, path, true);
			}
		}
	}
#undef __nc_type_case

	status = nc_close(ncid);
	throw_if_error(ncid, status, path, false);
	
	
	return nc_cntr{dimensions, variables, attributes};
}

struct nc_type_visitor : public boost::static_visitor<std::optional<nc_type>> {
	#define __nc_type_case(__nc_type, __type)                                                                          \
		std::optional<nc_type>                                                                                         \
		operator()(std::vector<__type> const&) const {                                                                 \
			return {__nc_type};                                                                                        \
		}
	
	__nc_type_case(NC_DOUBLE, double)
	__nc_type_case(NC_FLOAT, float)
	__nc_type_case(NC_INT, int)
	__nc_type_case(NC_LONG, long)
	__nc_type_case(NC_SHORT, short)
	__nc_type_case(NC_CHAR, char)
	__nc_type_case(NC_UINT, unsigned int)
	__nc_type_case(NC_USHORT, unsigned short)
	
	#undef __nc_type_case
	
	template<typename T>
	std::optional<nc_type>
	operator()(T const&) const {
		return {};
	}
};

struct array_ptr_visitor : public boost::static_visitor<void const *> {
	template<typename T>
	void const *
	operator()(std::vector<T> const& v) const {
		return v.data();
	}
};

struct array_size_visitor : public boost::static_visitor<std::size_t> {
	template<typename T>
	std::size_t
	operator()(std::vector<T> const& v) const {
		return v.size();
	}
};

void
write_nc_cntr(
	nc_cntr const& cntr,
	std::string const& path,
	bool overwrite
) {
	int ncid, status, ndims = 0, nvars = 0;
	
	status = nc_create(path.data(), overwrite ? NC_CLOBBER : NC_NOCLOBBER, &ncid);
	throw_if_error(ncid, status, path, false);
	
	for(nc_dimension const& dim : cntr.dimensions()) {
		int dimid;
		status = nc_def_dim(ncid, dim.name().data(), dim.length(), &dimid);
		throw_if_error(ncid, status, path, true);
		
		++ndims;
		BOOST_ASSERT(dimid == ndims-1);
	}
	
	for(nc_variable const& var : cntr.variables()) {
		int varid;
		std::vector<int> dimids;
		for(nc_dimension const& dim : var.dimensions()) {
			auto dimid_it = std::find_if(
				cntr.dimensions().begin(),
				cntr.dimensions().end(),
				[&dim](nc_dimension const& test_dim) {
					return dim == test_dim;
				}
			);
			BOOST_ASSERT(dimid_it != cntr.dimensions().end());
			dimids.push_back(static_cast<int>(
				std::distance(cntr.dimensions().begin(), dimid_it)
			));
		}
		BOOST_ASSERT(dimids.size() == var.dimensions().size());
		
		std::optional<nc_type> xtype = boost::apply_visitor(nc_type_visitor{}, var.data());
		if (!xtype) {
			status = NC_EBADTYPID;
			throw_if_error(ncid, status, path, true);
		}
		
		status = nc_def_var(ncid, var.name().data(), *xtype, dimids.size(), dimids.data(), &varid);
		throw_if_error(ncid, status, path, true);
		++nvars;
		BOOST_ASSERT(varid == nvars-1);
	}
	
	for(nc_attribute const& attr : cntr.attributes()) {
		int gattrid;
		
		std::optional<nc_type> xtype = boost::apply_visitor(nc_type_visitor{}, attr.value());
		if (!xtype) {
			status = NC_EBADTYPID;
			throw_if_error(ncid, status, path, true);
		}
		
		status = nc_put_att(ncid, NC_GLOBAL, attr.name().data(), *xtype,
			boost::apply_visitor(array_size_visitor{}, attr.value()),
			boost::apply_visitor(array_ptr_visitor{}, attr.value())
		);
		throw_if_error(ncid, status, path, true);
	}
	
	status = nc_enddef(ncid);
	throw_if_error(ncid, status, path, true);
	
	for(std::size_t i = 0; i < cntr.variables().size(); ++i) {
		status = nc_put_var(
			ncid,
			static_cast<int>(i),
			boost::apply_visitor(array_ptr_visitor{}, cntr.variables()[i].data())
		);
		throw_if_error(ncid, status, path, true);
	}
	
	status = nc_close(ncid);
	throw_if_error(ncid, status, path, false);
}

HBRS_THETA_UTILS_NAMESPACE_END