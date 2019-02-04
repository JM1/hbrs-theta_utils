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

#pragma once

#ifndef HBRS_THETA_UTILS_DT_NC_VARIABLE_HPP
#define HBRS_THETA_UTILS_DT_NC_VARIABLE_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/dt/nc_variable.hpp>
#include <hbrs/theta_utils/dt/nc_dimension.hpp>
#include <hbrs/theta_utils/preprocessor/core.hpp>
#include <boost/hana/core.hpp>
#include <vector>
#include <string>
#include <boost/variant.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;

struct nc_variable {
public:
	/* TODO: Add more data types here and in nc_ctnr.cpp!
	 * Ref.: 
	 *  https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_data_model.html#classic_model
	 *  https://www.unidata.ucar.edu/software/netcdf/docs/dvarget_8c.html
	 */
	
	typedef boost::variant<
		std::vector<double>, 
		std::vector<float>, 
		std::vector<int>,
		std::vector<long>, 
		std::vector<short>, 
		std::vector<char>, 
		std::vector<unsigned int>, 
		std::vector<unsigned short>
	> array;
	
	nc_variable(std::string name, std::vector<nc_dimension> dims, array data);
	
	nc_variable(nc_variable const&) = default;
	nc_variable(nc_variable &&) = default;
	
	nc_variable&
	operator=(nc_variable const&) = default;
	nc_variable&
	operator=(nc_variable &&) = default;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(name, std::string)
	HBRS_THETA_UTILS_DECLARE_ATTR(dimensions, std::vector<nc_dimension>)
	HBRS_THETA_UTILS_DECLARE_ATTR(data, array)
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::nc_variable > {
	using type = hbrs::theta_utils::nc_variable_tag;
};

template <>
struct make_impl<hbrs::theta_utils::nc_variable_tag> {
	static hbrs::theta_utils::nc_variable
	apply(std::string name, std::vector<hbrs::theta_utils::nc_dimension> dims, hbrs::theta_utils::nc_variable::array data) {
		return {name, dims, data};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_NC_VARIABLE_HPP
