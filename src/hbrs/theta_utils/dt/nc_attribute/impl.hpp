/* Copyright (c) 2019 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_DT_NC_ATTRIBUTE_IMPL_HPP
#define HBRS_THETA_UTILS_DT_NC_ATTRIBUTE_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/core/preprocessor.hpp>
#include <boost/hana/core.hpp>
#include <vector>
#include <string>
#include <boost/variant.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;

struct nc_attribute {
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
	
	nc_attribute(std::string name, array value);
	
	nc_attribute(nc_attribute const&) = default;
	nc_attribute(nc_attribute &&) = default;
	
	nc_attribute&
	operator=(nc_attribute const&) = default;
	nc_attribute&
	operator=(nc_attribute &&) = default;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(name, std::string)
	HBRS_THETA_UTILS_DECLARE_ATTR(value, array)
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::nc_attribute > {
	using type = hbrs::theta_utils::nc_attribute_tag;
};

template <>
struct make_impl<hbrs::theta_utils::nc_attribute_tag> {
	static hbrs::theta_utils::nc_attribute
	apply(std::string name, hbrs::theta_utils::nc_attribute::array value) {
		return {name, value};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_NC_ATTRIBUTE_IMPL_HPP
