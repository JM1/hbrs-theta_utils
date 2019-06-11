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

#ifndef HBRS_THETA_UTILS_DT_NC_CNTR_HPP
#define HBRS_THETA_UTILS_DT_NC_CNTR_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/dt/nc_cntr.hpp>
#include <hbrs/theta_utils/dt/nc_dimension.hpp>
#include <hbrs/theta_utils/dt/nc_variable.hpp>
#include <hbrs/theta_utils/dt/nc_attribute.hpp>
#include <hbrs/theta_utils/preprocessor/core.hpp>
#include <boost/optional.hpp>
#include <boost/hana/core.hpp>
#include <vector>
#include <string>


HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;

struct nc_cntr {
public:
	nc_cntr(
		std::vector<nc_dimension> dims,
		std::vector<nc_variable> vars,
		std::vector<nc_attribute> attrs
	);
	
	nc_cntr(nc_cntr const&) = default;
	nc_cntr(nc_cntr &&) = default;
	
	nc_cntr&
	operator=(nc_cntr const&) = default;
	nc_cntr&
	operator=(nc_cntr &&) = default;
	
	boost::optional<nc_dimension&>
	dimension(std::string const& name);
	
	boost::optional<nc_dimension const&>
	dimension(std::string const& name) const;
	
	boost::optional<nc_variable &>
	variable(std::string const& name);
	
	boost::optional<nc_variable const&>
	variable(std::string const& name) const;
	
	boost::optional<nc_attribute &>
	attribute(std::string const& name);
	
	boost::optional<nc_attribute const&>
	attribute(std::string const& name) const;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(dimensions, std::vector<nc_dimension>)
	HBRS_THETA_UTILS_DECLARE_ATTR(variables, std::vector<nc_variable>)
	HBRS_THETA_UTILS_DECLARE_ATTR(attributes, std::vector<nc_attribute>)
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::nc_cntr > {
	using type = hbrs::theta_utils::nc_cntr_tag;
};

template <>
struct make_impl<hbrs::theta_utils::nc_cntr_tag> {
	
	static hbrs::theta_utils::nc_cntr
	apply(
		std::vector<hbrs::theta_utils::nc_dimension> dims,
		std::vector<hbrs::theta_utils::nc_variable> vars,
		std::vector<hbrs::theta_utils::nc_attribute> attrs
	) {
		return {dims, vars, attrs};
	}
	
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_NC_CNTR_HPP
