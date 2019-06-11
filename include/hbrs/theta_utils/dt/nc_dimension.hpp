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

#ifndef HBRS_THETA_UTILS_DT_NC_DIMENSION_HPP
#define HBRS_THETA_UTILS_DT_NC_DIMENSION_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/dt/nc_dimension.hpp>
#include <hbrs/theta_utils/preprocessor/core.hpp>
#include <boost/hana/core.hpp>

#include <string>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;

struct nc_dimension {
public:
	nc_dimension(
		std::string name,
		std::size_t length
	);
	
	nc_dimension(nc_dimension const&) = default;
	nc_dimension(nc_dimension &&) = default;
	
	nc_dimension&
	operator=(nc_dimension const&) = default;
	nc_dimension&
	operator=(nc_dimension &&) = default;
	
	bool
	operator==(nc_dimension const& rhs) const;
    bool
    operator!=(nc_dimension const& rhs) const;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(name, std::string)
	HBRS_THETA_UTILS_DECLARE_ATTR(length, std::size_t)
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::nc_dimension > {
	using type = hbrs::theta_utils::nc_dimension_tag;
};

template <>
struct make_impl<hbrs::theta_utils::nc_dimension_tag> {
	
	static decltype(auto)
	apply(std::string name, std::size_t length) {
		return hbrs::theta_utils::nc_dimension{name, length};
	}
	
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_NC_DIMENSION_HPP
