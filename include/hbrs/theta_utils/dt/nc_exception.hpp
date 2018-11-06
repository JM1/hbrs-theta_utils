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

#pragma once

#ifndef HBRS_THETA_UTILS_DT_NC_EXCEPTION_HPP
#define HBRS_THETA_UTILS_DT_NC_EXCEPTION_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/dt/nc_exception.hpp>
#include <hbrs/theta_utils/dt/exception.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;

struct nc_exception : public virtual mpl::exception {};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::nc_exception > {
	using type = hbrs::theta_utils::nc_exception_tag;
};

template <>
struct make_impl<hbrs::theta_utils::nc_exception_tag> {
	
	static decltype(auto)
	apply() {
		return hbrs::theta_utils::nc_exception{};
	}
	
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_NC_EXCEPTION_HPP
