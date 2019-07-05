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

#ifndef HBRS_THETA_UTILS_DT_NC_EXCEPTION_FWD_HPP
#define HBRS_THETA_UTILS_DT_NC_EXCEPTION_FWD_HPP

#include <hbrs/theta_utils/config.hpp>
#include <boost/hana/fwd/core/make.hpp>
#include <boost/hana/fwd/core/to.hpp>

#include <boost/exception/error_info.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;

struct nc_exception;
struct nc_exception_tag {};
constexpr auto make_nc_exception = hana::make<nc_exception_tag>;
constexpr auto to_nc_exception = hana::to<nc_exception_tag>;

typedef boost::error_info<struct errinfo_nc_status_, int> errinfo_nc_status;

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_NC_EXCEPTION_FWD_HPP
