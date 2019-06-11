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

#ifndef HBRS_THETA_UTILS_FWD_DT_NC_CNTR_HPP
#define HBRS_THETA_UTILS_FWD_DT_NC_CNTR_HPP

#include <hbrs/theta_utils/config.hpp>
#include <boost/hana/fwd/core/make.hpp>
#include <boost/hana/fwd/core/to.hpp>

#include <vector>
#include <string>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;

struct nc_cntr;
struct nc_cntr_tag {};
constexpr auto make_nc_cntr = hana::make<nc_cntr_tag>;
constexpr auto to_nc_cntr = hana::to<nc_cntr_tag>;

nc_cntr
read_nc_cntr(
	std::string const& path,
	std::vector<std::string> const& includes = {} /*regex filter*/,
	std::vector<std::string> const& excludes = {} /*regex filter*/
);

void
write_nc_cntr(
	nc_cntr const& cntr,
	std::string const& path,
	bool overwrite = false
);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_DT_NC_CNTR_HPP
