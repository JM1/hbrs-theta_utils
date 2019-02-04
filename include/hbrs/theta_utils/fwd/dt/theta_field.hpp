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

#ifndef HBRS_THETA_UTILS_FWD_DT_THETA_FIELD_HPP
#define HBRS_THETA_UTILS_FWD_DT_THETA_FIELD_HPP

#include <hbrs/theta_utils/config.hpp>
#include <boost/hana/fwd/core/make.hpp>
#include <boost/hana/fwd/core/to.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <tuple>
#include <string>
#include <vector>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace fs = boost::filesystem;

struct theta_field_path;

struct theta_field;
struct theta_field_tag {};
constexpr auto make_theta_field = hana::make<theta_field_tag>;
constexpr auto to_theta_field = hana::to<theta_field_tag>;

std::vector<theta_field_path>
find_theta_fields(
	fs::path const& dir,
	std::string const& prefix
);

std::vector<theta_field_path>
filter_theta_fields_by_domain_num(
	std::vector<theta_field_path> const& fields,
	boost::optional<int> const& domain_num
);

theta_field
read_theta_field(
	std::string const& file_path,
	std::vector<std::string> const& includes = {} /*regex filter*/,
	std::vector<std::string> const& excludes = {} /*regex filter*/
);

std::vector<theta_field>
read_theta_fields(
	std::vector<theta_field_path> const& paths,
	std::vector<std::string> const& includes = {} /*regex filter*/,
	std::vector<std::string> const& excludes = {} /*regex filter*/
);

void
write_theta_field(
	theta_field field,
	std::string const& file_path,
	bool overwrite = false
);

void
write_theta_fields(
	std::vector< std::tuple<theta_field, theta_field_path> > fields,
	bool overwrite = false
);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_DT_THETA_FIELD_HPP
