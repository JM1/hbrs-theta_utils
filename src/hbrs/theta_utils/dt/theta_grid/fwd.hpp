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

#ifndef HBRS_THETA_UTILS_DT_THETA_GRID_FWD_HPP
#define HBRS_THETA_UTILS_DT_THETA_GRID_FWD_HPP

#include <hbrs/theta_utils/config.hpp>
#include <boost/hana/fwd/core/make.hpp>
#include <boost/hana/fwd/core/to.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <string>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace fs = boost::filesystem;

struct theta_grid_path;

struct theta_grid;
struct theta_grid_tag {};
constexpr auto make_theta_grid = hana::make<theta_grid_tag>;
constexpr auto to_theta_grid = hana::to<theta_grid_tag>;

boost::optional<theta_grid_path>
find_theta_grid(
	fs::path const& dir,
	std::string const& prefix
);

theta_grid
read_theta_grid(theta_grid_path const& path);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_THETA_GRID_FWD_HPP
