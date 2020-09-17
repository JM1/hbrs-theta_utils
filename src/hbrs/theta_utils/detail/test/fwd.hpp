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

#ifndef HBRS_THETA_UTILS_DETAIL_TEST_FWD_HPP
#define HBRS_THETA_UTILS_DETAIL_TEST_FWD_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/dt/theta_field_matrix.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <vector>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace fs = boost::filesystem;
namespace detail {

struct HBRS_THETA_UTILS_API mpi_world_size_condition;
struct HBRS_THETA_UTILS_API temp_test_directory;

HBRS_THETA_UTILS_API
std::vector<theta_field_path>
make_theta_field_paths(
	fs::path const& dir,
	std::string const& prefix,
	theta_field_matrix const& series,
	enum theta_field_path::naming_scheme scheme
);

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_TEST_FWD_HPP
