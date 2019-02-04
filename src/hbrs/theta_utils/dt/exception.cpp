/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
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

#include <hbrs/theta_utils/dt/exception.hpp>
#include <sstream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

domain_num_mismatch_error_info::domain_num_mismatch_error_info(fs::path path, int expected, boost::optional<int> got)
: path_{path}, expected_{expected}, got_{got} {}

HBRS_THETA_UTILS_DEFINE_ATTR(path, fs::path, domain_num_mismatch_error_info)
HBRS_THETA_UTILS_DEFINE_ATTR(expected, int, domain_num_mismatch_error_info)
HBRS_THETA_UTILS_DEFINE_ATTR(got, boost::optional<int>, domain_num_mismatch_error_info)

std::string
to_string(errinfo_ambiguous_domain_num e) {
	fs::path const& first = std::get<0>(e.value());
	fs::path const& second = std::get<1>(e.value());
	
	std::stringstream str;
	str 
		<< '[' << boost::error_info_name(e) << "] = "
		<< "[ first: " << first.string() << ", second: " << second.string() << "]"
		<< '\n';
	return str.str();
}

HBRS_THETA_UTILS_NAMESPACE_END