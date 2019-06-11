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

#ifndef HBRS_THETA_UTILS_DT_EXCEPTION_HPP
#define HBRS_THETA_UTILS_DT_EXCEPTION_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/dt/exception.hpp>
#include <hbrs/mpl/dt/exception.hpp>
#include <hbrs/theta_utils/preprocessor/core.hpp>
#include <boost/hana/core.hpp>
#include <boost/optional.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;

struct invalid_backend_exception : virtual mpl::exception {};
struct ambiguous_domain_num_exception : virtual mpl::exception {};
struct ambiguous_naming_scheme_exception : virtual mpl::exception {};
struct domain_num_mismatch_exception : virtual mpl::exception {};
struct mpi_not_initialized_exception : virtual mpl::exception {};
struct unsupported_format_exception : virtual mpl::exception {};
struct ambiguous_grid_exception : virtual mpl::exception {};
struct invalid_number_range_spec_exception : virtual mpl::exception {};
struct invalid_grid_exception : virtual mpl::exception {};
struct vtk_exception : virtual mpl::exception {};

struct domain_num_mismatch_error_info {
	domain_num_mismatch_error_info(fs::path path, int expected, boost::optional<int> got);
	
	HBRS_THETA_UTILS_DECLARE_ATTR(path, fs::path)
	HBRS_THETA_UTILS_DECLARE_ATTR(expected, int)
	HBRS_THETA_UTILS_DECLARE_ATTR(got, boost::optional<int>)
};

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_EXCEPTION_HPP
