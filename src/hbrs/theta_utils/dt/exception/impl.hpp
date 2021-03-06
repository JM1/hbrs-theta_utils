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

#ifndef HBRS_THETA_UTILS_DT_EXCEPTION_IMPL_HPP
#define HBRS_THETA_UTILS_DT_EXCEPTION_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/dt/exception.hpp>
#include <hbrs/theta_utils/core/preprocessor.hpp>
#include <boost/hana/core.hpp>
#include <boost/optional.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;

struct HBRS_THETA_UTILS_API invalid_backend_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API ambiguous_domain_num_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API ambiguous_naming_scheme_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API domain_num_mismatch_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API mpi_not_initialized_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API unsupported_format_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API ambiguous_grid_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API invalid_number_range_spec_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API invalid_grid_exception : virtual mpl::exception {};
struct HBRS_THETA_UTILS_API vtk_exception : virtual mpl::exception {};

struct HBRS_THETA_UTILS_API domain_num_mismatch_error_info {
	domain_num_mismatch_error_info(fs::path path, int expected, boost::optional<int> got);
	
	HBRS_THETA_UTILS_DECLARE_ATTR(path, fs::path)
	HBRS_THETA_UTILS_DECLARE_ATTR(expected, int)
	HBRS_THETA_UTILS_DECLARE_ATTR(got, boost::optional<int>)
};

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_EXCEPTION_IMPL_HPP
