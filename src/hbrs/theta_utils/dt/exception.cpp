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
#include <mpi.h>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

domain_num_mismatch_error_info::domain_num_mismatch_error_info(fs::path path, int expected, boost::optional<int> got)
: path_{path}, expected_{expected}, got_{got} {}

HBRS_THETA_UTILS_DEFINE_ATTR(path, fs::path, domain_num_mismatch_error_info)
HBRS_THETA_UTILS_DEFINE_ATTR(expected, int, domain_num_mismatch_error_info)
HBRS_THETA_UTILS_DEFINE_ATTR(got, boost::optional<int>, domain_num_mismatch_error_info)

mpi_error_info::mpi_error_info(int code) : code_{code}, class__{0}, string_{} {
	char c_str[MPI_MAX_ERROR_STRING];
    int c_str_len;
    
    int ec = MPI_Error_string(code, c_str, &c_str_len);
	if(ec == MPI_SUCCESS) {
		string_ = std::string(c_str, c_str_len);
	}
	
	int class_;
	ec = MPI_Error_class(code, &class_);
	if (ec == MPI_SUCCESS) {
		class__ = class_;
	}
}
mpi_error_info::mpi_error_info(int code, int class_, std::string string) 
: code_{code}, class__{class_}, string_{string} {}

HBRS_THETA_UTILS_DEFINE_ATTR(code, int, mpi_error_info)
HBRS_THETA_UTILS_DEFINE_ATTR(class_, int, mpi_error_info)
HBRS_THETA_UTILS_DEFINE_ATTR(string, std::string, mpi_error_info)

HBRS_THETA_UTILS_NAMESPACE_END