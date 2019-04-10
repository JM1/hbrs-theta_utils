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

#ifndef HBRS_THETA_UTILS_FWD_DT_EXCEPTION_HPP
#define HBRS_THETA_UTILS_FWD_DT_EXCEPTION_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/fwd/dt/exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/filesystem.hpp>
#include <hbrs/theta_utils/fwd/fn/vtk.hpp>
#include <hbrs/theta_utils/fwd/dt/command_option.hpp>
#include <tuple>
#include <string>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace fs = boost::filesystem;

struct invalid_backend_exception;
struct ambiguous_domain_num_exception;
struct domain_num_mismatch_exception;
struct mpi_not_initialized_exception;
struct unsupported_format_exception;
struct ambiguous_grid_exception;
struct invalid_number_range_spec_exception;
struct invalid_grid_exception;
struct vtk_exception;

typedef boost::error_info<struct errinfo_ambiguous_domain_num_, std::tuple<fs::path, fs::path> > errinfo_ambiguous_domain_num;
struct domain_num_mismatch_error_info;
typedef boost::error_info<struct errinfo_domain_num_mismatch_, domain_num_mismatch_error_info > errinfo_domain_num_mismatch;
typedef boost::error_info<struct errinfo_vtk_file_format_, vtk_file_format> errinfo_vtk_file_format;
typedef boost::error_info<struct errinfo_number_range_spec_, std::string> errinfo_number_range_spec;
typedef boost::error_info<struct errinfo_vtk_error_, std::string> errinfo_vtk_error;
typedef boost::error_info<struct errinfo_pca_backend_, pca_backend> errinfo_pca_backend;

std::string
to_string(errinfo_ambiguous_domain_num e);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_DT_EXCEPTION_HPP
