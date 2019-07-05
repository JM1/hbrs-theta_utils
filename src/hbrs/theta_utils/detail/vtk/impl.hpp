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

#ifndef HBRS_THETA_UTILS_DETAIL_VTK_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_VTK_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/core/preprocessor.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace fs = boost::filesystem;

struct vtk_path {
	vtk_path(fs::path folder, std::string basename, bool distributed, vtk_file_format format);
	vtk_path(vtk_path const&) = default;
	vtk_path(vtk_path &&) = default;
	
	vtk_path&
	operator=(vtk_path const&) = default;
	vtk_path&
	operator=(vtk_path &&) = default;
	
	std::string
	file_extension() const;
	
	fs::path
	filename() const;
	
	fs::path
	full_path() const;
	
	explicit
	operator fs::path() const;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(folder, fs::path)
	HBRS_THETA_UTILS_DECLARE_ATTR(basename, std::string)
	HBRS_THETA_UTILS_DECLARE_ATTR(distributed, bool)
	HBRS_THETA_UTILS_DECLARE_ATTR(format, vtk_file_format)
};

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_VTK_IMPL_HPP
