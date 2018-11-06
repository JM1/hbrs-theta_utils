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

#ifndef HBRS_THETA_UTILS_FWD_FN_VTK_HPP
#define HBRS_THETA_UTILS_FWD_FN_VTK_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/dt/theta_grid.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <boost/filesystem.hpp>
#include <ostream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace fs = boost::filesystem;

enum class vtk_file_format { legacy_ascii, xml_binary };

void
write_vtk_legacy_ascii(theta_grid const& grid, theta_field const& field, std::ostream & out);

void
write_vtk_legacy_ascii(theta_grid const& grid, theta_field const& field, std::string const& out);

void
write_vtk_xml_binary(theta_grid const& grid, theta_field const& field, std::string const& out);

void
write_vtk(theta_grid const& grid, theta_field const& field, std::string const& out, vtk_file_format frmt);

void
write_vtk(
	theta_grid const& grid, 
	std::vector<std::tuple<theta_field, fs::path>> const& field_series,
	vtk_file_format frmt
);

/* ParaView Data (PVD) File Format supports only XML-based VTK file format, but NOT legacy *.vtk file format files.
 * Ref.: https://www.paraview.org/Wiki/ParaView/Data_formats#PVD_File_Format
 */
void
write_pvd(std::vector<fs::path> const& series, fs::path filename);

std::string
vtk_file_extension(vtk_file_format const& frmt);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_FN_VTK_HPP
