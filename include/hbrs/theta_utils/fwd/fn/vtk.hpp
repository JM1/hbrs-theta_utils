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

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <ostream>
#include <vector>
#include <string>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace fs = boost::filesystem;

enum class vtk_file_format { legacy_ascii, xml_binary };

struct vtk_path;

vtkSmartPointer<vtkUnstructuredGrid>
make_vtk_unstructured_grid(theta_grid const& grid, theta_field const& field);

void
write_vtk_legacy_ascii(vtkSmartPointer<vtkUnstructuredGrid> grid, char const * file_path);

void
write_vtk_xml_binary(vtkSmartPointer<vtkUnstructuredGrid> grid, char const * file_path);

void
write_vtk_xml_binary_parallel(vtkSmartPointer<vtkUnstructuredGrid> grid, char const * file_path);

void
convert_to_vtk(
	theta_grid_path const& grid_path,
	std::vector<theta_field_path> const& field_paths,
	fs::path const& folder,
	std::string const& prefix,
	std::vector<std::string> const& includes,
	std::vector<std::string> const& excludes,
	bool simple_numbering,
	vtk_file_format format,
	bool overwrite);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_FN_VTK_HPP
