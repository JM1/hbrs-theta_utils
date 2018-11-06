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

#include <hbrs/theta_utils/fn/vtk.hpp>

#include <boost/format.hpp>
#include <vtkCellType.h>
#include <boost/throw_exception.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

void
write_vtk_legacy_ascii(theta_grid const& grid, theta_field const& field, std::ostream & out) {
	static const char * type = "float" /* TODO: change to double? */;
	out << "# vtk DataFile Version 3.2\n";
	out << "vtk grid and solution\n";
	out << "ASCII\n";
	out << "DATASET UNSTRUCTURED_GRID\n";
	out << boost::format("POINTS %i %s\n") % grid.no_of_points() % type;
	for(int i = 0; i < grid.no_of_points(); ++i) {
		out << boost::format("%8.8f %8.8f %8.8f\n") % grid.points_xc()[i] % grid.points_yc()[i] % grid.points_zc()[i];
	}
	
	int cells_n = grid.no_of_tetraeders() + grid.no_of_prisms() + grid.no_of_hexaeders() + grid.no_of_pyramids() +
		grid.no_of_surfacetriangles() + grid.no_of_surfacequadrilaterals();
	
	int cells_size = 
		(std::tuple_size<theta_grid::tetraeder>::value+1) * grid.no_of_tetraeders() +
		(std::tuple_size<theta_grid::prism>::value+1) * grid.no_of_prisms() +
		(std::tuple_size<theta_grid::hexaeder>::value+1) * grid.no_of_hexaeders() +
		(std::tuple_size<theta_grid::pyramid>::value+1) * grid.no_of_pyramids() +
		(std::tuple_size<theta_grid::surfacetriangle>::value+1) * grid.no_of_surfacetriangles() +
		(std::tuple_size<theta_grid::surfacequadrilateral>::value+1) * grid.no_of_surfacequadrilaterals();
	
	out << boost::format("\nCELLS %i %i\n") % cells_n % cells_size;
	
#define __print_vtk_value(__var)                                                                                       \
	for(int i = 0; i < grid.no_of_ ## __var ## s(); ++i) {                                                             \
		out << boost::format("%i") % std::tuple_size<theta_grid::__var>::value;                                        \
		for (std::size_t j = 0; j < std::tuple_size<theta_grid::__var>::value; ++j) {                                  \
			out << boost::format(" %i") % grid.points_of_ ## __var ## s()[i][j];                                       \
		}                                                                                                              \
		out << '\n';                                                                                                   \
	}
	
	__print_vtk_value(tetraeder)
	__print_vtk_value(prism)
	__print_vtk_value(hexaeder)
	__print_vtk_value(pyramid)
	__print_vtk_value(surfacetriangle)
	__print_vtk_value(surfacequadrilateral)
	
#undef __print_vtk_value
	
	out << boost::format("\nCELL_TYPES %i\n") % cells_n;
	
#define __print_vtk_type(__var, __vtk_type)                                                                            \
	for(int i = 0; i < grid.no_of_ ## __var ## s(); ++i) {                                                             \
		out << boost::format("%i\n") % __vtk_type;                                                                     \
	}
	
	__print_vtk_type(tetraeder, VTK_TETRA)
	__print_vtk_type(prism, VTK_WEDGE)
	__print_vtk_type(hexaeder, VTK_HEXAHEDRON)
	__print_vtk_type(pyramid, VTK_PYRAMID)
	__print_vtk_type(surfacetriangle, VTK_TRIANGLE)
	__print_vtk_type(surfacequadrilateral, VTK_QUAD)
	
#undef __print_vtk_type
	
	out << boost::format("\nPOINT_DATA %i\n") % grid.no_of_points();
	
	
#define __has_var(__var)                                                                                               \
	auto has_ ## __var = field.__var().size() > 0;                                                                     \
	if (has_ ## __var && field.__var().size() != (unsigned)grid.no_of_points()) {                                      \
		BOOST_THROW_EXCEPTION(std::runtime_error{                                                                      \
			std::string{"dimensions of variable "} + #__var + " do not match size of grid"                             \
		});                                                                                                            \
	}                                                                                                                  \
	
	__has_var(density)
	__has_var(x_velocity)
	__has_var(y_velocity)
	__has_var(z_velocity)
	__has_var(x_velocity_old)
	__has_var(y_velocity_old)
	__has_var(z_velocity_old)
	__has_var(pressure)
	__has_var(pressure_old)
	__has_var(residual)
	
#undef __has_var

	int field_n = 0;
#define __check_var(__var)                                                                                             \
	if (has_ ## __var) { field_n += 1; };

#define __check_var_vec(__var, __var1, __var2, __var3)                                                                 \
	auto has_ ## __var = has_ ## __var1 && has_ ## __var2 && has_ ## __var3;                                           \
	if (has_ ## __var) { field_n += 1; };
	
	__check_var(density)
	__check_var_vec(velocity, x_velocity, y_velocity, z_velocity)
	__check_var_vec(velocity_old, x_velocity_old, y_velocity_old, z_velocity_old)
	__check_var(pressure)
	__check_var(pressure_old)
	__check_var(residual)
#undef __check_var_vec
#undef __check_var
	
    out << boost::format("FIELD  FieldData %i\n") % field_n;
	
#define __print_vtk_field(__field)                                                                                     \
	if (has_ ## __field) {                                                                                             \
		out << boost::format("%s %i %i %s\n") % #__field % 1 % grid.no_of_points() % type;                             \
		for (int i = 0; i < grid.no_of_points(); ++i) {                                                                \
			out << boost::format("%8.6f ") % field.__field()[i];                                                       \
		}                                                                                                              \
		out << '\n';                                                                                                   \
	}

#define __print_vtk_field_vec(__name, __field1, __field2, __field3)                                                    \
	if (has_ ## __name) {                                                                                              \
		out << boost::format("%s %i %i %s\n") % #__name % 3 % grid.no_of_points() % type;                              \
		for (int i = 0; i < grid.no_of_points(); ++i) {                                                                \
			out << boost::format("%8.6f %8.6f %8.6f\n") %                                                              \
				field.__field1()[i] % field.__field2()[i] % field.__field3()[i];                                       \
		}                                                                                                              \
		out << '\n';                                                                                                   \
	}
	
	/* TODO: Make it possible to select fields that gonna be written */
	__print_vtk_field(density)
	__print_vtk_field_vec(velocity, x_velocity, y_velocity, z_velocity)
	__print_vtk_field_vec(velocity_old, x_velocity_old, y_velocity_old, z_velocity_old)
	__print_vtk_field(pressure)
	__print_vtk_field(pressure_old)
	__print_vtk_field(residual)

#undef __print_vtk_field_vec
#undef __print_vtk_field
}

HBRS_THETA_UTILS_NAMESPACE_END

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkWedge.h>
#include <vtkHexahedron.h>
#include <vtkPyramid.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

static vtkSmartPointer<vtkUnstructuredGrid> 
make_vtk_grid(theta_grid const& grid, theta_field const& field) {
	

	vtkSmartPointer<vtkUnstructuredGrid> vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for(int i = 0; i < grid.no_of_points(); ++i) {
		points->InsertNextPoint(grid.points_xc()[i], grid.points_yc()[i], grid.points_zc()[i]);
	}
	vtk_grid->SetPoints(points);
	
	
#define __insert_vtk_cell(__var, __cell_type)                                                                          \
	for(int i = 0; i < grid.no_of_ ## __var ## s(); ++i) {                                                             \
		vtkSmartPointer<__cell_type> cell = vtkSmartPointer<__cell_type>::New();                                       \
		for (std::size_t j = 0; j < std::tuple_size<theta_grid::__var>::value; ++j) {                                  \
			cell->GetPointIds()->SetId(j, grid.points_of_ ## __var ## s()[i][j]);                                      \
		}                                                                                                              \
		vtk_grid->InsertNextCell(cell->GetCellType(),cell->GetPointIds());                                             \
	}
	
	__insert_vtk_cell(tetraeder, vtkTetra)
	__insert_vtk_cell(prism, vtkWedge)
	__insert_vtk_cell(hexaeder, vtkHexahedron)
	__insert_vtk_cell(pyramid, vtkPyramid)
	__insert_vtk_cell(surfacetriangle, vtkTriangle)
	__insert_vtk_cell(surfacequadrilateral, vtkQuad)
	
#undef __insert_vtk_cell

#define __has_var(__var)                                                                                               \
	auto has_ ## __var = field.__var().size() > 0;                                                                     \
	if (has_ ## __var && field.__var().size() != (unsigned)grid.no_of_points()) {                                      \
		BOOST_THROW_EXCEPTION(std::runtime_error{                                                                      \
			std::string{"dimensions of variable "} + #__var + " do not match size of grid"                             \
		});                                                                                                            \
	}                                                                                                                  \
	
	__has_var(density)
	__has_var(x_velocity)
	__has_var(y_velocity)
	__has_var(z_velocity)
	__has_var(x_velocity_old)
	__has_var(y_velocity_old)
	__has_var(z_velocity_old)
	__has_var(pressure)
	__has_var(pressure_old)
	__has_var(residual)
	
#undef __has_var

#define __insert_vtk_pointdata(__field)                                                                                \
	{                                                                                                                  \
		if (has_ ## __field) {                                                                                         \
			vtkSmartPointer<vtkDoubleArray> __field = vtkSmartPointer<vtkDoubleArray>::New();                          \
			__field->SetNumberOfValues(grid.no_of_points());                                                           \
			for (int i = 0; i < grid.no_of_points(); ++i) {                                                            \
				__field->SetValue(i, field.__field()[i]);                                                              \
			}                                                                                                          \
			__field->SetName(#__field);                                                                                \
			vtk_grid->GetPointData()->AddArray(__field);                                                               \
		}                                                                                                              \
	}

#define __insert_vtk_pointdata_vec(__name, __field1, __field2, __field3)                                               \
	{                                                                                                                  \
		auto has_ ## __name = has_ ## __field1 && has_ ## __field2 && has_ ## __field3;                                \
		if (has_ ## __name) {                                                                                          \
			vtkSmartPointer<vtkDoubleArray> __name = vtkSmartPointer<vtkDoubleArray>::New();                           \
			__name->SetNumberOfComponents(3);                                                                          \
			__name->SetName(#__name);                                                                                  \
			for (int i = 0; i < grid.no_of_points(); ++i) {                                                            \
				__name->InsertNextTuple3(field.__field1()[i], field.__field2()[i], field.__field3()[i]);               \
			}                                                                                                          \
			vtk_grid->GetPointData()->AddArray(__name);                                                                \
		}                                                                                                              \
	}
	
	/* TODO: Make it possible to select fields that gonna be written */
	__insert_vtk_pointdata(density)
	__insert_vtk_pointdata_vec(velocity, x_velocity, y_velocity, z_velocity)
	__insert_vtk_pointdata_vec(velocity_old, x_velocity_old, y_velocity_old, z_velocity_old)
	__insert_vtk_pointdata(pressure)
	__insert_vtk_pointdata(pressure_old)
	__insert_vtk_pointdata(residual)
	
#undef __insert_vtk_pointdata_vec
#undef __insert_vtk_pointdata
	
	return vtk_grid;
}

void
write_vtk_legacy_ascii(theta_grid const& grid, theta_field const& field, std::string const& out) {
	auto vtk_grid = make_vtk_grid(grid, field);
	
	vtkSmartPointer<vtkUnstructuredGridWriter> wtr = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	wtr->SetFileName( out.data() );
	wtr->SetInputData(vtk_grid);
	wtr->SetHeader("vtk grid and solution");
	wtr->Write();
}

void
write_vtk_xml_binary(theta_grid const& grid, theta_field const& field, std::string const& out) {
	auto vtk_grid = make_vtk_grid(grid, field);
	
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> wtr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	wtr->SetFileName( out.data() );
	//wtr->SetDataModeToAscii();
	wtr->SetDataModeToBinary();
	wtr->SetCompressorTypeToZLib();
	wtr->SetInputData(vtk_grid);
	wtr->Write();
}


HBRS_THETA_UTILS_NAMESPACE_END

#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

void
write_vtk(theta_grid const& grid, theta_field const& field, std::string const& out, vtk_file_format frmt) {
	if (frmt == vtk_file_format::legacy_ascii) {
		write_vtk_legacy_ascii(grid, field, out);
	} else if (frmt == vtk_file_format::xml_binary) {
		write_vtk_xml_binary(grid, field, out);
	} else {
		BOOST_THROW_EXCEPTION(unsupported_format_exception{} << errinfo_vtk_file_format{frmt});
	}
}

void
write_vtk(
	theta_grid const& grid,
	std::vector<std::tuple<theta_field, fs::path>> const& field_series,
	vtk_file_format frmt
) {
	for(auto && fnp : field_series) {
		auto && field = std::get<0>(fnp);
		auto && path = std::get<1>(fnp);
		
		write_vtk(grid, field, path.string(), frmt);
	}
}

std::string
vtk_file_extension(vtk_file_format const& frmt) {
	if (frmt == vtk_file_format::legacy_ascii) {
		return "vtk";
	} else if (frmt == vtk_file_format::xml_binary) {
		return "vtu";
	} else {
		BOOST_THROW_EXCEPTION(unsupported_format_exception{} << errinfo_vtk_file_format{frmt});
	}
}

HBRS_THETA_UTILS_NAMESPACE_END

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/detail/endian.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace bpt = boost::property_tree;

void
write_pvd(std::vector<fs::path> const& series, fs::path filename) {
	/* Example:
	 * <VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">
	 *   <Collection>
	 *     <DataSet part="0" file="asd/asd_0.vtu"/>
	 *   </Collection>
	 * </VTKFile>
	 * 
	 * Ref.: https://www.paraview.org/Wiki/ParaView/Data_formats#PVD_File_Format
	 */
	
	#if defined(BOOST_LITTLE_ENDIAN)
		auto && endianness = "LittleEndian";
	#elif defined(BOOST_BIG_ENDIAN)
		auto && endianness = "BigEndian";
	#else
		#error "unsupported endianness"
	#endif
	
	using bpt::ptree;
	ptree pt;
	ptree & nd_vf = pt.add("VTKFile", "");
	nd_vf.put("<xmlattr>.type", "Collection");
	
	nd_vf.put("<xmlattr>.byte_order", endianness);
	ptree & nd_col = nd_vf.add("Collection", "");
	std::size_t i = 1;
	for(auto && path : series) {
		ptree & nd_ds = nd_col.add("DataSet", "");
		nd_ds.put("<xmlattr>.timestep", i);
		nd_ds.put("<xmlattr>.part", "0");
		
		if (fs::equivalent(path.parent_path(), filename.parent_path())) {
			nd_ds.put("<xmlattr>.file", path.filename().string());
		} else {
			nd_ds.put("<xmlattr>.file", path.string());
		}
		
		++i;
	}
	
	bpt::write_xml(filename.string(), pt);
}

HBRS_THETA_UTILS_NAMESPACE_END