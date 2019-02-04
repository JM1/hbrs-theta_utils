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

/* unused */
static void
write_vtk_legacy_ascii(theta_grid const& grid, theta_field const& field, std::ostream & out) {
	static const char * type = "float" /* TODO: change to double? */;
	out << "# vtk DataFile Version 3.2\n";
	out << "vtk grid and solution\n";
	out << "ASCII\n";
	out << "DATASET UNSTRUCTURED_GRID\n";
	out << boost::format("POINTS %i %s\n") % grid.no_of_points() % type;
	for(int i = 0; i < grid.no_of_points(); ++i) {
		out << boost::format("%8.8f %8.8f %8.8f\n") % grid.points()[i].x % grid.points()[i].y % grid.points()[i].z;
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
	__has_var(pressure)
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
	__check_var(pressure)
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
	__print_vtk_field(pressure)
	__print_vtk_field(residual)

#undef __print_vtk_field_vec
#undef __print_vtk_field
}

HBRS_THETA_UTILS_NAMESPACE_END

#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
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
#include <vtkCommand.h>
#include <vtkObject.h>
#include <vtkMPIController.h>
#include <vtkMultiProcessController.h>
#include <boost/lexical_cast.hpp>
#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpi = hbrs::mpl::detail::mpi;

namespace detail {

template<vtkCommand::EventIds EventId, typename CallData>
struct Observer : public vtkCommand {
	typedef std::function<void(vtkObject*, CallData)> Callback;
	
	Observer(Callback callback) : callback_{callback} {}
	
	virtual void
	Execute(vtkObject * caller, unsigned long event, void *calldata) {
		if (event == EventId) {
			callback_(caller, static_cast<CallData>(calldata));
		}
	}
	
private:
	Callback callback_;
};

typedef Observer<vtkCommand::ErrorEvent, const char *> ErrorObserver;
typedef Observer<vtkCommand::WarningEvent, const char *> WarningObserver;
/* namespace detail */ }

vtkSmartPointer<vtkUnstructuredGrid>
make_vtk_unstructured_grid(theta_grid const& grid, theta_field const& field) {
	static constexpr auto INVALID_ID = std::numeric_limits<std::size_t>::max();
	
	auto const mpi_size = mpi::size();
	auto const mpi_rank = mpi::rank();
	bool distributed = !field.global_id().empty();
	BOOST_ASSERT(!distributed ? mpi_size == 1 : true);
	BOOST_ASSERT(mpi_size > 1 ? distributed : true);
	
	std::size_t grid_no_of_points = boost::lexical_cast<std::size_t>(grid.no_of_points());
	std::size_t no_of_points = field.global_id().empty() ? grid_no_of_points : field.global_id().size();
	
	// no_of_points is smaller than grid.no_of_points() if grid was distributed among several processes
	BOOST_ASSERT(no_of_points <= grid_no_of_points);
	
#define __has_var(__var)                                                                                               \
	auto const has_ ## __var = field.__var().size() > 0;                                                               \
	if (has_ ## __var && field.__var().size() > no_of_points) {                                                        \
		BOOST_THROW_EXCEPTION(std::runtime_error{                                                                      \
			std::string{"dimensions of variable "} + #__var + " do not match size of grid"                             \
		});                                                                                                            \
	}                                                                                                                  \
	
	__has_var(density)
	__has_var(x_velocity)
	__has_var(y_velocity)
	__has_var(z_velocity)
	__has_var(pressure)
	__has_var(residual)
	
#undef __has_var
	
	vtkSmartPointer<vtkUnstructuredGrid> vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	
	vtkSmartPointer<detail::ErrorObserver> throw_error{new detail::ErrorObserver{
		[](auto caller, auto calldata){
			BOOST_THROW_EXCEPTION(
				vtk_exception{} << errinfo_vtk_error{std::string{calldata}}
			);
		}
	}};
	vtk_grid->AddObserver(vtkCommand::ErrorEvent, throw_error);
	
	vtkSmartPointer<detail::WarningObserver> print_warning {new detail::WarningObserver{
		[](auto caller, auto calldata){
			std::cerr << calldata << std::endl;
		}
	}};
	vtk_grid->AddObserver(vtkCommand::WarningEvent,print_warning);
	
	std::function<std::size_t(std::size_t)> get_id;
	if (distributed) {
		get_id = [&field](std::size_t i) {
			BOOST_ASSERT(i < field.global_id().size());
			return field.global_id()[i];
		};
	} else {
		get_id = [](std::size_t i) {
			return i;
		};
	}
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtk_grid->SetPoints(points);
	
	// add points of local grid
	for(std::size_t i = 0; i < no_of_points; ++i) {
		std::size_t global_id = get_id(i);
		BOOST_ASSERT(global_id < grid_no_of_points);
		points->InsertNextPoint(
			grid.points()[global_id].x,
			grid.points()[global_id].y,
			grid.points()[global_id].z
		);
	}
	
	std::vector<std::size_t> global_to_local_id(grid_no_of_points, INVALID_ID);
	for(std::size_t i = 0; i < no_of_points; ++i) {
		global_to_local_id[get_id(i)] = i;
	}
	
	std::vector<std::size_t> missing_global_ids;
	
	auto insert_vtk_cell = [&](
		auto const& no_of_objects,
		auto const& points_of_objects,
		auto const& object_size,
		auto cell_type
	) -> std::vector<std::size_t> {
		typedef typename decltype(cell_type)::type CellType;
		std::vector<std::size_t> bdry_objects;
		
		if (distributed) {
			for(int i = 0; i < no_of_objects; ++i) {
				std::size_t no_in_grid = 0;
				std::vector<std::size_t> points_not_in_grid;
				for (std::size_t j = 0; j < object_size; ++j) {
					auto global_id = points_of_objects[i][j];
					auto local_id = global_to_local_id[global_id];
					if (local_id == INVALID_ID) {
						points_not_in_grid.push_back(global_id);
					} else {
						++no_in_grid;
					}
				}
				
				bool partly_in_grid = no_in_grid > 0;
				bool partly_not_in_grid = !points_not_in_grid.empty();
				
				if (!partly_in_grid) {
					continue;
				}
				
				if (partly_not_in_grid) {
					static auto const insert_global_ids = [](
						std::vector<size_t> & seq, 
						std::vector<size_t> & more_seq
					) {
						seq.insert(seq.end(), more_seq.begin(), more_seq.end());
						std::sort(seq.begin(), seq.end());
						auto last = std::unique(seq.begin(), seq.end());
						seq.erase(last, seq.end());
					};
					
					insert_global_ids(missing_global_ids, points_not_in_grid);
					bdry_objects.push_back(i);
					continue;
				}
				
				vtkSmartPointer<CellType> cell = vtkSmartPointer<CellType>::New();
				for (std::size_t j = 0; j < object_size; ++j) {
					std::size_t global_id = points_of_objects[i][j];
					std::size_t local_id = global_to_local_id[global_id];
					BOOST_ASSERT(local_id != INVALID_ID);
					cell->GetPointIds()->SetId(j, local_id);
				}
				
				vtk_grid->InsertNextCell(cell->GetCellType(),cell->GetPointIds());
			}
		} else {
			for(int i = 0; i < no_of_objects; ++i) {
				vtkSmartPointer<CellType> cell = vtkSmartPointer<CellType>::New();
				for (std::size_t j = 0; j < object_size; ++j) {
					cell->GetPointIds()->SetId(j, points_of_objects[i][j]);
				}
				vtk_grid->InsertNextCell(cell->GetCellType(),cell->GetPointIds());
			}
		}
		
		return bdry_objects;
	};
	
#define __insert_vtk_cell(__var, __cell_type)                                                                          \
	std::vector<std::size_t> missing_ ## __var ## s = insert_vtk_cell(                                                 \
		grid.no_of_ ## __var ## s(),                                                                                   \
		grid.points_of_ ## __var ## s(),                                                                               \
		std::tuple_size<theta_grid::__var>::value,                                                                     \
		hana::type_c<__cell_type>                                                                                      \
	);
	__insert_vtk_cell(tetraeder, vtkTetra)
	__insert_vtk_cell(prism, vtkWedge)
	__insert_vtk_cell(hexaeder, vtkHexahedron)
	__insert_vtk_cell(pyramid, vtkPyramid)
	__insert_vtk_cell(surfacetriangle, vtkTriangle)
	__insert_vtk_cell(surfacequadrilateral, vtkQuad)
	
#undef __insert_vtk_cell
	
	typedef std::vector<std::size_t> provided_global_ids;
	std::vector<provided_global_ids> pro_gbl_ids_for_rank;
	std::vector<provided_global_ids> pro_gbl_ids_from_rank;
	std::size_t no_of_provided_global_ids = 0;
	
	// exchange boundary points
	if (distributed) {
		typedef std::vector<std::size_t> required_global_ids;
		std::vector<required_global_ids> req_gbl_ids_by_rank(mpi_size);
		req_gbl_ids_by_rank[mpi_rank] = missing_global_ids;
		
		for(std::size_t i = 0; i < missing_global_ids.size(); ++i) {
			BOOST_ASSERT(global_to_local_id[missing_global_ids[i]] == INVALID_ID);
		}
		
		for(std::size_t i = 1; i < missing_global_ids.size(); ++i) {
			BOOST_ASSERT(missing_global_ids[i-1] < missing_global_ids[i]);
		}
		
		// exchange sizes of req_gbl_ids_by_rank across all nodes using bcasts
		{
			std::vector<MPI_Request> reqs;
			reqs.reserve(mpi_size);
			
			std::vector<std::size_t> sizes;
			sizes.resize(mpi_size, 0);
			sizes[mpi_rank] = missing_global_ids.size();
			
			for(int i = 0; i < mpi_size; ++i) {
				reqs.push_back(
					mpi::ibcast(&sizes[i], 1, i, MPI_COMM_WORLD)
				);
			}
			for(int i = 0; i < mpi_size; ++i) {
				auto stat = mpi::wait(reqs[i]);
				if (i != mpi_rank) {
					req_gbl_ids_by_rank[i].resize(sizes[i], 0);
				}
			}
		}
		
		// exchange req_gbl_ids_by_rank across all nodes using bcasts
		{
			std::vector<MPI_Request> reqs;
			reqs.reserve(mpi_size);
			for(int i = 0; i < mpi_size; ++i) {
				reqs.push_back(
					mpi::ibcast(req_gbl_ids_by_rank[i].data(), req_gbl_ids_by_rank[i].size(), i, MPI_COMM_WORLD)
				);
			}
			for(int i = 0; i < mpi_size; ++i) {
				auto stat = mpi::wait(reqs[i]);
			}
		}
		
		// prepare global ids that this node can provide for other nodes
		pro_gbl_ids_for_rank.resize(mpi_size);
		for(std::size_t i = 0; i < mpi_size; ++i) {
			if (i != mpi_rank) {
				auto & missing = req_gbl_ids_by_rank[i];
				
				for(std::size_t x = 1; x < missing.size(); ++x) {
					BOOST_ASSERT(missing[x-1] < missing[x]);
				}
				
				auto & provided = pro_gbl_ids_for_rank[i];
				for(std::size_t j = 0; j < missing.size(); ++j) {
					auto global_id = missing[j];
					auto local_id = global_to_local_id[global_id];
					if (local_id != INVALID_ID) {
						provided.push_back(global_id);
					}
				}
				
				for(std::size_t x = 1; x < provided.size(); ++x) {
					BOOST_ASSERT(provided[x-1] < provided[x]);
				}
			}
		}
		
		pro_gbl_ids_from_rank.resize(mpi_size);
		{
			// send global_id of missing points between nodes using point-to-point communications
			std::vector<MPI_Request> send_reqs;
			send_reqs.reserve(mpi_size-1);
			for(int i = 0; i < mpi_size; ++i) {
				if (i != mpi_rank) {
					
					for(std::size_t x = 1; x < pro_gbl_ids_for_rank[i].size(); ++x) {
						BOOST_ASSERT(pro_gbl_ids_for_rank[i][x-1] < pro_gbl_ids_for_rank[i][x]);
					}
					
					send_reqs.push_back(
						mpi::isend(
							pro_gbl_ids_for_rank[i].data(),
							pro_gbl_ids_for_rank[i].size(),
							i/*dest*/,
							mpi_rank/*tag*/,
							MPI_COMM_WORLD
						)
					);
				}
			}
			
			//MPI_Barrier is not required here because MPI_Probe is blocking
			
			// probe no of provided points
			for(int i = 0; i < mpi_size; ++i) {
				if (i != mpi_rank) {
					MPI_Status stat = mpi::probe(i, i, MPI_COMM_WORLD);
					int count = mpi::get_count(stat, mpi::datatype(hana::type_c<std::size_t>));
					pro_gbl_ids_from_rank[i].resize(boost::numeric_cast<std::size_t>(count), 0);
				}
			}
			
			// receive global_id of missing points between nodes using point-to-point communications
			std::vector<MPI_Request> recv_reqs;
			recv_reqs.reserve(mpi_size-1);
			for(int i = 0; i < mpi_size; ++i) {
				if (i != mpi_rank) {
					recv_reqs.push_back(
						mpi::irecv(
							pro_gbl_ids_from_rank[i].data(),
							pro_gbl_ids_from_rank[i].size(),
							i /*source*/, 
							i /*tag*/, 
							MPI_COMM_WORLD
						)
					);
				}
			}
			
			for(int i = 0; i < send_reqs.size(); ++i) {
				auto stat = mpi::wait(send_reqs[i]);
			}
			for(int i = 0; i < recv_reqs.size(); ++i) {
				auto stat = mpi::wait(recv_reqs[i]);
			}
		}
		
		for(int i = 0; i < mpi_size; ++i) {
			if (i != mpi_rank) {
				auto & provided_ids = pro_gbl_ids_from_rank[i];
				
				for(std::size_t x = 1; x < provided_ids.size(); ++x) {
					BOOST_ASSERT(provided_ids[x-1] < provided_ids[x]);
				}
				
				for(std::size_t gi = 0; gi < provided_ids.size(); ++gi) {
					auto global_id = provided_ids[gi];
					BOOST_ASSERT(global_to_local_id[global_id] == INVALID_ID);
					auto local_id = no_of_points+no_of_provided_global_ids;
					BOOST_ASSERT(local_id != INVALID_ID);
					global_to_local_id[global_id] = local_id;
					++no_of_provided_global_ids;
					
					points->InsertNextPoint(
						grid.points()[global_id].x,
						grid.points()[global_id].y,
						grid.points()[global_id].z
					);
				}
			}
		}
	}
	
	// add boundary geometry
	if (distributed) {
		auto insert_missing_vtk_cell = [&](
			auto const& missing_objects,
			auto const& points_of_objects,
			auto const& object_size,
			auto cell_type
		) {
			typedef typename decltype(cell_type)::type CellType;
			for(auto i : missing_objects) {
				bool partly_not_in_grid = false;
				for (std::size_t j = 0; j < object_size; ++j) {
					auto global_id = points_of_objects[i][j];
					auto local_id = global_to_local_id[global_id];
					if (local_id == INVALID_ID) {
						//no one has this global id so we dont add this object
						partly_not_in_grid = true;
						break;
					}
				}
				if (partly_not_in_grid) {
					continue;
				}
				
				vtkSmartPointer<CellType> cell = vtkSmartPointer<CellType>::New();
				for (std::size_t j = 0; j < object_size; ++j) {
					auto global_id = points_of_objects[i][j];
					auto local_id = global_to_local_id[global_id];
					BOOST_ASSERT(local_id != INVALID_ID);
					cell->GetPointIds()->SetId(j, local_id);
				}
				
				vtk_grid->InsertNextCell(cell->GetCellType(),cell->GetPointIds());
			}
		};
		
#define __insert_missing_vtk_cell(__var, __cell_type)                                                                  \
	insert_missing_vtk_cell(                                                                                           \
		missing_ ## __var ## s,                                                                                        \
		grid.points_of_ ## __var ## s(),                                                                               \
		std::tuple_size<theta_grid::__var>::value,                                                                     \
		hana::type_c<__cell_type>                                                                                      \
	);
	
		__insert_missing_vtk_cell(tetraeder, vtkTetra)
		__insert_missing_vtk_cell(prism, vtkWedge)
		__insert_missing_vtk_cell(hexaeder, vtkHexahedron)
		__insert_missing_vtk_cell(pyramid, vtkPyramid)
		__insert_missing_vtk_cell(surfacetriangle, vtkTriangle)
		__insert_missing_vtk_cell(surfacequadrilateral, vtkQuad)
	
#undef __insert_missing_vtk_cell
	}
	
	auto exchange_point_data = [&](
		std::vector<double> const& local_data,
		std::vector<std::vector<double>> & data_by_rank
	) -> void {
		data_by_rank.resize(mpi_size, std::vector<double>{});
		
		std::vector<std::vector<double>> local_data_for_rank(mpi_size, std::vector<double>{});
		for(int i = 0; i < mpi_size; ++i) {
			if (i != mpi_rank) {
				auto & data_for_remote = local_data_for_rank[i];
				auto & global_ids_for_remote = pro_gbl_ids_for_rank[i];
				data_for_remote.resize(global_ids_for_remote.size(), 0);
				for(std::size_t g = 0; g < global_ids_for_remote.size(); ++g) {
					auto global_id = global_ids_for_remote[g];
					auto local_id = global_to_local_id[global_id];
					BOOST_ASSERT(local_id != INVALID_ID);
					data_for_remote[g] = local_data[local_id];
				}
			}
		}
		
		std::vector<MPI_Request> send_reqs;
		send_reqs.reserve(mpi_size-1);
		
		for(int i = 0; i < mpi_size; ++i) {
			if (i != mpi_rank) {
				//allocate recv buffer
				data_by_rank[i].resize(pro_gbl_ids_from_rank[i].size(), 0);
				
				//send data
				send_reqs.push_back(
					mpi::isend(
						local_data_for_rank[i].data(),
						local_data_for_rank[i].size(),
						i/*dest*/,
						mpi_rank/*tag*/,
						MPI_COMM_WORLD
					)
				);
			}
		}
		
		for(int i = 0; i < mpi_size; ++i) {
			if (i != mpi_rank) {
				MPI_Status stat = mpi::probe(i, i, MPI_COMM_WORLD);
				int count = mpi::get_count(stat, mpi::datatype(hana::type_c<double>));
				data_by_rank[i].resize(boost::numeric_cast<std::size_t>(count), 0);
			}
		}
		
		std::vector<MPI_Request> recv_reqs;
		recv_reqs.reserve(mpi_size-1);
		for(int i = 0; i < mpi_size; ++i) {
			if (i != mpi_rank) {
				recv_reqs.push_back(
					mpi::irecv(
						data_by_rank[i].data(),
						data_by_rank[i].size(),
						i /*source*/,
						i /*tag*/,
						MPI_COMM_WORLD
					)
				);
			}
		}
		
		for(int i = 0; i < send_reqs.size(); ++i) {
			auto stat = mpi::wait(send_reqs[i]);
		}
		for(int i = 0; i < recv_reqs.size(); ++i) {
			auto stat = mpi::wait(recv_reqs[i]);
		}
	};
	
#define __exchange_var(__var)                                                                                          \
	std::vector<std::vector<double>> bdry_ ## __var ## _from_rank;                                                     \
	if (has_ ## __var && distributed) {                                                                                \
		exchange_point_data(field.__var(), bdry_ ## __var ## _from_rank);                                              \
	}
		
	__exchange_var(density)
	__exchange_var(x_velocity)
	__exchange_var(y_velocity)
	__exchange_var(z_velocity)
	__exchange_var(pressure)
	__exchange_var(residual)
	
#undef __exchange_var
	
	auto insert_vtk_pointdata = [&](
		const char * name,
		std::vector<double> const& f,
		std::vector<std::vector<double>> & more_f
	) {
		vtkSmartPointer<vtkDoubleArray> pd = vtkSmartPointer<vtkDoubleArray>::New();
		
		pd->SetNumberOfValues(no_of_points+no_of_provided_global_ids);
		pd->SetName(name);
		
		for (std::size_t i = 0; i < no_of_points; ++i) {
			pd->SetValue(i, f[i]);
		}
		
		for(int i = 0; i < mpi_size; ++i) {
			if (i != mpi_rank) {
				auto & more_f_from_rank = more_f[i];
				auto & gbl_ids_from_rank = pro_gbl_ids_from_rank[i];
				BOOST_ASSERT(more_f_from_rank.size() == gbl_ids_from_rank.size());
				
				for(std::size_t d = 0; d < more_f_from_rank.size(); ++d) {
					auto data_f = more_f_from_rank[d];
					auto global_id = gbl_ids_from_rank[d];
					auto local_id = global_to_local_id[global_id];
					BOOST_ASSERT(local_id != INVALID_ID);
					BOOST_ASSERT(local_id < (no_of_points + no_of_provided_global_ids));
					pd->SetValue(local_id, data_f);
				}
				
			}
		}
		
		vtk_grid->GetPointData()->AddArray(pd);
	};
	
#define __insert_vtk_pointdata(__field)                                                                                \
	{                                                                                                                  \
		if (has_ ## __field) {                                                                                         \
			insert_vtk_pointdata(                                                                                  \
				#__field,                                                                                               \
				field.__field(),                                                                                      \
				bdry_ ## __field ## _from_rank                                                                       \
			);                                                                                                         \
		}                                                                                                              \
	}

	auto insert_vtk_pointdata_vec = [&](
		const char * name,
		std::vector<double> const& f1,
		std::vector<double> const& f2,
		std::vector<double> const& f3,
		std::vector<std::vector<double>> & more_f1,
		std::vector<std::vector<double>> & more_f2,
		std::vector<std::vector<double>> & more_f3
	) {
		vtkSmartPointer<vtkDoubleArray> f = vtkSmartPointer<vtkDoubleArray>::New();
			f->SetNumberOfComponents(3);
			f->SetNumberOfTuples(no_of_points+no_of_provided_global_ids);
			f->SetName(name);
			
			for (std::size_t i = 0; i < no_of_points; ++i) {
				f->SetTuple3(i, f1[i], f2[i], f3[i]);
			}
			
			for(int i = 0; i < mpi_size; ++i) {
				if (i != mpi_rank) {
					auto & more_f1_from_rank = more_f1[i];
					auto & more_f2_from_rank = more_f2[i];
					auto & more_f3_from_rank = more_f3[i];
					auto & gbl_ids_from_rank = pro_gbl_ids_from_rank[i];
					BOOST_ASSERT(more_f1_from_rank.size() == gbl_ids_from_rank.size());
					BOOST_ASSERT(more_f2_from_rank.size() == gbl_ids_from_rank.size());
					BOOST_ASSERT(more_f3_from_rank.size() == gbl_ids_from_rank.size());
					
					for(std::size_t d = 0; d < more_f1_from_rank.size(); ++d) {
						auto data_f1 = more_f1_from_rank[d];
						auto data_f2 = more_f2_from_rank[d];
						auto data_f3 = more_f3_from_rank[d];
						auto global_id = gbl_ids_from_rank[d];
						auto local_id = global_to_local_id[global_id];
						BOOST_ASSERT(local_id != INVALID_ID);
						BOOST_ASSERT(local_id < (no_of_points + no_of_provided_global_ids));
						f->SetTuple3(local_id, data_f1, data_f2, data_f3);
					}
					
				}
			}
			
			vtk_grid->GetPointData()->AddArray(f);
	};
	
#define __insert_vtk_pointdata_vec(__name, __field1, __field2, __field3)                                               \
	{                                                                                                                  \
		auto has_ ## __name = has_ ## __field1 && has_ ## __field2 && has_ ## __field3;                                \
		if (has_ ## __name) {                                                                                          \
			insert_vtk_pointdata_vec(                                                                                  \
				#__name,                                                                                               \
				field.__field1(),                                                                                      \
				field.__field2(),                                                                                      \
				field.__field2(),                                                                                      \
				bdry_ ## __field1 ## _from_rank,                                                                       \
				bdry_ ## __field2 ## _from_rank,                                                                       \
				bdry_ ## __field3 ## _from_rank                                                                        \
			);                                                                                                         \
		}                                                                                                              \
	}
	
	// TODO: Make it possible to select fields that gonna be written
	__insert_vtk_pointdata(density)
	__insert_vtk_pointdata_vec(velocity, x_velocity, y_velocity, z_velocity)
	__insert_vtk_pointdata(pressure)
	__insert_vtk_pointdata(residual)
	
#undef __insert_vtk_pointdata_vec
#undef __insert_vtk_pointdata
	
	return vtk_grid;
}

void
write_vtk_legacy_ascii(vtkSmartPointer<vtkUnstructuredGrid> grid, char const * file_path) {
	vtkNew<vtkUnstructuredGridWriter> wtr;
	vtkSmartPointer<detail::ErrorObserver> throw_error{new detail::ErrorObserver{
		[](auto caller, auto calldata){
			BOOST_THROW_EXCEPTION(
				vtk_exception{} << errinfo_vtk_error{std::string{calldata}}
			);
		}
	}};
	wtr->AddObserver(vtkCommand::ErrorEvent, throw_error);
	wtr->SetFileName(file_path);
	wtr->SetInputData(grid);
	wtr->SetHeader("vtk grid and solution");
	wtr->Write();
}

void
write_vtk_xml_binary(vtkSmartPointer<vtkUnstructuredGrid> grid, char const * file_path) {
	vtkNew<vtkXMLUnstructuredGridWriter> wtr;
	vtkSmartPointer<detail::ErrorObserver> throw_error{new detail::ErrorObserver{
		[](auto caller, auto calldata){
			BOOST_THROW_EXCEPTION(
				vtk_exception{} << errinfo_vtk_error{std::string{calldata}}
			);
		}
	}};
	wtr->AddObserver(vtkCommand::ErrorEvent, throw_error);
	wtr->SetFileName(file_path);
// 	wtr->SetDataModeToAscii();
	wtr->SetDataModeToBinary();
	wtr->SetCompressorTypeToZLib();
	wtr->SetInputData(grid);
	wtr->Write();
}

void
write_vtk_xml_binary_parallel(vtkSmartPointer<vtkUnstructuredGrid> grid, char const * file_path) {
	vtkSmartPointer<vtkMPIController> ctrl = vtkSmartPointer<vtkMPIController>::New();
	ctrl->Initialize(0,0,true);
	vtkMultiProcessController::SetGlobalController(ctrl);

	vtkNew<vtkXMLPUnstructuredGridWriter> wtr;
	
	wtr->SetController(ctrl);
	wtr->SetNumberOfPieces(mpi::size());
	wtr->SetStartPiece(mpi::rank());
	wtr->SetEndPiece(mpi::rank());
	
	vtkSmartPointer<detail::ErrorObserver> throw_error{new detail::ErrorObserver{
		[](auto caller, auto calldata){
			BOOST_THROW_EXCEPTION(
				vtk_exception{} << errinfo_vtk_error{std::string{calldata}}
			);
		}
	}};
	wtr->AddObserver(vtkCommand::ErrorEvent, throw_error);
	wtr->SetFileName(file_path);
// 	wtr->SetDataModeToAscii();
	wtr->SetDataModeToBinary();
	wtr->SetCompressorTypeToZLib();
	wtr->SetInputData(grid);
	wtr->Write();
	
// 	wtr->Update();
	
	ctrl->Finalize(true);
}

HBRS_THETA_UTILS_NAMESPACE_END

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/detail/endian.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace bpt = boost::property_tree;

struct pvd_path {
	fs::path folder;
	std::string prefix;
	
	fs::path
	full_path() const {
		return folder / (prefix + ".pvd");
	}
};

/* ParaView Data (PVD) File Format supports only XML-based VTK file format, but NOT legacy *.vtk file format files.
 * Ref.: https://www.paraview.org/Wiki/ParaView/Data_formats#PVD_File_Format
 */
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

#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/algorithm/string/replace.hpp>
HBRS_THETA_UTILS_NAMESPACE_BEGIN

vtk_path::vtk_path(fs::path folder, std::string basename, bool distributed, vtk_file_format format)
: folder_{folder}, basename_{basename}, distributed_{distributed}, format_{format} {};

std::string
vtk_path::file_extension() const {
	if (format_ == vtk_file_format::legacy_ascii && !distributed_) {
		return "vtk";
	} else if (format_ == vtk_file_format::xml_binary) {
		if (distributed_) {
			return "pvtu";
		} else {
			return "vtu";
		}
	} else {
		BOOST_THROW_EXCEPTION(unsupported_format_exception{} << errinfo_vtk_file_format{format_});
	}
}

fs::path
vtk_path::filename() const {
	std::string stripped = basename_;
	boost::replace_all(stripped, ".", "_"); // only dot allowed is for file extension
	return { stripped + '.' + file_extension() };
}

fs::path
vtk_path::full_path() const {
	return folder_ / filename();
}

vtk_path::operator fs::path() const {
	return full_path();
}

HBRS_THETA_UTILS_DEFINE_ATTR(folder, fs::path, vtk_path)
HBRS_THETA_UTILS_DEFINE_ATTR(basename, std::string, vtk_path)
HBRS_THETA_UTILS_DEFINE_ATTR(distributed, bool, vtk_path)
HBRS_THETA_UTILS_DEFINE_ATTR(format, vtk_file_format, vtk_path)

HBRS_THETA_UTILS_NAMESPACE_END

#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/mpl/fn/transform.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpi = hbrs::mpl::detail::mpi;

static
void
safe_write(fs::path path, bool overwrite) {
	if (fs::exists(path) && !overwrite) {
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("output file %s already exists") % path.string()).str(),
				make_error_code(boost::system::errc::file_exists)
			}
		));
	}
}


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
	bool overwrite
) {
	bool distributed = mpi::size() > 1;
	theta_grid const grid = read_theta_grid(grid_path);

	// create vtk filenames
	std::vector<vtk_path> vtk_paths;
	vtk_paths.reserve(field_paths.size());
	for(std::size_t i = 0; i < field_paths.size(); ++i) {
		theta_field_path field_path = field_paths[i];
		std::string basename;
		if (simple_numbering) {
			basename = prefix + ".nr_" + boost::lexical_cast<std::string>(i);
		} else {
			theta_field_path copy = field_path;
			copy.prefix() = prefix;
			copy.domain_num() = boost::none; // domain num must not be set because it is set by vtk
			basename = copy.filename().string();
		}
		
		vtk_paths.emplace_back(folder, basename, distributed, format);
	};
	
	for(auto vtk_path : vtk_paths) {
		safe_write(vtk_path.full_path(), overwrite);
	}
	
	// create pvd filename
	fs::path pvd_path_ = pvd_path{folder, prefix}.full_path();
	safe_write(pvd_path_, overwrite);
	
	// write vtk files
	for(std::size_t i = 0; i < field_paths.size(); ++i) {
		theta_field_path field_path = field_paths[i];
		
		theta_field field = read_theta_field(
			field_path.full_path().string(),
			includes /* TODO: Or hardcode includes? {".*_velocity", "global_id"} */,
			excludes
		);
		
		vtk_path vtk_path = vtk_paths[i];
		auto vtk_grid = make_vtk_unstructured_grid(grid, field);
		
		if (format == vtk_file_format::legacy_ascii && !distributed) {
			write_vtk_legacy_ascii(vtk_grid, vtk_path.full_path().string().data());
		} else if (format == vtk_file_format::xml_binary) {
			if (distributed) {
				write_vtk_xml_binary_parallel(vtk_grid, vtk_path.full_path().string().data());
			} else {
				write_vtk_xml_binary(vtk_grid, vtk_path.full_path().string().data());
			}
		} else {
			BOOST_THROW_EXCEPTION(unsupported_format_exception{} << errinfo_vtk_file_format{format});
		}
	}
	
	// write pvd file for easier ParaView usage if output was in xml
	if (format == vtk_file_format::xml_binary) {
		std::vector<fs::path> vtk_plain_paths = (*mpl::transform)(
			vtk_paths, 
			[](vtk_path path) { return path.full_path(); }
		);
		write_pvd(vtk_plain_paths, pvd_path_.string());
	}
}
HBRS_THETA_UTILS_NAMESPACE_END
