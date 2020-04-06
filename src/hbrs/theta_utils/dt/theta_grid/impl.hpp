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

#ifndef HBRS_THETA_UTILS_DT_THETA_GRID_IMPL_HPP
#define HBRS_THETA_UTILS_DT_THETA_GRID_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/theta_utils/dt/nc_cntr.hpp>
#include <hbrs/theta_utils/core/preprocessor.hpp>
#include <boost/hana/core.hpp>
#include <boost/optional.hpp>
#include <vector>
#include <array>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;

struct theta_grid_path {
	theta_grid_path(fs::path folder, std::string prefix);
	theta_grid_path(theta_grid_path const&) = default;
	theta_grid_path(theta_grid_path &&) = default;
	
	theta_grid_path&
	operator=(theta_grid_path const&) = default;
	theta_grid_path&
	operator=(theta_grid_path &&) = default;
	
	fs::path
	filename() const;
	
	fs::path
	full_path() const;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(folder, fs::path)
	HBRS_THETA_UTILS_DECLARE_ATTR(prefix, std::string)
};

struct theta_grid {
public:
    
    //TODO: Support other cell types? e.g. https://vtk.org/doc/nightly/html/vtkCellType_8h.html
	typedef std::array<int,4> tetraeder;
	typedef std::array<int,6> prism;
	typedef std::array<int,8> hexaeder;
	typedef std::array<int,5> pyramid;
	typedef std::array<int,3> surfacetriangle;
	typedef std::array<int,4> surfacequadrilateral;
	struct point {
		double x,y,z;
	};
	
	theta_grid(
		boost::optional<int> points_per_tetraeder,
		boost::optional<int> points_per_prism,
		boost::optional<int> points_per_hexaeder,
		boost::optional<int> points_per_pyramid,
		boost::optional<int> points_per_surfacetriangle,
		boost::optional<int> points_per_surfacequadrilateral,

		std::vector<tetraeder> points_of_tetraeders,
		std::vector<prism> points_of_prisms,
		std::vector<hexaeder> points_of_hexaeders,
		std::vector<pyramid> points_of_pyramids,
		std::vector<surfacetriangle> points_of_surfacetriangles,
		std::vector<surfacequadrilateral> points_of_surfacequadrilaterals,
		std::vector<int> boundarymarker_of_surfaces,
		
		std::vector<point> points
	);
	
	theta_grid(nc_cntr cntr);
	
	theta_grid(theta_grid const&) = default;
	theta_grid(theta_grid &&) = default;
	
	theta_grid&
	operator=(theta_grid const&) = default;
	theta_grid&
	operator=(theta_grid &&) = default;
	
	int no_of_points() const;
	/* no_of_elements = no_of_tetraeders + no_of_prisms + no_of_hexaeders + no_of_pyramids */ 
	int no_of_elements() const;
	/* no_of_surfaceelements = no_of_surfacetriangles + no_of_surfacequadrilaterals */
	int no_of_surfaceelements() const;
	
	int no_of_tetraeders() const;
	int no_of_prisms() const;
	int no_of_hexaeders() const;
	int no_of_pyramids() const;
	int no_of_surfacetriangles() const;
	int no_of_surfacequadrilaterals() const;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(points_per_tetraeder, boost::optional<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_per_prism, boost::optional<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_per_hexaeder, boost::optional<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_per_pyramid, boost::optional<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_per_surfacetriangle, boost::optional<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_per_surfacequadrilateral, boost::optional<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_of_tetraeders, std::vector<tetraeder>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_of_prisms, std::vector<prism>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_of_hexaeders, std::vector<hexaeder>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_of_pyramids, std::vector<pyramid>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_of_surfacetriangles, std::vector<surfacetriangle>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points_of_surfacequadrilaterals, std::vector<surfacequadrilateral>)
	HBRS_THETA_UTILS_DECLARE_ATTR(boundarymarker_of_surfaces, std::vector<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(points, std::vector<point>)
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::theta_grid > {
	using type = hbrs::theta_utils::theta_grid_tag;
};

template <>
struct make_impl<hbrs::theta_utils::theta_grid_tag> {
	template<typename... Args>
	static hbrs::theta_utils::theta_grid
	apply(Args&&... args) {
		return {HBRS_MPL_FWD(args)...};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_THETA_GRID_IMPL_HPP
