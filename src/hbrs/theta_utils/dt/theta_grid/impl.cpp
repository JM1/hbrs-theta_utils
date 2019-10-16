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

#include "impl.hpp"

#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <netcdf.h>
#include <stdexcept>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

theta_grid_path::theta_grid_path(
	fs::path folder,
	std::string prefix)
: folder_{folder}, prefix_{prefix} {}

fs::path
theta_grid_path::filename() const {
	// filename example: karman.grid
	std::stringstream fn;
	fn  << prefix_ << ".grid";
	return { fn.str() };
}

fs::path
theta_grid_path::full_path() const {
	return folder_ / filename();
}

HBRS_THETA_UTILS_DEFINE_ATTR(folder, fs::path, theta_grid_path)
HBRS_THETA_UTILS_DEFINE_ATTR(prefix, std::string, theta_grid_path)

theta_grid::theta_grid(
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
) :
	points_per_tetraeder_{points_per_tetraeder},
	points_per_prism_{points_per_prism},
	points_per_hexaeder_{points_per_hexaeder},
	points_per_pyramid_{points_per_pyramid},
	points_per_surfacetriangle_{points_per_surfacetriangle},
	points_per_surfacequadrilateral_{points_per_surfacequadrilateral},

	points_of_tetraeders_{points_of_tetraeders},
	points_of_prisms_{points_of_prisms},
	points_of_hexaeders_{points_of_hexaeders},
	points_of_pyramids_{points_of_pyramids},
	points_of_surfacetriangles_{points_of_surfacetriangles},
	points_of_surfacequadrilaterals_{points_of_surfacequadrilaterals},
	boundarymarker_of_surfaces_{boundarymarker_of_surfaces},
	
	points_{points}
	{}

theta_grid::theta_grid(nc_cntr cntr) {
#define __check(x)                                                                                                     \
	if (!(x)) { BOOST_THROW_EXCEPTION(invalid_grid_exception{}); }

#define __get_opt_dim(__name)                                                                                          \
	{                                                                                                                  \
		 auto opt_dim = cntr.dimension(#__name);                                                                       \
		 if (opt_dim) { __name ## _ = opt_dim->length(); }                                                             \
	}
	
	__get_opt_dim(points_per_tetraeder)
	__get_opt_dim(points_per_prism)
	__get_opt_dim(points_per_hexaeder)
	__get_opt_dim(points_per_pyramid)
	__get_opt_dim(points_per_surfacetriangle)
	__get_opt_dim(points_per_surfacequadrilateral)
	
#define __get_var(__name, __type)                                                                                      \
	__name ## _ = boost::get< __type >(cntr.variable(#__name)->data() );

	__get_var(boundarymarker_of_surfaces, std::vector<int>)
	std::vector<double> points_xc_;
	std::vector<double> points_yc_;
	std::vector<double> points_zc_;
	__get_var(points_xc, std::vector<double>)
	__get_var(points_yc, std::vector<double>)
	__get_var(points_zc, std::vector<double>)
	__check(points_xc_.size() == points_yc_.size());
	__check(points_xc_.size() == points_zc_.size());
	
	points_.reserve(points_xc_.size());
	for(std::size_t i = 0; i < points_xc_.size(); ++i) {
		points_.push_back({points_xc_[i], points_yc_[i], points_zc_[i]});
	}

#define __copy_var(__name, __raw_type, __cast_type, __type)                                                            \
	{                                                                                                                  \
		 auto opt_var = cntr.variable(#__name);                                                                        \
		 if (opt_var) {                                                                                                \
			 auto raw_data = boost::get< __raw_type >(opt_var->data());                                                \
			 auto begin = reinterpret_cast<__cast_type*>(raw_data.data());                                             \
			 auto end = begin + raw_data.size()/std::tuple_size<__cast_type>::value;                                                   \
			 __name ## _ = __type{begin, end};                                                                         \
		 }                                                                                                             \
	}
	
	__copy_var(points_of_tetraeders, std::vector<int>, tetraeder, std::vector<tetraeder>)
	__copy_var(points_of_prisms, std::vector<int>, prism, std::vector<prism>)
	__copy_var(points_of_hexaeders, std::vector<int>, hexaeder, std::vector<hexaeder>)
	__copy_var(points_of_pyramids, std::vector<int>, pyramid, std::vector<pyramid>)
	__copy_var(points_of_surfacetriangles, std::vector<int>, surfacetriangle, std::vector<surfacetriangle>)
	__copy_var(points_of_surfacequadrilaterals, std::vector<int>, surfacequadrilateral, std::vector<surfacequadrilateral>)
	
	if (points_per_tetraeder()) {
		__check(*points_per_tetraeder() == 4);
	}
	if (points_per_prism()) {
		__check(*points_per_prism() == 6);
	}
	if (points_per_hexaeder()) {
		__check(*points_per_hexaeder() == 8);
	}
	if (points_per_pyramid()) {
		__check(*points_per_pyramid() == 5);
	}
	if (points_per_surfacetriangle()) {
		__check(*points_per_surfacetriangle() == 3);
	}
	if (points_per_surfacequadrilateral()) {
		__check(*points_per_surfacequadrilateral() == 4);
	}
	
	int no_of_surfaceelements = boost::numeric_cast<int>( cntr.dimension("no_of_surfaceelements")->length() );
	int no_of_elements = boost::numeric_cast<int>( cntr.dimension("no_of_elements")->length() );
	int no_of_points = boost::numeric_cast<int>( cntr.dimension("no_of_points")->length() );
	
#define __check_array(__array, __dimension)                                                                            \
	int __dimension = boost::numeric_cast<int>(__array().size());                                                      \
	{                                                                                                                  \
		auto opt_dim = cntr.dimension(#__dimension);                                                                   \
		if (opt_dim) {                                                                                                 \
			__check(opt_dim->length() == (unsigned)__dimension);                                                       \
		}                                                                                                              \
	}
	
	__check_array(points_of_tetraeders, no_of_tetraeders)
	__check_array(points_of_prisms, no_of_prisms)
	__check_array(points_of_hexaeders, no_of_hexaeders)
	__check_array(points_of_pyramids, no_of_pyramids)
	__check_array(points_of_surfacetriangles, no_of_surfacetriangles)
	__check_array(points_of_surfacequadrilaterals, no_of_surfacequadrilaterals)
	
	__check(no_of_surfaceelements == (no_of_surfacetriangles + no_of_surfacequadrilaterals));
	__check(no_of_elements == (no_of_tetraeders + no_of_prisms + no_of_hexaeders + no_of_pyramids));
	__check(points_.size() == (unsigned)no_of_points);
}

HBRS_THETA_UTILS_DEFINE_ATTR(points_per_tetraeder, boost::optional<int>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_per_prism, boost::optional<int>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_per_hexaeder, boost::optional<int>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_per_pyramid, boost::optional<int>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_per_surfacetriangle, boost::optional<int>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_per_surfacequadrilateral, boost::optional<int>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_of_tetraeders, std::vector<theta_grid::tetraeder>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_of_prisms, std::vector<theta_grid::prism>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_of_hexaeders, std::vector<theta_grid::hexaeder>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_of_pyramids, std::vector<theta_grid::pyramid>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_of_surfacetriangles, std::vector<theta_grid::surfacetriangle>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_of_surfacequadrilaterals, std::vector<theta_grid::surfacequadrilateral>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(boundarymarker_of_surfaces, std::vector<int>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points, std::vector<theta_grid::point>, theta_grid)

int theta_grid::no_of_points() const { return points_.size(); }
int theta_grid::no_of_tetraeders() const { return points_of_tetraeders_.size(); }
int theta_grid::no_of_prisms() const { return points_of_prisms_.size(); }
int theta_grid::no_of_hexaeders() const { return points_of_hexaeders_.size(); }
int theta_grid::no_of_pyramids() const { return points_of_pyramids_.size(); }
int theta_grid::no_of_surfacetriangles() const { return points_of_surfacetriangles_.size(); }
int theta_grid::no_of_surfacequadrilaterals() const { return points_of_surfacequadrilaterals_.size(); }

int theta_grid::no_of_elements() const { 
	return points_of_tetraeders_.size() + points_of_prisms_.size() 
		+ points_of_hexaeders_.size() + points_of_pyramids_.size();
}

int theta_grid::no_of_surfaceelements() const { 
	return points_of_surfacetriangles_.size() + points_of_surfacequadrilaterals_.size(); 
}



boost::optional<theta_grid_path>
find_theta_grid(
	fs::path const& dir,
	std::string const& prefix
) {
	theta_grid_path path{dir, prefix};
	if (fs::exists(path.full_path())) {
		return path;
	} else {
		return {boost::none};
	}
}

theta_grid
read_theta_grid(theta_grid_path const& path) {
	return { read_nc_cntr(path.full_path().string()) };
}

HBRS_THETA_UTILS_NAMESPACE_END
