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

#include <hbrs/theta_utils/dt/theta_grid.hpp>
#include <netcdf.h>
#include <stdexcept>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

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
	
	std::vector<double> points_xc,
	std::vector<double> points_yc,
	std::vector<double> points_zc
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
	
	points_xc_{points_xc},
	points_yc_{points_yc},
	points_zc_{points_zc}
	{}

theta_grid::theta_grid(nc_cntr cntr) {
	
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
	__get_var(points_xc, std::vector<double>)
	__get_var(points_yc, std::vector<double>)
	__get_var(points_zc, std::vector<double>)

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
		BOOST_ASSERT(*points_per_tetraeder() == 4);
	}
	if (points_per_prism()) {
		BOOST_ASSERT(*points_per_prism() == 6);
	}
	if (points_per_hexaeder()) {
		BOOST_ASSERT(*points_per_hexaeder() == 8);
	}
	if (points_per_pyramid()) {
		BOOST_ASSERT(*points_per_pyramid() == 5);
	}
	if (points_per_surfacetriangle()) {
		BOOST_ASSERT(*points_per_surfacetriangle() == 3);
	}
	if (points_per_surfacequadrilateral()) {
		BOOST_ASSERT(*points_per_surfacequadrilateral() == 4);
	}
	
	int no_of_surfaceelements = cntr.dimension("no_of_surfaceelements")->length();
	int no_of_elements = cntr.dimension("no_of_elements")->length();
	int no_of_points = cntr.dimension("no_of_points")->length();
	
#define __assert(__array, __dimension)                                                                                 \
	int __dimension = __array().size();                                                                                \
	{                                                                                                                  \
		auto opt_dim = cntr.dimension(#__dimension);                                                                   \
		if (opt_dim) {                                                                                                 \
			BOOST_ASSERT(opt_dim->length() == (unsigned)__dimension);                                                            \
		}                                                                                                              \
	}
	
	__assert(points_of_tetraeders, no_of_tetraeders)
	__assert(points_of_prisms, no_of_prisms)
	__assert(points_of_hexaeders, no_of_hexaeders)
	__assert(points_of_pyramids, no_of_pyramids)
	__assert(points_of_surfacetriangles, no_of_surfacetriangles)
	__assert(points_of_surfacequadrilaterals, no_of_surfacequadrilaterals)
	
	BOOST_ASSERT(no_of_surfaceelements == (no_of_surfacetriangles + no_of_surfacequadrilaterals));
	BOOST_ASSERT(no_of_elements == (no_of_tetraeders + no_of_prisms + no_of_hexaeders + no_of_pyramids));
	BOOST_ASSERT(points_xc_.size() == points_yc_.size());
	BOOST_ASSERT(points_xc_.size() == points_zc_.size());
	BOOST_ASSERT(points_xc_.size() == (unsigned)no_of_points);
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
HBRS_THETA_UTILS_DEFINE_ATTR(points_xc, std::vector<double>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_yc, std::vector<double>, theta_grid)
HBRS_THETA_UTILS_DEFINE_ATTR(points_zc, std::vector<double>, theta_grid)

int theta_grid::no_of_points() const { return points_xc_.size(); }
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

HBRS_THETA_UTILS_NAMESPACE_END

#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/format.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/system/error_code.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

static
boost::optional<fs::path>
parse_theta_grid_path(fs::path path, std::string const& prefix) {
	namespace qi = boost::spirit::qi;
	namespace spirit = boost::spirit;
	namespace phoenix = boost::phoenix;
	
	auto filename = path.filename().string();
	// filename example: karman.para
	
	auto begin = filename.begin();
	auto end = filename.end();
	bool ok = qi::parse(begin, end,
		qi::lexeme[
			qi::lit(prefix) >> qi::lit(".grid")
		]
	);
	if ((begin != end) || !ok) { 
		return {};
	};
	
	return path;
}

static
std::vector<fs::path> 
find_theta_grid_files(
	std::string const& folder_path,
	std::string const& prefix
) {
	fs::path dir{folder_path};
	std::vector<fs::path> all_files;
	std::copy(fs::directory_iterator(dir), fs::directory_iterator(), std::back_inserter(all_files));
	
	std::vector<fs::path> grid_files;
	for(auto && path : all_files) {
		auto field_path = parse_theta_grid_path(path, prefix);
		if (!field_path) {
			continue;
		}
		grid_files.push_back(*field_path);
	}
	
	return grid_files;
}


static
void
sanity_check(std::vector<fs::path> const& paths, std::string const& folder_path, std::string const& prefix) {
	if (paths.empty()) { 
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("no grid file with prefix %s found in folder %s") % prefix % folder_path).str(),
				boost::system::errc::make_error_code(boost::system::errc::no_such_file_or_directory)
			}
		));
	}
	
	if (paths.size() > 1) {
		BOOST_THROW_EXCEPTION(ambiguous_grid_exception{});
	}
}


theta_grid
read_theta_grid(std::string const& folder_path, std::string const& prefix) {
	std::vector<fs::path> grid_files = find_theta_grid_files(folder_path, prefix);
	sanity_check(grid_files, folder_path, prefix);
	
	return read_theta_grid(grid_files[0]);
}

theta_grid
read_theta_grid(fs::path const& file_path) {
	return { read_nc_cntr(file_path.string()) };
}

HBRS_THETA_UTILS_NAMESPACE_END