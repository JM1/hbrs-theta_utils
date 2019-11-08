/* Copyright (c) 2019 Jakob Meng, <jakobmeng@web.de>
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

#include <hbrs/mpl/core/preprocessor.hpp>

#include <hbrs/mpl/detail/mpi.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/test/tree/test_unit.hpp>
#include <boost/assert.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace mpi = hbrs::mpl::detail::mpi;
namespace detail {

std::vector<theta_field_path>
make_theta_field_paths(
	fs::path const& dir,
	std::string const& prefix,
	std::vector<theta_field> const& field,
	enum theta_field_path::naming_scheme scheme
) {
	auto sz = local_size(field);
	
	std::vector<theta_field_path> fields;
	fields.reserve(sz.m());
	
	for(std::size_t j = 0; j < sz.n(); ++j) {
		boost::optional<int> domain_num;
		if (mpi::size() > 1) {
			domain_num = mpi::rank();
		}
		
		theta_field_path path{
			dir,
			prefix,
			{
				{boost::lexical_cast<std::string>(j), "000"} /* significand */,
				"00" /* exponent */
			} /* timestamp */,
			static_cast<int>(j) * 10 /* step */,
			domain_num,
			scheme
		};
		
		fields.push_back(path);
	}
	
	return fields;
}

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END
