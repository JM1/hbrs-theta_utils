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

#include <boost/test/tree/test_unit.hpp>
#include <hbrs/mpl/detail/mpi.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/lexical_cast.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace mpi = hbrs::mpl::detail::mpi;
namespace detail {

mpi_world_size_condition::mpi_world_size_condition(
	boost::integer_range<std::size_t> supported
) : supported_{supported} {}

boost::test_tools::assertion_result
mpi_world_size_condition::operator()(boost::unit_test::test_unit_id) {
	auto sz = mpi::comm_size();
	if (boost::range::find(supported_, sz) != supported_.end()) {
		return true;
	}
	
	boost::test_tools::assertion_result ans(false);
	ans.message() << "world size " << sz << " is not supported";
	return ans;
}

fs::path
temp_test_path(fs::path base /*= fs::temp_directory_path()*/) {
	return base / "theta_field_test" / fs::unique_path();
}

temp_test_directory::temp_test_directory(
	fs::path path /*= temp_test_path()*/
) : path_{path} {
	fs::create_directories(path_);
}

temp_test_directory::~temp_test_directory() {
	fs::remove_all(path_);
}

fs::path const&
temp_test_directory::path() const { return (path_); }

io_fixture::io_fixture(std::string tag) : wd_{}, prefix_{} {
	prefix_ = tag + "test_wsz" + boost::lexical_cast<std::string>(mpi::comm_size());
}

temp_test_directory const&
io_fixture::wd() const { return (wd_); }

std::string const&
io_fixture::prefix() const { return (prefix_); }
	

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END
