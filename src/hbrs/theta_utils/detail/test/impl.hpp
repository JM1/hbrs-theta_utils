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

#ifndef HBRS_THETA_UTILS_DETAIL_TEST_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_TEST_IMPL_HPP

#include "fwd.hpp"

#include <boost/test/tools/assertion_result.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include <boost/test/tree/test_unit.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <fstream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace fs = boost::filesystem;
namespace detail {

struct HBRS_THETA_UTILS_API mpi_world_size_condition {
	mpi_world_size_condition(boost::integer_range<std::size_t> supported);
	
	boost::test_tools::assertion_result
	operator()(boost::unit_test::test_unit_id);
	
private:
	boost::integer_range<std::size_t> supported_;
};

fs::path
temp_test_path(fs::path base = fs::temp_directory_path());

struct HBRS_THETA_UTILS_API temp_test_directory {
	temp_test_directory(fs::path path = temp_test_path());
	
	virtual ~temp_test_directory();
	
	fs::path const&
	path() const;
	
private:
	fs::path path_;
};

struct HBRS_THETA_UTILS_API io_fixture {
	io_fixture(std::string tag);
	
	temp_test_directory const&
	wd() const;
	
	std::string const&
	prefix() const;
	
private:
	temp_test_directory wd_ /* working directory */;
	std::string prefix_;
};

template<typename Path>
void
write_binary(Path path, char const * data, std::size_t sz) {
	std::ofstream file;
	file.exceptions(std::ofstream::badbit | std::ofstream::failbit);
	file.open(path, std::ios::binary);
	file.write(data, boost::numeric_cast<long>(sz));
}

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_TEST_IMPL_HPP
