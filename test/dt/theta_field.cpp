/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
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

#define BOOST_TEST_MODULE theta_field_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/mpl/detail/test.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/storage_order.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/detail/make_theta_fields.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <array>
#include <cstdio>
#include <fstream>
#include <iterator>

namespace mpl = hbrs::mpl;
namespace mpi = hbrs::mpl::detail::mpi;
namespace utf = boost::unit_test;
namespace tt = boost::test_tools;
namespace fs = boost::filesystem;
using namespace hbrs::theta_utils;

namespace {

struct mpi_world_size {
	mpi_world_size(boost::integer_range<std::size_t> rng) : rng{rng} {};
	
	tt::assertion_result
	operator()(utf::test_unit_id) {
		auto sz = mpi::size();
		if (boost::range::find(rng, sz) != rng.end()) {
			return true;
		}
		
		tt::assertion_result ans(false);
		ans.message() << "world size " << sz << " is not supported";
		return ans;
	}
	
	boost::integer_range<std::size_t> rng;
};

fs::path
tmp_path(fs::path base = fs::temp_directory_path()) {
	return base / "theta_field_test" / fs::unique_path();
}

struct tmp_directory {
	tmp_directory(fs::path path = tmp_path()) : path_{path} {
		fs::create_directories(path_);
	}
	
	virtual ~tmp_directory() {
		 fs::remove_all(path_);
	}
	
	decltype(auto)
	path() const { return (path_); }
private:
	fs::path const path_;
};

auto
make_theta_field_paths(
	fs::path const& dir,
	std::string const& prefix,
	std::vector<theta_field> const& field
) {
	auto sz = detail::size(field);
	
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
			{{static_cast<int>(j), 0}, 0} /* timestamp */,
			static_cast<int>(j) * 10 /* step */,
			domain_num
		};
		
		fields.push_back(path);
	}
	
	return fields;
}

template<typename Path>
void
write_binary(Path path, char const * data, std::size_t sz) {
	std::ofstream file;
	file.exceptions(std::ofstream::badbit | std::ofstream::failbit);
	file.open(path, std::ios::binary);
	file.write(data, sz);
}

struct io_fixture {
	io_fixture(std::string tag) : wd_{}, prefix_{} {
		prefix_ = tag + "test_wsz" + boost::lexical_cast<std::string>(mpi::size());
	}
	
	auto const&
	wd() const { return (wd_); }
	
	auto const&
	prefix() const { return (prefix_); }
	
	tmp_directory wd_ /* working directory */;
	std::string prefix_;
};

/* test data */
inline static auto const
fields0 = hbrs::theta_utils::detail::make_theta_fields(
	{
		{
			1,2,3,
			4,5,6,
			7,8,9,
			10,11,12,
			13,14,15,
			16,17,18
		},
		mpl::matrix_size<std::size_t, std::size_t>{6u,3u}
	}
);

inline static constexpr std::size_t
field0_t0_size = 232;

inline static constexpr unsigned char
field0_t0[field0_t0_size] = {
  0x43, 0x44, 0x46, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0c, 0x6e, 0x6f, 0x5f, 0x6f,
  0x66, 0x5f, 0x70, 0x6f, 0x69, 0x6e, 0x74, 0x73, 0x00, 0x00, 0x00, 0x02,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0b,
  0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x0a, 0x78, 0x5f, 0x76, 0x65,
  0x6c, 0x6f, 0x63, 0x69, 0x74, 0x79, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0xb8,
  0x00, 0x00, 0x00, 0x0a, 0x79, 0x5f, 0x76, 0x65, 0x6c, 0x6f, 0x63, 0x69,
  0x74, 0x79, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0xc8, 0x00, 0x00, 0x00, 0x0a,
  0x7a, 0x5f, 0x76, 0x65, 0x6c, 0x6f, 0x63, 0x69, 0x74, 0x79, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x10,
  0x00, 0x00, 0x00, 0xd8, 0x3f, 0xf0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x40, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x1c, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x40, 0x24, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x40, 0x2a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x30, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
};

inline static constexpr std::size_t
field0_t1_size = 232;

inline static constexpr unsigned char
field0_t1[field0_t1_size] = {
  0x43, 0x44, 0x46, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0c, 0x6e, 0x6f, 0x5f, 0x6f,
  0x66, 0x5f, 0x70, 0x6f, 0x69, 0x6e, 0x74, 0x73, 0x00, 0x00, 0x00, 0x02,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0b,
  0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x0a, 0x78, 0x5f, 0x76, 0x65,
  0x6c, 0x6f, 0x63, 0x69, 0x74, 0x79, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0xb8,
  0x00, 0x00, 0x00, 0x0a, 0x79, 0x5f, 0x76, 0x65, 0x6c, 0x6f, 0x63, 0x69,
  0x74, 0x79, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0xc8, 0x00, 0x00, 0x00, 0x0a,
  0x7a, 0x5f, 0x76, 0x65, 0x6c, 0x6f, 0x63, 0x69, 0x74, 0x79, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x10,
  0x00, 0x00, 0x00, 0xd8, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x40, 0x14, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x20, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x40, 0x26, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x40, 0x2c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x31, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
};

inline static constexpr std::size_t
field0_t2_size = 232;

inline static constexpr unsigned char
field0_t2[field0_t2_size] = {
  0x43, 0x44, 0x46, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0c, 0x6e, 0x6f, 0x5f, 0x6f,
  0x66, 0x5f, 0x70, 0x6f, 0x69, 0x6e, 0x74, 0x73, 0x00, 0x00, 0x00, 0x02,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0b,
  0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x0a, 0x78, 0x5f, 0x76, 0x65,
  0x6c, 0x6f, 0x63, 0x69, 0x74, 0x79, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0xb8,
  0x00, 0x00, 0x00, 0x0a, 0x79, 0x5f, 0x76, 0x65, 0x6c, 0x6f, 0x63, 0x69,
  0x74, 0x79, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0xc8, 0x00, 0x00, 0x00, 0x0a,
  0x7a, 0x5f, 0x76, 0x65, 0x6c, 0x6f, 0x63, 0x69, 0x74, 0x79, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x10,
  0x00, 0x00, 0x00, 0xd8, 0x40, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x40, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x22, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x40, 0x28, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x40, 0x2e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x32, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
};

typedef std::tuple<unsigned char const *, unsigned char const *> ptr_range;
std::array<ptr_range, 3> fields0_bytes{{
	{field0_t0+0, field0_t0+field0_t0_size},
	{field0_t1+0, field0_t1+field0_t1_size},
	{field0_t2+0, field0_t2+field0_t2_size}
}};

/* unnamed namespace */ }

BOOST_AUTO_TEST_SUITE(theta_field_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(read, * utf::precondition(mpi_world_size{{0,4}}) ) {
	io_fixture fx{"read"};
	
	BOOST_TEST_MESSAGE("Working directory: " << fx.wd().path().string());
	
	//prepare test data
	{
		auto paths = make_theta_field_paths(
			fx.wd().path(), fx.prefix(), fields0
		);
		
		for(std::size_t j = 0; j < paths.size(); ++j) {
			auto path = paths.at(j);
			auto field = fields0_bytes.at(j);
			auto fsz = std::distance(std::get<0>(field), std::get<1>(field));
			BOOST_TEST_MESSAGE("Writing " << fsz << " bytes of test data to file " << path.full_path().string());
			
			write_binary(path.full_path().string(), reinterpret_cast<char const*>(std::get<0>(field)), fsz);
		}
	}
	
	// read test data
	auto paths = find_theta_fields(fx.wd().path(), fx.prefix());
	auto got = read_theta_fields(paths);
	
	//compare to reference data
	auto ref_sz = detail::size(fields0);
	auto ref = detail::copy_matrix(fields0, El::Matrix<double>{(El::Int)ref_sz.m(),(El::Int)ref_sz.n()});
	auto got_sz = detail::size(got);
	auto got_ = detail::copy_matrix(got, El::Matrix<double>{(El::Int)got_sz.m(),(El::Int)got_sz.n()});
	HBRS_MPL_TEST_MMEQ(ref, got_, false);
}

BOOST_AUTO_TEST_CASE(write, * utf::precondition(mpi_world_size{{0,4}}) ) {
	io_fixture fx{"write"};
	BOOST_TEST_MESSAGE("Working directory: " << fx.wd().path().string());
	
	auto paths = make_theta_field_paths(
		fx.wd().path(), fx.prefix(), fields0
	);
	
	// write test data
	std::vector< std::tuple<theta_field, theta_field_path> > fields;
	fields.reserve(paths.size());
	
	for(std::size_t j = 0; j < paths.size(); ++j) {
		fields.push_back({fields0.at(j), paths.at(j)});
	}
	write_theta_fields(fields, false);
	
	//compare to reference data
	auto read_and_compare = [](ptr_range const& ref, fs::path const& got) {
		std::ifstream got_ifs{got.string(), std::ios::binary};
		got_ifs >> std::noskipws;
		
		unsigned char const * ref_begin = std::get<0>(ref);
		unsigned char const * ref_end = std::get<1>(ref);
		std::istream_iterator<unsigned char> got_begin{got_ifs}, got_end;
		
		BOOST_CHECK_EQUAL_COLLECTIONS(ref_begin, ref_end, got_begin, got_end);
	};
	
	for(std::size_t j = 0; j < paths.size(); ++j) {
		BOOST_TEST_MESSAGE("timestep := " << j);
		read_and_compare(fields0_bytes.at(j), std::get<1>(fields.at(j)).full_path());
	}
}

BOOST_AUTO_TEST_SUITE_END()
