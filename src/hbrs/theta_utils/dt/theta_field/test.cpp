/* Copyright (c) 2018-2019 Jakob Meng, <jakobmeng@web.de>
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

#define BOOST_TEST_MODULE dt_theta_field_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <hbrs/mpl/config.hpp>
#include <hbrs/mpl/detail/test.hpp>
#include <hbrs/theta_utils/detail/test.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>

#include <array>

namespace utf = boost::unit_test;
using namespace hbrs::theta_utils;

namespace {

/* test data */
inline static auto const
fields0 = hbrs::theta_utils::make_theta_field_matrix(
	mpl::rtsam<double, mpl::storage_order::row_major>{
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

BOOST_AUTO_TEST_SUITE(dt_theta_field_test)

using hbrs::mpl::detail::environment_fixture;
BOOST_TEST_GLOBAL_FIXTURE(environment_fixture);

BOOST_AUTO_TEST_CASE(read, * utf::precondition(detail::mpi_world_size_condition{{0,4}}) ) {
	for(auto scheme: { theta_field_path::naming_scheme::theta, theta_field_path::naming_scheme::tau_unsteady }) {
		detail::io_fixture fx{"read"};
		
		BOOST_TEST_MESSAGE("Working directory: " << fx.wd().path().string());
		
		//prepare test data
		{
			auto paths = detail::make_theta_field_paths(
				fx.wd().path(), fx.prefix(), fields0, scheme
			);
			
			for(std::size_t j = 0; j < paths.size(); ++j) {
				auto path = paths.at(j);
				auto field = fields0_bytes.at(j);
				auto fsz = std::distance(std::get<0>(field), std::get<1>(field));
				BOOST_TEST_MESSAGE("Writing " << fsz << " bytes of test data to file " << path.full_path().string());
				
				detail::write_binary(path.full_path().string(), reinterpret_cast<char const*>(std::get<0>(field)), boost::numeric_cast<std::size_t>(fsz));
			}
		}
		
		// read test data
		auto paths = find_theta_fields(fx.wd().path(), fx.prefix());
		auto got = theta_field_matrix{read_theta_fields(paths)};
		
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
			//compare to reference data
			auto ref_sz = fields0.size();
			auto ref = detail::copy_matrix(fields0, mpl::el_matrix<double>{(El::Int)ref_sz.m(),(El::Int)ref_sz.n()});
			auto got_sz = got.size();
			auto got_ = detail::copy_matrix(got, mpl::el_matrix<double>{(El::Int)got_sz.m(),(El::Int)got_sz.n()});
			HBRS_MPL_TEST_MMEQ(ref, got_, false);
		#else
			//TODO: Implement test without dependency to Elemental
		#endif
	}
}


BOOST_AUTO_TEST_CASE(write, * utf::precondition(detail::mpi_world_size_condition{{0,4}}) ) {
	for(auto scheme: { theta_field_path::naming_scheme::theta, theta_field_path::naming_scheme::tau_unsteady }) {
		detail::io_fixture fx{"write"};
		BOOST_TEST_MESSAGE("Working directory: " << fx.wd().path().string());
	
		auto paths = detail::make_theta_field_paths(
			fx.wd().path(), fx.prefix(), fields0, scheme
		);
		
		// write test data
		std::vector< std::tuple<theta_field, theta_field_path> > fields;
		fields.reserve(paths.size());
		
		for(std::size_t j = 0; j < paths.size(); ++j) {
			fields.push_back({fields0.data().at(j), paths.at(j)});
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
}

BOOST_AUTO_TEST_SUITE_END()
