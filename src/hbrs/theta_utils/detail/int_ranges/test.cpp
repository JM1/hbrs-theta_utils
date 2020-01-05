/* Copyright (c) 2018-2020 Jakob Meng, <jakobmeng@web.de>
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

#define BOOST_TEST_MODULE int_ranges_test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <hbrs/theta_utils/detail/int_ranges.hpp>
#include <limits>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;
BOOST_AUTO_TEST_SUITE(int_ranges_test)

BOOST_AUTO_TEST_CASE(read) {
	using namespace hbrs::theta_utils;
	using namespace hbrs::theta_utils::detail;
	
	BOOST_TEST((int_range<std::size_t>{0, std::numeric_limits<std::size_t>::max()}.empty() == false));
	
	std::vector<std::string> seqs {
		{"1-3"},
		{"1-3,7-9"}
	};
	
	std::vector<std::vector<std::size_t>> ref_seqs {
		{1,2,3},
		{1,2,3,7,8,9}
	};
	
	std::vector<std::vector<std::size_t>> got_seqs;
	for(auto seq : seqs) {
		int_ranges<std::size_t> parsed;
		bool pass = parse_int_ranges(seq.begin(), seq.end(), parsed);
		BOOST_TEST(pass);
		
		std::vector<std::size_t> got_seq;
		for(std::size_t i = 0; i < 11; ++i) {
			if (in_int_ranges(parsed, i)) {
				got_seq.push_back(i);
			}
		}
		got_seqs.push_back(got_seq);
	}
	
	BOOST_TEST(ref_seqs == got_seqs);
	
	BOOST_TEST(in_int_ranges(int_ranges<std::size_t>{{0, std::numeric_limits<std::size_t>::max()}}, 0ul));
	BOOST_TEST(in_int_ranges(int_ranges<std::size_t>{{0, std::numeric_limits<std::size_t>::max()}}, 1ul));
	
	BOOST_TEST(in_int_ranges(
		int_ranges<std::size_t>{{0, std::numeric_limits<std::size_t>::max()}},
		std::numeric_limits<std::size_t>::max()-1
	));
	
	BOOST_TEST(in_int_ranges(
		int_ranges<std::size_t>{{0, std::numeric_limits<std::size_t>::max()}},
		std::numeric_limits<std::size_t>::max()
	) == false);
	
	BOOST_TEST(in_int_ranges(
		int_ranges<std::size_t>{{0, std::numeric_limits<std::size_t>::max()-1}},
		std::numeric_limits<std::size_t>::max()
	) == false);
	
	// test empty intervals
	for(std::size_t i : {std::numeric_limits<std::size_t>::min(), std::size_t(0), std::numeric_limits<std::size_t>::max()}) {
		BOOST_TEST(in_int_ranges(
			int_ranges<std::size_t>{{0, 0}},
			i
		) == false, i << " is in int_ranges{{0, 0}}");
		
		BOOST_TEST(in_int_ranges(
			int_ranges<std::size_t>{{std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()}},
			i
		) == false, i << " is in int_ranges{{numeric_limits<size_t>::max(), numeric_limits<size_t>::max(}}");
	}
}

BOOST_AUTO_TEST_SUITE_END()
