/* Copyright (c) 2016-2020 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_DETAIL_INT_RANGES_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_INT_RANGES_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>
#include <vector>
#include <algorithm>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace detail {

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

template<typename Integer>
using int_range = boost::integer_range<Integer>;

template<typename Integer>
using int_ranges = std::vector<int_range<Integer>> ;

template<typename Integer, typename Enable = void>
struct int_parser;

template<typename Integer>
struct int_parser<Integer, std::enable_if_t<std::is_unsigned<Integer>::value>>
: qi::uint_parser<Integer>{};

template<typename Integer>
struct int_parser<Integer, std::enable_if_t<std::is_signed<Integer>::value>>
: qi::int_parser<Integer>{};

template <typename Iterator, typename Integer>
struct int_ranges_grammar : qi::grammar<Iterator, int_ranges<Integer>()> {

	int_ranges_grammar() : int_ranges_grammar::base_type{query} {
		auto parse_range = [](
			fusion::vector<Integer, Integer> const& attr,
			auto & ctx, bool & pass
		) {
			int_ranges<Integer> & seqs = fusion::at_c<0>(ctx.attributes);
			
			auto low = fusion::at_c<0>(attr);
			auto high = fusion::at_c<1>(attr);
			
			if (low > high) {
				pass = false;
				return;
			}
			
			seqs.push_back(int_range<Integer>{low, high+1});
			pass = true;
		};
		auto parse_no = [](
			auto && attr,
			auto & ctx, bool & pass
		) {
			int_ranges<Integer> & seqs = fusion::at_c<0>(ctx.attributes);
			
			Integer no;
			if constexpr(std::is_arithmetic<std::decay_t<decltype(attr)>>::value) {
				no = attr;
			} else {
				no = fusion::at_c<0>(attr);
			}
			seqs.push_back(int_range<Integer>{no,no+1});
			pass = true;
		};
		
		query = (numbers >> (*((qi::lit(',') | ';') >> numbers)));
		numbers =  (int_ >> '-' >> int_)[parse_range] | int_[parse_no];
	}
	
	int_parser<Integer> const int_{};
	qi::rule<Iterator, int_ranges<Integer>()> query;
	qi::rule<Iterator, int_ranges<Integer>()> numbers;
};

template <typename Iterator, typename Integer>
bool
parse_int_ranges(
	Iterator const& first,
	Iterator const& last,
	int_ranges<Integer> & s
) {
	return qi::parse(first, last, detail::int_ranges_grammar<Iterator, Integer>{}, s);
}

template <typename Integer>
bool
in_int_ranges(int_ranges<Integer> const& s, Integer const& val) {
	for(auto const& rng : s) {
		if (
			!rng.empty() /* nothing is in an empty interval */ &&
			(rng.front() <= val && val <= rng.back())
		) {
			return true;
		}
	}
	return false;
}

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_INT_RANGES_IMPL_HPP
