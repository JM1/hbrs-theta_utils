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

#pragma once

#ifndef HBRS_THETA_UTILS_DT_NUMBER_SEQUENCE_HPP
#define HBRS_THETA_UTILS_DT_NUMBER_SEQUENCE_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/dt/number_sequence.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <vector>
#include <algorithm>
#include <boost/assert.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

template <typename Iterator>
struct number_sequence : qi::grammar<Iterator, std::vector<unsigned long>()> {
	
	number_sequence() : number_sequence::base_type(query) {
		auto parse_range = [](
			fusion::vector<unsigned long, unsigned long> const& attr, 
			auto & ctx, bool & pass
		) { 
			std::vector<unsigned long> & val = fusion::at_c<0>(ctx.attributes);
			BOOST_ASSERT(val.size() == 0);
			
			auto low = fusion::at_c<0>(attr);
			auto high = fusion::at_c<1>(attr);
			
			if (low > high) {
				pass = false;
				return;
			}
			
			val.reserve(high-low+1);
			
			for(;low != high+1; ++low) {
				val.push_back(low);
			}
			
			pass = true;
		};
		
		auto parse_no = [](fusion::vector<unsigned long> const& attr, auto & ctx, bool & pass) {
			std::vector<unsigned long> & val = fusion::at_c<0>(ctx.attributes);
			BOOST_ASSERT(val.size() == 0);
			
			unsigned long no = fusion::at_c<0>(attr);
			val.push_back(no);
			pass = true;
		};
		
		auto append = [](
			fusion::vector< std::vector<unsigned long>, std::vector<std::vector<unsigned long>> > const& attr, 
			auto & ctx, bool & pass
		) {
			std::vector<unsigned long> & val = fusion::at_c<0>(ctx.attributes);
			
			std::vector<unsigned long> numbers = fusion::at_c<0>(attr);
			std::vector<std::vector<unsigned long>> more_numbers = fusion::at_c<1>(attr);
			
			val = std::accumulate(more_numbers.begin(), more_numbers.end(), numbers, 
				[](std::vector<unsigned long> v1, std::vector<unsigned long> v2){
					return std::accumulate(v2.begin(), v2.end(), v1, 
						[](std::vector<unsigned long> v, unsigned long e) { v.push_back(e); return v; }
					);
				}
			);
			
			pass = true;
		};
		
		auto sort_uniq = [](fusion::vector<std::vector<unsigned long>> const& attr, auto & ctx, bool & pass) {
			std::vector<unsigned long> & val = fusion::at_c<0>(ctx.attributes);
			std::vector<unsigned long> numbers = fusion::at_c<0>(attr);
			
			std::sort(numbers.begin(), numbers.end());
			std::unique(numbers.begin(), numbers.end());
			
			val = numbers;
			pass = true;
		};
		
		query =  numbers0[sort_uniq];
		numbers0 = (numbers >> (*((qi::lit(',') | ';') >> numbers)))[append];
		numbers =  (qi::ulong_ >> '-' >> qi::ulong_)[parse_range] | qi::ulong_[parse_no];
		
	}

	qi::rule<Iterator, std::vector<unsigned long>()> query;
	qi::rule<Iterator, std::vector<unsigned long>()> numbers;
	qi::rule<Iterator, std::vector<unsigned long>()> numbers0;
};


template <typename Iterator>
bool
parse_number_sequence(
	Iterator const& first,
	Iterator const& last,
	std::vector<unsigned long> & v
) {
	return qi::parse(first, last, number_sequence<Iterator>{}, v);
}

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DT_NUMBER_SEQUENCE_HPP
