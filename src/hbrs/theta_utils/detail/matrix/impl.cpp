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

#include <hbrs/theta_utils/dt/theta_field_matrix.hpp>

namespace boost { namespace hana {

#ifdef HBRS_MPL_ENABLE_MATLAB

HBRS_THETA_UTILS_API
hbrs::mpl::ml_matrix<double>
to_impl<
	hbrs::mpl::ml_matrix_tag,
	hbrs::theta_utils::theta_field_matrix_tag
>::apply(hbrs::theta_utils::theta_field_matrix const& rhs) {
	namespace mpl = hbrs::mpl;
	
	auto sz = rhs.size();
	int m_ = boost::numeric_cast<int>(mpl::m(sz));
	int n_ = boost::numeric_cast<int>(mpl::n(sz));
	mpl::ml_matrix<double> lhs{ m_, n_ };
	hbrs::theta_utils::detail::copy_matrix(rhs, lhs);
	return lhs;
}
#endif // !HBRS_MPL_ENABLE_MATLAB

#ifdef HBRS_MPL_ENABLE_ELEMENTAL

HBRS_THETA_UTILS_API
hbrs::mpl::el_matrix<double>
to_impl<
	hbrs::mpl::el_matrix_tag,
	hbrs::theta_utils::theta_field_matrix_tag
>::apply(hbrs::theta_utils::theta_field_matrix const& rhs) {
	namespace mpl = hbrs::mpl;
	
	auto sz = rhs.size();
	El::Int m_ = boost::numeric_cast<El::Int>(mpl::m(sz));
	El::Int n_ = boost::numeric_cast<El::Int>(mpl::n(sz));
	mpl::el_matrix<double> lhs{ m_, n_ };
	hbrs::theta_utils::detail::copy_matrix(rhs, lhs);
	return lhs;
}
#endif // !HBRS_MPL_ENABLE_ELEMENTAL

/* namespace hana */ } /* namespace boost */ }
