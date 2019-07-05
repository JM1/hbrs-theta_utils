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

#ifndef HBRS_THETA_UTILS_CORE_PREPROCESSOR_FWD_HPP
#define HBRS_THETA_UTILS_CORE_PREPROCESSOR_FWD_HPP

#define HBRS_THETA_UTILS_DECLARE_ATTR(__name, __type)                                                                  \
public:                                                                                                                \
	__type &      __name() &      ;                                                                                    \
	__type const& __name() const& ;                                                                                    \
	__type &&     __name() &&     ;                                                                                    \
private:                                                                                                               \
	__type __name ## _;

#define HBRS_THETA_UTILS_DEFINE_ATTR(__name, __type, __class)                                                          \
	__type &      __class::__name() &       { return                     (__class::__name ## _) ; }                    \
	__type const& __class::__name() const & { return                     (__class::__name ## _) ; }                    \
	__type &&     __class::__name() &&      { return std::forward<__type>(__class::__name ## _); }


#endif // !HBRS_THETA_UTILS_CORE_PREPROCESSOR_FWD_HPP
