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

#ifndef HBRS_THETA_UTILS_FN_PCA_FILTER_FWD_HPP
#define HBRS_THETA_UTILS_FN_PCA_FILTER_FWD_HPP

#include <hbrs/theta_utils/config.hpp>

#include <hbrs/mpl/dt/pca_filter_result/fwd.hpp>
#include <hbrs/theta_utils/detail/int_ranges.hpp>
#include <hbrs/theta_utils/dt/command_option/fwd.hpp>
#include <hbrs/theta_utils/dt/theta_field/fwd.hpp>
#include <vector>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace mpl = hbrs::mpl;

//TODO: Make pca_filter() a generic function using HBRS_MPL_DEC_F1 macros!

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, detail::int_ranges<std::size_t> const& keep, matlab_lapack_backend);

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, detail::int_ranges<std::size_t> const& keep, elemental_openmp_backend);

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, detail::int_ranges<std::size_t> const& keep, elemental_mpi_backend);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FN_PCA_FILTER_FWD_HPP
