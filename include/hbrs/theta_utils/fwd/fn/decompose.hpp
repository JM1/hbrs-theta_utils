/* Copyright (c) 2016-2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_FWD_FN_DECOMPOSE_HPP
#define HBRS_THETA_UTILS_FWD_FN_DECOMPOSE_HPP

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fwd/fn/visualize.hpp>
#include <boost/hana/integral_constant.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

enum class pca_backend { matlab_lapack, elemental_openmp, elemental_mpi };
template <pca_backend backend>
using pca_backend_ = hana::integral_constant<pca_backend, backend>;
template <pca_backend backend>
constexpr pca_backend_<backend> pca_backend_c{};

using matlab_lapack_backend = pca_backend_<pca_backend::matlab_lapack>;
constexpr auto matlab_lapack_backend_c = matlab_lapack_backend{};
using elemental_openmp_backend = pca_backend_<pca_backend::elemental_openmp>;
constexpr auto elemental_openmp_backend_c = elemental_openmp_backend{};
using elemental_mpi_backend = pca_backend_<pca_backend::elemental_mpi>;
constexpr auto elemental_mpi_backend_c = elemental_mpi_backend{};

struct pca_decomposition_options {
	std::vector<std::string> pc_nr_seqs;
	pca_backend backend;
};

void
decompose_with_pca(
	visualize_options v_opts,
	pca_decomposition_options pca_opts,
	bool verbose
);

HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_FWD_FN_DECOMPOSE_HPP
