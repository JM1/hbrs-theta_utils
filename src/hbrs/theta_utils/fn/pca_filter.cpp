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

#include <hbrs/theta_utils/fn/pca_filter.hpp>
#include <hbrs/theta_utils/detail/copy_matrix.hpp>
#include <hbrs/mpl/fn/pca_filter.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/dt/exception.hpp>

#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_ADDON_MATLAB
    #include <matlab/dt/matrix.hpp>
#endif
#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
    #include <elemental/dt/matrix.hpp>
    #include <elemental/dt/dist_matrix.hpp>
#endif

#include <vector>
#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/throw_exception.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace mpi = hbrs::mpl::detail::mpi;

namespace {

#if defined(HBRS_MPL_ENABLE_ADDON_MATLAB) || defined(HBRS_MPL_ENABLE_ADDON_ELEMENTAL)
	template<
		typename From,
		typename std::enable_if_t<
			#ifdef HBRS_MPL_ENABLE_ADDON_MATLAB
				std::is_same< hana::tag_of_t<From>, matlab::column_vector_tag >::value ||
			#endif
			#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
				std::is_same< hana::tag_of_t<From>, elemental::column_vector_tag >::value ||
			#endif
			false
		>* = nullptr
	>
	auto
	to_vector(From const& from) {
		auto from_sz = (*mpl::size)(from);
		std::vector<double> to(from_sz);
		BOOST_ASSERT(to.size() == from_sz);
		
		for (std::size_t i = 0; i < from_sz; ++i) {
			to[i] = (*mpl::at)(from, i);
		}
		
		return to;
	}
#endif

#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
	template<
		typename Matrix,
		typename std::enable_if_t<
			std::is_same< hana::tag_of_t<Matrix>, hana::ext::El::DistMatrix_tag>::value
		>* = nullptr
	>
	auto
	to_vector(elemental::dist_column_vector<Matrix> const& from) {
		auto from_local = from.data().LockedMatrix();
		auto from_sz = from_local.Height();
		BOOST_ASSERT(from_local.Width() == 1);
		
		std::vector<double> to(from_sz);
		BOOST_ASSERT(to.size() == from_sz);
		
		for (std::size_t i = 0; i < from_sz; ++i) {
			to[i] = from_local.Get(i, 0);
		}
		
		return to;
	}
#endif

#ifdef HBRS_MPL_ENABLE_ADDON_MATLAB
	auto
	make_matrix(hana::basic_type<matlab::matrix_tag>, mpl::matrix_size<int, int> sz) {
		return matlab::make_matrix(hana::type_c<double>, sz);
	}
#endif

#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
	auto
	make_matrix(hana::basic_type<hana::ext::El::Matrix_tag>, mpl::matrix_size<El::Int, El::Int> sz) {
		return elemental::make_matrix(hana::type_c<double>, sz);
	}

	auto
	make_matrix(hana::basic_type<hana::ext::El::DistMatrix_tag>, mpl::matrix_size<El::Int, El::Int> sz) {
		static El::Grid const grid{El::mpi::COMM_WORLD};
		

		auto min_m = El::mpi::AllReduce(sz.m(), El::mpi::MIN, grid.Comm());
		auto max_m = El::mpi::AllReduce(sz.m(), El::mpi::MAX, grid.Comm());
		BOOST_ASSERT(min_m <= max_m || min_m > max_m);
		auto full_m = max_m * mpi::size(); //TODO: Is this correct?

	// 	auto full_m = El::mpi::AllReduce(sz.m(), El::mpi::SUM, grid.Comm());
	// 	BOOST_ASSERT(min_m * mpi::size() == full_m);	
		
		auto min_n = El::mpi::AllReduce(sz.n(), El::mpi::MIN, grid.Comm());
		auto max_n = El::mpi::AllReduce(sz.n(), El::mpi::MAX, grid.Comm());
		BOOST_ASSERT(min_n == max_n);
		
		El::DistMatrix<double, El::VC, El::STAR, El::ELEMENT> a{grid};
		a.Resize(full_m, sz.n());
		
		BOOST_ASSERT(a.Matrix().Width() == sz.n());
		return a;
	}
#endif

/* unnamed namespace */ }

template<typename ToTag>
mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, detail::int_ranges<std::size_t> const& keep, hana::basic_type<ToTag>) {
	auto sz = detail::size(series);
	auto m = (*mpl::m)(sz);
	auto n = (*mpl::n)(sz);
	
	auto mat = detail::copy_matrix(series, make_matrix(hana::type_c<ToTag>, {m, n}));
	auto mat_sz = (*mpl::size)(mat);
	auto mat_m = (*mpl::m)(mat_sz);
	auto mat_n = (*mpl::n)(mat_sz);
	auto latent_sz = (mat_m-1)<mat_n ? mat_m-1 : std::min(mat_m,mat_n);
	
	auto reduced = (*mpl::pca_filter)(
		mat,
		[&keep](std::size_t i) {
			return detail::in_int_ranges(keep, i); 
		}
	);
	
	series = detail::copy_matrix(reduced.data(), series);

	BOOST_ASSERT((*mpl::equal)(mpl::size(reduced.latent()), latent_sz));
	auto latent = to_vector(reduced.latent());
	
	return {series, latent};
}

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, detail::int_ranges<std::size_t> const& keep, matlab_lapack_backend) {
	#ifdef HBRS_MPL_ENABLE_ADDON_MATLAB
		return pca_filter(series, keep, hana::type_c<matlab::matrix_tag>);
	#else
		BOOST_THROW_EXCEPTION(invalid_backend_exception{} << errinfo_pca_backend{matlab_lapack_backend_c});
	#endif
}

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, detail::int_ranges<std::size_t> const& keep, elemental_openmp_backend) {
	#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
		return pca_filter(series, keep, hana::type_c<hana::ext::El::Matrix_tag>);
	#else
		BOOST_THROW_EXCEPTION(invalid_backend_exception{} << errinfo_pca_backend{elemental_openmp_backend_c});
	#endif
}

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, detail::int_ranges<std::size_t> const& keep, elemental_mpi_backend) {
	#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
		return pca_filter(series, keep, hana::type_c<hana::ext::El::DistMatrix_tag>);
	#else
		BOOST_THROW_EXCEPTION(invalid_backend_exception{} << errinfo_pca_backend{elemental_mpi_backend_c});
	#endif
}

HBRS_THETA_UTILS_NAMESPACE_END