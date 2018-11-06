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
#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <matlab/dt/matrix.hpp>
#include <elemental/dt/matrix.hpp>
#include <vector>
#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;

//TODO: Uncomment and implement!
// mpl::pca_filter_result<
// 	matrix::Matrix<double> /* data */,
// 	vector::Vector<double> /* latent*/
// >
// pca_filter_impl1::operator()(matrix::Matrix<double> a, mpl::rtsav<bool> const& keep) const {
// 	auto && dec = pca2::pca(a, true);
// 	auto && coeff = dec.coeff();
// 	auto && score = dec.score();
// 	auto && latent = dec.latent();
// 	auto && mu = dec.mean();
// 	
// 	auto && size = matrix::size_of(coeff);
// 	auto && m = std::get<0>(size);
// 	auto && n = std::get<1>(size);
// 	
// 	for(std::size_t i = 0; i < n; ++i) {
// 		if(!keep[i]) {
// 			for(std::size_t j = 0; j < m; ++j) {
// 				coeff[j][i] = 0;
// 			}
// 		}
// 	}
// 	
// 	typedef double Real;
// 	auto && reduced = matrix::plus<Real>(
// 		matrix::product<Real>(
// 			score, matrix::transpose<Real>(coeff)
// 		), 
// 		mu,
// 		matrix::VectorType::RowVector
// 	);
// 	
// 	return { reduced, latent };
// }
// 
// mpl::pca_filter_result<
// 	std::vector<theta_field> /* data */,
// 	std::vector<double> /* latent*/
// >
// pca_filter_impl2::operator()(std::vector<theta_field> series, mpl::rtsav<bool> const& keep) const {
// 	
// 	BOOST_ASSERT(series.size() > 0);
// 	std::size_t const m = series[0].x_velocity().size() + series[0].y_velocity().size() + series[0].z_velocity().size();
// 	std::size_t const n = series.size();
// 	BOOST_ASSERT(m > 0);
// 	BOOST_ASSERT(n > 0);
// 	
// 	matrix::Matrix<double> mat { boost::extents[m][n] };
// 	typedef matrix::Matrix<double>::index _index;
// 	typedef boost::multi_array_types::index_range _range;
// 	
// 	matrix::Matrix<double>::index_gen indices;
// 	for (_index i = 0; i < (signed)n; ++i) {
// 		auto column = mat[ indices[_range()][i] ];
// 		
// 		std::copy(series[i].x_velocity().begin(), series[i].x_velocity().end(), column.begin());
// 		std::copy(series[i].y_velocity().begin(), series[i].y_velocity().end(), column.begin()+series[i].x_velocity().size());
// 		std::copy(series[i].z_velocity().begin(), series[i].z_velocity().end(), column.begin()+series[i].x_velocity().size()+series[i].y_velocity().size());
// 	}
// 	
// 	auto && reduced = pca_filter_impl1{}(mat, keep);
// 	auto && data_ = reduced.data();
// 	auto && latent_ = reduced.latent();
// 	
// 	std::vector<double> latent__;
// 	latent__.resize(latent_.shape()[0]);
// 	
// 	std::copy(latent_.begin(), latent_.end(), latent__.data());
// 	
// 	for (_index i = 0; i < (signed)n; ++i) {
// 		auto && column = data_[ indices[_range()][i] ];
// 		
// 		std::copy(column.begin()+0*series[i].x_velocity().size(), column.begin()+1*series[i].x_velocity().size(), series[i].x_velocity().begin());
// 		std::copy(column.begin()+1*series[i].x_velocity().size(), column.begin()+2*series[i].x_velocity().size(), series[i].y_velocity().begin());
// 		std::copy(column.begin()+2*series[i].x_velocity().size(), column.begin()+3*series[i].x_velocity().size(), series[i].z_velocity().begin());
// 	}
// 	
// 	return { series, latent__ };
// }

template<typename From>
decltype(auto)
copy_vector(From const& from, std::vector<double> && to) {
	auto from_sz = (*mpl::size)(from);
	auto to_sz = (*mpl::size)(to);
	
	BOOST_ASSERT((*mpl::equal)(from_sz, to_sz));
	
	for (std::size_t i = 0; i < from_sz; ++i) {
		to[i] = (*mpl::at)(from, i);
	}
	
	return HBRS_MPL_FWD(to);
}

template<typename To>
mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, std::vector<bool> const& keep, hana::basic_type<To>) {
	auto sz = detail::size(series);
	auto m = (*mpl::m)(sz);
	auto n = (*mpl::n)(sz);
	auto latent_sz = (m-1)<n ? m-1 : std::min(m,n);
	BOOST_ASSERT((*mpl::equal)(keep.size(), latent_sz));
	
	using M = decltype( (*mpl::m)(mpl::size(std::declval<To>())) );
	using N = decltype( (*mpl::n)(mpl::size(std::declval<To>())) );
	auto to_m = boost::numeric_cast<M>(m);
	auto to_n = boost::numeric_cast<N>(n);
	
	auto mat = detail::copy_matrix(series, To{to_m, to_n});
	auto reduced = (*mpl::pca_filter)(mat, keep);
	
	series = detail::copy_matrix(reduced.data(), series);
	auto latent = copy_vector(reduced.latent(), std::vector<double>(latent_sz));
	
	return {series, latent};
}

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, std::vector<bool> const& keep, matlab_lapack_backend) {
	return pca_filter(series, keep, hana::type_c<matlab::matrix<real_T>>);
}

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, std::vector<bool> const& keep, elemental_openmp_backend) {
	return pca_filter(series, keep, hana::type_c<El::Matrix<double>>);
}

mpl::pca_filter_result<
	std::vector<theta_field> /* data */,
	std::vector<double> /* latent*/
>
pca_filter(std::vector<theta_field> series, std::vector<bool> const& keep, elemental_mpi_backend) {
	BOOST_ASSERT_MSG(false, "TODO");
	return {}; //TODO
}

HBRS_THETA_UTILS_NAMESPACE_END