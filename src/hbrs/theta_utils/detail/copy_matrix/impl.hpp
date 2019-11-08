/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
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

#ifndef HBRS_THETA_UTILS_DETAIL_COPY_MATRIX_IMPL_HPP
#define HBRS_THETA_UTILS_DETAIL_COPY_MATRIX_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/mpl/core/preprocessor.hpp>

#include <hbrs/mpl/dt/matrix_size.hpp>
#include <hbrs/mpl/dt/rtsam.hpp>
#include <hbrs/mpl/dt/ctsam.hpp>
#include <hbrs/mpl/dt/sm.hpp>

#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/less_equal.hpp>
#include <hbrs/mpl/fn/greater_equal.hpp>


#include <hbrs/mpl/config.hpp>
#ifdef HBRS_MPL_ENABLE_MATLAB
    #include <hbrs/mpl/dt/ml_matrix.hpp>
#endif
#ifdef HBRS_MPL_ENABLE_ELEMENTAL
    #include <hbrs/mpl/dt/el_matrix.hpp>
    #include <hbrs/mpl/dt/el_dist_matrix.hpp>
#endif

#include <hbrs/mpl/detail/mpi.hpp>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/assert.hpp>
#include <algorithm>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpl = hbrs::mpl;
namespace mpi = hbrs::mpl::detail::mpi;
namespace detail {

template<
	typename To,
	typename std::enable_if_t<
		#ifdef HBRS_MPL_ENABLE_MATLAB
			std::is_same< hana::tag_of_t<To>, mpl::ml_matrix_tag >::value ||
		#endif
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
			std::is_same< hana::tag_of_t<To>, mpl::el_matrix_tag >::value ||
		#endif
		std::is_same< hana::tag_of_t<To>, mpl::rtsam_tag >::value ||
		std::is_same< hana::tag_of_t<To>, mpl::ctsam_tag >::value ||
		std::is_same< hana::tag_of_t<To>, mpl::sm_tag >::value
	>* = nullptr
>
decltype(auto)
copy_matrix(std::vector<theta_field> const& from, To && to) {
	auto from_sz = local_size(from);
	auto from_m = (*mpl::m)(from_sz);
	auto from_n = (*mpl::n)(from_sz);
	
	auto to_sz = (*mpl::size)(to);
	auto to_m = (*mpl::m)(to_sz);
	auto to_n = (*mpl::n)(to_sz);
	
	BOOST_ASSERT((*mpl::less_equal)(from_m, to_m));
	BOOST_ASSERT((*mpl::equal)(from_n, to_n));
	
	for (std::size_t j = 0; j < from_n; ++j) {
		BOOST_ASSERT(from.at(j).density().empty());
		BOOST_ASSERT(from.at(j).pressure().empty());
		BOOST_ASSERT(from.at(j).residual().empty());
	
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i, j)) = from.at(j).x_velocity().at(i);
		}
		
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i+from_m/3, j)) = from.at(j).y_velocity().at(i);
		}
		
		for(std::size_t i = 0; i < from_m/3; ++i) {
			(*mpl::at)(to, mpl::make_matrix_index(i+from_m/3*2, j)) = from.at(j).z_velocity().at(i);
		}
	}
	
	return HBRS_MPL_FWD(to);
}

template<typename T>
std::ostream &
operator<<(std::ostream & o, std::vector<T> const& v) {
	o << '[';
	if (!v.empty()) {
		o << "0:" << v[0];
		for(std::size_t i = 1; i < v.size(); ++i) {
			o << ", " << i << ':' << v[i];
		}
	}
	o << ']';
	return o;
}

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
template<typename T>
std::ostream &
operator<<(std::ostream & o, El::Matrix<T> const& m) {
	El::Print(m, "Matrix", o);
	return o;
}

template<typename T>
std::ostream &
operator<<(std::ostream & o, El::AbstractDistMatrix<T> const& m) {
	El::Print(m, "DistMatrix", o);
	return o;
}
#endif //!HBRS_MPL_ENABLE_ELEMENTAL

template<typename T, typename... Ts>
void
debug_output(T && t, Ts && ... ts) {
	#ifdef HBRS_THETA_UTILS_DEBUG_OUTPUT
		std::cerr << HBRS_MPL_FWD(t);
		if constexpr(sizeof...(ts) > 0) {
			debug_output(HBRS_MPL_FWD(ts)...);
		}
		std::cerr << std::flush;
	#endif // !HBRS_THETA_UTILS_DEBUG_OUTPUT
}

bool
iff(bool lhs, bool rhs);

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	inline decltype(auto)
	copy_matrix(
		std::vector<theta_field> const& from,
		hbrs::mpl::el_dist_matrix<double, El::VC, El::STAR, El::ELEMENT> && to
	) {
		using std::size_t;
		
		#if !defined(NDEBUG)
		{
			mpl::matrix_size<size_t, size_t> from_lcl_sz = local_size(from);
			mpl::matrix_size<size_t, size_t> from_gbl_sz = global_size(from);
			mpl::matrix_size<size_t, size_t> to_sz   = {to.size()};
			BOOST_ASSERT(from_lcl_sz.n() == to_sz.n());
			BOOST_ASSERT(from_gbl_sz.n() == to.data().Matrix().Width());
			BOOST_ASSERT(from_gbl_sz.m() == to.data().Height());
		}
		#endif
		
		size_t mpi_sz = boost::numeric_cast<size_t>(mpi::size());
		size_t mpi_rank = boost::numeric_cast<size_t>(mpi::rank());
		mpl::matrix_size<size_t, size_t> lcl_sz = local_size(from);
		
		/* Example for 3 processes:
		 * 
		 * mpi_sz = 3
		 * 
		 * lcl_to@proc0=  [  1  2  3
		 *                   4  5  6
		 *                   7  8  9
		 *                  10 11 12 ]
		 * lcl_to@proc1 = [ 13 14 15 ]
		 * lcl_to@proc2 = [ 16 17 18
		 *                  19 20 21 ]
		 * 
		 * lcl_m@proc0 = 4
		 * lcl_m@proc1 = 1
		 * lcl_m@proc2 = 2
		 *
		 * lcl_n@proc0 = 3
		 * lcl_n@proc1 = 3
		 * lcl_n@proc2 = 3
		 * 
		 * gbl_m = 7
		 * 
		 * balanced_chunk_size = 2
		 * nr_of_bigger_chunks = 1
		 * nr_of_normal_chunks = 2
		 * 
		 * lcl_ms_sums     = [ 4 5 7 ]
		 * bal_lcl_ms_sums = [ 3 5 7 ]
		 * 
		 * block distribution:
		 * | gbl_row | send_proc | recv_proc | lcl_send_row | lcl_recv_row |
		 * |---------|-----------|-----------|--------------|--------------|
		 * |       0 |         0 |         0 |            0 |            0 |
		 * |       1 |         0 |         0 |            1 |            1 |
		 * |       2 |         0 |         0 |            2 |            2 |
		 * |       3 |         0 |         1 |            3 |            0 |
		 * |       4 |         1 |         1 |            0 |            1 |
		 * |       5 |         2 |         2 |            0 |            0 |
		 * |       6 |         2 |         2 |            1 |            1 |
		 * 
		 * bal_lcl_to@proc0=  [  1  2  3
		 *                       4  5  6
		 *                       7  8  9
		 * bal_lcl_to@proc1 = [ 10 11 12
		 *                      13 14 15 ]
		 * bal_lcl_to@proc2 = [ 16 17 18
		 *                      19 20 21 ]
		 * 
		 * round robin distribution:
		 * | gbl_row | send_proc | recv_proc | lcl_send_row | lcl_recv_row |
		 * |---------|-----------|-----------|--------------|--------------|
		 * |       0 |         0 |         0 |            0 |            0 |
		 * |       1 |         0 |         1 |            1 |            0 |
		 * |       2 |         0 |         2 |            2 |            0 |
		 * |       3 |         0 |         0 |            3 |            1 |
		 * |       4 |         1 |         1 |            0 |            1 |
		 * |       5 |         2 |         2 |            0 |            1 |
		 * |       6 |         2 |         0 |            1 |            2 |
		 * 
		 * bal_lcl_to@proc0 = [ 1  2  3
		 *                      10 11 12
		 *                      19 20 21 ]
		 * bal_lcl_to@proc1 = [ 4  5  6
		 *                      13 14 15 ]
		 * bal_lcl_to@proc2 = [ 7  8  9
		 *                      16 17 18 ]
		 */
		
		size_t lcl_m = lcl_sz.m();
		size_t lcl_n = lcl_sz.n();
		std::vector<size_t> lcl_ms(mpi_sz, 0u);
		mpi::allgather(&lcl_m, 1, lcl_ms.data(), 1, to.data().Grid().Comm().comm);
		debug_output("lcl_ms:", lcl_ms, "@mpi_rank:", mpi_rank, '\n');
		
		std::vector<size_t> lcl_ms_sums(mpi_sz, 0u);
		std::partial_sum(lcl_ms.begin(), lcl_ms.end(), lcl_ms_sums.begin());
		debug_output("lcl_ms_sums:", lcl_ms_sums, "@mpi_rank:", mpi_rank, '\n');
		
		size_t gbl_m = lcl_ms_sums.back();
		debug_output("gbl_m:", gbl_m, "@mpi_rank:", mpi_rank, '\n');
		
		size_t balanced_chunk_size = gbl_m / mpi_sz; // integer division
		size_t nr_of_bigger_chunks = gbl_m % mpi_sz; // number of processes that get chunks of size=chunk_size+1
		size_t nr_of_normal_chunks = mpi_sz - nr_of_bigger_chunks;
		BOOST_ASSERT((nr_of_bigger_chunks+nr_of_normal_chunks) == mpi_sz);
		
		std::vector<size_t> bal_lcl_ms(mpi_sz, 0u);
		for(size_t i = 0; i < nr_of_bigger_chunks; ++i) {
			bal_lcl_ms.at(i) = balanced_chunk_size+1;
		}
		for(size_t i = nr_of_bigger_chunks; i < mpi_sz; ++i) {
			bal_lcl_ms.at(i) = balanced_chunk_size;
		}
		std::vector<size_t> bal_lcl_ms_sums(mpi_sz, 0u);
		std::partial_sum(bal_lcl_ms.begin(), bal_lcl_ms.end(), bal_lcl_ms_sums.begin());
		
		debug_output("bal_lcl_ms:", bal_lcl_ms, "@mpi_rank:", mpi_rank, '\n');
		debug_output("bal_lcl_ms_sums:", bal_lcl_ms_sums, "@mpi_rank:", mpi_rank, '\n');
		
		size_t bal_lcl_m = bal_lcl_ms.at(mpi_rank);
		mpl::rtsam<double, mpl::storage_order::row_major> bal_lcl_to =
			{ mpl::make_matrix_size(bal_lcl_m, lcl_n) };
		
		{
			mpl::rtsam<double, mpl::storage_order::row_major> lcl_to =
				copy_matrix(from, mpl::rtsam<double, mpl::storage_order::row_major>{lcl_sz});
			std::vector<MPI_Request> reqs;
			for(size_t gbl_row = 0; gbl_row<gbl_m; ++gbl_row) {
				size_t send_proc = std::distance(
					lcl_ms_sums.begin(),
					std::lower_bound(lcl_ms_sums.begin(), lcl_ms_sums.end(), gbl_row+1)
				);
				
				// block distribution
// 				size_t recv_proc = std::distance(
// 					bal_lcl_ms_sums.begin(),
// 					std::lower_bound(bal_lcl_ms_sums.begin(), bal_lcl_ms_sums.end(), gbl_row)
// 				);
				// round robin distribution
				size_t recv_proc = gbl_row % mpi_sz;
				
				
				size_t lcl_send_row = (send_proc == 0) ? gbl_row : (gbl_row - lcl_ms_sums.at(send_proc-1));
				// block distribution
// 				size_t lcl_recv_row = (recv_proc == 0) ? gbl_row : (gbl_row - bal_lcl_ms_sums.at(recv_proc-1));
				// round robin distribution
				size_t lcl_recv_row = gbl_row / mpi_sz;
				
				debug_output("gbl_row:", gbl_row, ",send_proc:",send_proc, ",recv_proc:", recv_proc, "@mpi_rank:", mpi_rank, '\n');
				debug_output("gbl_row:", gbl_row, ",lcl_send_row:",lcl_send_row, ",lcl_recv_row:", lcl_recv_row, "@mpi_rank:", mpi_rank, '\n');
				
				if ((send_proc == recv_proc) && (send_proc == mpi_rank)) {
					debug_output("gbl_row:", gbl_row, ", COPY", "@mpi_rank:", mpi_rank, '\n');
					
					double const& lcl_send_ref = lcl_to.at(mpl::make_matrix_index(lcl_send_row, 0));
					double      & lcl_recv_ref = bal_lcl_to.at(mpl::make_matrix_index(lcl_recv_row, 0));
					std::copy_n(&lcl_send_ref, lcl_n, &lcl_recv_ref);
					
				} else if (send_proc == mpi_rank) {
					debug_output("gbl_row:", gbl_row, ", SEND_ROW", "@mpi_rank:", mpi_rank, '\n');
					
					double const& lcl_send_ref = lcl_to.at(mpl::make_matrix_index(lcl_send_row, 0));
					reqs.push_back(
						mpi::isend(
							&lcl_send_ref,
							lcl_n,
							recv_proc /*dest*/,
							send_proc /*tag*/,
							MPI_COMM_WORLD
						)
					);
				} else if (recv_proc == mpi_rank) {
					debug_output("gbl_row:", gbl_row, ", RECV_ROW", "@mpi_rank:", mpi_rank, '\n');
					
					double & lcl_recv_ref = bal_lcl_to.at(mpl::make_matrix_index(lcl_recv_row, 0));
					reqs.push_back(
						mpi::irecv(
							&lcl_recv_ref,
							lcl_n,
							send_proc /*source*/,
							send_proc /*tag*/,
							MPI_COMM_WORLD
						)
					);
				} else {
					debug_output("gbl_row:", gbl_row, ", SKIP", "@mpi_rank:", mpi_rank, '\n');
				}
			}
			
			for(size_t i = 0; i < reqs.size(); ++i) {
				debug_output("reqs[", i, "], WAIT_BEGIN", "@mpi_rank:", mpi_rank, '\n');
				
				[[maybe_unused]] auto stat = mpi::wait(reqs[i]);
				//TODO: Do anything with stat?
				debug_output("reqs[", i, "], WAIT_END", "@mpi_rank:", mpi_rank, '\n');
			}
		}
		
		El::DistMatrix<double, El::VC, El::STAR, El::ELEMENT> & gbl_to = to.data();
		El::Matrix<double> & lcl_to = gbl_to.Matrix();
		
		BOOST_ASSERT(gbl_to.Matrix().Height() == bal_lcl_to.size().m());
		BOOST_ASSERT(gbl_to.Matrix().Width() == bal_lcl_to.size().n());
		
		for(size_t gbl_row = 0; gbl_row<gbl_m; ++gbl_row) {
			size_t lcl_row = gbl_row / mpi_sz;
			bool is_lcl_row = ((gbl_row % mpi_sz) == mpi_rank);
			BOOST_ASSERT(iff(is_lcl_row, gbl_to.IsLocalRow(gbl_row)));
			
			if (is_lcl_row) {
				debug_output("lcl_row:", lcl_row, "@mpi_rank:", mpi_rank, '\n');
				
				for(size_t col = 0; col < lcl_n; ++col) {
					BOOST_ASSERT(gbl_to.IsLocal(gbl_row, col));
					
					gbl_to.SetLocal(
						gbl_to.LocalRow((El::Int)gbl_row),
						gbl_to.LocalCol((El::Int)col),
						bal_lcl_to.at(mpl::make_matrix_index(lcl_row, col))
					);
				}
			}
			
			debug_output("lcl_to:", lcl_to, "@mpi_rank:", mpi_rank, '\n');
			debug_output("gbl_to:", gbl_to, "@mpi_rank:", mpi_rank, '\n');
		}
		
		debug_output("DONE@mpi_rank:", mpi_rank, '\n');
		return HBRS_MPL_FWD(to);
	}
#endif

template<
	typename From,
	typename std::enable_if_t<
		#ifdef HBRS_MPL_ENABLE_MATLAB
			std::is_same< hana::tag_of_t<From>, mpl::ml_matrix_tag >::value ||
		#endif
		#ifdef HBRS_MPL_ENABLE_ELEMENTAL
			std::is_same< hana::tag_of_t<From>, mpl::el_matrix_tag >::value ||
		#endif
		std::is_same< hana::tag_of_t<From>, mpl::rtsam_tag >::value ||
		std::is_same< hana::tag_of_t<From>, mpl::ctsam_tag >::value ||
		std::is_same< hana::tag_of_t<From>, mpl::sm_tag >::value
	>* = nullptr
>
decltype(auto)
copy_matrix(From const& from, std::vector<theta_field> & to) {
	auto from_sz = (*mpl::size)(from);
	auto from_m = (*mpl::m)(from_sz);
	auto from_n = (*mpl::n)(from_sz);
	
	auto to_sz = local_size(to);
	auto to_m = (*mpl::m)(to_sz);
	auto to_n = (*mpl::n)(to_sz);
	
	BOOST_ASSERT((*mpl::greater_equal)(from_m, to_m));
	BOOST_ASSERT((*mpl::equal)(from_n, to_n));
	
	for (std::size_t j = 0; j < to_n; ++j) {
		to.at(j).density().clear();
		to.at(j).pressure().clear();
		to.at(j).residual().clear();
		BOOST_ASSERT(to.at(j).x_velocity().size() * 3 == to_m);
		BOOST_ASSERT(to.at(j).y_velocity().size() * 3 == to_m);
		BOOST_ASSERT(to.at(j).z_velocity().size() * 3 == to_m);
	
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to.at(j).x_velocity().at(i) = (*mpl::at)(from, mpl::make_matrix_index(i, j));
		}
		
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to.at(j).y_velocity().at(i) = (*mpl::at)(from, mpl::make_matrix_index(i+to_m/3, j));
		}
		
		for(std::size_t i = 0; i < to_m/3; ++i) {
			to.at(j).z_velocity().at(i) = (*mpl::at)(from, mpl::make_matrix_index(i+to_m/3*2, j));
		}
	}
	
	return to;
}

#ifdef HBRS_MPL_ENABLE_ELEMENTAL
	inline decltype(auto)
	copy_matrix(
		hbrs::mpl::el_dist_matrix<double, El::VC, El::STAR, El::ELEMENT> const& from,
		std::vector<theta_field> & to
	) {
		#if !defined(NDEBUG)
		{
			mpl::matrix_size<std::size_t, std::size_t> from_sz = {from.size()};
			mpl::matrix_size<std::size_t, std::size_t> to_lcl_sz   = local_size(to);
			mpl::matrix_size<std::size_t, std::size_t> to_gbl_sz   = global_size(to);
			BOOST_ASSERT(from_sz.n() == to_lcl_sz.n());
			BOOST_ASSERT(from.data().Height() == to_gbl_sz.m());
			BOOST_ASSERT(from.data().LockedMatrix().Width() == to_gbl_sz.n());
		}
		#endif
		
		size_t mpi_sz = boost::numeric_cast<size_t>(mpi::size());
		size_t mpi_rank = boost::numeric_cast<size_t>(mpi::rank());
		mpl::matrix_size<size_t, size_t> lcl_sz = local_size(to);
		
		size_t lcl_m = lcl_sz.m();
		size_t lcl_n = lcl_sz.n();
		std::vector<size_t> lcl_ms(mpi_sz, 0u);
		mpi::allgather(&lcl_m, 1, lcl_ms.data(), 1, from.data().Grid().Comm().comm);
		debug_output("lcl_ms:", lcl_ms, "@mpi_rank:", mpi_rank, '\n');
		
		std::vector<size_t> lcl_ms_sums(mpi_sz, 0u);
		std::partial_sum(lcl_ms.begin(), lcl_ms.end(), lcl_ms_sums.begin());
		debug_output("lcl_ms_sums:", lcl_ms_sums, "@mpi_rank:", mpi_rank, '\n');
		
		size_t gbl_m = lcl_ms_sums.back();
		debug_output("gbl_m:", gbl_m, "@mpi_rank:", mpi_rank, '\n');
		
		size_t balanced_chunk_size = gbl_m / mpi_sz; // integer division
		size_t nr_of_bigger_chunks = gbl_m % mpi_sz; // number of processes that get chunks of size=chunk_size+1
		size_t nr_of_normal_chunks = mpi_sz - nr_of_bigger_chunks;
		BOOST_ASSERT((nr_of_bigger_chunks+nr_of_normal_chunks) == mpi_sz);
		
		std::vector<size_t> bal_lcl_ms(mpi_sz, 0u);
		for(size_t i = 0; i < nr_of_bigger_chunks; ++i) {
			bal_lcl_ms.at(i) = balanced_chunk_size+1;
		}
		for(size_t i = nr_of_bigger_chunks; i < mpi_sz; ++i) {
			bal_lcl_ms.at(i) = balanced_chunk_size;
		}
		std::vector<size_t> bal_lcl_ms_sums(mpi_sz, 0u);
		std::partial_sum(bal_lcl_ms.begin(), bal_lcl_ms.end(), bal_lcl_ms_sums.begin());
		
		debug_output("bal_lcl_ms:", bal_lcl_ms, "@mpi_rank:", mpi_rank, '\n');
		debug_output("bal_lcl_ms_sums:", bal_lcl_ms_sums, "@mpi_rank:", mpi_rank, '\n');
		
		size_t bal_lcl_m = bal_lcl_ms.at(mpi_rank);
		mpl::rtsam<double, mpl::storage_order::row_major> bal_lcl_from =
			{ mpl::make_matrix_size(bal_lcl_m, lcl_n) };
		
		El::DistMatrix<double, El::VC, El::STAR, El::ELEMENT> const& gbl_from = from.data();
		
		BOOST_ASSERT(gbl_from.LockedMatrix().Height() == bal_lcl_from.size().m());
		BOOST_ASSERT(gbl_from.LockedMatrix().Width() == bal_lcl_from.size().n());
		
		for(size_t gbl_row = 0; gbl_row<gbl_m; ++gbl_row) {
			size_t lcl_row = gbl_row / mpi_sz;
			bool is_lcl_row = ((gbl_row % mpi_sz) == mpi_rank);
			BOOST_ASSERT(iff(is_lcl_row, gbl_from.IsLocalRow(gbl_row)));
			
			if (is_lcl_row) {
				debug_output("lcl_row:", lcl_row, "@mpi_rank:", mpi_rank, '\n');
				
				for(size_t col = 0; col < lcl_n; ++col) {
					BOOST_ASSERT(gbl_from.IsLocal(gbl_row, col));
					
					bal_lcl_from.at(mpl::make_matrix_index(lcl_row, col)) = 
						gbl_from.GetLocal(
							gbl_from.LocalRow((El::Int)gbl_row),
							gbl_from.LocalCol((El::Int)col)
						);
				}
			}
		}
		
	
		{
			mpl::rtsam<double, mpl::storage_order::row_major> lcl_to {lcl_sz};
			
			std::vector<MPI_Request> reqs;
			for(size_t gbl_row = 0; gbl_row<gbl_m; ++gbl_row) {
				size_t recv_proc = std::distance(
					lcl_ms_sums.begin(),
					std::lower_bound(lcl_ms_sums.begin(), lcl_ms_sums.end(), gbl_row+1)
				);
				
				// block distribution
// 				size_t recv_proc = std::distance(
// 					bal_lcl_ms_sums.begin(),
// 					std::lower_bound(bal_lcl_ms_sums.begin(), bal_lcl_ms_sums.end(), gbl_row)
// 				);
				// round robin distribution
				size_t send_proc = gbl_row % mpi_sz;
				
				
				size_t lcl_recv_row = (recv_proc == 0) ? gbl_row : (gbl_row - lcl_ms_sums.at(recv_proc-1));
				// block distribution
// 				size_t lcl_send_row = (send_proc == 0) ? gbl_row : (gbl_row - bal_lcl_ms_sums.at(send_proc-1));
				// round robin distribution
				size_t lcl_send_row = gbl_row / mpi_sz;
				
				debug_output("gbl_row:", gbl_row, ",send_proc:",send_proc, ",recv_proc:", recv_proc, "@mpi_rank:", mpi_rank, '\n');
				debug_output("gbl_row:", gbl_row, ",lcl_send_row:",lcl_send_row, ",lcl_recv_row:", lcl_recv_row, "@mpi_rank:", mpi_rank, '\n');
				
				if ((send_proc == recv_proc) && (send_proc == mpi_rank)) {
					debug_output("gbl_row:", gbl_row, ", COPY", "@mpi_rank:", mpi_rank, '\n');
					
					double       & lcl_recv_ref = lcl_to.at(mpl::make_matrix_index(lcl_recv_row, 0));
					double const & lcl_send_ref = bal_lcl_from.at(mpl::make_matrix_index(lcl_send_row, 0));
					std::copy_n(&lcl_send_ref, lcl_n, &lcl_recv_ref);
					
				} else if (send_proc == mpi_rank) {
					debug_output("gbl_row:", gbl_row, ", SEND_ROW", "@mpi_rank:", mpi_rank, '\n');
					
					double const& lcl_send_ref = bal_lcl_from.at(mpl::make_matrix_index(lcl_send_row, 0));
					reqs.push_back(
						mpi::isend(
							&lcl_send_ref,
							lcl_n,
							recv_proc /*dest*/,
							send_proc /*tag*/,
							MPI_COMM_WORLD
						)
					);
				} else if (recv_proc == mpi_rank) {
					debug_output("gbl_row:", gbl_row, ", RECV_ROW", "@mpi_rank:", mpi_rank, '\n');
					
					double & lcl_recv_ref = lcl_to.at(mpl::make_matrix_index(lcl_recv_row, 0));
					
					reqs.push_back(
						mpi::irecv(
							&lcl_recv_ref,
							lcl_n,
							send_proc /*source*/,
							send_proc /*tag*/,
							MPI_COMM_WORLD
						)
					);
				} else {
					debug_output("gbl_row:", gbl_row, ", SKIP", "@mpi_rank:", mpi_rank, '\n');
				}
			}
			
			for(size_t i = 0; i < reqs.size(); ++i) {
				debug_output("reqs[", i, "], WAIT_BEGIN", "@mpi_rank:", mpi_rank, '\n');
				
				[[maybe_unused]] auto stat = mpi::wait(reqs[i]);
				//TODO: Do anything with stat?
				debug_output("reqs[", i, "], WAIT_END", "@mpi_rank:", mpi_rank, '\n');
			}
			
			to = copy_matrix(lcl_to, to);
		}
		
		debug_output("DONE@mpi_rank:", mpi_rank, '\n');
		return HBRS_MPL_FWD(to);
		
		
	}
#endif

/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END

#endif // !HBRS_THETA_UTILS_DETAIL_COPY_MATRIX_IMPL_HPP
