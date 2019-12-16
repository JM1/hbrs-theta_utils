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

#include <hbrs/theta_utils/dt/exception.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>

#include <hbrs/mpl/fn/size.hpp>
#include <hbrs/mpl/fn/at.hpp>
#include <hbrs/mpl/fn/m.hpp>
#include <hbrs/mpl/fn/n.hpp>
#include <hbrs/mpl/fn/equal.hpp>
#include <hbrs/mpl/fn/less_equal.hpp>
#include <hbrs/mpl/fn/greater_equal.hpp>

#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/mpl/detail/log.hpp>
#include <hbrs/theta_utils/detail/iff.hpp>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/assert.hpp>
#include <algorithm>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace mpi = hbrs::mpl::detail::mpi;
namespace detail {

mpl::matrix_size<std::size_t, std::size_t>
distributed_size(
	theta_field_matrix const& series,
	theta_field_distribution_1
) {
	auto lcl_sz = series.size();
	
	if (mpi::comm_size() == 1) {
		return lcl_sz;
	}
	
	std::size_t lcl_m = lcl_sz.m();
	std::size_t lcl_n = lcl_sz.n();
	std::size_t gbl_m;
	mpi::allreduce(&lcl_m, &gbl_m, 1, MPI_SUM, MPI_COMM_WORLD);
	
	std::size_t gbl_min_n, gbl_max_n;
	mpi::allreduce(&lcl_n, &gbl_min_n, 1, MPI_MIN, MPI_COMM_WORLD);
	mpi::allreduce(&lcl_n, &gbl_max_n, 1, MPI_MAX, MPI_COMM_WORLD);
	
	if(gbl_min_n != gbl_max_n) {
		BOOST_THROW_EXCEPTION((mpl::incompatible_matrix_exception{} << mpl::errinfo_matrix_size{lcl_sz}));
	}
	
	return { gbl_m, gbl_min_n };
}

#ifdef HBRS_MPL_ENABLE_ELEMENTAL

mpl::el_dist_matrix<double, El::VC, El::STAR, El::ELEMENT>
scatter(
	theta_field_matrix const& from,
	scatter_control<theta_field_distribution_1>
) {
	using std::size_t;
	using hbrs::mpl::detail::loggable;
	
	mpl::matrix_size<size_t, size_t> lcl_sz = from.size();
	mpl::matrix_size<size_t, size_t> gbl_sz = distributed_size(from, theta_field_distribution_1{});
	
	static El::Grid const grid{El::mpi::COMM_WORLD};
	mpl::el_dist_matrix<double, El::VC, El::STAR, El::ELEMENT> to{
		grid,
		boost::numeric_cast<El::Int>(gbl_sz.m()),
		boost::numeric_cast<El::Int>(gbl_sz.n())
	};
	
	BOOST_ASSERT(lcl_sz.n() == to.size().n());
	BOOST_ASSERT(gbl_sz.n() == to.data().Matrix().Width());
	BOOST_ASSERT(gbl_sz.m() == to.data().Height());
	
	size_t mpi_sz = boost::numeric_cast<size_t>(mpi::comm_size());
	size_t mpi_rank = boost::numeric_cast<size_t>(mpi::comm_rank());
	
	
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
	HBRS_MPL_LOG_TRIVIAL(trace) << "lcl_ms:" << loggable{lcl_ms} << "@mpi_rank:" << mpi_rank;
	
	std::vector<size_t> lcl_ms_sums(mpi_sz, 0u);
	std::partial_sum(lcl_ms.begin(), lcl_ms.end(), lcl_ms_sums.begin());
	HBRS_MPL_LOG_TRIVIAL(trace) << "lcl_ms_sums:" << loggable{lcl_ms_sums} << "@mpi_rank:" << mpi_rank;
	
	size_t gbl_m = lcl_ms_sums.back();
	HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_m:" << gbl_m << "@mpi_rank:" << mpi_rank;
	
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
	
	HBRS_MPL_LOG_TRIVIAL(trace) << "bal_lcl_ms:" << loggable{bal_lcl_ms} << "@mpi_rank:" << mpi_rank;
	HBRS_MPL_LOG_TRIVIAL(trace) << "bal_lcl_ms_sums:" << loggable{bal_lcl_ms_sums} << "@mpi_rank:" << mpi_rank;
	
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
// 			size_t recv_proc = std::distance(
// 				bal_lcl_ms_sums.begin(),
// 				std::lower_bound(bal_lcl_ms_sums.begin(), bal_lcl_ms_sums.end(), gbl_row)
// 			);
			// round robin distribution
			size_t recv_proc = gbl_row % mpi_sz;
			
			
			size_t lcl_send_row = (send_proc == 0) ? gbl_row : (gbl_row - lcl_ms_sums.at(send_proc-1));
			// block distribution
// 			size_t lcl_recv_row = (recv_proc == 0) ? gbl_row : (gbl_row - bal_lcl_ms_sums.at(recv_proc-1));
			// round robin distribution
			size_t lcl_recv_row = gbl_row / mpi_sz;
			
			HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_row:" << gbl_row << ",send_proc:" << send_proc << ",recv_proc:" << recv_proc << "@mpi_rank:" << mpi_rank;
			HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_row:" << gbl_row << ",lcl_send_row:" << lcl_send_row << ",lcl_recv_row:" << lcl_recv_row << "@mpi_rank:" << mpi_rank;
			
			if ((send_proc == recv_proc) && (send_proc == mpi_rank)) {
				HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_row:" << gbl_row << ", COPY" << "@mpi_rank:" << mpi_rank;
				
				double const& lcl_send_ref = lcl_to.at(mpl::make_matrix_index(lcl_send_row, 0));
				double      & lcl_recv_ref = bal_lcl_to.at(mpl::make_matrix_index(lcl_recv_row, 0));
				std::copy_n(&lcl_send_ref, lcl_n, &lcl_recv_ref);
				
			} else if (send_proc == mpi_rank) {
				HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_row:" << gbl_row << ", SEND_ROW" << "@mpi_rank:" << mpi_rank;
				
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
				HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_row:" << gbl_row << ", RECV_ROW" << "@mpi_rank:" << mpi_rank;
				
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
				HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_row:" << gbl_row << ", SKIP" << "@mpi_rank:" << mpi_rank;
			}
		}
		
		for(size_t i = 0; i < reqs.size(); ++i) {
			HBRS_MPL_LOG_TRIVIAL(trace) << "reqs[" << i << "], WAIT_BEGIN" << "@mpi_rank:" << mpi_rank;
			
			[[maybe_unused]] auto stat = mpi::wait(reqs[i]);
			//TODO: Do anything with stat?
			HBRS_MPL_LOG_TRIVIAL(trace) << "reqs[" << i << "], WAIT_END" << "@mpi_rank:" << mpi_rank;
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
			HBRS_MPL_LOG_TRIVIAL(trace) << "lcl_row:" << lcl_row << "@mpi_rank:" << mpi_rank;
			
			for(size_t col = 0; col < lcl_n; ++col) {
				BOOST_ASSERT(gbl_to.IsLocal(gbl_row, col));
				
				gbl_to.SetLocal(
					gbl_to.LocalRow((El::Int)gbl_row),
					gbl_to.LocalCol((El::Int)col),
					bal_lcl_to.at(mpl::make_matrix_index(lcl_row, col))
				);
			}
		}
		
		HBRS_MPL_LOG_TRIVIAL(trace) << "lcl_to:" << loggable{lcl_to} << "@mpi_rank:" << mpi_rank;
		HBRS_MPL_LOG_TRIVIAL(trace) << "gbl_to:" << loggable{gbl_to} << "@mpi_rank:" << mpi_rank;
	}
	
	HBRS_MPL_LOG_TRIVIAL(trace) << "DONE@mpi_rank:" << mpi_rank;
	return HBRS_MPL_FWD(to);
}
#endif // !HBRS_MPL_ENABLE_ELEMENTAL
/* namespace detail */ }
HBRS_THETA_UTILS_NAMESPACE_END
