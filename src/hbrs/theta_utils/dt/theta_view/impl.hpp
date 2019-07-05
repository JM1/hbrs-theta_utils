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

#ifndef HBRS_THETA_UTILS_DT_THETA_VIEW_IMPL_HPP
#define HBRS_THETA_UTILS_DT_THETA_VIEW_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/mpl/dt/rtsav.hpp>

#include <boost/assert.hpp>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace mpl = hbrs::mpl;

//rho := density
//p := pressure
template<typename T>
struct theta_view {
public:
	
	theta_view(mpl::rtsav<T> rho, mpl::rtsav<T> vx, mpl::rtsav<T> vy, mpl::rtsav<T> vz, mpl::rtsav<T> p) noexcept
	: rho_{rho}, vx_{vx}, vy_{vy}, vz_{vz}, p_{p}
	{}
	
	theta_view(theta_view const&) noexcept = default;
	theta_view(theta_view &&) noexcept = default;
	
	theta_view&
	operator=(theta_view const&) noexcept = delete;
	theta_view&
	operator=(theta_view &&) noexcept = delete;
	
	decltype(auto)
	rho() noexcept {
		return (rho_);
	}
	
	decltype(auto)
	rho() const noexcept {
		return (rho_);
	}
	
	decltype(auto)
	vx() noexcept {
		return (vx_);
	}
	
	decltype(auto)
	vx() const noexcept {
		return (vx_);
	}
	
	decltype(auto)
	vy() noexcept {
		return (vy_);
	}
	
	decltype(auto)
	vy() const noexcept {
		return (vy_);
	}
	
	decltype(auto)
	vz() noexcept {
		return (vz_);
	}
	
	decltype(auto)
	vz() const noexcept {
		return (vz_);
	}
	
	decltype(auto)
	p() noexcept {
		return (p_);
	}
	
	decltype(auto)
	p() const noexcept {
		return (p_);
	}
	
private:
	mpl::rtsav<T> rho_;
	mpl::rtsav<T> vx_;
	mpl::rtsav<T> vy_;
	mpl::rtsav<T> vz_;
	mpl::rtsav<T> p_;
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<typename T>
struct tag_of< hbrs::theta_utils::theta_view<T> > {
	using type = hbrs::theta_utils::theta_view_tag;
};

template <>
struct make_impl<hbrs::theta_utils::theta_view_tag> {
	
	template<typename T>
	static constexpr decltype(auto)
	reinterpret(T raw) {
		std::size_t pos{0};
		typedef std::conditional_t<
			std::is_const<std::remove_reference_t<T>>::value,
			double const,
			double
		> theta_data;
		//typedef std::remove_pointer_t<decltype(raw.data())> theta_data;
		typedef int const theta_size;
		
		//TODO: Replace ...data()+pos... with fold()
		auto rho_size = (unsigned) *reinterpret_cast<theta_size*>( raw.data()+pos );
		pos += sizeof(theta_size);
		auto vx_size  = (unsigned) *reinterpret_cast<theta_size*>( raw.data()+pos );
		pos += sizeof(theta_size);
		auto vy_size  = (unsigned) *reinterpret_cast<theta_size*>( raw.data()+pos );
		pos += sizeof(theta_size);
		auto vz_size  = (unsigned) *reinterpret_cast<theta_size*>( raw.data()+pos );
		pos += sizeof(theta_size);
		auto p_size   = (unsigned) *reinterpret_cast<theta_size*>( raw.data()+pos );
		pos += sizeof(theta_size);
		
		auto rho  = reinterpret_cast<theta_data*>( raw.data()+pos );
		pos += rho_size * sizeof(theta_data);
		auto vx   = reinterpret_cast<theta_data*>( raw.data()+pos );
		pos += vx_size  * sizeof(theta_data);
		auto vy   = reinterpret_cast<theta_data*>( raw.data()+pos );
		pos += vy_size  * sizeof(theta_data);
		auto vz   = reinterpret_cast<theta_data*>( raw.data()+pos );
		pos += vz_size  * sizeof(theta_data);
		auto p    = reinterpret_cast<theta_data*>( raw.data()+pos );
		pos += p_size   * sizeof(theta_data);

		BOOST_ASSERT(raw.length() == pos);
		
		return hbrs::theta_utils::theta_view<theta_data>{ 
			hbrs::mpl::rtsav<theta_data>{rho, rho_size}, 
			hbrs::mpl::rtsav<theta_data>{vx, vx_size}, 
			hbrs::mpl::rtsav<theta_data>{vy, vy_size}, 
			hbrs::mpl::rtsav<theta_data>{vz, vz_size}, 
			hbrs::mpl::rtsav<theta_data>{p, p_size}
		};
	}
	
	static decltype(auto)
	apply(hbrs::mpl::rtsav<char> const& raw) {
		return reinterpret(raw);
	}
	
	static decltype(auto)
	apply(hbrs::mpl::rtsav<char> & raw) {
		return reinterpret(raw);
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_THETA_VIEW_IMPL_HPP
