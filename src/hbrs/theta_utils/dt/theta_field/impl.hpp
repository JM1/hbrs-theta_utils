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

#ifndef HBRS_THETA_UTILS_DT_THETA_FIELD_IMPL_HPP
#define HBRS_THETA_UTILS_DT_THETA_FIELD_IMPL_HPP

#include "fwd.hpp"

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/dt/nc_cntr.hpp>
#include <hbrs/theta_utils/core/preprocessor.hpp>
#include <boost/hana/core.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include <array>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace hana = boost::hana;
namespace fs = boost::filesystem;

struct HBRS_THETA_UTILS_API theta_field_path {
	struct timestamp {
		typedef std::array<std::string, 2> significand_t;
		timestamp(significand_t significand, std::string exponent);
		timestamp(timestamp const&) = default;
		timestamp(timestamp &&) = default;
		
		timestamp&
		operator=(timestamp const&) = default;
		timestamp&
		operator=(timestamp &&) = default;
		
		bool
		operator<(timestamp const& o) const;
		
		std::string
		string() const;
		
		HBRS_THETA_UTILS_DECLARE_ATTR(significand, significand_t)
		HBRS_THETA_UTILS_DECLARE_ATTR(exponent, std::string)
	};
	
	
	enum class naming_scheme { theta, tau_unsteady };
	
	theta_field_path(
		fs::path folder,
		std::string prefix,
		struct timestamp timestamp,
		int step,
		boost::optional<int> domain_num,
		enum naming_scheme naming_scheme
	);
	
	theta_field_path(theta_field_path const&) = default;
	theta_field_path(theta_field_path &&) = default;
	
	theta_field_path&
	operator=(theta_field_path const&) = default;
	theta_field_path&
	operator=(theta_field_path &&) = default;
	
	fs::path
	filename() const;
	
	fs::path
	full_path() const;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(folder, fs::path)
	HBRS_THETA_UTILS_DECLARE_ATTR(prefix, std::string)
	HBRS_THETA_UTILS_DECLARE_ATTR(timestamp, struct timestamp)
	HBRS_THETA_UTILS_DECLARE_ATTR(step, int)
	HBRS_THETA_UTILS_DECLARE_ATTR(domain_num, boost::optional<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(naming_scheme, enum naming_scheme)
};

struct HBRS_THETA_UTILS_API theta_field {
public:
	theta_field(
		std::vector<double> density,
		std::vector<double> x_velocity,
		std::vector<double> y_velocity,
		std::vector<double> z_velocity,
		std::vector<double> pressure,
		std::vector<double> residual,
		std::vector<int> global_id,
		boost::optional<int> ndomains
	);
	
	theta_field(nc_cntr cntr);
	
	theta_field(theta_field const&) = default;
	theta_field(theta_field &&) = default;
	
	theta_field&
	operator=(theta_field const&) = default;
	theta_field&
	operator=(theta_field &&) = default;
	
	HBRS_THETA_UTILS_DECLARE_ATTR(density, std::vector<double>)
	HBRS_THETA_UTILS_DECLARE_ATTR(x_velocity, std::vector<double>)
	HBRS_THETA_UTILS_DECLARE_ATTR(y_velocity, std::vector<double>)
	HBRS_THETA_UTILS_DECLARE_ATTR(z_velocity, std::vector<double>)
	HBRS_THETA_UTILS_DECLARE_ATTR(pressure, std::vector<double>)
	HBRS_THETA_UTILS_DECLARE_ATTR(residual, std::vector<double>)
	HBRS_THETA_UTILS_DECLARE_ATTR(global_id, std::vector<int>)
	HBRS_THETA_UTILS_DECLARE_ATTR(ndomains, boost::optional<int>)
};

HBRS_THETA_UTILS_NAMESPACE_END

namespace boost { namespace hana {

template<>
struct tag_of< hbrs::theta_utils::theta_field > {
	using type = hbrs::theta_utils::theta_field_tag;
};

template <>
struct make_impl<hbrs::theta_utils::theta_field_tag> {
	static hbrs::theta_utils::theta_field
	apply(hbrs::theta_utils::nc_cntr cntr) {
		return {cntr};
	}
};

/* namespace hana */ } /* namespace boost */ }

#endif // !HBRS_THETA_UTILS_DT_THETA_FIELD_IMPL_HPP
