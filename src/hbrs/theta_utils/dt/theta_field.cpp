/* Copyright (c) 2016 Jakob Meng, <jakobmeng@web.de>
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

#include <hbrs/theta_utils/dt/theta_field.hpp>
#include <hbrs/theta_utils/dt/nc_exception.hpp>
#include <hbrs/theta_utils/dt/exception.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/system/error_code.hpp>
#include <mpi.h>
#include <iterator>
#include <sstream>

HBRS_THETA_UTILS_NAMESPACE_BEGIN
namespace ba = boost::algorithm;

theta_field::theta_field(
	std::vector<double> density,
	std::vector<double> x_velocity,
	std::vector<double> x_velocity_old,
	std::vector<double> y_velocity,
	std::vector<double> y_velocity_old,
	std::vector<double> z_velocity,
	std::vector<double> z_velocity_old,
	std::vector<double> pressure,
	std::vector<double> pressure_old,
	std::vector<double> residual
) : density_{density},
	x_velocity_{x_velocity},
	x_velocity_old_{x_velocity_old},
	y_velocity_{y_velocity},
	y_velocity_old_{y_velocity_old},
	z_velocity_{z_velocity},
	z_velocity_old_{z_velocity_old},
	pressure_{pressure},
	pressure_old_{pressure_old},
	residual_{residual}
	{}

theta_field::theta_field(nc_cntr cntr) {
#define __get(__name)                                                                                                  \
	{                                                                                                                  \
		auto opt = cntr.variable(#__name);                                                                             \
		if (opt) {                                                                                                     \
			__name ## _ = boost::get< std::vector<double> >( opt->data() );                                            \
		}                                                                                                              \
	}
	
	__get(density)
	__get(x_velocity)
	__get(x_velocity_old)
	__get(y_velocity)
	__get(y_velocity_old)
	__get(z_velocity)
	__get(z_velocity_old)
	__get(pressure)
	__get(pressure_old)
	__get(residual)
}

HBRS_THETA_UTILS_DEFINE_ATTR(density, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(x_velocity, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(x_velocity_old, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(y_velocity, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(y_velocity_old, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(z_velocity, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(z_velocity_old, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(pressure, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(pressure_old, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(residual, std::vector<double>, theta_field)

theta_field_path::timestamp::timestamp(std::array<int, 2> significand, int exponent)
: significand_{significand}, exponent_{exponent} {}

bool
theta_field_path::timestamp::operator<(timestamp const& o) const {
	return boost::lexical_cast<double>(string()) < boost::lexical_cast<double>(o.string());
}

std::string
theta_field_path::timestamp::string() const {
	std::stringstream ss;
	ss << significand_[0] << '.' << significand_[1] << 'e' << exponent_;
	return ss.str();
}

HBRS_THETA_UTILS_DEFINE_ATTR(significand, theta_field_path::timestamp::significand_t, theta_field_path::timestamp)
HBRS_THETA_UTILS_DEFINE_ATTR(exponent, int, theta_field_path::timestamp)

theta_field_path::theta_field_path(fs::path folder, std::string prefix, struct timestamp timestamp, int step, boost::optional<int> domain_num)
: folder_{folder}, prefix_{prefix}, timestamp_{timestamp}, step_{step}, domain_num_{domain_num} {}

fs::path
theta_field_path::filename() {
	// filename example: karman.pval.t5_000e-02.100.domain_1
	std::stringstream fn;
	fn  << prefix_
		<< ".pval."
		<< (boost::format("t%d_%03de%+03d") % timestamp_.significand()[0] % timestamp_.significand()[1] % timestamp_.exponent())
		<< "." << step_;
	if (domain_num_) {
		fn << "domain_" << *domain_num_;
	}
	
	return { fn.str() };
}

fs::path
theta_field_path::full_path() {
	return folder_ / filename();
}

HBRS_THETA_UTILS_DEFINE_ATTR(folder, fs::path, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(prefix, std::string, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(timestamp, struct theta_field_path::timestamp, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(step, int, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(domain_num, boost::optional<int>, theta_field_path)

theta_field
read_theta_field(
	std::string const& file_path,
	std::vector<std::string> const& includes,
	std::vector<std::string> const& excludes
) {
	if (includes.empty() && excludes.empty()) {	
		return { read_nc_cntr(file_path, { "density", ".*_velocity", "pressure", "residual" ".*_old" }) };
	} else {
		return { read_nc_cntr(file_path, includes, excludes) };
	}
}

static
boost::optional<theta_field_path>
parse_theta_field_path(fs::path file_path, std::string const& input_prefix) {
	namespace qi = boost::spirit::qi;
	namespace spirit = boost::spirit;
	namespace phoenix = boost::phoenix;
	
	auto filename = file_path.filename().string();
	// filename example: karman.pval.t5_000e-02.100.domain_1
	
	
	int timestamp[3];
	int step;
	std::vector<int> domain_num;
	
	auto begin = filename.begin();
	auto end = filename.end();
	bool ok = qi::parse(begin, end,
		qi::lexeme[
			qi::lit(input_prefix) 
			>> qi::lit(".pval.") 
			>> qi::char_('t') 
			>> (
				qi::int_[phoenix::ref(timestamp[0]) = qi::_1]
				>> qi::char_('_') 
				>> qi::int_[phoenix::ref(timestamp[1]) = qi::_1] 
				>> qi::char_('e') 
				>> qi::int_[phoenix::ref(timestamp[2]) = qi::_1]
			)
			>> qi::char_('.') 
			>> qi::int_[phoenix::ref(step) = qi::_1]
			>> qi::repeat(0, 1)[
				qi::lit("domain_") >> qi::int_[phoenix::push_back(phoenix::ref(domain_num), qi::_1)]
			]
		]
	);
	if ((begin != end) || !ok) { 
		return {}; 
	};
	
	return {
		{ 
			file_path.parent_path(),
			input_prefix, 
			{ {timestamp[0], timestamp[1]}, timestamp[2] },
			step,
			domain_num.empty() ? boost::optional<int>{} : boost::optional<int>{ domain_num[0] }
		}
	};
}

static
std::vector<theta_field_path> 
find_theta_field_files(
	std::string const& folder_path,
	std::string const& input_prefix
) {
	fs::path dir{folder_path};
	std::vector<fs::path> all_files;
	std::copy(fs::directory_iterator(dir), fs::directory_iterator(), std::back_inserter(all_files));
	
	std::vector<theta_field_path> field_files;
	for(auto && path : all_files) {
		auto field_path = parse_theta_field_path(path, input_prefix);
		if (!field_path) {
			continue;
		}
		field_files.push_back(*field_path);
	}
	
	std::sort(
		field_files.begin(), 
		field_files.end(), 
		[](theta_field_path const& l, theta_field_path const& r) {
			BOOST_ASSERT((
				(l.timestamp() < r.timestamp()) == (l.step() < r.step())
			));
			return l.timestamp() < r.timestamp();
		}
	);
	
	return field_files;
}


static
void
sanity_check(std::vector<theta_field_path> paths, std::string const& folder_path, std::string const& prefix) {
	if (paths.empty()) { 
		BOOST_THROW_EXCEPTION((
			fs::filesystem_error{
				(boost::format("no theta field files (*.pval.*) with prefix %s found in folder %s") % prefix % folder_path).str(),
				boost::system::errc::make_error_code(boost::system::errc::no_such_file_or_directory)
			}
		));
	}
	
	boost::optional<int> domain_num = paths.front().domain_num();
	for (auto && path : paths) {
		if (path.domain_num() != domain_num) {
			BOOST_THROW_EXCEPTION((
				ambiguous_domain_num_exception{}
				<< errinfo_ambiguous_domain_num{{paths.front().full_path(), path.full_path()}}
			));
		}
	}
	
	int init = false;
	
	int ec = MPI_Initialized(&init);
	if (ec != MPI_SUCCESS) { 
		BOOST_THROW_EXCEPTION(mpi_exception{} << errinfo_mpi_error_info{mpi_error_info{ec}});
	}
	
	if (!init && !domain_num) {
		return;
	}
	
	if (!init && domain_num) {
		BOOST_THROW_EXCEPTION(mpi_not_initialized_exception{});
	}
	
	int size = 0;
	ec = MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (ec != MPI_SUCCESS) { 
		BOOST_THROW_EXCEPTION(mpi_exception{} << errinfo_mpi_error_info{mpi_error_info{ec}});
	}
	
	int rank = -1;
	ec = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (ec != MPI_SUCCESS) { 
		BOOST_THROW_EXCEPTION(mpi_exception{} << errinfo_mpi_error_info{mpi_error_info{ec}});
	}
	
	if (size == 1 && !domain_num) {
		return;
	}
	
	if (!domain_num) {
		BOOST_THROW_EXCEPTION((
			domain_num_mismatch_exception{} 
			<< errinfo_domain_num_mismatch{domain_num_mismatch_error_info{paths.front().full_path(), rank, domain_num}}
		));
	}
	
	if (*domain_num == rank) {
		return;
	}
	
	BOOST_THROW_EXCEPTION((
		domain_num_mismatch_exception{} 
		<< errinfo_domain_num_mismatch{domain_num_mismatch_error_info{paths.front().full_path(), rank, domain_num}}
	));
}

std::vector<
	std::tuple<theta_field, theta_field_path>
>
read_theta_fields(
	std::string const& folder_path,
	std::string const& input_prefix,
	std::vector<std::string> const& includes,
	std::vector<std::string> const& excludes
) {
	std::vector<theta_field_path> field_files = find_theta_field_files(folder_path, input_prefix);
	sanity_check(field_files, folder_path, input_prefix);
	
	std::vector<
		std::tuple<theta_field, theta_field_path>
	> zipped;
	zipped.reserve(field_files.size());
	
	for(auto && file : field_files) {
		zipped.push_back({
			read_theta_field(file.full_path().string(), includes, excludes),
			file
		});
	}
	
	return zipped;
}

HBRS_THETA_UTILS_NAMESPACE_END
