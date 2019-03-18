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
#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/mpl/dt/exception.hpp>
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
namespace mpi = hbrs::mpl::detail::mpi;

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

theta_field_path::theta_field_path(
	fs::path folder,
	std::string prefix,
	struct timestamp timestamp,
	int step,
	boost::optional<int> domain_num)
: folder_{folder}, prefix_{prefix}, timestamp_{timestamp}, step_{step}, domain_num_{domain_num} {}

fs::path
theta_field_path::filename() const {
	// filename example: karman.pval.t5_000e-02.100.domain_1
	std::stringstream fn;
	fn  << prefix_
		<< ".pval."
		<< (boost::format("t%d_%03de%+03d") % timestamp_.significand()[0] % timestamp_.significand()[1] % timestamp_.exponent())
		<< "." << step_;
	if (domain_num_) {
		fn << ".domain_" << *domain_num_;
	}
	
	return { fn.str() };
}

fs::path
theta_field_path::full_path() const {
	return folder_ / filename();
}

HBRS_THETA_UTILS_DEFINE_ATTR(folder, fs::path, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(prefix, std::string, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(timestamp, struct theta_field_path::timestamp, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(step, int, theta_field_path)
HBRS_THETA_UTILS_DEFINE_ATTR(domain_num, boost::optional<int>, theta_field_path)

theta_field::theta_field(
	std::vector<double> density,
	std::vector<double> x_velocity,
	std::vector<double> y_velocity,
	std::vector<double> z_velocity,
	std::vector<double> pressure,
	std::vector<double> residual,
	std::vector<int> global_id,
	boost::optional<int> ndomains
) : density_{density},
	x_velocity_{x_velocity},
	y_velocity_{y_velocity},
	z_velocity_{z_velocity},
	pressure_{pressure},
	residual_{residual},
	global_id_{global_id},
	ndomains_{ndomains}
	{}

namespace {
bool
operator==(std::vector<nc_dimension> const& dims, std::vector<std::string> const& dim_names) {
	if (dims.size() != dim_names.size()) {
		return false;
	}
	
	for(std::size_t i = 0; i < dims.size(); ++i) {
		if (dims[i].name() != dim_names[i]) {
			return false;
		}
	}
	
	return true;
}

bool
operator!=(std::vector<nc_dimension> const& dims, std::vector<std::string> const& dim_names) {
	return !(dims == dim_names);
}

/* unnamed namespace */ }

theta_field::theta_field(nc_cntr cntr) {
	#define __get(__name, __type)                                                                                      \
		{                                                                                                              \
			auto opt = cntr.variable(#__name);                                                                         \
			if (opt) {                                                                                                 \
				if (opt->dimensions() != std::vector<std::string>{"no_of_points"}) {                                   \
					BOOST_THROW_EXCEPTION(unsupported_format_exception{});                                             \
				}                                                                                                      \
				__name ## _ = boost::get< std::vector<__type> >( opt->data() );                                        \
			}                                                                                                          \
		}
	
	__get(density, double)
	__get(x_velocity, double)
	__get(y_velocity, double)
	__get(z_velocity, double)
	__get(pressure, double)
	__get(residual, double)
	__get(global_id, int)
	
	#undef __get
	
	{
		auto opt = cntr.attribute("ndomains");
		if (opt) {
			std::vector<int> ndomains = boost::get<std::vector<int>>(opt->value());
			BOOST_ASSERT(ndomains.size() == 1);
			ndomains_ = ndomains[0];
		}
	}
}

HBRS_THETA_UTILS_DEFINE_ATTR(density, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(x_velocity, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(y_velocity, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(z_velocity, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(pressure, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(residual, std::vector<double>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(global_id, std::vector<int>, theta_field)
HBRS_THETA_UTILS_DEFINE_ATTR(ndomains, boost::optional<int>, theta_field)

theta_field
read_theta_field(
	std::string const& file_path,
	std::vector<std::string> const& includes,
	std::vector<std::string> const& excludes
) {
	if (includes.empty() && excludes.empty()) {
		// only include currently supported:
		return { read_nc_cntr(file_path, { "density", ".*_velocity", "pressure", "residual", /*".*_old",*/ "global_id" }) };
	} else {
		return { read_nc_cntr(file_path, includes, excludes) };
	}
}

namespace {

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
				qi::lit(".domain_") >> qi::int_[phoenix::push_back(phoenix::ref(domain_num), qi::_1)]
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

/* unnamed namespace */ }

std::vector<theta_field_path>
find_theta_fields(
	fs::path const& dir,
	std::string const& prefix
) {
	std::vector<fs::path> all_files;
	std::copy(fs::directory_iterator(dir), fs::directory_iterator(), std::back_inserter(all_files));
	
	std::vector<theta_field_path> field_files;
	for(auto && path : all_files) {
		auto field_path = parse_theta_field_path(path, prefix);
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

std::vector<theta_field_path>
filter_theta_fields_by_domain_num(
	std::vector<theta_field_path> const& fields,
	boost::optional<int> const& domain_num
) {
	BOOST_ASSERT(domain_num ? *domain_num == mpi::rank() : true);
	
	std::vector<theta_field_path> kept;
	if (!domain_num) {
		BOOST_ASSERT(mpi::size() == 1);
		// keep files without domain suffix
		std::copy_if(
			fields.begin(),
			fields.end(),
			std::back_inserter(kept),
			[](auto path){return !path.domain_num();}
		);
		
		if (!kept.empty()) {
			boost::optional<int> first_num = kept.front().domain_num();
			for (auto && path : kept) {
				if (path.domain_num() != first_num) {
					BOOST_THROW_EXCEPTION((
						ambiguous_domain_num_exception{}
						<< errinfo_ambiguous_domain_num{{kept.front().full_path(), path.full_path()}}
					));
				}
			}
		}
	} else {
		// keep files with a domain suffix
		std::copy_if(
			fields.begin(),
			fields.end(),
			std::back_inserter(kept),
			[&domain_num](auto path){return path.domain_num() == domain_num;}
		);
	}
	return kept;
}

std::vector<theta_field>
read_theta_fields(
	std::vector<theta_field_path> const& paths,
	std::vector<std::string> const& includes,
	std::vector<std::string> const& excludes
) {
	std::vector<theta_field> fields;
	fields.reserve(paths.size());
	
	for(auto && path : paths) {
		fields.push_back(
			read_theta_field(path.full_path().string(), includes, excludes)
		);
	}
	
	return fields;
}

namespace {

nc_cntr
gen_nc_cntr(theta_field field) {
	std::vector<nc_dimension> dims;
	std::vector<nc_variable> vars;
	std::vector<nc_attribute> atts;
	
	#define __add(__name)                                                                                              \
		{                                                                                                              \
			auto vec = std::move(field.__name());                                                                      \
			if (!vec.empty()) {                                                                                        \
				auto dim_it = std::find_if(                                                                            \
					dims.begin(),                                                                                      \
					dims.end(),                                                                                        \
					[](auto const& dim) { return dim.name() == "no_of_points"; }                                       \
				);                                                                                                     \
				if (dim_it != dims.end()) {                                                                            \
					if (vec.size() != dim_it->length()) {                                                              \
						BOOST_THROW_EXCEPTION(unsupported_format_exception{});                                         \
						/* TODO: Add more exception info, e.g. expected and received sizes */                          \
					}                                                                                                  \
				} else {                                                                                               \
					dims.push_back({"no_of_points", vec.size()});                                                      \
					dim_it = --dims.end();                                                                             \
				}                                                                                                      \
				vars.push_back({#__name, {*dim_it}, {std::move(vec)}});                                                \
			}                                                                                                          \
		}
	
	__add(density)
	__add(x_velocity)
	__add(y_velocity)
	__add(z_velocity)
	__add(pressure)
	__add(residual)
	__add(global_id)
	
	#undef __add
	
	{
		if (field.ndomains()) {
			atts.push_back({
				"ndomains",
				std::vector<int>{{ *field.ndomains() }}
			});
		}
	}
	
	return {dims, vars, atts};
}

/* unnamed namespace */ }

void
write_theta_field(
	theta_field field,
	std::string const& file_path,
	bool overwrite
) {
	write_nc_cntr(gen_nc_cntr(std::move(field)), file_path, overwrite);
}

void
write_theta_fields(
	std::vector< std::tuple<theta_field, theta_field_path> > fields,
	bool overwrite
) {
	for(auto & pack : fields) {
		auto & [field, path] = pack;
		auto file_path = (path.folder() / path.filename()).string();
		write_theta_field(std::move(field), file_path, overwrite);
	}
}

HBRS_THETA_UTILS_NAMESPACE_END
