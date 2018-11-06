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

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fn/visualize.hpp>
#include <hbrs/theta_utils/fn/decompose.hpp>
#include <hbrs/theta_utils/dt/exception.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include <boost/program_options.hpp>
#include <boost/throw_exception.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/static_visitor.hpp>
#include <mpi.h>

HBRS_THETA_UTILS_NAMESPACE_BEGIN

struct generic_options {
	bool verbose = false;
	bool debug = false;
};

struct help_cmd {
	generic_options g_opts;
	std::string help;
};

struct version_cmd {
	generic_options g_opts;
};

struct visualize_cmd {
	generic_options g_opts;
	visualize_options v_opts;
};

struct pca_cmd {
	generic_options g_opts;
	visualize_options v_opts;
	pca_decomposition_options pca_opts;
};

static void
execute(visualize_cmd const& cmd) {
	visualize(cmd.v_opts, cmd.g_opts.verbose);
}

static void
execute(pca_cmd & cmd) {
	decompose_with_pca(cmd.v_opts, cmd.pca_opts, cmd.g_opts.verbose);
}


static void
execute(version_cmd const&) {
	std::cout 
		<< "Copyright (c) 2016-2018 Jakob Meng, <jakobmeng@web.de>" << std::endl
		<< "hbrs::theta_utils " << HBRS_THETA_UTILS_VERSION_STRING << std::endl;
}

static void
execute(help_cmd const& cmd) {
	execute(version_cmd{});
	std::cout 
		<< "Analyse THETA simulations." << std::endl
		<< cmd.help << std::endl;
}

boost::variant<
	help_cmd,
	version_cmd,
	visualize_cmd,
	pca_cmd
>
parse_options(int argc, char *argv[]) {
	namespace bpo = boost::program_options;
	
	/* NOTE:
	 * If you want to use bpo::notify(...) then do not use option required(). bpo::notify(...) will throw exceptions in 
	 * case a required option is missing. This means you will not be able to show the progam's help without also 
	 * specifying those required options at the same time, too!
	 */
	
	bpo::options_description generic;
	generic.add_options()
		("verbose,v", "verbose output");
	
	bpo::options_description hidden;
	hidden.add_options()
		("debug", "enable debug mode, e.g. no exceptions get caught in this mode")
		;
	
	bpo::options_description misc;
	misc.add_options()
		(
			"options-file",
			bpo::value<std::string>()->value_name("FILENAME"),
			"load options from FILENAME"
		)
		("help,h", "show help message")
		("version,V", "show version info");
	
	bpo::options_description subcommands;
	subcommands.add_options()
		(
			"command,C",
			bpo::value<std::string>(),
			"command to execute, one of: passthrough, pca"
		)
		(
			"subargs,A", 
			bpo::value<std::vector<std::string> >(), 
			"arguments for command"
		);
	
	bpo::options_description cmdline_options;
	cmdline_options.add(generic).add(hidden).add(subcommands).add(misc);
	
	bpo::options_description config_file_options;
	config_file_options.add(generic).add(hidden);
	
	bpo::options_description visible;
	visible.add(generic).add(subcommands).add(misc);
	
	bpo::positional_options_description positional_options;
    positional_options.add("command", 1).add("subargs", -1);

	/* Parse command line arguments */
	bpo::variables_map vm;
	
	bpo::parsed_options parsed = bpo::command_line_parser(argc, argv).
		options(cmdline_options).
		positional(positional_options).
		allow_unregistered().
		run();
	
	bpo::store(parsed, vm);
	
	if (vm.count("options-file")) {
		const std::string config_filename = vm["options-file"].as<std::string>();
		
		std::ifstream config_file(config_filename);
		if (!config_file.is_open()) {
			BOOST_THROW_EXCEPTION(bpo::reading_file{config_filename.c_str()});
		}
		bpo::store(bpo::parse_config_file(config_file, config_file_options, false /* allow_unregistered */), vm);
	}
	
	/* Parse generic options */
	generic_options g_opts;
	g_opts.verbose = (vm.count("verbose") > 0);
	g_opts.debug = (vm.count("debug") > 0);
	
	
	/* Parse commands */
	
	if (vm.count("help")) {
		std::stringstream help;
		help << visible;
		
		return help_cmd{g_opts, help.str()};
	}

	if (vm.count("version")) {
		return version_cmd{g_opts};
	}
	
	std::string cmd = vm["command"].as<std::string>();
	
	auto make_visualize_options = []() {
		bpo::options_description opts;
		opts.add_options()
			(
				"path",
				bpo::value< std::string >()->value_name("PREFIX"),
				"directory path to grid and *.pval.* files, defaults to current working directory"
			)
			(
				"input-prefix",
				bpo::value< std::string >()->value_name("PREFIX"),
				"load variables from PREFIX.pval.* files"
			)
			(
				"grid-prefix",
				bpo::value< std::string >()->value_name("PREFIX"),
				"load grid from file PREFIX.grid, defaults to input-prefix"
			)
			(
				"output-prefix",
				bpo::value< std::string >()->value_name("PREFIX"),
				"output filenames start with PREFIX, defaults to input-prefix"
			)
			(
				"output-format",
				bpo::value<std::string>()->value_name("FORMAT"),
				"output format to use, currently only VTK_LEGACY_ASCII and VTK_XML_BINARY are supported"
			)
			(
				"overwrite",
				"overwrite existing files"
			)
			(
				"simple-numbering",
				"drop timestamps from filenames and use ascending numbers (1, 2, 3...) instead, e.g. to animate time series in ParaView"
			)
			(
				"include",
				bpo::value< std::vector<std::string> >()->multitoken()->composing()->value_name("PATTERN"),
				"include (whitelist) variables from *.pval.* files matching a PATTERN, multiple listings are possible"
			)
			(
				"exclude",
				bpo::value< std::vector<std::string> >()->multitoken()->composing()->value_name("PATTERN"),
				"exclude (blacklist) variables from *.pval.* files matching a PATTERN, multiple listings are possible"
			);
		return opts;
	};
	
	auto parse_visualize_options = [](bpo::variables_map & vm) {
		visualize_options v_opts;
		
		if (vm.count("path")) {
			v_opts.path = vm["path"].as<std::string>();
		} else {
			v_opts.path = fs::current_path().string();
		}
		
		if (vm.count("input-prefix")) {
			v_opts.input_prefix = vm["input-prefix"].as< std::string >();
		} else {
			BOOST_THROW_EXCEPTION(bpo::required_option{"input-prefix"});
		}
		
		if (vm.count("grid-prefix")) {
			v_opts.grid_prefix = vm["grid-prefix"].as<std::string>();
		} else {
			v_opts.grid_prefix = v_opts.input_prefix;
		}
		
		if (vm.count("output-prefix")) {
			v_opts.output_prefix = vm["output-prefix"].as< std::string >();
		} else {
			v_opts.output_prefix = v_opts.input_prefix;
		}
		
		if (vm.count("output-format")) {
			std::string frmt = vm["output-format"].as<std::string>();
			
			BOOST_ASSERT(!frmt.empty());
			
			if (boost::iequals(frmt, "VTK_LEGACY_ASCII")) {
				v_opts.output_format = vtk_file_format::legacy_ascii;
			} else if (boost::iequals(frmt, "VTK_XML_BINARY")) {
				v_opts.output_format = vtk_file_format::xml_binary;
			} else {
				BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
					(boost::format("output format %s is unknown / not supported") % frmt).str()
				});
			}
		} else {
			v_opts.output_format = vtk_file_format::xml_binary;
		}
		
		v_opts.overwrite = (vm.count("overwrite") > 0);
		v_opts.simple_numbering = (vm.count("simple-numbering") > 0);
		
		if (vm.count("include")) {
			v_opts.includes = vm["include"].as< std::vector<std::string> >();
		}
		
		if (vm.count("exclude")) {
			v_opts.excludes = vm["exclude"].as< std::vector<std::string> >();
		}
		
		return v_opts;
	};
	
	/* Collect all the unrecognized options from the first pass. 
	 * This will include the (positional) command name, so we need to erase that.
	 * Ref.: https://stackoverflow.com/a/23098581/6490710
	 */
	std::vector<std::string> unreg_opts = bpo::collect_unrecognized(parsed.options, bpo::include_positional);
	unreg_opts.erase(unreg_opts.begin());
	
	if (cmd == "visualize") {
		bpo::options_description cmd_options("visualize options");
		cmd_options.add(make_visualize_options());
		bpo::store(bpo::command_line_parser(unreg_opts).options(cmd_options).run(), vm);
		
		visualize_cmd cmd;
		cmd.g_opts = g_opts;
		cmd.v_opts = parse_visualize_options(vm);
		
		return cmd;
	} else if (cmd == "pca") {
		bpo::options_description cmd_options("pca options");
		cmd_options.add(make_visualize_options()).add_options()
			(
				"backend",
				bpo::value<std::string>()->value_name("NAME"),
				"pca implementation to use, currently MATLAB_LAPACK, ELEMENTAL_OPENMP and ELEMENTAL_MPI are supported"
			)
			(
				"pcs",
				bpo::value< std::vector<std::string> >()->multitoken()->value_name("SELECTIONS"),
				"include only principal components within SELECTIONS. Multiple listings are possible, e.g. 1-3 or 4-last"
			);
		
		bpo::store(bpo::command_line_parser(unreg_opts).options(cmd_options).run(), vm);
		
		pca_cmd cmd;
		cmd.g_opts = g_opts;
		cmd.v_opts = parse_visualize_options(vm);
		
		if (vm.count("backend")) {
			std::string backend = vm["backend"].as<std::string>();
			
			BOOST_ASSERT(!backend.empty());
			
			if (boost::iequals(backend, "MATLAB_LAPACK")) {
				cmd.pca_opts.backend = pca_backend::matlab_lapack;
			} else if (boost::iequals(backend, "ELEMENTAL_OPENMP")) {
				cmd.pca_opts.backend = pca_backend::elemental_openmp;
			} else if (boost::iequals(backend, "ELEMENTAL_MPI")) {
				cmd.pca_opts.backend = pca_backend::elemental_mpi;
			} else {
				BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
					(boost::format("pca backend %s is unknown / not supported") % backend).str()
				});
			}
		} else {
			cmd.pca_opts.backend = pca_backend::elemental_mpi;
		}
		
		if (vm.count("pcs")) {
			cmd.pca_opts.pc_nr_seqs = vm["pcs"].as< std::vector<std::string> >();
		}
		
		return cmd;
	}
	
	BOOST_THROW_EXCEPTION(bpo::invalid_option_value{cmd});
}

HBRS_THETA_UTILS_NAMESPACE_END

int
main(int argc, char *argv[]) {
	namespace bpo = boost::program_options;
	using namespace hbrs::theta_utils;
	
	struct mpi {
		mpi(int & argc, char** & argv) {
			int ec = MPI_Init(&argc, &argv);
			if (ec != MPI_SUCCESS) {
				BOOST_THROW_EXCEPTION((
					mpi_exception{} 
					<< errinfo_mpi_error_info{mpi_error_info{ec}}
				));
			}
		}
		~mpi() {
			MPI_Finalize();
		}
	};
	mpi mpi{argc, argv};
	
	typedef decltype(parse_options(argc, argv)) cmd_t;
	
	cmd_t cmd;
	try {
		cmd = parse_options(argc, argv);
	} catch (bpo::error & e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	return boost::apply_visitor(
		[](auto && cmd) {
			if (cmd.g_opts.debug) {
				execute(cmd);
			} else {
				try {
					execute(cmd);
				} catch(boost::exception & ex) {
					std::cerr << boost::diagnostic_information(ex, cmd.g_opts.verbose) << std::endl;
					return EXIT_FAILURE;
				} catch(std::exception & ex) {
					std::cerr << boost::diagnostic_information(ex, cmd.g_opts.verbose) << std::endl;
					return EXIT_FAILURE;
				}
			}
			return EXIT_SUCCESS;
		}, cmd);
}
