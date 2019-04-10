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

#include <hbrs/theta_utils/config.hpp>
#include <hbrs/theta_utils/fn/execute.hpp>
#include <hbrs/theta_utils/dt/exception.hpp>
#include <hbrs/mpl/detail/environment.hpp>
#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/mpl/preprocessor/core.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <hbrs/mpl/config.hpp>

#include <boost/program_options.hpp>
#include <boost/throw_exception.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/static_visitor.hpp>
#include <sstream>

using namespace hbrs::theta_utils;
namespace {
namespace mpi = hbrs::mpl::detail::mpi;

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
	
	bpo::options_description generic{"generic"};
	generic.add_options()
		("verbose,v", "verbose output");
	
	bpo::options_description hidden;
	hidden.add_options()
		("debug", "enable debug mode, e.g. no exceptions get caught in this mode")
		;
	
	bpo::options_description misc{"misc"};
	misc.add_options()
		(
			"options-file",
			bpo::value<std::string>()->value_name("FILENAME"),
			"load options from FILENAME"
		)
		("help,h", "show help message")
		("version,V", "show version info");
	
	bpo::options_description commands{"commands"};
	commands.add_options()
		(
			"command",
			bpo::value<std::string>(),
			"command to execute, one of: visualize, pca"
		)
		(
			"command-options",
			bpo::value<std::vector<std::string> >(), 
			"command options"
		);
	
	bpo::options_description cmdline_options;
	cmdline_options.add(generic).add(hidden).add(misc).add(commands);
	
	bpo::options_description config_file_options;
	config_file_options.add(generic).add(hidden);
	
	bpo::options_description visible;
	visible.add(generic).add(misc).add(commands);
	
	bpo::positional_options_description positional_options;
    positional_options.add("command", 1).add("command-options", -1);

	/* Parse command line arguments */
	bpo::variables_map vm;
	
	bpo::parsed_options parsed = bpo::command_line_parser(argc, argv).
		options(cmdline_options).
		positional(positional_options).
		allow_unregistered().
		run();
	
	bpo::store(parsed, vm);
	
	if (vm.count("options-file")) {
		std::string config_filename = vm["options-file"].as<std::string>();
		
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
	
	fs::path exe{argv[0]};
	
	if (vm.count("version")) {
		std::stringstream version;
		version 
			<< exe.filename().string() 
			<< HBRS_THETA_UTILS_VERSION_STRING 
			<< std::endl;
		return version_cmd{g_opts, version.str()};
	}
	
	if (!vm.count("command")) {
		std::stringstream help;
		help
			<< "Usage: " << exe.filename().string() << " [generic/misc-options] command [command-options]" << std::endl
			<< visible;
		return help_cmd{g_opts, help.str()};
	}
	
	std::string cmd = vm["command"].as<std::string>();
	
	auto make_theta_input_options = []() {
		bpo::options_description opts;
		opts.add_options()
			(
				"path",
				bpo::value< std::string >()->value_name("PATH"),
				"directory path to grid and *.pval.* files, defaults to current working directory"
			)
			(
				"pval-prefix",
				bpo::value< std::string >()->value_name("PREFIX"),
				"load variables from PREFIX.pval.* files"
			)
			(
				"grid-prefix",
				bpo::value< std::string >()->value_name("PREFIX"),
				"load grid from file PREFIX.grid, defaults to value of --pval-prefix"
			);
		return opts;
	};
	
	auto parse_theta_input_options = [](bpo::variables_map & vm) {
		theta_input_options opts;
		
		if (vm.count("path")) {
			opts.path = vm["path"].as<std::string>();
		} else {
			opts.path = fs::current_path().string();
		}
		
		if (vm.count("pval-prefix")) {
			opts.pval_prefix = vm["pval-prefix"].as< std::string >();
		} else {
			BOOST_THROW_EXCEPTION(bpo::required_option{"pval-prefix"});
		}
		
		if (vm.count("grid-prefix")) {
			opts.grid_prefix = vm["grid-prefix"].as<std::string>();
		} else {
			opts.grid_prefix = opts.pval_prefix;
		}
		
		return opts;
	};
	
	auto make_theta_output_options = []() {
		bpo::options_description opts;
		opts.add_options()
			(
				"output-path",
				bpo::value< std::string >()->value_name("PATH"),
				"directory path to write output files to, defaults to value of --path"
			)
			(
				"output-prefix",
				bpo::value< std::string >()->value_name("PREFIX"),
				"output filenames start with PREFIX, defaults to value of --pval-prefix"
			)
			(
				"overwrite",
				"overwrite existing files"
			);
		return opts;
	};
	
	auto parse_theta_output_options = [](theta_input_options & i_opts, bpo::variables_map & vm) {
		theta_output_options opts;
		
		if (vm.count("output-path")) {
			opts.path = vm["output-path"].as<std::string>();
		} else {
			opts.path = i_opts.path;
		}
		
		if (vm.count("output-prefix")) {
			opts.prefix = vm["output-prefix"].as< std::string >();
		} else {
			opts.prefix = i_opts.pval_prefix;
		}
		
		opts.overwrite = (vm.count("overwrite") > 0);
		
		return opts;
	};
	
	/* Collect all the unrecognized options from the first pass. 
	 * This will include the (positional) command name, so we need to erase that.
	 * Ref.: https://stackoverflow.com/a/23098581/6490710
	 */
	std::vector<std::string> unreg_opts = bpo::collect_unrecognized(parsed.options, bpo::include_positional);
	if (!unreg_opts.empty()) {
		unreg_opts.erase(unreg_opts.begin());
	}
	
	/* NOTE: Global options (from e.g. cmdline_options, positional_options and config_file_options) can be parsed by 
	 *       commands, but it not possible to redefine those options.
	 */
	
	if (cmd == "visualize") {
		bpo::options_description cmd_options("visualize options");
		cmd_options.add(make_theta_input_options()).add(make_theta_output_options()).add_options()
			(
				"include",
				bpo::value< std::vector<std::string> >()->multitoken()->composing()->value_name("PATTERN"),
				"include (whitelist) variables from *.pval.* files matching a PATTERN, multiple listings are possible"
			)
			(
				"exclude",
				bpo::value< std::vector<std::string> >()->multitoken()->composing()->value_name("PATTERN"),
				"exclude (blacklist) variables from *.pval.* files matching a PATTERN, multiple listings are possible"
			)
			(
				"output-format",
				bpo::value<std::string>()->value_name("FORMAT"),
				"output format to use, currently only VTK_LEGACY_ASCII and VTK_XML_BINARY are supported"
			)
			(
				"simple-numbering",
				"drop timestamps from filenames and use ascending numbers (1, 2, 3...) instead, e.g. to animate time series in ParaView"
			)
		;
		
		bpo::parsed_options unreg_parsed = bpo::command_line_parser(unreg_opts).options(cmd_options).run();
		bpo::store(unreg_parsed, vm);
		
		unreg_opts = bpo::collect_unrecognized(unreg_parsed.options, bpo::include_positional);
		if (!unreg_opts.empty()) {
			BOOST_THROW_EXCEPTION(bpo::unknown_option{unreg_opts.front()});
		}
		
		if (vm.count("help")) {
			bpo::options_description visible;
			visible.add(generic).add(misc).add(cmd_options);
			
			std::stringstream help;
			help
				<< "Usage: " << exe.filename().string() << " [generic/misc-options] visualize [visualize-options]" << std::endl
				<< visible;
			return help_cmd{g_opts, help.str()};
		}
		
		visualize_cmd cmd;
		cmd.g_opts = g_opts;
		cmd.i_opts = parse_theta_input_options(vm);
		cmd.o_opts = parse_theta_output_options(cmd.i_opts, vm);


		if (vm.count("include")) {
			cmd.v_opts.includes = vm["include"].as< std::vector<std::string> >();
		}
		
		if (vm.count("exclude")) {
			cmd.v_opts.excludes = vm["exclude"].as< std::vector<std::string> >();
		}		

		if (vm.count("output-format")) {
			std::string frmt = vm["output-format"].as<std::string>();
			
			BOOST_ASSERT(!frmt.empty());
			
			if (boost::iequals(frmt, "VTK_LEGACY_ASCII")) {
				cmd.v_opts.format = vtk_file_format::legacy_ascii;
			} else if (boost::iequals(frmt, "VTK_XML_BINARY")) {
				cmd.v_opts.format = vtk_file_format::xml_binary;
			} else {
				BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
					(boost::format("output format %s is unknown / not supported") % frmt).str()
				});
			}
		} else {
			cmd.v_opts.format = vtk_file_format::xml_binary;
		}
		
		cmd.v_opts.simple_numbering = (vm.count("simple-numbering") > 0);
		
		return cmd;
	} else if (cmd == "pca") {
		bpo::options_description cmd_options("pca options");
		cmd_options.add(make_theta_input_options()).add(make_theta_output_options()).add_options()
			(
				"backend",
				bpo::value<std::string>()->value_name("NAME"),
				"pca implementation to use, currently MATLAB_LAPACK, ELEMENTAL_OPENMP and ELEMENTAL_MPI are supported"
			)
			(
				"pcs",
				bpo::value< std::vector<std::string> >()->multitoken()->value_name("SELECTIONS"),
				"include only principal components within SELECTIONS, e.g. \"0\", \"0,1,2\", \"0-2,6-8\", \"first\" (equal to \"0\") or \"last\". Multiple listings are possible."
			);
		
		bpo::parsed_options unreg_parsed = bpo::command_line_parser(unreg_opts).options(cmd_options).run();
		bpo::store(unreg_parsed, vm);
		
		unreg_opts = bpo::collect_unrecognized(unreg_parsed.options, bpo::include_positional);
		if (!unreg_opts.empty()) {
			BOOST_THROW_EXCEPTION(bpo::unknown_option{unreg_opts.front()});
		}
		
		if (vm.count("help")) {
			bpo::options_description visible;
			visible.add(generic).add(misc).add(cmd_options);
			
			std::stringstream help;
			help
				<< "Usage: " << exe.filename().string() << " [generic/misc-options] pca [pca-options]" << std::endl
				<< visible;
			return help_cmd{g_opts, help.str()};
		}
		
		pca_cmd cmd;
		cmd.g_opts = g_opts;
		cmd.i_opts = parse_theta_input_options(vm);
		cmd.o_opts = parse_theta_output_options(cmd.i_opts, vm);
		
		if (vm.count("backend")) {
			std::string backend = vm["backend"].as<std::string>();
			
			BOOST_ASSERT(!backend.empty());
			
			if (boost::iequals(backend, "MATLAB_LAPACK")) {
				#ifdef HBRS_MPL_ENABLE_ADDON_MATLAB
					cmd.pca_opts.backend = pca_backend::matlab_lapack;
				#else
					BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
						(boost::format("pca backend %s was not enabled during build") % backend).str()
					});
				#endif
			} else if (boost::iequals(backend, "ELEMENTAL_OPENMP")) {
				#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
					cmd.pca_opts.backend = pca_backend::elemental_openmp;
				#else
					BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
						(boost::format("pca backend %s was not enabled during build") % backend).str()
					});
				#endif
			} else if (boost::iequals(backend, "ELEMENTAL_MPI")) {
				#ifdef HBRS_MPL_ENABLE_ADDON_ELEMENTAL
					cmd.pca_opts.backend = pca_backend::elemental_mpi;
				#else
					BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
						(boost::format("pca backend %s was not enabled during build") % backend).str()
					});
				#endif
			} else {
				BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
					(boost::format("pca backend %s is unknown / not supported") % backend).str()
				});
			}
			
			if (mpi::size() > 1 && cmd.pca_opts.backend != pca_backend::elemental_mpi) {
				BOOST_THROW_EXCEPTION(bpo::invalid_option_value{
					(boost::format("Running with %d processes, but pca backend %s is single process only. Use ELEMENTAL_MPI instead!") % mpi::size() % backend).str()
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

/* unnamed namespace */ }

int
main(int argc, char *argv[]) {
	namespace bpo = boost::program_options;
	using namespace hbrs::theta_utils;
	
	mpl::detail::environment env{argc, argv};
	
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
				} catch(mpl::mpi_exception & ex) {
					mpi::abort();
					std::cerr << boost::diagnostic_information(ex, true) << std::endl;
					return EXIT_FAILURE;
				} catch(boost::exception & ex) {
					std::cerr << boost::diagnostic_information(ex, true) << std::endl;
					return EXIT_FAILURE;
				} catch(std::exception & ex) {
					std::cerr << boost::diagnostic_information(ex, true) << std::endl;
					return EXIT_FAILURE;
				}
			}
			return EXIT_SUCCESS;
		}, cmd);
}
