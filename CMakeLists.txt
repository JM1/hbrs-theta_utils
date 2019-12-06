# Copyright (c) 2018-2019 Jakob Meng, <jakobmeng@web.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

#################### project settings ####################

project(hbrs-theta_utils VERSION 2019.0.0.0)

include(FeatureSummary)

#################### options ####################

option(HBRS_THETA_UTILS_ENABLE_TESTS "Build unit tests." OFF)

#################### find all used packages ####################

# NOTE: Keep list of packages in sync with dependencies in *-cmake.in!

find_package(hbrs-cmake)
set_package_properties(hbrs-cmake PROPERTIES
    PURPOSE "Required for CMake modules."
    TYPE REQUIRED)

find_package(hbrs-mpl)
set_package_properties(hbrs-mpl PROPERTIES
    PURPOSE "Required for linear algebra algorithms."
    TYPE REQUIRED)

find_package(Boost 1.62 COMPONENTS unit_test_framework system thread serialization chrono program_options filesystem iostreams regex exception signals log)
set_package_properties(Boost PROPERTIES
    PURPOSE "Required for meta programming, unit tests and others."
    TYPE REQUIRED)

find_package(MPI)
set_package_properties(MPI PROPERTIES
    PURPOSE "Required for distributed linear algebra algortihms."
    TYPE REQUIRED)

find_package(netcdf)
set_package_properties(netcdf PROPERTIES
    PURPOSE "Required for reading theta result files."
    TYPE REQUIRED)

find_package(VTK)
set_package_properties(VTK PROPERTIES
    PURPOSE "Required for writing visualization output files."
    TYPE REQUIRED)

feature_summary(WHAT REQUIRED_PACKAGES_NOT_FOUND FATAL_ON_MISSING_REQUIRED_PACKAGES)

#################### source settings ####################

# put the include dirs which are in the source or build tree before all other include dirs, so the headers in the sources are preferred over the already installed ones
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
# Src: kdelibs/cmake/modules/KDE4Defaults.cmake

if(VTK_FOUND)
    include(${VTK_USE_FILE})
endif()

set(CXX_TEMPLATE_DEPTH 250            CACHE STRING "Set the maximum instantiation depth for template classes to n (-ftemplate-depth=n).") # increased depth is required by operators.hpp
set(CXX_TEMPLATE_BACKTRACE_LIMIT 250  CACHE STRING "Set the maximum number of template instantiation notes for a single warning or error to n (-ftemplate-backtrace-limit=n).")
set(CXX_ERROR_LIMIT 0                 CACHE STRING "Limits the maximum number of error messages to n (-fmax-errors=n / -ferror-limit=n).")
set(CXX_CONSTEXPR_BACKTRACE_LIMIT 250 CACHE STRING "-fconstexpr-backtrace-limit=n")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# "CMake will generate tests only if the enable_testing() command has been invoked." [1]
# "Note that ctest expects to find a test file in the build directory root. Therefore, this command should be in the 
#  source directory root." [2]
# References:
# [1] https://cmake.org/cmake/help/latest/command/add_test.html#command:add_test
# [2] https://cmake.org/cmake/help/latest/command/enable_testing.html
if(HBRS_THETA_UTILS_ENABLE_TESTS)
    enable_testing()
endif()

#################### build settings ####################

include(MaybeAddFlags)
maybe_add_c_flag  (HBRS_THETA_UTILS_HAS_C_WALL                         "-Wall")
maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_WALL                       "-Wall")
#maybe_add_c_flag  (HBRS_THETA_UTILS_HAS_C_WEXTRA                       "-Wextra")
#maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_WEXTRA                     "-Wextra")
maybe_add_c_flag  (HBRS_THETA_UTILS_HAS_C_WPEDANTIC                    "-Wpedantic")
maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_WPEDANTIC                  "-Wpedantic")
maybe_add_c_flag  (HBRS_THETA_UTILS_HAS_C_WERROR_RETURN_TYPE           "-Werror=return-type")
maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_WERROR_RETURN_TYPE         "-Werror=return-type")
maybe_add_c_flag  (HBRS_THETA_UTILS_HAS_C_PEDANTIC_ERRORS              "-pedantic-errors")
maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_PEDANTIC_ERRORS            "-pedantic-errors")

maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_FTEMPLATE_DEPTH            "-ftemplate-depth=${CXX_TEMPLATE_DEPTH}")
maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_FTEMPLATE_BACKTRACE_LIMIT  "-ftemplate-backtrace-limit=${CXX_TEMPLATE_BACKTRACE_LIMIT}")
maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_FCONSTEXPR_BACKTRACE_LIMIT "-fconstexpr-backtrace-limit=${CXX_CONSTEXPR_BACKTRACE_LIMIT}")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_FERROR_LIMIT           "-ferror-limit=${CXX_ERROR_LIMIT}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    maybe_add_cxx_flag(HBRS_THETA_UTILS_HAS_CXX_FMAX_ERRORS            "-fmax-errors=${CXX_ERROR_LIMIT}")
endif()

if(NOT HBRS_THETA_UTILS_HAS_CXX_FTEMPLATE_DEPTH)
    message(FATAL_ERROR "Your C++ compiler ${CMAKE_CXX_COMPILER} does not support '-ftemplate-depth=${CXX_TEMPLATE_DEPTH}', use it at your own risk.")
endif()

if(NOT HBRS_THETA_UTILS_HAS_CXX_FTEMPLATE_BACKTRACE_LIMIT)
    message(WARNING "Your C++ compiler ${CMAKE_CXX_COMPILER} does not support '-ftemplate-backtrace-limit=${CXX_TEMPLATE_BACKTRACE_LIMIT}', use it at your own risk.")
endif()

if(NOT HBRS_THETA_UTILS_HAS_CXX_FCONSTEXPR_BACKTRACE_LIMIT)
    message(WARNING "Your C++ compiler ${CMAKE_CXX_COMPILER} does not support '-fconstexpr-backtrace-limit=${CXX_CONSTEXPR_BACKTRACE_LIMIT}', use it at your own risk.")
endif()

#################### install settings ####################

include(GNUInstallDirs)
# Src.: http://www.cmake.org/cmake/help/v2.8.8/cmake.html#module:GNUInstallDirs

# Offer a choice of overriding the installation directories
set(HBRS_THETA_UTILS_INSTALL_BINDIR      ${CMAKE_INSTALL_BINDIR}                               CACHE PATH "install dir for executables")
set(HBRS_THETA_UTILS_INSTALL_LIBDIR      ${CMAKE_INSTALL_LIBDIR}                               CACHE PATH "install dir for libraries")
set(HBRS_THETA_UTILS_INSTALL_INCLUDEDIR  ${CMAKE_INSTALL_INCLUDEDIR}/hbrs/theta_utils          CACHE PATH "install dir for headers")
set(HBRS_THETA_UTILS_INSTALL_CMAKEDIR    ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}   CACHE PATH "install dir for CMake files")

set(HBRS_THETA_UTILS_INSTALL_TARGETS_DEFAULT_DESTINATIONS
    RUNTIME DESTINATION "${HBRS_THETA_UTILS_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${HBRS_THETA_UTILS_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${HBRS_THETA_UTILS_INSTALL_LIBDIR}"
    PUBLIC_HEADER DESTINATION "${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}"
    PRIVATE_HEADER DESTINATION "${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}"
    INCLUDES DESTINATION "${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}")


#################### tests ####################

# CMake's target test does not build but only executes tests
# Ref.: https://stackoverflow.com/a/736838/6490710
if(HBRS_THETA_UTILS_ENABLE_TESTS)
   set(_TESTS_MAYBE_ALL ALL)
endif()

add_custom_target(tests ${_TESTS_MAYBE_ALL} COMMENT "Build all tests.")

include(DevTools)
function(hbrs_theta_utils_add_test target)
    add_executable(${target} EXCLUDE_FROM_ALL "${ARGN}")
    target_link_libraries(${target} hbrs_theta_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_include_directories(${target}
        SYSTEM
        PUBLIC ${Boost_INCLUDE_DIRS})

    target_compile_definitions(${target} PUBLIC ${MPI_CXX_COMPILE_DEFINITIONS} BOOST_LOG_DYN_LINK)

    add_test(${target} ${CMAKE_CURRENT_BINARY_DIR}/${target})
    add_dependencies(tests ${target})
endfunction()

#################### list the subdirectories ####################

add_subdirectory(src)

#################### install files ####################

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/hbrs-theta_utils-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/hbrs-theta_utils-config.cmake"
    INSTALL_DESTINATION "${HBRS_THETA_UTILS_INSTALL_CMAKEDIR}"
    PATH_VARS
        HBRS_THETA_UTILS_INSTALL_BINDIR
        HBRS_THETA_UTILS_INSTALL_LIBDIR
        HBRS_THETA_UTILS_INSTALL_INCLUDEDIR
        HBRS_THETA_UTILS_INSTALL_CMAKEDIR)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/hbrs-theta_utils-config-version.cmake"
    VERSION ${hbrs-theta_utils_VERSION}
    COMPATIBILITY AnyNewerVersion)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/hbrs-theta_utils-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/hbrs-theta_utils-config-version.cmake"
    DESTINATION ${HBRS_THETA_UTILS_INSTALL_CMAKEDIR}
    COMPONENT development)

install(EXPORT hbrs-theta_utils-targets
  FILE hbrs-theta_utils-targets.cmake
  NAMESPACE hbrs-theta_utils::
  DESTINATION ${HBRS_THETA_UTILS_INSTALL_CMAKEDIR})

#################### summary ####################

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)
