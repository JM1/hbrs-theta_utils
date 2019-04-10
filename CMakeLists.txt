# Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

#################### project settings ####################

project(hbrs-theta_utils VERSION 2019.0.0.0)

include(FeatureSummary)

#################### options ####################

option(HBRS_THETA_UTILS_ENABLE_TESTS "Build unit tests." ON)

#################### find all used packages ####################

# NOTE: Keep list of packages in sync with dependencies in *-cmake.in!

find_package(hbrs-cmake REQUIRED)
set_package_properties(hbrs-cmake PROPERTIES
    PURPOSE "Required for CMake modules."
    TYPE REQUIRED
)

find_package(hbrs-mpl REQUIRED)
set_package_properties(hbrs-mpl PROPERTIES
    PURPOSE "Required for linear algebra algorithms."
    TYPE REQUIRED
)

find_package(Boost 1.62 REQUIRED COMPONENTS unit_test_framework system thread serialization chrono program_options filesystem iostreams regex exception signals)
set_package_properties(Boost PROPERTIES
    PURPOSE "Required for meta programming, unit tests and others."
    TYPE REQUIRED
)

find_package(MPI REQUIRED)
set_package_properties(MPI PROPERTIES
    PURPOSE "Required for distributed linear algebra algortihms."
    TYPE REQUIRED
)

find_package(netcdf REQUIRED)
set_package_properties(netcdf PROPERTIES
    PURPOSE "Required for reading theta result files."
    TYPE REQUIRED
)

find_package(VTK REQUIRED)
set_package_properties(VTK PROPERTIES
    PURPOSE "Required for writing visualization output files."
    TYPE REQUIRED
)

#################### source settings ####################

# put the include dirs which are in the source or build tree before all other include dirs, so the headers in the sources are preferred over the already installed ones
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
# Src: kdelibs/cmake/modules/KDE4Defaults.cmake

if(VTK_FOUND)
    include(${VTK_USE_FILE})
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# C++14 is required e.g. because of decltype(auto) as return type and constexpr std::forward_as_tuple()!
# C++17 is required for constexpr lambdas

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

if (netcdf_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${netcdf_CFLAGS}")
endif()

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

#################### list the subdirectories ####################

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(test)

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
