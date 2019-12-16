# Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

#################### configure files ####################

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/config.hpp"
    @ONLY)

#################### build ####################

add_library(hbrs_theta_utils STATIC)

get_target_property(hbrs_theta_utils_TARGET_TYPE hbrs_theta_utils TYPE)
if (hbrs_theta_utils_TARGET_TYPE STREQUAL "STATIC_LIBRARY")
    target_compile_definitions(hbrs_theta_utils INTERFACE HBRS_THETA_UTILS_STATIC_DEFINE)
endif ()

set_target_properties(hbrs_theta_utils PROPERTIES
    VERSION ${hbrs-theta_utils_VERSION}
    SOVERSION 0
    INTERFACE_hbrs_theta_utils_MAJOR_VERSION ${hbrs-theta_utils_VERSION})

set_property(TARGET hbrs_theta_utils
    APPEND PROPERTY COMPATIBLE_INTERFACE_STRING hbrs_theta_utils_MAJOR_VERSION)

target_include_directories(hbrs_theta_utils
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src> # for export.hpp
        $<INSTALL_INTERFACE:${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}>)

target_include_directories(hbrs_theta_utils
    SYSTEM
    PUBLIC ${MPI_CXX_INCLUDE_DIRS}
    PUBLIC ${Boost_INCLUDE_DIRS}
    PUBLIC ${VTK_INCLUDE_DIRS})

target_compile_definitions(hbrs_theta_utils
    PUBLIC ${MPI_CXX_COMPILE_DEFINITIONS} BOOST_LOG_DYN_LINK)

set_target_properties(hbrs_theta_utils PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_OPTIONS}"
    COMPILE_FLAGS "${netcdf_CFLAGS}"
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

target_link_libraries(hbrs_theta_utils PUBLIC
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${netcdf_LIBRARIES}
    ${VTK_LIBRARIES}
    hbrs-mpl::hbrs_mpl)

include(GenerateExportHeader)
generate_export_header(hbrs_theta_utils
    BASE_NAME HBRS_THETA_UTILS
    EXPORT_FILE_NAME export.hpp)

add_executable(hbrs_theta_utils_cli main.cpp)
target_link_libraries(hbrs_theta_utils_cli hbrs_theta_utils)

#################### list the subdirectories ####################

add_subdirectory(core)
add_subdirectory(detail)
add_subdirectory(dt)
add_subdirectory(fn)

#################### install files ####################

install(TARGETS hbrs_theta_utils hbrs_theta_utils_cli
    EXPORT hbrs-theta_utils-targets
    ${HBRS_THETA_UTILS_INSTALL_TARGETS_DEFAULT_DESTINATIONS})

file(RELATIVE_PATH hbrs_theta_utils_includedir "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/export.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/config.hpp"
    DESTINATION "${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}/${hbrs_theta_utils_includedir}"
    COMPONENT development)

install(DIRECTORY ./
    DESTINATION "${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}/${hbrs_theta_utils_includedir}"
    COMPONENT development
    FILES_MATCHING PATTERN "*.hpp")
