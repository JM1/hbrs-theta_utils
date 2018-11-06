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

#################### install files ####################

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/config.hpp
    DESTINATION "${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}"
    COMPONENT development)

install(DIRECTORY ./
    DESTINATION "${HBRS_THETA_UTILS_INSTALL_INCLUDEDIR}"
    COMPONENT development
    FILES_MATCHING PATTERN "*.hpp")
