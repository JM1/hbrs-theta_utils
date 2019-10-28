# Copyright (c) 2019 Jakob Meng, <jakobmeng@web.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

#################### tests ####################

foreach(name
    pca)
    hbrs_theta_utils_add_test(fn_execute_${name} "${name}.cpp")
endforeach()
