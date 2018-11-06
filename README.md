# hbrs-theta_utils
[![Build Status](https://travis-ci.com/JM1/hbrs-theta_utils.svg?branch=master)](https://travis-ci.com/JM1/hbrs-theta_utils)

# How to build this code using Docker

```sh
# install Docker CE for Debian or derivatives
# please follow guide at https://docs.docker.com/install/linux/docker-ce/debian/

# docker version 18.06.0-ce or later is recommended
docker --version

# fetch and setup dockerctl
git clone --depth 1 https://gist.github.com/ab8c3beea108ea14a6b8955050f36357.git
export PATH="$PWD/ab8c3beea108ea14a6b8955050f36357/:$PATH"

# fetch docker image
docker pull jm1337/debian-testing-cpp-dev

# log into docker container
dockerctl run --login --systemd -- jm1337/debian-testing-cpp-dev



# the following commands are executed from within the docker container

# fetch, compile and install prerequisites
git clone --depth 1 https://github.com/JM1/hbrs-cmake.git
cd hbrs-cmake
mkdir build && cd build/
cmake -DHBRS_CMAKE_ENABLE_TESTS=OFF ..
make -j$(nproc)
sudo make install
cd ../../

git clone --depth 1 https://github.com/JM1/hbrs-mpl.git
cd hbrs-mpl
mkdir build && cd build/
cmake -DHBRS_MPL_ENABLE_TESTS=OFF -DHBRS_MPL_ENABLE_ADDON_ELEMENTAL=ON -DHBRS_MPL_ENABLE_ADDON_MATLAB=OFF ..
make -j$(nproc)
sudo make install
cd ../../

# fetch, compile and install hbrs-theta_utils
git clone --depth 1 https://github.com/JM1/hbrs-theta_utils.git
cd hbrs-theta_utils
mkdir build && cd build/
cmake -DHBRS_THETA_UTILS_ENABLE_TESTS=ON ..
make -j$(nproc)
ctest --output-on-failure
sudo make install
```

For more examples how to build and test code see [`.travis.yml`](https://github.com/JM1/hbrs-theta_utils/blob/master/.travis.yml).
