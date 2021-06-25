# Postprocessing tool `hbrs-theta_utils` to CFD solvers TAU and THETA for data-driven modal decompositions

`hbrs-theta_utils` ([GitHub.com][hbrs-theta-utils], [H-BRS GitLab][hbrs-gitlab-hbrs-theta-utils]) is a tool for coherent
structure analysis in fluid dynamics. It helps scientists and engineers with characterizing energetic structures in CFD
data and with making fluids more tractable to analysis and engineering design. Its development started in 2015 as a
research project at [Bonn-Rhein-Sieg University of Applied Sciences][hbrs], from 2016-2019 it was funded by BMBF
project [AErOmAt][aeromat].

`hbrs-theta_utils` is a postprocessing tool to CFD solvers [TAU and THETA][tau] for data-driven modal decompositions.
TAU and THETA do numerical flow simulations, e.g. predicting flows around airfoils of wind turbines or aircrafts.
`hbrs-theta_utils` allows us to decompose these (in)compressible flows using various methods, e.g. [PCA/POD][wiki-pca]
or [DMD][wiki-dmd], to extract relevant features. It also allows to export these flow fields into a [VTK][vtk] files for
visualization, e.g. with [ParaView][paraview]. It is designed operate in parallel on distributed large-scale datasets at
HPC clusters like the [Platform for Scientific Computing at BRSU][hbrs-wr].

## Under the hood

`hbrs-theta_utils` is a CLI application written in [C++17][wiki-cpp17] and using [MPI][wiki-mpi] for distributed
computations at HPC clusters.

Its `pca` command first reads time series of 3d velocity fields from distributed [netCDF][netcdf] files, which were e.g.
created during flow simulations in TAU or THETA. These snapshots of x-, y- and z-velocities are then decomposed using
[Principal Component Analysis (PCA) / Proper Orthogonal Decomposition (POD)][wiki-pca] from [hbrs-mpl][hbrs-mpl], a
generic C++ library for math and statistics. PCA is a dimensionality reduction technique that transforms data in
high-dimensional space to a space of fewer dimensions. With command line argument `--pcs` the user selects which
principal components to keep and drop. The reassembled and possibly reduced dataset is then written to disk using the
same distributed file format as the input.

The `visualize` command reads an unstructured 3d grid and a time series of 3d velocity fields, both from distributed
[netCDF][netcdf] files. The grid contains all geometries (tetraeders, prisms, surfacetriangles, ...) that are used
within flow simulations, stored in a proprietary TAU format. Before simulation, this grid is split and distributed
across MPI processes, using TAU's preprocessing tool. For visualization, each process has to exchange x-, y- and
z-velocities at its local grid boundaries with neighbours. Then, for each simulated time step, the 3d velocity field and
the grid are written as a set of `*.pvtu` files to disk. For that, the [Visualization Toolkit (VTK)][vtk] is used to
export the 3d geometry objects, their surfaces colored with 3d velocities, to
[parallel unstructured grid (`*.pvtu`) files][vtk-file-formats]. These `*.pvtu` files can then be opened and viewed in
[ParaView][paraview] or using [pvserver][pvserver-setup] for distributed visualization on a cluster.

All functionality is heavily being tested using [automated and extensive unit tests][hbrs-gitlab-hbrs-theta-utils-ci].
It has been applied to a real-world dataset, the airflows around a side mirror of a [car][drivaer], utilizing 19 Mio.
grid points, 1.000 time steps, 330GB simulation files, 600 MPI processes and 1.6TB of RAM. Its decompositions have been
verified by scientists from DLR.

`hbrs-theta_utils` builds heavily upon C++ libraries [`hbrs-mpl`][hbrs-mpl] and [`Elemental`][elemental] which provide
HPC-ready data structures and algorithms for linear algebra and dimension reduction.

The full tech stack consists of:
* [`C++17`][cpp-ref] for generic and efficient library code
* C++ library [`hbrs-mpl`][hbrs-mpl] ([GitHub.com][hbrs-mpl], [H-BRS GitLab][hbrs-gitlab-hbrs-mpl])
* C++ library [`Elemental`][elemental]
* C++ metaprogramming library [`Boost.Hana`][boost-hana-ref]
* [MPI][wiki-mpi] for [distributed][hbrs-theta-utils-detail-vtk] [computations][hbrs-theta-utils-detail-scatter]
* [netCDF][netcdf] for reading and writing unstructured 3d grids and time series of 3d velocity fields
* [Visualization Toolkit (VTK)][vtk] for exporting 3d geometry objects to
  [parallel unstructured grid (`*.pvtu`) files][vtk-file-formats] for [ParaView][paraview]
* [`Boost.Test`][boost-test] for unit tests, e.g.
  [`hbrs::theta_utils::fn::execute`][hbrs-theta-utils-fn-execute-test-pca] and
  [`hbrs::theta_utils::dt::theta_field`][hbrs-theta-utils-dt-theta-field-test]
* [CMake 3][cmake3-tut] and [`hbrs-cmake`][hbrs-cmake] to [build, export and install our library](CMakeLists.txt)
* [GitLab CI][hbrs-gitlab-hbrs-theta-utils-ci] to continuously build and test our code

## About TAU & THETA

[TAU][tau] is a
> software system for the prediction of viscous and inviscid flows about complex geometries from the low subsonic to 
> the hypersonic flow regime, employing hybrid unstructured grids.

It is developed and distributed by German Aerospace Center aka Deutsches Zentrum für Luft- und Raumfahrt (DLR).

> TAU\
> The DLR-TAU-Code (TAU=Triangular Adaptive Upwind) is a software for the numerical flow simulation based on the (U)RANS
> or the hybrid RANS/LES approach using a finite-volume discretization for adaptable unstructured grids. TAU allows for
> flow predictions around complex moving geometries over a wide range of Mach numbers and has been established as a
> production code in the European aircraft industry, as well as a research tool for new aerospace technologies.
>
> THETA\
> The DLR-THETA-Code (THETA=Turbulent Heat Release Extension of the TAU-Code) was developed for the simulation of
> incompressible combustion chamber flows. Further areas of application include two-phase flow using the Volume-of-Fluid
> method, which is used to simulate fuel sloshing in tanks of upper stages of rockets, as well as the high-fidelity
> simulation of flows around wind turbines and in the surrounding terrain taking important atmospheric parameters into
> account.

For details about TAU and THETA navigate to:
* [`DLR > Institute of Aerodynamics and Flow Technology > Departments > C²A²S²E (AS-CAS)`][dlr-as-case]
* [`DLR > Institute of Combustion Technology > Computer Simulation > THETA Code`][dlr-vt-theta]

## How to build, install and run code using `Docker` or `Podman`

For a quick and easy start into developing with C++, a set of ready-to-use `Docker`/`Podman` images
`jm1337/debian-dev-hbrs` and `jm1337/debian-dev-full` (supports more languages) has been created. They contain a full
development system including all tools and libraries necessary to hack on distributed decomposition algorithms and more
([Docker Hub][docker-hub-jm1337], [source files for Docker images][docker-artifacts]).

### Install `Docker` or `Podman`

* On `Debian 10 (Buster)` or `Debian 11 (Bullseye)` just run `sudo apt install docker.io`
  or follow the [official install guide][docker-install-debian] for Docker Engine on Debian
* On `Ubuntu 18.04 LTS (Bionic Beaver)` and `Ubuntu 20.04 LTS (Focal Fossa)` just run `sudo apt install docker.io`
  (from `bionic/universe` and `focal/universe` repositories)
  or follow the [official install guide][docker-install-ubuntu] for Docker Engine on Ubuntu
* On `Windows 10` follow the [official install guide][docker-install-windows] for Docker Desktop on Windows
* On `Mac` follow the [official install guide][docker-install-mac] for Docker Desktop on Mac
* On `Fedora`, `Red Hat Enterprise Linux (RHEL)` and `CentOS` follow the [official install guide][podman-install] for
  Podman

### Setup and run container

```sh
# docker version 18.06.0-ce or later is recommended
docker --version

# fetch docker image
docker pull jm1337/debian-dev-hbrs:bullseye

# log into docker container
docker run -ti jm1337/debian-dev-hbrs:bullseye
# or using a persistent home directory, e.g.
docker run -ti -v /HOST_DIR:/home/devil/ jm1337/debian-dev-hbrs:bullseye
# or using a persistent home directory on Windows hosts, e.g.
docker run -ti -v C:\YOUR_DIR:/home/devil/ jm1337/debian-dev-hbrs:bullseye
```

Podman strives for complete CLI compatibility with Docker, hence
[you may use the `alias` command to create a `docker` alias for Podman][docker-to-podman-transition]:
```sh
alias docker=podman
```

### Build and run code inside container

Execute the following commands within the `Docker`/`Podman` container:

```sh
# choose a compiler
export CC=clang-10
export CXX=clang++-10
# or
export CC=gcc-10
export CXX=g++-10

# fetch, compile and install prerequisites
git clone --depth 1 https://github.com/JM1/hbrs-cmake.git
cd hbrs-cmake
mkdir build && cd build/
# install to non-system directory because sudo is not allowed in this docker container
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/.local \
    ..
make -j$(nproc)
make install
cd ../../

git clone --depth 1 https://github.com/JM1/hbrs-mpl.git
cd hbrs-mpl
mkdir build && cd build/
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/.local \
    -DHBRS_MPL_ENABLE_ELEMENTAL=ON \
    -DHBRS_MPL_ENABLE_MATLAB=OFF \
    -DHBRS_MPL_ENABLE_TESTS=OFF \
    -DHBRS_MPL_ENABLE_BENCHMARKS=OFF \
    ..
make -j$(nproc)
make install
cd ../../

# fetch, compile and install hbrs-theta_utils
git clone --depth 1 https://github.com/JM1/hbrs-theta_utils.git
cd hbrs-theta_utils
mkdir build && cd build/
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/.local \
    -DHBRS_THETA_UTILS_ENABLE_TESTS=ON \
    ..
make -j$(nproc)
ctest --verbose --output-on-failure
make install
```

For more examples on how to build and test this code see [`.gitlab-ci.yml`](.gitlab-ci.yml).

## License

GNU General Public License v3.0 or later

See [LICENSE.md](LICENSE.md) to see the full text.

## Author

Jakob Meng
@jm1 ([GitHub.com][github-jm1], [Web][jm])

[//]: # (References)

[aeromat]: https://www.h-brs.de/de/aeromat
[boost-hana-ref]: https://boostorg.github.io/hana/
[boost-test]: https://www.boost.org/doc/libs/release/libs/test/
[cmake3-tut]: https://cmake.org/cmake/help/latest/guide/tutorial/index.html
[cpp-ref]: https://en.cppreference.com/w/cpp
[dlr-as-case]: https://www.dlr.de/as/en/desktopdefault.aspx/tabid-4083/6455_read-9239/
[dlr-vt-theta]: https://www.dlr.de/vt/de/desktopdefault.aspx/tabid-3082/4659_read-15475/
[docker-artifacts]: https://github.com/JM1/docker-artifacts
[docker-hub-jm1337]: https://hub.docker.com/r/jm1337/
[docker-install-debian]: https://docs.docker.com/engine/install/debian/
[docker-install-mac]: https://docs.docker.com/docker-for-mac/install/
[docker-install-ubuntu]: https://docs.docker.com/engine/install/ubuntu/
[docker-install-windows]: https://docs.docker.com/docker-for-windows/install/
[docker-to-podman-transition]: https://developers.redhat.com/blog/2020/11/19/transitioning-from-docker-to-podman/
[drivaer]: https://www.aer.mw.tum.de/en/research-groups/automotive/drivaer/
[elemental]: https://github.com/elemental/Elemental
[github-jm1]: https://github.com/jm1
[hbrs]: https://www.h-brs.de
[hbrs-gitlab-hbrs-mpl]: https://git.inf.h-brs.de/jmeng2m/hbrs-mpl/
[hbrs-gitlab-hbrs-theta-utils]: https://git.inf.h-brs.de/jmeng2m/hbrs-theta_utils/
[hbrs-gitlab-hbrs-theta-utils-ci]: https://git.inf.h-brs.de/jmeng2m/hbrs-theta_utils/-/pipelines
[hbrs-cmake]: https://github.com/JM1/hbrs-cmake/
[hbrs-mpl]: https://github.com/JM1/hbrs-mpl/
[hbrs-theta-utils]: https://github.com/JM1/hbrs-theta_utils/
[hbrs-theta-utils-detail-scatter]: https://github.com/JM1/hbrs-theta_utils/blob/master/src/hbrs/theta_utils/detail/scatter/impl.cpp
[hbrs-theta-utils-detail-vtk]: https://github.com/JM1/hbrs-theta_utils/blob/master/src/hbrs/theta_utils/detail/vtk/impl.cpp
[hbrs-theta-utils-dt-theta-field-test]: https://github.com/JM1/hbrs-theta_utils/blob/master/src/hbrs/theta_utils/dt/theta_field/test.cpp
[hbrs-theta-utils-fn-execute-test-pca]: https://github.com/JM1/hbrs-theta_utils/blob/master/src/hbrs/theta_utils/fn/execute/test/pca.cpp
[hbrs-wr]: https://wr0.wr.inf.h-brs.de/wr/index.html
[jm]: http://www.jakobmeng.de
[netcdf]: https://www.unidata.ucar.edu/software/netcdf/
[paraview]: https://www.paraview.org/
[podman-install]: https://podman.io/getting-started/installation
[pvserver-setup]: https://www.paraview.org/Wiki/Setting_up_a_ParaView_Server
[tau]: http://tau.dlr.de/
[vtk]: https://vtk.org/
[vtk-file-formats]: https://www.vtk.org/VTK/img/file-formats.pdf
[wiki-cpp17]: https://en.wikipedia.org/wiki/C++17
[wiki-dmd]: https://en.wikipedia.org/wiki/Dynamic_mode_decomposition
[wiki-mpi]: https://en.wikipedia.org/wiki/Message_Passing_Interface
[wiki-pca]: https://en.wikipedia.org/wiki/Principal_component_analysis
