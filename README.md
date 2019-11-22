# `hbrs-theta_utils`
[![Build Status](https://travis-ci.com/JM1/hbrs-theta_utils.svg?branch=master)](https://travis-ci.com/JM1/hbrs-theta_utils)

`hbrs-theta_utils` is a postprocessing tool to CFD solver [TAU (and THETA)](http://tau.dlr.de/) for data-driven modal decompositions.
TAU is for numerical flow simulations, e.g. predicting flows around airfoils of wind turbines or aircrafts.
`hbrs-theta_utils` allows is to decompose these (in)compressible flows by various methods, e.g. [PCA/POD](https://en.wikipedia.org/wiki/Principal_component_analysis) or [DMD](https://en.wikipedia.org/wiki/Dynamic_mode_decomposition), to extract relevant features.
It also allows to export these flow fields into a [VTK](https://vtk.org/) files for visualization, e.g. with [ParaView](https://www.paraview.org/).
It is designed operate in parallel on distributed large-scale datasets at [HPC clusters (e.g. Platform for Scientific Computing at BRSU)](https://wr0.wr.inf.h-brs.de/wr/index.html).

`hbrs-theta_utils` is a tool for coherent structure analysis in fluid dynamics.
It helps scientists and engineers with characterizing energetic structures from data and thus making fluids e.g. more tractable to analysis and engineering design.

## History

Its development started in 2015 as a research project at Bonn-Rhein-Sieg University of Applied Sciences, from 2016-2019 it was funded by BMBF project [AErOmAt](https://www.h-brs.de/de/aeromat) and till today its actively developed by its initial author, Jakob Meng.

## Under the hood

`hbrs-theta_utils` is a CLI application written in [C++17](https://en.wikipedia.org/wiki/C++17) and using [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) for distributed computations at HPC clusters.
Its `pca` command first reads time series of 3d velocity fields from distributed [netCDF](https://www.unidata.ucar.edu/software/netcdf/) files, which were e.g. created during flow simulations in TAU or THETA.
These snapshots of x-, y- and z-velocities are then decomposed using [Principal Component Analysis (PCA) / Proper Orthogonal Decomposition (POD)](https://en.wikipedia.org/wiki/Principal_component_analysis) from [hbrs-mpl](https://github.com/JM1/hbrs-mpl/), a generic C++ library for math and statistics.
PCA is a dimensionality reduction technique that transforms data in high-dimensional space to a space of fewer dimensions.
With command line argument `--pcs` the user selects which principal components to keep and drop.
The reassembled and possibly reduced dataset is then written to disk using the same distributed file format as the input.

The `visualize` command reads an unstructured 3d grid and a time series of 3d velocity fields, both from distributed [netCDF](https://www.unidata.ucar.edu/software/netcdf/) files.
The grid contains all geometries (tetraeders, prisms, surfacetriangles, ...) that are used within flow simulations, stored in a proprietary TAU format.
Before simulation, this grid is split and distributed across MPI processes, using TAU's preprocessing tool.
For visualization, each process has to exchange x-, y- and z-velocities at its local grid boundaries with neighbours.
Then, for each simulated time step, the 3d velocity field and the grid are written as a set of `*.pvtu` files to disk.
For that, the [Visualization Toolkit (VTK)](https://vtk.org/) is used to export the 3d geometry objects, their surfaces colored with 3d velocities, to [parallel unstructured grid (`*.pvtu`) files](https://www.vtk.org/VTK/img/file-formats.pdf).
These `*.pvtu` files can then be opened and viewed in [ParaView](https://www.paraview.org/) or using [pvserver](https://www.paraview.org/Wiki/Setting_up_a_ParaView_Server) for distributed visualization on a cluster.

All functionality is thoroughly being checked using [automated and extensive unit tests](https://travis-ci.com/JM1/hbrs-theta_utils/).
It has been applied to a real-world dataset, the airflows around a side mirror of a [car](https://www.aer.mw.tum.de/en/research-groups/automotive/drivaer/), utilizing 19 Mio. grid points, 1.000 time steps, 330GB simulation files, 600 MPI processes and 1.6TB of RAM.
Its decompositions have been verified by scientists from DLR.

## About TAU & THETA

TAU is a ["software system for the prediction of viscous and inviscid flows about complex geometries from the low subsonic to the hypersonic flow regime, employing hybrid unstructured grids"](http://tau.dlr.de/).
It is developed and distributed by German Aerospace Center a.k.a. Deutsches Zentrum für Luft- und Raumfahrt (DLR).

> TAU\
> The DLR-TAU-Code (TAU=Triangular Adaptive Upwind) is a software for the numerical flow simulation based on the (U)RANS or the hybrid RANS/LES approach using a finite-volume discretization for adaptable unstructured grids. TAU allows for flow predictions around complex moving geometries over a wide range of Mach numbers and has been established as a production code in the European aircraft industry, as well as a research tool for new aerospace technologies.
>
> THETA\
> The DLR-THETA-Code (THETA=Turbulent Heat Release Extension of the TAU-Code) was developed for the simulation of incompressible combustion chamber flows. Further areas of application include two-phase flow using the Volume-of-Fluid method, which is used to simulate fuel sloshing in tanks of upper stages of rockets, as well as the high-fidelity simulation of flows around wind turbines and in the surrounding terrain taking important atmospheric parameters into account.

References:
- [tau.dlr.de](http://tau.dlr.de/)
- [DLR : Institute of Aerodynamics and Flow Technology : C²A²S²E (AS-CAS) : Software Products / Research Tools](https://www.dlr.de/as/en/desktopdefault.aspx/tabid-4083/6455_read-9239/)
- [DLR : Institute of Combustion Technology: Computer Simulation : THETA Code](https://www.dlr.de/vt/de/desktopdefault.aspx/tabid-3082/4659_read-15475/)

# How to build this code using Docker

```sh
# install Docker CE for Debian or derivatives
# please follow guide at https://docs.docker.com/install/linux/docker-ce/debian/

# docker version 18.06.0-ce or later is recommended
docker --version

# fetch docker image
docker pull jm1337/debian-dev-hbrs:buster

# log into docker container
docker run -ti jm1337/debian-dev-hbrs:buster



# the following commands are executed from within the docker container

# fetch, compile and install prerequisites
git clone --depth 1 https://github.com/JM1/hbrs-cmake.git
cd hbrs-cmake
mkdir build && cd build/
cmake ..
make -j$(nproc)
sudo make install
cd ../../

git clone --depth 1 https://github.com/JM1/hbrs-mpl.git
cd hbrs-mpl
mkdir build && cd build/
cmake -DHBRS_MPL_ENABLE_TESTS=OFF -DHBRS_MPL_ENABLE_BENCHMARKS=OFF -DHBRS_MPL_ENABLE_ADDON_ELEMENTAL=ON -DHBRS_MPL_ENABLE_ADDON_MATLAB=OFF ..
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
