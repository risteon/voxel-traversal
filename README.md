# Voxel Traversal

[![Build Status](https://github.com/risteon/voxel-traversal/actions/workflows/test.yml/badge.svg)](https://github.com/risteon/voxel-traversal/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![voxel traversal cover image](voxel_traversal.png?raw=true)

A small library to compute voxel traversal on the CPU.
The implemented algorithm is J. Amanatides, A. Woo: 
["A Fast Voxel Traversal Algorithm"](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf).
It is based on the code prototype of Chris Gyurgyik
[Fast-Voxel-Traversal-Algorithm](https://github.com/cgyurgyik/fast-voxel-traversal-algorithm).

The contributions of this repository are:
* **Tests!**
* Python bindings (in progress...)
* Installation and cmake packaging
* Use Eigen for readability, vectorization, and grid counting.
* Execution and timing on real LiDAR data from the [nuScenes dataset](https://www.nuscenes.org/). The demonstration data files are bundled with git lfs.

## Requirements and Dependencies
* Eigen3
* C++17 compiler

## Python bindings

```bash
# setup python environment, make cmake >= 3.21 available, e.g. with conda install cmake=3.22
$ python setup.py install
# or
$ python setup.py develop
$ pytest tests
```

## Run the Tests and Install
```bash
$ git clone https://github.com/risteon/voxel-traversal.git
$ mkdir build && cd build
# set CMAKE_INSTALL_PREFIX to your liking
# specify python version with -DPYTHON_EXECUTABLE=<path/to/python> or -DPYBIND11_PYTHON_VERSION=3.XX
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ..
# build and install
$ cmake --build . --target install -- -j 4
# run tests
$ ctest
```

## Find with cmake & Usage 

Your project's CMakeLists.txt:
```cmake
cmake_minimum_required(VERSION 3.21)
project(your_project)
set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VoxelTraversal REQUIRED)

add_executable(your_executable main.cpp)
target_link_libraries(your_executable VoxelTraversal::VoxelTraversal)
```
When configuring, set `CMAKE_PREFIX_PATH` to this project's install directory.

Code example:
```c++
#include <VoxelTraversal/voxel_traversal.h>

using namespace voxel_traversal;

// types for vectors and indices
using V3 = typename Grid3DSpatialDef<double>::Vector3d;
using C3 = typename Grid3DSpatialDef<double>::Index3d;
using R = Ray<double>;

const V3 bound_min(0.0, 0.0, 0.0);
const V3 bound_max(2.0, 2.0, 2.0);
const C3 voxel_count(2, 2, 2);
Grid3DSpatialDef<double> grid(bound_min, bound_max, voxel_count);
// use this subclass to count traversed voxels
Grid3DTraversalCounter<double> grid_counter(bound_min, bound_max, voxel_count);

// Ray
const auto ray = R::fromOriginDir({.5, .5, .5}, {1., 0., 0.});

// determine which voxels are traversed (in order from origin to end)
TraversedVoxels<TypeParam> traversed{};
const auto does_intersect = traverseVoxelGrid(ray, grid, traversed);
// count traversed voxels
traverseVoxelGrid(ray, grid_counter);
```
