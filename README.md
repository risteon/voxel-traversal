# Voxel Traversal

A small library to compute voxel traversal on the CPU.
The implemented algorithm is J. Amanatides, A. Woo: 
["A Fast Voxel Traversal Algorithm"](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf).
It is based on the code prototype of Chris Gyurgyik
[Fast-Voxel-Traversal-Algorithm](https://github.com/cgyurgyik/fast-voxel-traversal-algorithm).

The contributions of this repository are:
* **Tests!**
* Installation and cmake packaging
* Use Eigen for readability, vectorization, and grid counting.
* Execution and timing on real LiDAR data from the [nuScenes dataset](https://www.nuscenes.org/). The demonstration data files are bundled with git lfs.

## Requirements and Dependencies
* Eigen3
* C++17 compiler

## Run the Tests and Install
```bash
$ git clone https://github.com/risteon/voxel-traversal.git
$ mkdir build && cd build
# set CMAKE_INSTALL_PREFIX to your liking
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ..
# build and install
$ cmake --build . --target install -- -j 4
# run tests
$ ctest
```

## Include & Usage Example

Your projects' CMakeLists.txt:
```cmake
cmake_minimum_required(VERSION 3.21)
project(your_project)
set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VoxelTraversal REQUIRED)

add_executable(your_executable main.cpp)
target_link_libraries(your_executable VoxelTraversal::VoxelTraversal)
```
When configuring, set `CMAKE_PREFIX_PATH` to this projects install directory.
You will need this include:
```c++
#include <VoxelTraversal/voxel_traversal.h>
```



## ToDo's
* More tests
* Algorithm and structures as templates for single and double precision version.