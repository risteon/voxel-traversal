# Voxel Traversal

A small library to compute voxel traversal on the CPU.
The implemented algorithm is J. Amanatides, A. Woo: 
["A Fast Voxel Traversal Algorithm"](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf).
It is based on the code prototype of Chris Gyurgyik
[Fast-Voxel-Traversal-Algorithm](https://github.com/cgyurgyik/fast-voxel-traversal-algorithm).

The contributions of this repository are:
* **Tests!**
* C++17 compatibility
* Installation and cmake packaging
* Use Eigen for readability, vectorization, and grid counting.
* Execution and timing on real LiDAR data from the [nuScenes dataset](https://www.nuscenes.org/). The demonstration data files are bundled with git lfs.

## Run the tests and install
```bash
$ git clone https://github.com/risteon/voxel-traversal.git
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ..
$ cmake --build . --target all
$ cmake --build . --target install
$ ctest
```

## Missing parts
* More tests
* Algorithm and structures as templates for single and double precision version.