#ifndef VOXEL_TRAVERSAL_PYTHON_BINDINGS_H
#define VOXEL_TRAVERSAL_PYTHON_BINDINGS_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "grid.h"

namespace pytraversal {

using grid_type = voxel_traversal::Grid3DSpatialDef<double>;

pybind11::array_t<int64_t, pybind11::array::c_style> traverse(
    const grid_type& grid, const grid_type::Vector3d& ray_origin,
    const grid_type::Vector3d& ray_end);

}  // namespace pytraversal

#endif  // VOXEL_TRAVERSAL_PYTHON_BINDINGS_H
