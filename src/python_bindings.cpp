#include "python_bindings.h"

#include <pybind11/stl.h>

#include "voxel_traversal.h"

namespace py = pybind11;
using namespace voxel_traversal;

namespace pytraversal {

py::array_t<int64_t, py::array::c_style> traverse(
    const grid_type& grid, const grid_type::Vector3d& ray_origin,
    const grid_type::Vector3d& ray_end) {
  TraversedVoxels<double> traversed_voxels{};
  const auto ray = Ray<double>::fromOriginEnd(ray_origin, ray_end);
  const bool hit = traverseVoxelGrid(ray, grid, traversed_voxels);

  if (hit) {
    return py::cast(traversed_voxels);
  } else {
    // make sure to return empty array with correct shape [0, 3]
    const py::array::ShapeContainer shape({0, 3});
    return py::array_t<int64_t, py::array::c_style>(shape);
  }
}

}  // namespace pytraversal

PYBIND11_MODULE(pytraversal, m) {
  using namespace pytraversal;
  m.doc() = "pytraversal bindings";

  py::class_<grid_type>(m, "Grid3D")
      .def(py::init<>())
      .def(py::init<const grid_type::Vector3d&, const grid_type::Vector3d&,
                    const grid_type::Index3d&>())
      .def("traverse", traverse);
}
