#ifndef VOXEL_TRAVERSAL_H
#define VOXEL_TRAVERSAL_H

#include "grid.h"
#include "ray.h"

namespace voxel_traversal {
// Implements the algorithm presented in Amanatides & Woo's "A Fast Voxel
// Traversal Algorithm for Ray Tracing." See:
// https://www.researchgate.net/publication/2611491_A_Fast_Voxel_Traversal_Algorithm_for_Ray_Tracing
// If the ray origin is outside the voxel grid, uses a safer version of Smit's
// ray box intersection algorithm to determine intersection. The bounds [t0, t1]
// determine the begin and end parameter for which the ray travels. The
// voxel_traversal occurs in two phases, initialization and traversal. Requires:
//     t1 > t0
//     To encapsulate entire ray traversal, set t0 = 0.0, t1 = 1.0
//     'grid' encapsulates a valid voxel grid system.
template <typename float_type>
using TraversedVoxels =
    std::vector<typename Grid3DSpatialDef<float_type>::Index3d>;

template <typename float_type = double>
bool traverseVoxelGrid(const Ray<float_type>& ray,
                       const Grid3DSpatialDef<float_type>& grid,
                       TraversedVoxels<float_type>& traversed_voxels,
                       float_type t0 = float_type{0.0},
                       float_type t1 = float_type{1.0}) noexcept;

template <typename float_type = double>
bool traverseVoxelGrid(const Ray<float_type>& ray,
                       Grid3DTraversalCounter<float_type>& grid,
                       float_type t0 = float_type{0.0},
                       float_type t1 = float_type{1.0}) noexcept;

template <typename float_type = double>
[[nodiscard]] bool rayBoxIntersection(const Ray<float_type>& ray,
                                      const Grid3DSpatialDef<float_type>& grid,
                                      float_type& tMin, float_type& tMax,
                                      float_type t0 = 0.0,
                                      float_type t1 = 1.0) noexcept;
}  // namespace voxel_traversal

#endif  // VOXEL_TRAVERSAL_H
