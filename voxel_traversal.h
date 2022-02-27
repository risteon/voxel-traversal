#ifndef VOXEL_TRAVERSAL_H
#define VOXEL_TRAVERSAL_H

#include "Grid3DSpatialDef.h"
#include "Ray.h"

namespace algorithm {
// Implements the algorithm presented in Amanatides & Woo's "A Fast Voxel
// Traversal Algorithm for Ray Tracing." See:
// https://www.researchgate.net/publication/2611491_A_Fast_Voxel_Traversal_Algorithm_for_Ray_Tracing
// If the ray origin is outside the voxel grid, uses a safer version of Smit's
// ray box intersection algorithm to determine intersection. The bounds [t0, t1]
// determine the begin and end parameter for which the ray travels. The
// algorithm occurs in two phases, initialization and traversal. Requires:
//     t1 > t0
//     0.0 <= t0 <= 1.0
//     0.0 <= t1 <= 1.0
//     To encapsulate entire ray traversal, set t0 = 0.0, t1 = 1.0
//     'grid' encapsulates a valid voxel grid system.
//
// Notes:
//     Assumes that indices for voxel coordinates begin at 1.

using float_type = double;

bool traverseVoxelGrid(
    const Ray& ray, const Grid3DSpatialDef& grid,
    std::vector<Grid3DSpatialDef::Index3d>& traversed_voxels,
                       float_type t0 = float_type{0.0},
                       float_type t1 = float_type{1.0}) noexcept;

[[nodiscard]] bool rayBoxIntersection(const Ray& ray,
                                      const Grid3DSpatialDef& grid,
                                      float_type& tMin, float_type& tMax,
                                      float_type t0 = 0.0,
                                      float_type t1 = 1.0) noexcept;
}  // namespace algorithm

#endif  // VOXEL_TRAVERSAL_H
