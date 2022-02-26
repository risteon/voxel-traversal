#ifndef VOXEL_TRAVERSAL_AMANATIDESWOOALGORITHM_H
#define VOXEL_TRAVERSAL_AMANATIDESWOOALGORITHM_H

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

using value_type = double;

void amanatidesWooAlgorithm(const Ray& ray, const Grid3DSpatialDef& grid,
                            value_type t0, value_type t1) noexcept;

[[nodiscard]] bool rayBoxIntersection(const Ray& ray,
                                      const Grid3DSpatialDef& grid,
                                      value_type& tMin, value_type& tMax,
                                      value_type t0 = 0.0, value_type t1 = 1.0) noexcept;
}  // namespace algorithm

#endif  // VOXEL_TRAVERSAL_AMANATIDESWOOALGORITHM_H
