#include "voxel_traversal.h"

#include <numeric>

namespace algorithm {

// Macro defined to avoid unnecessary checks with NaNs when using std::max
#define MAX(a, b) ((a > b ? a : b))
#define MIN(a, b) ((a < b ? a : b))

// Uses the improved version of Smit's algorithm to determine if the given ray
// will intersect the grid between tMin and tMax. This version causes an
// additional efficiency penalty, but takes into account the negative zero case.
// tMin and tMax are then updated to incorporate the new intersection values.
// Returns true if the ray intersects the grid, and false otherwise.
// See: http://www.cs.utah.edu/~awilliam/box/box.pdf
[[nodiscard]] bool rayBoxIntersection(const Ray& ray,
                                      const Grid3DSpatialDef& grid,
                                      float_type& tMin, float_type& tMax,
                                      float_type t0, float_type t1) noexcept {
  float_type tYMin, tYMax, tZMin, tZMax;
  const auto inverse_direction = ray.direction().cwiseInverse();

  if (inverse_direction.x() >= 0) {
    tMin = (grid.minBound().x() - ray.origin().x()) * inverse_direction.x();
    tMax = (grid.maxBound().x() - ray.origin().x()) * inverse_direction.x();
  } else {
    tMin = (grid.maxBound().x() - ray.origin().x()) * inverse_direction.x();
    tMax = (grid.minBound().x() - ray.origin().x()) * inverse_direction.x();
  }

  if (inverse_direction.y() >= 0) {
    tYMin = (grid.minBound().y() - ray.origin().y()) * inverse_direction.y();
    tYMax = (grid.maxBound().y() - ray.origin().y()) * inverse_direction.y();
  } else {
    tYMin = (grid.maxBound().y() - ray.origin().y()) * inverse_direction.y();
    tYMax = (grid.minBound().y() - ray.origin().y()) * inverse_direction.y();
  }

  if (tMin > tYMax || tYMin > tMax) return false;
  if (tYMin > tMin) tMin = tYMin;
  if (tYMax < tMax) tMax = tYMax;

  if (inverse_direction.z() >= 0) {
    tZMin = (grid.minBound().z() - ray.origin().z()) * inverse_direction.z();
    tZMax = (grid.maxBound().z() - ray.origin().z()) * inverse_direction.z();
  } else {
    tZMin = (grid.maxBound().z() - ray.origin().z()) * inverse_direction.z();
    tZMax = (grid.minBound().z() - ray.origin().z()) * inverse_direction.z();
  }

  if (tMin > tZMax || tZMin > tMax) return false;
  if (tZMin > tMin) tMin = tZMin;
  if (tZMax < tMax) tMax = tZMax;
  return (tMin < t1 && tMax > t0);
}

bool traverseVoxelGrid(const Ray& ray, const Grid3DSpatialDef& grid,
                       std::vector<Grid3DSpatialDef::Index3d>& traversed_voxels,
                       float_type t0, float_type t1) noexcept {
  traversed_voxels.clear();

  float_type tMin{};
  float_type tMax{};
  const bool ray_intersects_grid =
      rayBoxIntersection(ray, grid, tMin, tMax, t0, t1);
  if (!ray_intersects_grid) return false;

  tMin = MAX(tMin, t0);
  tMax = MIN(tMax, t1);
  const auto ray_start = ray.at(tMin);
  const auto ray_end = ray.at(tMax);

  // get voxel index of start and end within grid

  size_t current_X_index = MAX(
      1, std::ceil((ray_start.x() - grid.minBound().x()) / grid.voxelSize().x()));
  current_X_index = MIN(current_X_index, grid.numVoxels().x());

  size_t end_X_index = MAX(
      1, std::ceil((ray_end.x() - grid.minBound().x()) / grid.voxelSize().x()));
  end_X_index = MIN(end_X_index, grid.numVoxels().x());

  int stepX;
  float_type tDeltaX;
  float_type tMaxX;
  if (ray.direction().x() > 0.0) {
    stepX = 1;
    tDeltaX = grid.voxelSize().x() / ray.direction().x();
    tMaxX = tMin + (grid.minBound().x() +
                    current_X_index * grid.voxelSize().x() - ray_start.x()) /
                       ray.direction().x();
  } else if (ray.direction().x() < 0.0) {
    stepX = -1;
    tDeltaX = grid.voxelSize().x() / -ray.direction().x();
    const size_t previous_X_index = current_X_index - 1;
    tMaxX = tMin + (grid.minBound().x() +
                    previous_X_index * grid.voxelSize().x() - ray_start.x()) /
                       ray.direction().x();
  } else {
    stepX = 0;
    tDeltaX = tMax;
    tMaxX = tMax;
  }

  size_t current_Y_index = MAX(
      1, std::ceil((ray_start.y() - grid.minBound().y()) / grid.voxelSize().y()));
  current_Y_index = MIN(current_Y_index, grid.numVoxels().y());

  size_t end_Y_index = MAX(
      1, std::ceil((ray_end.y() - grid.minBound().y()) / grid.voxelSize().y()));
  end_Y_index = MIN(end_Y_index, grid.numVoxels().y());

  int stepY;
  float_type tDeltaY;
  float_type tMaxY;
  if (ray.direction().y() > 0.0) {
    stepY = 1;
    tDeltaY = grid.voxelSize().y() / ray.direction().y();
    tMaxY = tMin + (grid.minBound().y() +
                    current_Y_index * grid.voxelSize().y() - ray_start.y()) /
                       ray.direction().y();
  } else if (ray.direction().y() < 0.0) {
    stepY = -1;
    tDeltaY = grid.voxelSize().y() / -ray.direction().y();
    const size_t previous_Y_index = current_Y_index - 1;
    tMaxY = tMin + (grid.minBound().y() +
                    previous_Y_index * grid.voxelSize().y() - ray_start.y()) /
                       ray.direction().y();
  } else {
    stepY = 0;
    tDeltaY = tMax;
    tMaxY = tMax;
  }

  size_t current_Z_index = MAX(
      1, std::ceil((ray_start.z() - grid.minBound().z()) / grid.voxelSize().z()));
  current_Z_index = MIN(current_Z_index, grid.numVoxels().z());

  size_t end_Z_index = MAX(
      1, std::ceil((ray_end.z() - grid.minBound().z()) / grid.voxelSize().z()));
  end_Z_index = MIN(end_Z_index, grid.numVoxels().z());

  int stepZ;
  float_type tDeltaZ;
  float_type tMaxZ;
  if (ray.direction().z() > 0.0) {
    stepZ = 1;
    tDeltaZ = grid.voxelSize().z() / ray.direction().z();
    tMaxZ = tMin + (grid.minBound().z() +
                    current_Z_index * grid.voxelSize().z() - ray_start.z()) /
                       ray.direction().z();
  } else if (ray.direction().z() < 0.0) {
    stepZ = -1;
    tDeltaZ = grid.voxelSize().z() / -ray.direction().z();
    const size_t previous_Z_index = current_Z_index - 1;
    tMaxZ = tMin + (grid.minBound().z() +
                    previous_Z_index * grid.voxelSize().z() - ray_start.z()) /
                       ray.direction().z();
  } else {
    stepZ = 0;
    tDeltaZ = tMax;
    tMaxZ = tMax;
  }

  traversed_voxels.emplace_back(current_X_index - 1, current_Y_index - 1,
                                current_Z_index - 1);

  while (current_X_index != end_X_index || current_Y_index != end_Y_index ||
         current_Z_index != end_Z_index) {
    if (tMaxX < tMaxY && tMaxX < tMaxZ) {
      // X-axis traversal.
      current_X_index += stepX;
      tMaxX += tDeltaX;
    } else if (tMaxY < tMaxZ) {
      // Y-axis traversal.
      current_Y_index += stepY;
      tMaxY += tDeltaY;
    } else {
      // Z-axis traversal.
      current_Z_index += stepZ;
      tMaxZ += tDeltaZ;
    }

    traversed_voxels.emplace_back(current_X_index - 1, current_Y_index - 1,
                                  current_Z_index - 1);
  }

  return true;
}

}  // namespace algorithm