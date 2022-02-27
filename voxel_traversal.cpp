#include "voxel_traversal.h"

#include <numeric>

namespace algorithm {

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
  using int_type = Grid3DSpatialDef::int_type;
  traversed_voxels.clear();

  float_type tMin{};
  float_type tMax{};
  const bool ray_intersects_grid =
      rayBoxIntersection(ray, grid, tMin, tMax, t0, t1);
  if (!ray_intersects_grid) return false;

  tMin = std::max(tMin, t0);
  tMax = std::min(tMax, t1);
  const auto ray_start = ray.at(tMin);
  const auto ray_end = ray.at(tMax);

  // get voxel index of start and end within grid
  const auto voxelIndexStartUnlimited =
      ((ray_start - grid.minBound()).array() / grid.voxelSize().array())
          .floor()
          .cast<int_type>();
  const auto voxelIndexStart =
      voxelIndexStartUnlimited.cwiseMax(0).cwiseMin(grid.numVoxels() - 1);

  const auto voxelIndexEndUnlimited =
      ((ray_end - grid.minBound()).array() / grid.voxelSize().array())
          .floor()
          .cast<int_type>();
  const auto voxelIndexEnd =
      voxelIndexEndUnlimited.cwiseMax(0).cwiseMin(grid.numVoxels() - 1);

  size_t current_X_index = voxelIndexStart.x();
  size_t current_Y_index = voxelIndexStart.y();
  size_t current_Z_index = voxelIndexStart.z();
  size_t end_X_index = voxelIndexEnd.x();
  size_t end_Y_index = voxelIndexEnd.y();
  size_t end_Z_index = voxelIndexEnd.z();

  //  const auto t_max_xyz = grid.minBound() +
  const auto index_delta = (ray.direction().array() > 0.0).cast<int_type>();
  const auto start_index = current_X_index + index_delta;
  const auto t_max_xyz = ((grid.minBound().array() +
                          ((start_index.cast<float_type>() * grid.voxelSize()) -
                           ray_start.array())) /
                              ray.direction().array()) + tMin;

  auto tmax = (ray.direction().array() == 0.0).select(tMax, t_max_xyz).eval();

  int_type stepX;
  float_type tDeltaX;
  float_type tMaxX;
  if (ray.direction().x() > 0.0) {
    stepX = 1;
    tDeltaX = grid.voxelSize().x() / ray.direction().x();
    tMaxX =
        tMin + (grid.minBound().x() +
                (current_X_index + 1) * grid.voxelSize().x() - ray_start.x()) /
                   ray.direction().x();
  } else if (ray.direction().x() < 0.0) {
    stepX = -1;
    tDeltaX = grid.voxelSize().x() / -ray.direction().x();
    tMaxX = tMin + (grid.minBound().x() +
                    current_X_index * grid.voxelSize().x() - ray_start.x()) /
                       ray.direction().x();
  } else {
    stepX = 0;
    tDeltaX = tMax;
    tMaxX = tMax;
  }

  int_type stepY;
  float_type tDeltaY;
  float_type tMaxY;
  if (ray.direction().y() > 0.0) {
    stepY = 1;
    tDeltaY = grid.voxelSize().y() / ray.direction().y();
    tMaxY =
        tMin + (grid.minBound().y() +
                (current_Y_index + 1) * grid.voxelSize().y() - ray_start.y()) /
                   ray.direction().y();
  } else if (ray.direction().y() < 0.0) {
    stepY = -1;
    tDeltaY = grid.voxelSize().y() / -ray.direction().y();
    tMaxY = tMin + (grid.minBound().y() +
                    current_Y_index * grid.voxelSize().y() - ray_start.y()) /
                       ray.direction().y();
  } else {
    stepY = 0;
    tDeltaY = tMax;
    tMaxY = tMax;
  }

  int_type stepZ;
  float_type tDeltaZ;
  float_type tMaxZ;
  if (ray.direction().z() > 0.0) {
    stepZ = 1;
    tDeltaZ = grid.voxelSize().z() / ray.direction().z();
    tMaxZ =
        tMin + (grid.minBound().z() +
                (current_Z_index + 1) * grid.voxelSize().z() - ray_start.z()) /
                   ray.direction().z();
  } else if (ray.direction().z() < 0.0) {
    stepZ = -1;
    tDeltaZ = grid.voxelSize().z() / -ray.direction().z();
    tMaxZ = tMin + (grid.minBound().z() +
                    current_Z_index * grid.voxelSize().z() - ray_start.z()) /
                       ray.direction().z();
  } else {
    stepZ = 0;
    tDeltaZ = tMax;
    tMaxZ = tMax;
  }

  traversed_voxels.emplace_back(current_X_index, current_Y_index,
                                current_Z_index);

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

    traversed_voxels.emplace_back(current_X_index, current_Y_index,
                                  current_Z_index);
  }

  return true;
}

}  // namespace algorithm