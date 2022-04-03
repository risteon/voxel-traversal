#include "voxel_traversal.h"

#include <numeric>

namespace voxel_traversal {

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

template <typename Grid>
bool setupTraversal(const Ray& ray, const Grid& grid, float_type t0,
                    float_type t1, typename Grid::Size3d& delta,
                    typename Grid::Size3d& tmax,
                    typename Grid::Index3d& step_index,
                    typename Grid::Index3d& current_index,
                    typename Grid::Index3d& final_index) {
  using int_type = Grid3DSpatialDef::int_type;

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
          .template cast<int_type>();
  current_index =
      voxelIndexStartUnlimited.cwiseMax(0).cwiseMin(grid.numVoxels() - 1);

  const auto voxelIndexEndUnlimited =
      ((ray_end - grid.minBound()).array() / grid.voxelSize().array())
          .floor()
          .template cast<int_type>();
  final_index =
      voxelIndexEndUnlimited.cwiseMax(0).cwiseMin(grid.numVoxels() - 1);

  //  const auto t_max_xyz = grid.minBound() +
  const auto index_delta = (ray.direction().array() > 0.0).cast<int_type>();
  const auto start_index = current_index + index_delta;
  const auto t_max_xyz =
      ((grid.minBound().array() +
        ((start_index.template cast<float_type>() * grid.voxelSize()) -
         ray_start.array())) /
       ray.direction().array()) +
      tMin;

  tmax = (ray.direction().array() == 0.0).select(tMax, t_max_xyz);
  const auto step_float = ray.direction().array().sign().eval();
  step_index = step_float.cast<int_type>().eval();
  delta = (step_index == 0)
              .select(tMax,
                      grid.voxelSize() / ray.direction().array() * step_float);

  return true;
}

bool traverseVoxelGrid(const Ray& ray, Grid3DTraversalCounter& grid,
                       float_type t0, float_type t1) noexcept {
  Grid3DTraversalCounter::Size3d delta{};
  Grid3DTraversalCounter::Size3d tmax{};
  Grid3DTraversalCounter::Index3d step_index{};
  Grid3DTraversalCounter::Index3d current_index{};
  Grid3DTraversalCounter::Index3d final_index{};

  const auto intersect = setupTraversal(ray, grid, t0, t1, delta, tmax,
                                        step_index, current_index, final_index);

  if (!intersect) return false;

  grid.increaseAt(current_index);

  while ((current_index != final_index).any()) {
    if (tmax.x() < tmax.y() && tmax.x() < tmax.z()) {
      // X-axis traversal.
      current_index.x() += step_index.x();
      tmax.x() += delta.x();
    } else if (tmax.y() < tmax.z()) {
      // Y-axis traversal.
      current_index.y() += step_index.y();
      tmax.y() += delta.y();
    } else {
      // Z-axis traversal.
      current_index.z() += step_index.z();
      tmax.z() += delta.z();
    }
    grid.increaseAt(current_index);
  }

  return true;
}

bool traverseVoxelGrid(const Ray& ray, const Grid3DSpatialDef& grid,
                       std::vector<Grid3DSpatialDef::Index3d>& traversed_voxels,
                       float_type t0, float_type t1) noexcept {
  Grid3DSpatialDef::Size3d delta{};
  Grid3DSpatialDef::Size3d tmax{};
  Grid3DSpatialDef::Index3d step_index{};
  Grid3DSpatialDef::Index3d current_index{};
  Grid3DSpatialDef::Index3d final_index{};

  traversed_voxels.clear();

  const auto intersect = setupTraversal(ray, grid, t0, t1, delta, tmax,
                                        step_index, current_index, final_index);
  if (!intersect) return false;

  traversed_voxels.push_back(current_index);

  while ((current_index != final_index).any()) {
    if (tmax.x() < tmax.y() && tmax.x() < tmax.z()) {
      // X-axis traversal.
      current_index.x() += step_index.x();
      tmax.x() += delta.x();
    } else if (tmax.y() < tmax.z()) {
      // Y-axis traversal.
      current_index.y() += step_index.y();
      tmax.y() += delta.y();
    } else {
      // Z-axis traversal.
      current_index.z() += step_index.z();
      tmax.z() += delta.z();
    }

    traversed_voxels.push_back(current_index);
  }

  return true;
}

}  // namespace voxel_traversal