#include "voxel_traversal.h"

#include <numeric>

namespace voxel_traversal {

namespace detail {
template <typename Grid>
bool setupTraversal(const Ray<typename Grid::float_t>& ray, const Grid& grid,
                    typename Grid::float_t t0, typename Grid::float_t t1,
                    typename Grid::Size3d& delta, typename Grid::Size3d& tmax,
                    typename Grid::Index3d& step_index,
                    typename Grid::Index3d& current_index,
                    typename Grid::Index3d& final_index) {
  using float_type = typename Grid::float_t;
  using int_type = typename Grid::int_type;

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
  const auto index_delta =
      (ray.direction().array() > 0.0).template cast<int_type>();
  const auto start_index = current_index + index_delta;
  const auto t_max_xyz =
      ((grid.minBound().array() +
        ((start_index.template cast<float_type>() * grid.voxelSize()) -
         ray_start.array())) /
       ray.direction().array()) +
      tMin;

  tmax = (ray.direction().array() == 0.0).select(tMax, t_max_xyz);
  const auto step_float = ray.direction().template array().sign().eval();
  step_index = step_float.template cast<int_type>().eval();
  delta = (step_index == 0)
              .select(tMax,
                      grid.voxelSize() / ray.direction().array() * step_float);

  return true;
}
}  // namespace detail

// Uses the improved version of Smit's algorithm to determine if the given ray
// will intersect the grid between tMin and tMax. This version causes an
// additional efficiency penalty, but takes into account the negative zero case.
// tMin and tMax are then updated to incorporate the new intersection values.
// Returns true if the ray intersects the grid, and false otherwise.
// See: http://www.cs.utah.edu/~awilliam/box/box.pdf
template <typename float_type>
[[nodiscard]] bool rayBoxIntersection(const Ray<float_type>& ray,
                                      const Grid3DSpatialDef<float_type>& grid,
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

template <typename float_type>
bool traverseVoxelGrid(const Ray<float_type>& ray,
                       Grid3DTraversalCounter<float_type>& grid, float_type t0,
                       float_type t1) noexcept {
  using grid_type = Grid3DTraversalCounter<float_type>;
  typename grid_type::Size3d delta{};
  typename grid_type::Size3d tmax{};
  typename grid_type::Index3d step_index{};
  typename grid_type::Index3d current_index{};
  typename grid_type::Index3d final_index{};

  const auto intersect = detail::setupTraversal(
      ray, grid, t0, t1, delta, tmax, step_index, current_index, final_index);

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

template <typename float_type>
bool traverseVoxelGrid(const Ray<float_type>& ray,
                       const Grid3DSpatialDef<float_type>& grid,
                       TraversedVoxels<float_type>& traversed_voxels,
                       float_type t0, float_type t1) noexcept {
  using grid_type = Grid3DSpatialDef<float_type>;
  typename grid_type::Size3d delta{};
  typename grid_type::Size3d tmax{};
  typename grid_type::Index3d step_index{};
  typename grid_type::Index3d current_index{};
  typename grid_type::Index3d final_index{};

  traversed_voxels.clear();

  const auto intersect = detail::setupTraversal<Grid3DSpatialDef<float_type>>(
      ray, grid, t0, t1, delta, tmax, step_index, current_index, final_index);
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

// instantiations
template bool traverseVoxelGrid<float>(const Ray<float>&,
                                       const Grid3DSpatialDef<float>&,
                                       TraversedVoxels<float>&, float, float);
template bool traverseVoxelGrid<double>(const Ray<double>&,
                                        const Grid3DSpatialDef<double>&,
                                        TraversedVoxels<double>&, double,
                                        double);
template bool traverseVoxelGrid<long double>(const Ray<long double>&,
                                        const Grid3DSpatialDef<long double>&,
                                        TraversedVoxels<long double>&, long double,
                                        long double);

template bool traverseVoxelGrid<float>(const Ray<float>&,
                                       Grid3DTraversalCounter<float>&, float t0,
                                       float t1);
template bool traverseVoxelGrid<double>(const Ray<double>&,
                                        Grid3DTraversalCounter<double>&,
                                        double t0, double t1);
template bool traverseVoxelGrid<long double>(const Ray<long double>&,
                                        Grid3DTraversalCounter<long double>&,
                                        long double t0, long double t1);

}  // namespace voxel_traversal