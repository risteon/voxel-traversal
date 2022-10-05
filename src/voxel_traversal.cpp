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

  float_type tmin_temp{};
  float_type tmax_temp{};
  const bool ray_intersects_grid =
      rayBoxIntersection(ray, grid, tmin_temp, tmax_temp, t0, t1);
  if (!ray_intersects_grid) return false;

  tmin_temp = std::max(tmin_temp, t0);
  tmax_temp = std::min(tmax_temp, t1);
  const auto ray_start = ray.at(tmin_temp);
  const auto ray_end = ray.at(tmax_temp);

  // get voxel index of start and end within grid
  const auto voxelIndexStartUnlimited = grid.getIndex(ray_start);
  current_index =
      voxelIndexStartUnlimited.cwiseMax(0).cwiseMin(grid.numVoxels() - 1);

  const auto voxelIndexEndUnlimited = grid.getIndex(ray_end);
  final_index =
      voxelIndexEndUnlimited.cwiseMax(0).cwiseMin(grid.numVoxels() - 1);

  //  const auto tmax_xyz = grid.minBound() +
  const auto index_delta =
      (ray.direction().array() > 0.0).template cast<int_type>();
  const auto start_index = current_index + index_delta;
  const auto tmax_xyz =
      ((grid.minBound().array() +
        ((start_index.template cast<float_type>() * grid.voxelSize()) -
         ray_start.array())) /
       ray.direction().array()) +
      tmin_temp;

  tmax = (ray.direction().array() == 0.0).select(tmax_temp, tmax_xyz);
  const auto step_float = ray.direction().template array().sign().eval();
  step_index = step_float.template cast<int_type>().eval();
  delta = (step_index == 0)
              .select(tmax_temp,
                      grid.voxelSize() / ray.direction().array() * step_float);

  return true;
}

template <typename float_type>
bool traverseSingle(
    typename Grid3DSpatialDef<float_type>::Size3d& tmax,
    typename Grid3DSpatialDef<float_type>::Index3d& current_index,
    const typename Grid3DSpatialDef<float_type>::Index3d& overflow_index,
    const typename Grid3DSpatialDef<float_type>::Index3d& step_index,
    const typename Grid3DSpatialDef<float_type>::Size3d& delta) {
  if (tmax.x() < tmax.y() && tmax.x() < tmax.z()) {
    // X-axis traversal.
    current_index.x() += step_index.x();
    tmax.x() += delta.x();
    if (current_index.x() == overflow_index.x()) {
      return false;
    }
  } else if (tmax.y() < tmax.z()) {
    // Y-axis traversal.
    current_index.y() += step_index.y();
    tmax.y() += delta.y();
    if (current_index.y() == overflow_index.y()) {
      return false;
    }
  } else {
    // Z-axis traversal.
    current_index.z() += step_index.z();
    tmax.z() += delta.z();
    if (current_index.z() == overflow_index.z()) {
      return false;
    }
  }
  return true;
}
}  // namespace detail

// Uses the improved version of Smit's algorithm to determine if the given ray
// will intersect the grid between tmin and tmax. This version causes an
// additional efficiency penalty, but takes into account the negative zero case.
// tMin and tmax are then updated to incorporate the new intersection values.
// Returns true if the ray intersects the grid, and false otherwise.
// See: http://www.cs.utah.edu/~awilliam/box/box.pdf
template <typename float_type>
[[nodiscard]] bool rayBoxIntersection(const Ray<float_type>& ray,
                                      const Grid3DSpatialDef<float_type>& grid,
                                      float_type& tmin, float_type& tmax,
                                      float_type t0, float_type t1) noexcept {
  float_type tmin_temp{};
  float_type tmax_temp{};
  const auto inverse_direction = ray.direction().cwiseInverse();

  if (inverse_direction.x() >= 0) {
    tmin = (grid.minBound().x() - ray.origin().x()) * inverse_direction.x();
    tmax = (grid.maxBound().x() - ray.origin().x()) * inverse_direction.x();
  } else {
    tmin = (grid.maxBound().x() - ray.origin().x()) * inverse_direction.x();
    tmax = (grid.minBound().x() - ray.origin().x()) * inverse_direction.x();
  }

  if (inverse_direction.y() >= 0) {
    tmin_temp = (grid.minBound().y() - ray.origin().y()) * inverse_direction.y();
    tmax_temp = (grid.maxBound().y() - ray.origin().y()) * inverse_direction.y();
  } else {
    tmin_temp = (grid.maxBound().y() - ray.origin().y()) * inverse_direction.y();
    tmax_temp = (grid.minBound().y() - ray.origin().y()) * inverse_direction.y();
  }

  if (tmin > tmax_temp || tmin_temp > tmax) return false;
  if (tmin_temp > tmin) tmin = tmin_temp;
  if (tmax_temp < tmax) tmax = tmax_temp;

  if (inverse_direction.z() >= 0) {
    tmin_temp = (grid.minBound().z() - ray.origin().z()) * inverse_direction.z();
    tmax_temp = (grid.maxBound().z() - ray.origin().z()) * inverse_direction.z();
  } else {
    tmin_temp = (grid.maxBound().z() - ray.origin().z()) * inverse_direction.z();
    tmax_temp = (grid.minBound().z() - ray.origin().z()) * inverse_direction.z();
  }

  if (tmin > tmax_temp || tmin_temp > tmax) return false;
  if (tmin_temp > tmin) tmin = tmin_temp;
  if (tmax_temp < tmax) tmax = tmax_temp;
  return (tmin < t1 && tmax > t0);
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

  // one too far in every direction. Stop as soon as we hit any of these "walls"
  // It can happen, that the final_index is not exacty hit (float errors)
  const typename grid_type::Index3d overflow_index = final_index + step_index;

  while (true) {
    if (!detail::traverseSingle<float_type>(tmax, current_index, overflow_index,
                                            step_index, delta)) {
      break;
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

  // one too far in every direction. Stop as soon as we hit any of these "walls"
  // It can happen, that the final_index is not exacty hit (float errors)
  const typename grid_type::Index3d overflow_index = final_index + step_index;

  while (true) {
    if (!detail::traverseSingle<float_type>(tmax, current_index, overflow_index,
                                            step_index, delta)) {
      break;
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
template bool traverseVoxelGrid<long double>(
    const Ray<long double>&, const Grid3DSpatialDef<long double>&,
    TraversedVoxels<long double>&, long double, long double);

template bool traverseVoxelGrid<float>(const Ray<float>&,
                                       Grid3DTraversalCounter<float>&, float t0,
                                       float t1);
template bool traverseVoxelGrid<double>(const Ray<double>&,
                                        Grid3DTraversalCounter<double>&,
                                        double t0, double t1);
template bool traverseVoxelGrid<long double>(
    const Ray<long double>&, Grid3DTraversalCounter<long double>&,
    long double t0, long double t1);

}  // namespace voxel_traversal