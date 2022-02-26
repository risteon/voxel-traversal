#ifndef VOXEL_TRAVERSAL_GRID3DSPATIALDEF_H
#define VOXEL_TRAVERSAL_GRID3DSPATIALDEF_H

#include <Eigen/Dense>
#include <cstddef>

// TODO(risteon): make templates
// template<typename float_type = double>
class Grid3DSpatialDef {
 public:
  using float_type = double;
  using int_type = int32_t;

  using Vector3d = Eigen::Matrix<float_type, 3, 1>;
  using Size3d = Eigen::Array<float_type, 3, 1>;
  using Count3d = Eigen::Array<int_type, 3, 1>;

  Grid3DSpatialDef() = default;
  Grid3DSpatialDef& operator=(Grid3DSpatialDef other) {
    swap(*this, other);
    return *this;
  }
  Grid3DSpatialDef(Grid3DSpatialDef&& other) noexcept : Grid3DSpatialDef() {
    swap(*this, other);
  }

  Grid3DSpatialDef(const Vector3d& min_bound, const Vector3d& max_bound,
                   const Count3d& num_voxels)
      : min_bound_{min_bound},
        max_bound_{max_bound},
        grid_size_{max_bound - min_bound},
        num_voxels_{num_voxels},
        voxel_size_{grid_size_ / num_voxels.cast<float_type>()} {
    assert((num_voxels_ > 0).all());
    assert((min_bound_.array() < max_bound_.array()).all());
  }

  [[nodiscard]] const Count3d& numVoxels() const { return num_voxels_; }

  [[nodiscard]] const Vector3d& minBound() const { return min_bound_; }
  [[nodiscard]] const Vector3d&  maxBound() const { return max_bound_; }
  [[nodiscard]] const Size3d& gridSize() const { return grid_size_; }
  [[nodiscard]] const Size3d& voxelSize() const { return voxel_size_; }

  friend void swap(Grid3DSpatialDef& first, Grid3DSpatialDef& second) noexcept {
    // enable ADL
    using std::swap;
    swap(first.min_bound_, second.min_bound_);
    swap(first.max_bound_, second.max_bound_);
    swap(first.grid_size_, second.grid_size_);
    swap(first.num_voxels_, second.num_voxels_);
    swap(first.voxel_size_, second.voxel_size_);
  }

 private:
  // The minimum bound vector of the voxel grid.
  Vector3d min_bound_;
  // The maximum bound vector of the voxel grid.
  Vector3d max_bound_;
  // The grid size, determined by (max_bound_ - min_bound_).
  Size3d grid_size_;
  // The number of voxels in each of the x, y, z directions.
  Count3d num_voxels_;
  // The size of the voxel's x dimension.
  Size3d voxel_size_;
};

class Grid3DTraversalCounter : public Grid3DSpatialDef {};

#endif  // VOXEL_TRAVERSAL_GRID3DSPATIALDEF_H
