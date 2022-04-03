#ifndef VOXEL_TRAVERSAL_GRID_H
#define VOXEL_TRAVERSAL_GRID_H

#include <Eigen/Geometry>
#include <cstddef>
#include <unsupported/Eigen/CXX11/Tensor>

namespace voxel_traversal {

template <typename float_type = double>
class Grid3DSpatialDef {
 public:
  using int_type = int32_t;
  using float_t = float_type;

  using Vector3d = Eigen::Matrix<float_type, 3, 1>;
  using Size3d = Eigen::Array<float_type, 3, 1>;
  using Index3d = Eigen::Array<int_type, 3, 1>;

  Grid3DSpatialDef() = default;
  virtual ~Grid3DSpatialDef() = default;

  Grid3DSpatialDef& operator=(Grid3DSpatialDef other) {
    swap(*this, other);
    return *this;
  }
  Grid3DSpatialDef(Grid3DSpatialDef&& other) noexcept : Grid3DSpatialDef() {
    swap(*this, other);
  }

  Grid3DSpatialDef(const Vector3d& min_bound, const Vector3d& max_bound,
                   const Index3d& num_voxels)
      : min_bound_{min_bound},
        max_bound_{max_bound},
        grid_size_{max_bound - min_bound},
        num_voxels_{num_voxels},
        voxel_size_{grid_size_ / num_voxels.cast<float_type>()} {
    assert((num_voxels_ > 0).all());
    assert((min_bound_.array() < max_bound_.array()).all());
  }

  [[nodiscard]] const Index3d& numVoxels() const { return num_voxels_; }

  [[nodiscard]] const Vector3d& minBound() const { return min_bound_; }
  [[nodiscard]] const Vector3d& maxBound() const { return max_bound_; }
  [[nodiscard]] const Size3d& gridSize() const { return grid_size_; }
  [[nodiscard]] const Size3d& voxelSize() const { return voxel_size_; }

  friend void swap(Grid3DSpatialDef<float_type>& first,
                   Grid3DSpatialDef<float_type>& second) noexcept {
    // enable ADL
    using std::swap;
    swap(first.min_bound_, second.min_bound_);
    swap(first.max_bound_, second.max_bound_);
    swap(first.grid_size_, second.grid_size_);
    swap(first.num_voxels_, second.num_voxels_);
    swap(first.voxel_size_, second.voxel_size_);
  }

 protected:
  // The minimum bound vector of the voxel grid.
  Vector3d min_bound_;
  // The maximum bound vector of the voxel grid.
  Vector3d max_bound_;
  // The grid size, determined by (max_bound_ - min_bound_).
  Size3d grid_size_;
  // The number of voxels in each of the x, y, z directions.
  Index3d num_voxels_;
  // The size of the voxel's x dimension.
  Size3d voxel_size_;
};

template <typename float_type = double>
class Grid3DTraversalCounter : public Grid3DSpatialDef<float_type> {
 public:
  // declarations for independent base class types
  using typename Grid3DSpatialDef<float_type>::Vector3d;
  using typename Grid3DSpatialDef<float_type>::Index3d;

  using counter_type = uint64_t;
  using tensor_type = Eigen::Tensor<counter_type, 3>;

  Grid3DTraversalCounter() = default;
  virtual ~Grid3DTraversalCounter() override = default;

  Grid3DTraversalCounter& operator=(Grid3DTraversalCounter other) {
    swap(*this, other);
    return *this;
  }
  Grid3DTraversalCounter(Grid3DTraversalCounter&& other) noexcept
      : Grid3DTraversalCounter() {
    swap(*this, other);
  }

  Grid3DTraversalCounter(const Vector3d& min_bound, const Vector3d& max_bound,
                         const Index3d& num_voxels)
      : Grid3DSpatialDef<float_type>(min_bound, max_bound, num_voxels),
        counter_(tensor_type(num_voxels[0], num_voxels[1], num_voxels[2])) {
    assert((this->num_voxels_ > 0).all());
    assert((this->min_bound_.array() < this->max_bound_.array()).all());
    counter_.setZero();
  }

  [[nodiscard]] const tensor_type& getCounter() const { return counter_; }

  void clear() { counter_.setZero(); }
  void increaseAt(const Index3d& index) { counter_(index)++; }

  friend void swap(Grid3DTraversalCounter<float_type>& first,
                   Grid3DTraversalCounter<float_type>& second) noexcept {
    using std::swap;
    // TOOD(risteon) good option to call base class swap?
    swap(first.min_bound_, second.min_bound_);
    swap(first.max_bound_, second.max_bound_);
    swap(first.grid_size_, second.grid_size_);
    swap(first.num_voxels_, second.num_voxels_);
    swap(first.voxel_size_, second.voxel_size_);
    swap(first.counter_, second.counter_);
  }

 private:
  tensor_type counter_;
};

}  // namespace voxel_traversal
#endif  // VOXEL_TRAVERSAL_GRID_H
