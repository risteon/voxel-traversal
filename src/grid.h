#ifndef VOXEL_TRAVERSAL_GRID_H
#define VOXEL_TRAVERSAL_GRID_H

#include <Eigen/Geometry>
#include <cstddef>
#include <unsupported/Eigen/CXX11/Tensor>

namespace voxel_traversal {

template <typename float_type = double>
class Grid3DSpatialDef {
 public:
  // need to express negative voxel indices when calculating voxel traversal
  using int_type = int32_t;
  using float_t = float_type;

  using Vector3d = Eigen::Matrix<float_type, 3, 1>;
  using Size3d = Eigen::Array<float_type, 3, 1>;
  using Index3d = Eigen::Array<int_type, 3, 1>;

  Grid3DSpatialDef() = default;
  Grid3DSpatialDef(const Vector3d& min_bound, const Vector3d& max_bound,
                   const Index3d& num_voxels)
      : min_bound_{min_bound},
        max_bound_{max_bound},
        grid_size_{max_bound - min_bound},
        num_voxels_{num_voxels},
        voxel_size_{grid_size_ / num_voxels.cast<float_type>()},
        voxel_count_{static_cast<decltype(voxel_count_)>(num_voxels.prod())} {
    assert((num_voxels_ > 0).all());
    assert((min_bound_.array() < max_bound_.array()).all());
  }
  virtual ~Grid3DSpatialDef() = default;

  Grid3DSpatialDef(const Grid3DSpatialDef& other) = default;
  Grid3DSpatialDef(Grid3DSpatialDef&& other) noexcept : Grid3DSpatialDef() {
    swap(*this, other);
  }
  //! Copy-and-swap. Handles both lvalues and rvalues
  Grid3DSpatialDef& operator=(Grid3DSpatialDef other) noexcept {
    swap(*this, other);
    return *this;
  }

  [[nodiscard]] const Index3d& numVoxels() const { return num_voxels_; }

  [[nodiscard]] const Vector3d& minBound() const { return min_bound_; }
  [[nodiscard]] const Vector3d& maxBound() const { return max_bound_; }
  [[nodiscard]] const Size3d& gridSize() const { return grid_size_; }
  [[nodiscard]] const Size3d& voxelSize() const { return voxel_size_; }
  //! Total number of voxels in this grid
  [[nodiscard]] uint64_t voxelCount() const noexcept { return voxel_count_; }

  //! Maximum number of voxels that can be traversed by a single ray
  [[nodiscard]] int_type upperBoundVoxelTraversal() const {
    return num_voxels_.sum();
  }

  //! Get voxel index of given position
  [[nodiscard]] Index3d getIndex(const Vector3d& v) const noexcept {
    return ((v - min_bound_).array() / voxel_size_)
        .floor()
        .template cast<int_type>();
  }

  //! Check if given index falls into the voxel grid
  [[nodiscard]] bool isWithinGrid(const Index3d& index) const noexcept {
    return (index >= 0).all() && (index < num_voxels_).all();
  }

  //! Check if given points falls into the voxel grid
  [[nodiscard]] bool isWithinGrid(const Vector3d& v) const noexcept {
    return (v.array() >= min_bound_.array()).all() &&
           (v.array() <= max_bound_.array()).all();
  }

  //! For copy-and-swap
  friend void swap(Grid3DSpatialDef<float_type>& first,
                   Grid3DSpatialDef<float_type>& second) noexcept {
    // enable ADL
    using std::swap;
    swap(first.min_bound_, second.min_bound_);
    swap(first.max_bound_, second.max_bound_);
    swap(first.grid_size_, second.grid_size_);
    swap(first.num_voxels_, second.num_voxels_);
    swap(first.voxel_size_, second.voxel_size_);
    swap(first.voxel_count_, second.voxel_count_);
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
  // Number of voxels in this grid
  uint64_t voxel_count_;
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

  Grid3DTraversalCounter(const Grid3DTraversalCounter& other) = default;

  //! Copy-and-swap. Handles both lvalues and rvalues
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
    swap(first.voxel_count_, second.voxel_count_);
    swap(first.counter_, second.counter_);
  }

 private:
  tensor_type counter_;
};

}  // namespace voxel_traversal
#endif  // VOXEL_TRAVERSAL_GRID_H
