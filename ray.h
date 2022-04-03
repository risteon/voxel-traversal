#ifndef VOXEL_TRAVERSAL_RAY_H
#define VOXEL_TRAVERSAL_RAY_H

#include <Eigen/Dense>

namespace voxel_traversal {
// TODO(risteon): make templates
// template<typename float_type = double>
class Ray {
 public:
  using float_type = double;
  using Vector3d = Eigen::Matrix<float_type, 3, 1>;

  static Ray fromOriginDir(const Vector3d& origin, const Vector3d& dir) {
    return Ray(origin, origin + dir, dir);
  }
  static Ray fromOriginEnd(const Vector3d& start, const Vector3d& end) {
    return Ray(start, end, end - start);
  }

  // Represents the function p(t) = origin + t * direction,
  // where p is a 3-dimensional position, and t is a scalar.
  [[nodiscard]] Vector3d at(float_type t) const {
    return (origin_ * (float_type{1.0} - t)) + (end_point_ * t);
  }

  [[nodiscard]] const Vector3d& origin() const { return origin_; }
  [[nodiscard]] const Vector3d& direction() const { return direction_; }

 private:
  explicit Ray(Vector3d origin, Vector3d end_point, Vector3d direction)
      : origin_{std::move(origin)},
        end_point_{std::move(end_point)},
        direction_{std::move(direction)} {}

  // origin and end point of the ray
  Vector3d origin_;
  Vector3d end_point_;
  // end - origin
  Vector3d direction_;
};

}  // namespace voxel_traversal

#endif  // VOXEL_TRAVERSAL_RAY_H
