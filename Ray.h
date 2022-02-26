#ifndef VOXEL_TRAVERSAL_RAY_H
#define VOXEL_TRAVERSAL_RAY_H

#include <Eigen/Dense>

// Encapsulates the functionality of a ray.
// This consists of two components, the origin of the ray,
// and the direction of the ray.

// TODO(risteon): make templates
// template<typename float_type = double>
class Ray {
 public:
  using float_type = double;
  using Vector3d = Eigen::Matrix<float_type, 3, 1>;

  explicit Ray(Vector3d origin, Vector3d direction)
      : origin_{std::move(origin)}, direction_{std::move(direction)} {}

  static Ray fromEndpoints(const Vector3d& start, const Vector3d& end) {
    return Ray(start, end - start);
  }

  // Represents the function p(t) = origin + t * direction,
  // where p is a 3-dimensional position, and t is a scalar.
  [[nodiscard]] Vector3d point_at_parameter(float_type t) const {
    return origin_ + (direction_ * t);
  }

  [[nodiscard]] const Vector3d& origin() const { return origin_; }
  [[nodiscard]] const Vector3d& direction() const { return direction_; }

 private:
  // The origin of the ray.
  Vector3d origin_;
  // The normalized direction of the ray.
  Vector3d direction_;
};

#endif  // VOXEL_TRAVERSAL_RAY_H
