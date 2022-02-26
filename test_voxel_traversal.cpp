//
//
#include <gtest/gtest.h>

#include "voxel_traversal.h"
#include <limits>

using namespace algorithm;

// Small voxel grid 2x2x2
class TestVoxel2x2x2Intersect : public ::testing::Test {
 protected:
  using float_type = Grid3DSpatialDef::float_type;
  using V3 = Grid3DSpatialDef::Vector3d;
  using C3 = Grid3DSpatialDef::Count3d;

  static const constexpr Grid3DSpatialDef::float_type SQRT2 = sqrt(2.0);
  static const constexpr Grid3DSpatialDef::float_type HALF_SQRT2 = 0.5 * SQRT2;
  static const constexpr Grid3DSpatialDef::float_type SQRT3 = sqrt(3.0);

  void SetUp() override {
    const V3 bound_min(0.0, 0.0, 0.0);
    const V3 bound_max(2.0, 2.0, 2.0);
    const C3 voxel_count(2, 2, 2);
    grid_ = Grid3DSpatialDef(bound_min, bound_max, voxel_count);
  }

  Grid3DSpatialDef grid_;
};

TEST_F(TestVoxel2x2x2Intersect, XYPlaneFixedT) {
  float_type t_min, t_max;
  {
    const Ray ray({.5, .5, .5}, {1.0, 0., 0.});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-0.5, t_min);
    EXPECT_FLOAT_EQ(1.5, t_max);
  }
  {
    const Ray ray(V3(-.5, -.5, -.5), V3(1.0, 0., 0.));
    const auto intersect =
        rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_FALSE(intersect);
  }
  {
    const Ray ray(V3(.5, -.5, .5), V3(1.0, 0., 0.));
    const auto intersect =
        rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_FALSE(intersect);
  }
  {
    // miss corner
    const Ray ray(V3(.5, -.55, 1.0), V3(-1.0, 1.0, 0.));
    const auto intersect =
        rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_FALSE(intersect);
  }
  {
    // Non-unit direction vector
    const Ray ray(V3(.5, -.4, 1.0), V3(-1.0, 1.0, 0.));
    const auto intersect =
        rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_TRUE(intersect);
  }
  {
    // Unit direction vector
    const Ray ray({1.0, 1.0, 1.0}, {-HALF_SQRT2, -HALF_SQRT2, 0.});
    const auto intersect =
        rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-SQRT2, t_min);
    EXPECT_FLOAT_EQ(SQRT2, t_max);
  }
  {
    // Non-unit direction vector
    const Ray ray({1.0, 1.0, 1.0}, {-1.0, -1.0, 0.});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-1.0, t_min);
    EXPECT_FLOAT_EQ(1.0, t_max);
  }
  {
    // positive direction
    const Ray ray({1.0, 1.0, 1.0}, {1.0, 1.0, 0.});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-1.0, t_min);
    EXPECT_FLOAT_EQ(1.0, t_max);
  }
  {
    // too short, stops before grid
    const Ray ray({-1.0, 3.0, 1.0}, {1.1, -0.9, 0.});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_FALSE(intersect);
    EXPECT_GE(t_min, 1.0);
    EXPECT_GE(t_max, 1.0);
  }
  {
    // long enough, reaches into grid and ends within
    const Ray ray({-1.0, 3.0, 1.0}, {1.1, -1.1, 0.});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    EXPECT_LE(t_min, 1.0);
    EXPECT_GE(t_max, 1.0);
  }
  {
    // parallel to grid, does never intersect
    const auto ray = Ray::fromEndpoints({-0.5, 0.5, 0.5}, {-0.5, 1.5, 0.5});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_FALSE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_max));
  }
  {
    // parallel to grid, does intersect
    const auto ray = Ray::fromEndpoints({-0.5, 0.5, 0.5}, {0.5, 0.5, 0.5});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(0.5, t_min);
    EXPECT_FLOAT_EQ(2.5, t_max);
  }
}

TEST_F(TestVoxel2x2x2Intersect, RayOutsideXYPlane) {
  double t_min, t_max;
  {
    //
    const auto ray = Ray::fromEndpoints({-0.5, -0.5, -0.5}, {2.0, 2.0, 2.0});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_GE(t_min, 0.1);
    EXPECT_FLOAT_EQ(1.0, t_max);
  }
  {
    // outside of grid
    const auto ray = Ray::fromEndpoints({2.1, 2.1, 2.1}, {2.2, 2.2, 2.2});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_FALSE(intersect);
    EXPECT_FLOAT_EQ(t_min, -21.0);
    EXPECT_FLOAT_EQ(t_max, -1.0);
  }
}

TEST_F(TestVoxel2x2x2Intersect, DegenerateCase) {
  float_type t_min, t_max;
  {
    // zero direction within grid -> does intersect with single point
    const Ray ray({.5, .5, .5}, V3::Zero());
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_max));
  }
  {
    // zero direction within grid -> does intersect with single point
    // different position near corner
    const Ray ray({1.99, 1.99, 0.01}, V3::Zero());
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_max));
  }
  {
    // zero direction outside of grid -> does not intersect
    const Ray ray({.5, .5, -0.1}, V3::Zero());
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_FALSE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_max));
  }
  {
    // zero direction outside of grid -> does not intersect
    // different position
    const Ray ray({.5, -0.01, .5}, V3::Zero());
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_FALSE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<float_type>::infinity(), std::fabs(t_max));
  }
}
