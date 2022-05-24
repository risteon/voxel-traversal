//
//
#include <gtest/gtest.h>

#include <limits>

#include "voxel_traversal.h"

using namespace voxel_traversal;

// define acceptable eps for t_min, t_max calculation
template <typename float_type>
class Epsilon {
 public:
  static constexpr float_type eps = 0.0;
};
template <>
class Epsilon<float> {
 public:
  static constexpr float eps = 5e-5;
};
template <>
class Epsilon<double> {
 public:
  static constexpr double eps = 1e-11;
};
template <>
class Epsilon<long double> {
 public:
  static constexpr double eps = 2e-14;
};

// Small voxel grid 2x2x2
template <typename float_type>
class TestVoxel2x2x2Intersect : public ::testing::Test {
 protected:
  using V3 = typename Grid3DSpatialDef<float_type>::Vector3d;
  using C3 = typename Grid3DSpatialDef<float_type>::Index3d;
  using R = Ray<float_type>;

  static const constexpr float_type SQRT2 = sqrt(2.0);
  static const constexpr float_type HALF_SQRT2 = 0.5 * SQRT2;
  static const constexpr float_type SQRT3 = sqrt(3.0);

  void SetUp() override {
    const V3 bound_min(0.0, 0.0, 0.0);
    const V3 bound_max(2.0, 2.0, 2.0);
    const C3 voxel_count(2, 2, 2);
    grid_ = Grid3DSpatialDef(bound_min, bound_max, voxel_count);
  }

  Grid3DSpatialDef<float_type> grid_;
};

using Implementations = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(TestVoxel2x2x2Intersect, Implementations);

TYPED_TEST(TestVoxel2x2x2Intersect, XYPlaneFixedT) {
  TypeParam t_min, t_max;
  {
    const auto ray = TestFixture::R::fromOriginDir({.5, .5, .5}, {1.0, 0., 0.});
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-0.5, t_min);
    EXPECT_FLOAT_EQ(1.5, t_max);
  }
  {
    const auto ray =
        TestFixture::R::fromOriginDir(typename TestFixture::V3(-.5, -.5, -.5),
                                      typename TestFixture::V3(1.0, 0., 0.));
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_FALSE(intersect);
  }
  {
    const auto ray =
        TestFixture::R::fromOriginDir(typename TestFixture::V3(.5, -.5, .5),
                                      typename TestFixture::V3(1.0, 0., 0.));
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_FALSE(intersect);
  }
  {
    // miss corner
    const auto ray =
        TestFixture::R::fromOriginDir(typename TestFixture::V3(.5, -.55, 1.0),
                                      typename TestFixture::V3(-1.0, 1.0, 0.));
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_FALSE(intersect);
  }
  {
    // Non-unit direction vector
    const auto ray =
        TestFixture::R::fromOriginDir(typename TestFixture::V3(.5, -.4, 1.0),
                                      typename TestFixture::V3(-1.0, 1.0, 0.));
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
  }
  {
    // Unit direction vector
    const auto ray = TestFixture::R::fromOriginDir(
        {1.0, 1.0, 1.0},
        {-TestFixture::HALF_SQRT2, -TestFixture::HALF_SQRT2, 0.});
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-TestFixture::SQRT2, t_min);
    EXPECT_FLOAT_EQ(TestFixture::SQRT2, t_max);
  }
  {
    // Non-unit direction vector
    const auto ray =
        TestFixture::R::fromOriginDir({1.0, 1.0, 1.0}, {-1.0, -1.0, 0.});
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-1.0, t_min);
    EXPECT_FLOAT_EQ(1.0, t_max);
  }
  {
    // positive direction
    const auto ray =
        TestFixture::R::fromOriginDir({1.0, 1.0, 1.0}, {1.0, 1.0, 0.});
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(-1.0, t_min);
    EXPECT_FLOAT_EQ(1.0, t_max);
  }
  {
    // too short, stops before grid
    const auto ray =
        TestFixture::R::fromOriginDir({-1.0, 3.0, 1.0}, {1.1, -0.9, 0.});
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_FALSE(intersect);
    EXPECT_GE(t_min, 1.0);
    EXPECT_GE(t_max, 1.0);
  }
  {
    // long enough, reaches into grid and ends within
    const auto ray =
        TestFixture::R::fromOriginDir({-1.0, 3.0, 1.0}, {1.1, -1.1, 0.});
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    EXPECT_LE(t_min, 1.0);
    EXPECT_GE(t_max, 1.0);
  }
  {
    // parallel to grid, does never intersect
    const auto ray =
        TestFixture::R::fromOriginEnd({-0.5, 0.5, 0.5}, {-0.5, 1.5, 0.5});
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_FALSE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_max));
  }
  {
    // parallel to grid, does intersect
    const auto ray =
        TestFixture::R::fromOriginEnd({-0.5, 0.5, 0.5}, {0.5, 0.5, 0.5});
    const auto intersect = rayBoxIntersection(
        ray, TestFixture::grid_, t_min, t_max, TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(0.5, t_min);
    EXPECT_FLOAT_EQ(2.5, t_max);
  }
}

TYPED_TEST(TestVoxel2x2x2Intersect, RayOutsideXYPlane) {
  TypeParam t_min, t_max;
  {
    //
    const auto ray =
        TestFixture::R::fromOriginEnd({-0.5, -0.5, -0.5}, {2.0, 2.0, 2.0});
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_GE(t_min, TypeParam{0.1});
    EXPECT_NEAR(t_max, TypeParam{1.0}, Epsilon<TypeParam>::eps);
  }
  {
    // outside of grid
    const auto ray =
        TestFixture::R::fromOriginEnd({2.1, 2.1, 2.1}, {2.2, 2.2, 2.2});
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_FALSE(intersect);

    EXPECT_NEAR(t_min, TypeParam{-21.0}, Epsilon<TypeParam>::eps);
    EXPECT_NEAR(t_max, TypeParam{-1.0}, Epsilon<TypeParam>::eps);
  }
}

TYPED_TEST(TestVoxel2x2x2Intersect, DegenerateCase) {
  TypeParam t_min, t_max;
  {
    // zero direction within grid -> does intersect with single point
    const auto ray =
        TestFixture::R::fromOriginDir({.5, .5, .5}, TestFixture::V3::Zero());
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_max));
  }
  {
    // zero direction within grid -> does intersect with single point
    // different position near corner
    const auto ray = TestFixture::R::fromOriginDir({1.99, 1.99, 0.01},
                                                   TestFixture::V3::Zero());
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_TRUE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_max));
  }
  {
    // zero direction outside of grid -> does not intersect
    const auto ray =
        TestFixture::R::fromOriginDir({.5, .5, -0.1}, TestFixture::V3::Zero());
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_FALSE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_max));
  }
  {
    // zero direction outside of grid -> does not intersect
    // different position
    const auto ray =
        TestFixture::R::fromOriginDir({.5, -0.01, .5}, TestFixture::V3::Zero());
    const auto intersect =
        rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max);
    EXPECT_FALSE(intersect);
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_min));
    EXPECT_FLOAT_EQ(std::numeric_limits<TypeParam>::infinity(),
                    std::fabs(t_max));
  }
}
