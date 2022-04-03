//
//
#include <gtest/gtest.h>

#include <limits>

#include "voxel_traversal.h"

using namespace voxel_traversal;

// Small voxel grid 2x2x2
template <typename float_type>
class TestVoxelCounter : public ::testing::Test {
 protected:
  static const constexpr float_type SQRT2 = sqrt(2.0);
  static const constexpr float_type HALF_SQRT2 = 0.5 * SQRT2;
  static const constexpr float_type SQRT3 = sqrt(3.0);

  void expectCounted(const TraversedVoxels<float_type>& expected,
                     bool more_allowed = false) const {
    // copy
    typename Grid3DTraversalCounter<float_type>::tensor_type counts =
        grid_.getCounter();
    for (const auto& index : expected) {
      EXPECT_GE(counts(index), 1);
      if (counts(index) > 0) counts(index)--;
    }
    if (!more_allowed) {
      const Eigen::Tensor<uint64_t, 0> max_value = counts.maximum();
      EXPECT_EQ(0, max_value());
    }
  }

  Grid3DTraversalCounter<float_type> grid_;
};

// Small voxel grid 2x2x2
template <typename float_type>
class TestVoxel2x2x2Counter : public TestVoxelCounter<float_type> {
 protected:
  using V3 = typename Grid3DTraversalCounter<float_type>::Vector3d;
  using C3 = typename Grid3DTraversalCounter<float_type>::Index3d;
  using R = Ray<float_type>;

  void SetUp() override {
    const V3 bound_min(0.0, 0.0, 0.0);
    const V3 bound_max(2.0, 2.0, 2.0);
    const C3 voxel_count(2, 2, 2);
    this->grid_ =
        Grid3DTraversalCounter<float_type>(bound_min, bound_max, voxel_count);
  }
};

// Slightly bigger grid 5x5x5
template <typename float_type>
class TestVoxel5x5x5Counter : public TestVoxelCounter<float_type> {
 protected:
  using V3 = typename Grid3DTraversalCounter<float_type>::Vector3d;
  using C3 = typename Grid3DTraversalCounter<float_type>::Index3d;
  using R = Ray<float_type>;

  void SetUp() override {
    const V3 bound_min(-10.0, -10.0, -10.0);
    const V3 bound_max(10.0, 10.0, 10.0);
    const C3 voxel_count(5, 5, 5);
    this->grid_ =
        Grid3DTraversalCounter<float_type>(bound_min, bound_max, voxel_count);
  }
};

// cuboid grid
template <typename float_type>
class TestVoxel4x2x1Counter : public TestVoxelCounter<float_type> {
 protected:
  using V3 = typename Grid3DTraversalCounter<float_type>::Vector3d;
  using C3 = typename Grid3DTraversalCounter<float_type>::Index3d;
  using R = Ray<float_type>;

  void SetUp() override {
    const V3 bound_min(0.0, 0.0, -50.0);
    const V3 bound_max(4.0, 2.0, -15.0);
    const C3 voxel_count(4, 2, 1);
    this->grid_ =
        Grid3DTraversalCounter<float_type>(bound_min, bound_max, voxel_count);
  }
};

// setup typed test suites
using Implementations = ::testing::Types<float, double>;
TYPED_TEST_SUITE(TestVoxel2x2x2Counter, Implementations);
TYPED_TEST_SUITE(TestVoxel5x5x5Counter, Implementations);
TYPED_TEST_SUITE(TestVoxel4x2x1Counter, Implementations);

TYPED_TEST(TestVoxel2x2x2Counter, MultipleRays) {
  TestFixture::grid_.clear();
  TraversedVoxels<TypeParam> expected{
      {0, 0, 0}, {1, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 0, 0}, {1, 0, 1},
      {1, 0, 0}, {0, 0, 0}, {1, 1, 0}, {1, 0, 0}, {1, 0, 1}, {1, 0, 0}};
  {
    const auto ray = TestFixture::R::fromOriginDir({.5, .5, .5}, {1., 0., 0.});
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray = TestFixture::R::fromOriginDir({1.5, .5, .5}, {0., 1., 0.});
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray = TestFixture::R::fromOriginDir({1.5, .5, .5}, {0., 0., 1.});
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
  }
  // Negative directions
  {
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, .5, .5}, {-1., 0., 0.});
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, 1.5, .5}, {0., -1., 0.});
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, .5, 1.5}, {0., 0., -1.});
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
  }
  TestFixture::expectCounted(expected);
}

TYPED_TEST(TestVoxel2x2x2Counter, AllDirectionsWithinGrid) {
  {
    // should traverse two voxels in X dir. Ray completely within grid
    TestFixture::grid_.clear();
    const auto ray = TestFixture::R::fromOriginDir({.5, .5, .5}, {1., 0., 0.});
    TraversedVoxels<TypeParam> expected{{0, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    TestFixture::grid_.clear();
    const auto ray = TestFixture::R::fromOriginDir({1.5, .5, .5}, {0., 1., 0.});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {1, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    TestFixture::grid_.clear();
    const auto ray = TestFixture::R::fromOriginDir({1.5, .5, .5}, {0., 0., 1.});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {1, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  // Negative directions
  {
    // should traverse two voxels in X dir. Ray completely within grid
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, .5, .5}, {-1., 0., 0.});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {0, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, 1.5, .5}, {0., -1., 0.});
    TraversedVoxels<TypeParam> expected{{1, 1, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, .5, 1.5}, {0., 0., -1.});
    TraversedVoxels<TypeParam> expected{{1, 0, 1}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
}

TYPED_TEST(TestVoxel2x2x2Counter, SingleVoxel) {
  {
    // only single voxel, ray too short to reach second
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, 1.5, 1.5}, {0.4, 0., 0.});
    TraversedVoxels<TypeParam> expected{{1, 1, 1}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginEnd({-0.45, 0.5, 1.5}, {0.55, -0.5, 1.5});
    TraversedVoxels<TypeParam> expected{{0, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginEnd({-0.5, 1.5, 0.55}, {0.5, 1.5, -0.45});
    TraversedVoxels<TypeParam> expected{{0, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
}

TYPED_TEST(TestVoxel2x2x2Counter, NoVoxel) {
  {
    // only single voxel, ray too short to reach second
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, 1.5, 2.1}, {0., 1., 0.});
    TraversedVoxels<TypeParam> expected{};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_FALSE(intersect);
    TestFixture::expectCounted(expected);
  }
}

TYPED_TEST(TestVoxel5x5x5Counter, Diagonal) {
  TypeParam t_min, t_max;
  {
    // full diagonal. We do not assert specific order of off-diagonal voxels
    TestFixture::grid_.clear();
    const auto ray = TestFixture::R::fromOriginDir({-20.0, -20.0, -20.0},
                                                   {40.0, 40.0, 40.0});
    TraversedVoxels<TypeParam> expected{
        {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}};
    //
    EXPECT_TRUE(rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max));
    EXPECT_FLOAT_EQ(0.25, t_min);
    EXPECT_FLOAT_EQ(0.75, t_max);
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_);
    EXPECT_TRUE(intersect);
    const Eigen::Tensor<uint64_t, 0> summed =
        TestFixture::grid_.getCounter().sum();
    EXPECT_GE(summed(), 5);
    TestFixture::expectCounted(expected, true);
  }
}

TYPED_TEST(TestVoxel4x2x1Counter, StoppingStartingRay) {
  {
    TestFixture::grid_.clear();
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {6.0, 1.7, -25.0});
    TraversedVoxels<TypeParam> expected{
        {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    TestFixture::grid_.clear();
    // other way around
    const auto ray =
        TestFixture::R::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels<TypeParam> expected{
        {3, 1, 0}, {2, 1, 0}, {2, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    TestFixture::grid_.clear();
    // shorter
    const auto ray =
        TestFixture::R::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels<TypeParam> expected{{3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TypeParam{0.0}, TypeParam{0.5});
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    TestFixture::grid_.clear();
    // shorter ray, use t1
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {2, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TypeParam{0.0}, TypeParam{0.9});
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    TestFixture::grid_.clear();
    // shorter ray, use t1
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TypeParam{0.0}, TypeParam{1.25});
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
  {
    TestFixture::grid_.clear();
    // long ray, use t1
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels<TypeParam> expected{
        {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TypeParam{0.0}, TypeParam{15.0});
    EXPECT_TRUE(intersect);
    TestFixture::expectCounted(expected);
  }
}
