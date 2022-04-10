//
//
#include <gtest/gtest.h>

#include <limits>

#include "voxel_traversal.h"

using namespace voxel_traversal;

// Small voxel grid 2x2x2
template <typename float_type>
class TestVoxelTraversal : public ::testing::Test {
 protected:
  static const constexpr float_type SQRT2 = sqrt(2.0);
  static const constexpr float_type HALF_SQRT2 = 0.5 * SQRT2;
  static const constexpr float_type SQRT3 = sqrt(3.0);

  void expectTraversed(const TraversedVoxels<float_type>& expected,
                       const TraversedVoxels<float_type>& actual) {
    EXPECT_EQ(expected.size(), actual.size());
    EXPECT_TRUE(std::equal(
        expected.cbegin(), expected.cend(), actual.cbegin(),
        [](const auto& a, const auto& b) { return (a == b).all(); }));
  }
  void expectTraversedInOrderWithGaps(
      const TraversedVoxels<float_type>& expected,
      const TraversedVoxels<float_type>& actual) {
    auto it_traversed = actual.cbegin();

    for (const auto& exp : expected) {
      const auto fp =
          std::find_if(it_traversed, actual.cend(),
                       [&exp](const auto& a) { return (a == exp).all(); });
      if (fp == actual.cend()) {
        ADD_FAILURE();
        break;
      }
      it_traversed = fp;
    }
  }

  Grid3DSpatialDef<float_type> grid_;
  TraversedVoxels<float_type> traversed_voxels_;
};

// Small voxel grid 2x2x2
template <typename float_type>
class TestVoxel2x2x2Traversal : public TestVoxelTraversal<float_type> {
 protected:
  using V3 = typename Grid3DSpatialDef<float_type>::Vector3d;
  using C3 = typename Grid3DSpatialDef<float_type>::Index3d;
  using R = Ray<float_type>;

  void SetUp() override {
    const V3 bound_min(0.0, 0.0, 0.0);
    const V3 bound_max(2.0, 2.0, 2.0);
    const C3 voxel_count(2, 2, 2);
    this->grid_ =
        Grid3DSpatialDef<float_type>(bound_min, bound_max, voxel_count);
    this->traversed_voxels_.reserve(1000);
  }
};

// Slightly bigger grid 5x5x5
template <typename float_type>
class TestVoxel5x5x5Traversal : public TestVoxelTraversal<float_type> {
 protected:
  using V3 = typename Grid3DSpatialDef<float_type>::Vector3d;
  using C3 = typename Grid3DSpatialDef<float_type>::Index3d;
  using R = Ray<float_type>;

  void SetUp() override {
    const V3 bound_min(-10.0, -10.0, -10.0);
    const V3 bound_max(10.0, 10.0, 10.0);
    const C3 voxel_count(5, 5, 5);
    this->grid_ =
        Grid3DSpatialDef<float_type>(bound_min, bound_max, voxel_count);
    this->traversed_voxels_.reserve(1000);
  }
};

// cuboid grid
template <typename float_type>
class TestVoxel4x2x1Traversal : public TestVoxelTraversal<float_type> {
 protected:
  using V3 = typename Grid3DSpatialDef<float_type>::Vector3d;
  using C3 = typename Grid3DSpatialDef<float_type>::Index3d;
  using R = Ray<float_type>;

  void SetUp() override {
    const V3 bound_min(0.0, 0.0, -50.0);
    const V3 bound_max(4.0, 2.0, -15.0);
    const C3 voxel_count(4, 2, 1);
    this->grid_ =
        Grid3DSpatialDef<float_type>(bound_min, bound_max, voxel_count);
    this->traversed_voxels_.reserve(1000);
  }
};

// setup typed test suites
using Implementations = ::testing::Types<float, double>;
TYPED_TEST_SUITE(TestVoxel2x2x2Traversal, Implementations);
TYPED_TEST_SUITE(TestVoxel5x5x5Traversal, Implementations);
TYPED_TEST_SUITE(TestVoxel4x2x1Traversal, Implementations);

TYPED_TEST(TestVoxel2x2x2Traversal, GridProperties) {
  EXPECT_EQ(8ul, TestFixture::grid_.voxelCount());
  const Grid3DSpatialDef<TypeParam> grid_copy = TestFixture::grid_;
  EXPECT_EQ(8ul, grid_copy.voxelCount());
  EXPECT_TRUE((TestFixture::grid_.gridSize() == grid_copy.gridSize()).all());
  EXPECT_TRUE((TestFixture::grid_.voxelSize() == grid_copy.voxelSize()).all());
  EXPECT_TRUE((TestFixture::grid_.numVoxels() == grid_copy.numVoxels()).all());
  EXPECT_EQ(TestFixture::grid_.minBound(), grid_copy.minBound());
  EXPECT_EQ(TestFixture::grid_.maxBound(), grid_copy.maxBound());
  EXPECT_TRUE((TestFixture::grid_.gridSize() == grid_copy.gridSize()).all());
}

TYPED_TEST(TestVoxel2x2x2Traversal, AllDirectionsWithinGrid) {
  {
    // should traverse two voxels in X dir. Ray completely within grid
    const auto ray = TestFixture::R::fromOriginDir({.5, .5, .5}, {1., 0., 0.});
    TraversedVoxels<TypeParam> expected{{0, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    const auto ray = TestFixture::R::fromOriginDir({1.5, .5, .5}, {0., 1., 0.});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {1, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    const auto ray = TestFixture::R::fromOriginDir({1.5, .5, .5}, {0., 0., 1.});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {1, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  // Negative directions
  {
    // should traverse two voxels in X dir. Ray completely within grid
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, .5, .5}, {-1., 0., 0.});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {0, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, 1.5, .5}, {0., -1., 0.});
    TraversedVoxels<TypeParam> expected{{1, 1, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, .5, 1.5}, {0., 0., -1.});
    TraversedVoxels<TypeParam> expected{{1, 0, 1}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
}

TYPED_TEST(TestVoxel2x2x2Traversal, SingleVoxel) {
  {
    // only single voxel, ray too short to reach second
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, 1.5, 1.5}, {0.4, 0., 0.});
    TraversedVoxels<TypeParam> expected{{1, 1, 1}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    const auto ray =
        TestFixture::R::fromOriginEnd({-0.45, 0.5, 1.5}, {0.55, -0.5, 1.5});
    TraversedVoxels<TypeParam> expected{{0, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    const auto ray =
        TestFixture::R::fromOriginEnd({-0.5, 1.5, 0.55}, {0.5, 1.5, -0.45});
    TraversedVoxels<TypeParam> expected{{0, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
}

TYPED_TEST(TestVoxel2x2x2Traversal, NoVoxel) {
  {
    // only single voxel, ray too short to reach second
    const auto ray =
        TestFixture::R::fromOriginDir({1.5, 1.5, 2.1}, {0., 1., 0.});
    TraversedVoxels<TypeParam> expected{};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_FALSE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
}

TYPED_TEST(TestVoxel5x5x5Traversal, GridProperties) {
  EXPECT_EQ(125ul, TestFixture::grid_.voxelCount());
  const Grid3DSpatialDef<TypeParam> grid_copy = TestFixture::grid_;
  EXPECT_EQ(125ul, grid_copy.voxelCount());
  EXPECT_TRUE((TestFixture::grid_.gridSize() == grid_copy.gridSize()).all());
  EXPECT_TRUE((TestFixture::grid_.voxelSize() == grid_copy.voxelSize()).all());
  EXPECT_TRUE((TestFixture::grid_.numVoxels() == grid_copy.numVoxels()).all());
  EXPECT_EQ(TestFixture::grid_.minBound(), grid_copy.minBound());
  EXPECT_EQ(TestFixture::grid_.maxBound(), grid_copy.maxBound());
  EXPECT_TRUE((TestFixture::grid_.gridSize() == grid_copy.gridSize()).all());
}

TYPED_TEST(TestVoxel5x5x5Traversal, Diagonal) {
  TypeParam t_min, t_max;
  {
    // full diagonal. We do not assert specific order of off-diagonal voxels
    const auto ray = TestFixture::R::fromOriginDir({-20.0, -20.0, -20.0},
                                                   {40.0, 40.0, 40.0});
    TraversedVoxels<TypeParam> expected{
        {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}};
    //
    EXPECT_TRUE(rayBoxIntersection(ray, TestFixture::grid_, t_min, t_max));
    EXPECT_FLOAT_EQ(0.25, t_min);
    EXPECT_FLOAT_EQ(0.75, t_max);
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_);
    EXPECT_TRUE(intersect);
    EXPECT_GE(this->traversed_voxels_.size(), 5);
    TestFixture::expectTraversedInOrderWithGaps(expected,
                                                TestFixture::traversed_voxels_);
    EXPECT_TRUE((this->traversed_voxels_.front() == expected.front()).all());
    EXPECT_TRUE((this->traversed_voxels_.back() == expected.back()).all());
  }
}

TYPED_TEST(TestVoxel4x2x1Traversal, GridProperties) {
  EXPECT_EQ(8ul, TestFixture::grid_.voxelCount());
  const Grid3DSpatialDef<TypeParam> grid_copy = TestFixture::grid_;
  EXPECT_EQ(8ul, grid_copy.voxelCount());
  EXPECT_TRUE((TestFixture::grid_.gridSize() == grid_copy.gridSize()).all());
  EXPECT_TRUE((TestFixture::grid_.voxelSize() == grid_copy.voxelSize()).all());
  EXPECT_TRUE((TestFixture::grid_.numVoxels() == grid_copy.numVoxels()).all());
  EXPECT_EQ(TestFixture::grid_.minBound(), grid_copy.minBound());
  EXPECT_EQ(TestFixture::grid_.maxBound(), grid_copy.maxBound());
  EXPECT_TRUE((TestFixture::grid_.gridSize() == grid_copy.gridSize()).all());
}

TYPED_TEST(TestVoxel4x2x1Traversal, StoppingStartingRay) {
  {
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {6.0, 1.7, -25.0});
    TraversedVoxels<TypeParam> expected{
        {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_,
                                             TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // other way around
    const auto ray =
        TestFixture::R::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels<TypeParam> expected{
        {3, 1, 0}, {2, 1, 0}, {2, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_,
                                             TypeParam{0.0}, TypeParam{1.0});
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // shorter
    const auto ray =
        TestFixture::R::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels<TypeParam> expected{{3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_,
                                             TypeParam{0.0}, TypeParam{0.5});
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // shorter ray, use t1
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {2, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_,
                                             TypeParam{0.0}, TypeParam{0.9});
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // shorter ray, use t1
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels<TypeParam> expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_,
                                             TypeParam{0.0}, TypeParam{1.25});
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
  {
    // long ray, use t1
    const auto ray =
        TestFixture::R::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels<TypeParam> expected{
        {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, TestFixture::grid_,
                                             TestFixture::traversed_voxels_,
                                             TypeParam{0.0}, TypeParam{15.0});
    EXPECT_TRUE(intersect);
    TestFixture::expectTraversed(expected, TestFixture::traversed_voxels_);
  }
}
