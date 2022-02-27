//
//
#include <gtest/gtest.h>

#include <limits>

#include "voxel_traversal.h"

using namespace algorithm;

// Small voxel grid 2x2x2
class TestVoxelTraversal : public ::testing::Test {
 protected:
  using float_type = Grid3DSpatialDef::float_type;
  using V3 = Grid3DSpatialDef::Vector3d;
  using C3 = Grid3DSpatialDef::Index3d;
  using TraversedVoxels = std::vector<Grid3DSpatialDef::Index3d>;

  static const constexpr Grid3DSpatialDef::float_type SQRT2 = sqrt(2.0);
  static const constexpr Grid3DSpatialDef::float_type HALF_SQRT2 = 0.5 * SQRT2;
  static const constexpr Grid3DSpatialDef::float_type SQRT3 = sqrt(3.0);

  void expectTraversed(const TraversedVoxels& expected,
                       const TraversedVoxels& actual) {
    EXPECT_EQ(expected.size(), actual.size());
    EXPECT_TRUE(std::equal(
        expected.cbegin(), expected.cend(), actual.cbegin(),
        [](const auto& a, const auto& b) { return (a == b).all(); }));
  }
  void expectTraversedInOrderWithGaps(const TraversedVoxels& expected,
                                      const TraversedVoxels& actual) {
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

  Grid3DSpatialDef grid_;
  TraversedVoxels traversed_voxels_;
};

// Small voxel grid 2x2x2
class TestVoxel2x2x2Traversal : public TestVoxelTraversal {
 protected:
  void SetUp() override {
    const V3 bound_min(0.0, 0.0, 0.0);
    const V3 bound_max(2.0, 2.0, 2.0);
    const C3 voxel_count(2, 2, 2);
    grid_ = Grid3DSpatialDef(bound_min, bound_max, voxel_count);
    traversed_voxels_.reserve(1000);
  }
};

// Slightly bigger grid 5x5x5
class TestVoxel5x5x5Traversal : public TestVoxelTraversal {
 protected:
  void SetUp() override {
    const V3 bound_min(-10.0, -10.0, -10.0);
    const V3 bound_max(10.0, 10.0, 10.0);
    const C3 voxel_count(5, 5, 5);
    grid_ = Grid3DSpatialDef(bound_min, bound_max, voxel_count);
    traversed_voxels_.reserve(1000);
  }
};

// cuboid grid
class TestVoxel4x2x1Traversal : public TestVoxelTraversal {
 protected:
  void SetUp() override {
    const V3 bound_min(0.0, 0.0, -50.0);
    const V3 bound_max(4.0, 2.0, -15.0);
    const C3 voxel_count(4, 2, 1);
    grid_ = Grid3DSpatialDef(bound_min, bound_max, voxel_count);
    traversed_voxels_.reserve(1000);
  }
};

TEST_F(TestVoxel2x2x2Traversal, AllDirectionsWithinGrid) {
  {
    // should traverse two voxels in X dir. Ray completely within grid
    const auto ray = Ray::fromOriginDir({.5, .5, .5}, {1., 0., 0.});
    TraversedVoxels expected{{0, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {0., 1., 0.});
    TraversedVoxels expected{{1, 0, 0}, {1, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {0., 0., 1.});
    TraversedVoxels expected{{1, 0, 0}, {1, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  // Negative directions
  {
    // should traverse two voxels in X dir. Ray completely within grid
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {-1., 0., 0.});
    TraversedVoxels expected{{1, 0, 0}, {0, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    const auto ray = Ray::fromOriginDir({1.5, 1.5, .5}, {0., -1., 0.});
    TraversedVoxels expected{{1, 1, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    const auto ray = Ray::fromOriginDir({1.5, .5, 1.5}, {0., 0., -1.});
    TraversedVoxels expected{{1, 0, 1}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
}

TEST_F(TestVoxel2x2x2Traversal, SingleVoxel) {
  {
    // only single voxel, ray too short to reach second
    const auto ray = Ray::fromOriginDir({1.5, 1.5, 1.5}, {0.4, 0., 0.});
    TraversedVoxels expected{{1, 1, 1}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    const auto ray = Ray::fromOriginEnd({-0.45, 0.5, 1.5}, {0.55, -0.5, 1.5});
    TraversedVoxels expected{{0, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    const auto ray = Ray::fromOriginEnd({-0.5, 1.5, 0.55}, {0.5, 1.5, -0.45});
    TraversedVoxels expected{{0, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
}

TEST_F(TestVoxel2x2x2Traversal, NoVoxel) {
  {
    // only single voxel, ray too short to reach second
    const auto ray = Ray::fromOriginDir({1.5, 1.5, 2.1}, {0., 1., 0.});
    TraversedVoxels expected{};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_FALSE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
}

TEST_F(TestVoxel5x5x5Traversal, Diagonal) {
  float_type t_min, t_max;
  {
    // full diagonal. We do not assert specific order of off-diagonal voxels
    const auto ray =
        Ray::fromOriginDir({-20.0, -20.0, -20.0}, {40.0, 40.0, 40.0});
    TraversedVoxels expected{
        {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}};
    //
    EXPECT_TRUE(rayBoxIntersection(ray, grid_, t_min, t_max));
    EXPECT_FLOAT_EQ(0.25, t_min);
    EXPECT_FLOAT_EQ(0.75, t_max);
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_);
    EXPECT_TRUE(intersect);
    EXPECT_GE(traversed_voxels_.size(), 5);
    expectTraversedInOrderWithGaps(expected, traversed_voxels_);
    EXPECT_TRUE((traversed_voxels_.front() == expected.front()).all());
    EXPECT_TRUE((traversed_voxels_.back() == expected.back()).all());
  }
}

TEST_F(TestVoxel4x2x1Traversal, StoppingStartingRay) {
  {
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {6.0, 1.7, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // other way around
    const auto ray = Ray::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels expected{{3, 1, 0}, {2, 1, 0}, {2, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // shorter
    const auto ray = Ray::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels expected{{3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_, 0.0, 0.5);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // shorter ray, use t1
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_, 0.0, 0.9);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // shorter ray, use t1
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_, 0.0, 1.25);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
  {
    // long ray, use t1
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, traversed_voxels_, 0.0, 15.0);
    EXPECT_TRUE(intersect);
    expectTraversed(expected, traversed_voxels_);
  }
}
