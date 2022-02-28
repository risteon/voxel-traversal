//
//
#include <gtest/gtest.h>

#include <limits>

#include "voxel_traversal.h"

using namespace algorithm;

// Small voxel grid 2x2x2
class TestVoxelCounter : public ::testing::Test {
 protected:
  using float_type = Grid3DTraversalCounter::float_type;
  using V3 = Grid3DTraversalCounter::Vector3d;
  using C3 = Grid3DTraversalCounter::Index3d;
  using TraversedVoxels = std::vector<Grid3DSpatialDef::Index3d>;

  static const constexpr Grid3DTraversalCounter::float_type SQRT2 = sqrt(2.0);
  static const constexpr Grid3DTraversalCounter::float_type HALF_SQRT2 =
      0.5 * SQRT2;
  static const constexpr Grid3DTraversalCounter::float_type SQRT3 = sqrt(3.0);

  void expectCounted(const TraversedVoxels& expected,
                     bool more_allowed = false) const {
    // copy
    Grid3DTraversalCounter::tensor_type counts = grid_.getCounter();
    for (const auto& index : expected) {
      EXPECT_GE(counts(index), 1);
      if (counts(index) > 0) counts(index)--;
    }
    if (!more_allowed) {
      const Eigen::Tensor<uint64_t, 0> max_value = counts.maximum();
      EXPECT_EQ(0, max_value());
    }
  }

  Grid3DTraversalCounter grid_;
};

// Small voxel grid 2x2x2
class TestVoxel2x2x2Counter : public TestVoxelCounter {
 protected:
  void SetUp() override {
    const V3 bound_min(0.0, 0.0, 0.0);
    const V3 bound_max(2.0, 2.0, 2.0);
    const C3 voxel_count(2, 2, 2);
    grid_ = Grid3DTraversalCounter(bound_min, bound_max, voxel_count);
  }
};

// Slightly bigger grid 5x5x5
class TestVoxel5x5x5Counter : public TestVoxelCounter {
 protected:
  void SetUp() override {
    const V3 bound_min(-10.0, -10.0, -10.0);
    const V3 bound_max(10.0, 10.0, 10.0);
    const C3 voxel_count(5, 5, 5);
    grid_ = Grid3DTraversalCounter(bound_min, bound_max, voxel_count);
  }
};

// cuboid grid
class TestVoxel4x2x1Counter : public TestVoxelCounter {
 protected:
  void SetUp() override {
    const V3 bound_min(0.0, 0.0, -50.0);
    const V3 bound_max(4.0, 2.0, -15.0);
    const C3 voxel_count(4, 2, 1);
    grid_ = Grid3DTraversalCounter(bound_min, bound_max, voxel_count);
  }
};

TEST_F(TestVoxel2x2x2Counter, MultipleRays) {
  grid_.clear();
  TraversedVoxels expected{{0, 0, 0}, {1, 0, 0}, {1, 0, 0}, {1, 1, 0},
                           {1, 0, 0}, {1, 0, 1}, {1, 0, 0}, {0, 0, 0},
                           {1, 1, 0}, {1, 0, 0}, {1, 0, 1}, {1, 0, 0}};
  {
    const auto ray = Ray::fromOriginDir({.5, .5, .5}, {1., 0., 0.});
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {0., 1., 0.});
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {0., 0., 1.});
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
  }
  // Negative directions
  {
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {-1., 0., 0.});
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray = Ray::fromOriginDir({1.5, 1.5, .5}, {0., -1., 0.});
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
  }
  {
    const auto ray = Ray::fromOriginDir({1.5, .5, 1.5}, {0., 0., -1.});
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
  }
  expectCounted(expected);
}

TEST_F(TestVoxel2x2x2Counter, AllDirectionsWithinGrid) {
  {
    // should traverse two voxels in X dir. Ray completely within grid
    grid_.clear();
    const auto ray = Ray::fromOriginDir({.5, .5, .5}, {1., 0., 0.});
    TraversedVoxels expected{{0, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    grid_.clear();
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {0., 1., 0.});
    TraversedVoxels expected{{1, 0, 0}, {1, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    grid_.clear();
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {0., 0., 1.});
    TraversedVoxels expected{{1, 0, 0}, {1, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  // Negative directions
  {
    // should traverse two voxels in X dir. Ray completely within grid
    grid_.clear();
    const auto ray = Ray::fromOriginDir({1.5, .5, .5}, {-1., 0., 0.});
    TraversedVoxels expected{{1, 0, 0}, {0, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    // should traverse two voxels in Y dir. Ray completely within grid
    grid_.clear();
    const auto ray = Ray::fromOriginDir({1.5, 1.5, .5}, {0., -1., 0.});
    TraversedVoxels expected{{1, 1, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    // should traverse two voxels in Z dir. Ray completely within grid
    grid_.clear();
    const auto ray = Ray::fromOriginDir({1.5, .5, 1.5}, {0., 0., -1.});
    TraversedVoxels expected{{1, 0, 1}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
}

TEST_F(TestVoxel2x2x2Counter, SingleVoxel) {
  {
    // only single voxel, ray too short to reach second
    grid_.clear();
    const auto ray = Ray::fromOriginDir({1.5, 1.5, 1.5}, {0.4, 0., 0.});
    TraversedVoxels expected{{1, 1, 1}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    grid_.clear();
    const auto ray = Ray::fromOriginEnd({-0.45, 0.5, 1.5}, {0.55, -0.5, 1.5});
    TraversedVoxels expected{{0, 0, 1}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    // only single voxel, cut through corner
    // -> make sure that there is no infinite loop
    grid_.clear();
    const auto ray = Ray::fromOriginEnd({-0.5, 1.5, 0.55}, {0.5, 1.5, -0.45});
    TraversedVoxels expected{{0, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
}

TEST_F(TestVoxel2x2x2Counter, NoVoxel) {
  {
    // only single voxel, ray too short to reach second
    grid_.clear();
    const auto ray = Ray::fromOriginDir({1.5, 1.5, 2.1}, {0., 1., 0.});
    TraversedVoxels expected{};
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_FALSE(intersect);
    expectCounted(expected);
  }
}

TEST_F(TestVoxel5x5x5Counter, Diagonal) {
  float_type t_min, t_max;
  {
    // full diagonal. We do not assert specific order of off-diagonal voxels
    grid_.clear();
    const auto ray =
        Ray::fromOriginDir({-20.0, -20.0, -20.0}, {40.0, 40.0, 40.0});
    TraversedVoxels expected{
        {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}};
    //
    EXPECT_TRUE(rayBoxIntersection(ray, grid_, t_min, t_max));
    EXPECT_FLOAT_EQ(0.25, t_min);
    EXPECT_FLOAT_EQ(0.75, t_max);
    const auto intersect = traverseVoxelGrid(ray, grid_);
    EXPECT_TRUE(intersect);
    const Eigen::Tensor<uint64_t, 0> summed = grid_.getCounter().sum();
    EXPECT_GE(summed(), 5);
    expectCounted(expected, true);
  }
}

TEST_F(TestVoxel4x2x1Counter, StoppingStartingRay) {
  {
    grid_.clear();
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {6.0, 1.7, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    grid_.clear();
    // other way around
    const auto ray = Ray::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels expected{{3, 1, 0}, {2, 1, 0}, {2, 0, 0}, {1, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, 0.0, 1.0);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    grid_.clear();
    // shorter
    const auto ray = Ray::fromOriginEnd({6.0, 1.7, -25.0}, {1.5, 0.8, -25.0});
    TraversedVoxels expected{{3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, 0.0, 0.5);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    grid_.clear();
    // shorter ray, use t1
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, 0.0, 0.9);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    grid_.clear();
    // shorter ray, use t1
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, 0.0, 1.25);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
  {
    grid_.clear();
    // long ray, use t1
    const auto ray = Ray::fromOriginEnd({1.5, 0.8, -25.0}, {2.5, 1.0, -25.0});
    TraversedVoxels expected{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}};
    const auto intersect = traverseVoxelGrid(ray, grid_, 0.0, 15.0);
    EXPECT_TRUE(intersect);
    expectCounted(expected);
  }
}
