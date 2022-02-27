//
//
#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>

#include "voxel_traversal.h"

using namespace algorithm;

// Small voxel grid 2x2x2
class TestOnData : public ::testing::Test {
 protected:
  using float_type = Grid3DSpatialDef::float_type;
  using V3 = Grid3DSpatialDef::Vector3d;
  using C3 = Grid3DSpatialDef::Index3d;
  using V3f = Eigen::Vector3f;
  using TraversedVoxels = std::vector<Grid3DSpatialDef::Index3d>;

  void SetUp() override {
    const V3 bound_min(0.0, -25.6, -2.0);
    const V3 bound_max(51.2, 25.6, 4.4);
//    const C3 voxel_count(256, 256, 32);
    const C3 voxel_count(512, 512, 64);
    grid_ = Grid3DSpatialDef(bound_min, bound_max, voxel_count);

    ReadPoints(points_raw_, "points.bin");
    ReadPoints(origins_raw_, "ray_origins.bin");
    ReadPoints(frame_counts_, "counts.bin");

    TransformRawData();
    // free up memory
    points_raw_ = decltype(points_raw_){};
    origins_raw_ = decltype(origins_raw_){};
  }

  void TransformRawData() {
    // check data consistency
    ASSERT_EQ(points_raw_.size() % 3, 0);
    ASSERT_EQ(origins_raw_.size() % 3, 0);
    const auto num_frames = origins_raw_.size() / 3;
    const auto num_points = points_raw_.size() / 3;
    ASSERT_EQ(num_frames, frame_counts_.size());
    total_point_count_ =
        std::accumulate(frame_counts_.cbegin(), frame_counts_.cend(), 0);
    ASSERT_EQ(total_point_count_, num_points);

    origins_.resize(num_frames);
    points_.resize(num_frames);
    points_double_precision_.resize(num_frames);
    origins_double_precision_.resize(num_frames);

    for (std::size_t i = 0, j = 0; j < origins_.size(); ++j, i += 3) {
      origins_[j] = V3f(origins_raw_.data() + i);
      origins_double_precision_[j] = origins_[j].cast<double>();
    }
    std::size_t batch_count = 0;
    std::size_t current_frame = 0;
    for (const auto frame_count : frame_counts_) {
      auto& pp = points_[current_frame];
      auto& ppd = points_double_precision_[current_frame];
      pp.resize(frame_count);
      ppd.resize(frame_count);
      for (std::size_t i = 0, j = 0; j < frame_count; ++j, i += 3) {
        pp[j] = V3f(points_raw_.data() + batch_count * 3 + i);
        ppd[j] = pp[j].cast<double>();
      }
      batch_count += frame_count;
      ++current_frame;
    }

    ASSERT_EQ(points_.size(), origins_.size());
    ASSERT_EQ(points_double_precision_.size(), origins_double_precision_.size());
  }

  template <typename Scalar>
  static void ReadPoints(std::vector<Scalar>& into,
                         const std::string& filename) {
    const char* datadir = std::getenv("DATADIR");
    ASSERT_TRUE(datadir != nullptr) << "Env variable 'DATADIR' not available.";


    const auto filepath = std::string(datadir) + filename;
    std::ifstream f_points(filepath, std::ios::binary | std::ios::in);
    ASSERT_TRUE(f_points.good()) << "Cannot open " << filepath;
    f_points.seekg(0, std::ios::end);
    const auto num_bytes = f_points.tellg();
    ASSERT_TRUE(num_bytes % sizeof(Scalar) == 0);
    const auto num_values = num_bytes / sizeof(float);

    f_points.seekg(0);
    into.reserve(num_values);
    Scalar value{};
    while (f_points.read(reinterpret_cast<char*>(&value), sizeof(Scalar))) {
      into.push_back(value);
    }
    into.shrink_to_fit();
    ASSERT_EQ(into.size(), num_values);
  };

  Grid3DSpatialDef grid_;

  std::vector<float> points_raw_;
  std::vector<float> origins_raw_;
  std::vector<int32_t> frame_counts_;
  // multiple frames with their respective origin
  std::vector<std::vector<V3f>> points_;
  std::vector<V3f> origins_;
  std::size_t total_point_count_;
  // double precision
  std::vector<std::vector<V3>> points_double_precision_;
  std::vector<V3> origins_double_precision_;
};

TEST_F(TestOnData, timeIntersectionCalculation) {
  using std::chrono::duration;
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;

  double t_min, t_max;
  uint32_t intersection_count{0};
  TraversedVoxels traversedVoxels{};
  traversedVoxels.reserve(1000);

  auto t1 = high_resolution_clock::now();
  for (std::size_t frame = 0; frame < points_double_precision_.size(); ++frame) {
    const auto& pp = points_double_precision_[frame];
    const V3& origin = origins_double_precision_[frame];

    if (frame % 50 == 0) {
      std::cout << "...processing frame #" << frame << " with #" << pp.size()
                << " points..." << std::endl;
    }

    for (const auto& point : pp) {
      const auto ray =
          Ray::fromOriginEnd(origin, point);
      const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
      if (intersect) {
        ++intersection_count;
      }
    }
  }
  auto t2 = high_resolution_clock::now();
  const duration<double> ss = t2 - t1;

  std::cout << "Elapsed seconds: " << ss.count() << std::endl;
  std::cout << "Total number of rays: " << total_point_count_ << std::endl;
  std::cout << "Total number of intersections: " << intersection_count
            << std::endl;

  const std::size_t expected_hits{4097337};
  std::cout << "Expected #" << expected_hits << " hits." << std::endl;
  if (intersection_count != expected_hits) {
    std::cout << "!! Deviation of "
              << (static_cast<double>(intersection_count) /
                      static_cast<double>(expected_hits) -
                  1.0) *
                     100.0
              << "%" << std::endl;
  }
  else {
    std::cout << "... exactly as expected." << std::endl;
  }
}

TEST_F(TestOnData, timeTraversalCalculation) {
  using std::chrono::duration;
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;

  double t_min, t_max;
  uint32_t intersection_count{0};
  uint64_t traversal_count{0};

  TraversedVoxels traversedVoxels{};
  traversedVoxels.reserve(1000);

  auto t1 = high_resolution_clock::now();
  for (std::size_t frame = 0; frame < points_double_precision_.size(); ++frame) {
    const auto& pp = points_double_precision_[frame];
    const V3& origin = origins_double_precision_[frame];

    if (frame % 50 == 0) {
      std::cout << "...processing frame #" << frame << " with #" << pp.size()
                << " points..." << std::endl;
    }

    for (const auto& point : pp) {
      const auto ray =
          Ray::fromOriginEnd(origin, point);

      const auto intersect = traverseVoxelGrid(ray, grid_, traversedVoxels, 0.0, 1.0);
      if (intersect) {
        ++intersection_count;
      }
      traversal_count += traversedVoxels.size();
    }
  }
  auto t2 = high_resolution_clock::now();
  const duration<double> ss = t2 - t1;

  std::cout << "Elapsed seconds: " << ss.count() << std::endl;
  std::cout << "Total number of rays: " << total_point_count_ << std::endl;
  std::cout << "Total number of intersections: " << intersection_count
            << std::endl;
  // originally: 278 963 540
  std::cout << "Total number of traversed voxels: " << traversal_count
            << std::endl;
}

TEST_F(TestOnData, sanityCheckTraversalCalculation) {

  TraversedVoxels traversedVoxels{};
  traversedVoxels.reserve(1000);

  for (std::size_t frame = 0; frame < points_double_precision_.size(); ++frame) {
    const auto& pp = points_double_precision_[frame];
    const V3& origin = origins_double_precision_[frame];

    if (frame % 50 == 0) {
      std::cout << "...processing frame #" << frame << " with #" << pp.size()
                << " points..." << std::endl;
    }

    for (const auto& point : pp) {
      const auto ray =
          Ray::fromOriginEnd(origin, point);

      const auto intersect = traverseVoxelGrid(ray, grid_, traversedVoxels, 0.0, 1.0);
      for (const auto& tv : traversedVoxels){
        EXPECT_TRUE((tv >= 0).all());
        EXPECT_TRUE((tv < grid_.numVoxels()).all());
      }
      EXPECT_LE(traversedVoxels.size(), grid_.numVoxels().x() + grid_.numVoxels().y());
    }
  }
}
