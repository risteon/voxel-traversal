//
//
#include <gtest/gtest.h>

#include <fstream>
#include <limits>
#include <numeric>

#include "voxel_traversal.h"

using namespace algorithm;

// Small voxel grid 2x2x2
class TestOnData : public ::testing::Test {
 protected:
  using float_type = Grid3DSpatialDef::float_type;
  using V3 = Grid3DSpatialDef::Vector3d;
  using C3 = Grid3DSpatialDef::Count3d;
  using V3f = Eigen::Vector3f;

  static const constexpr Grid3DSpatialDef::float_type SQRT2 = sqrt(2.0);
  static const constexpr Grid3DSpatialDef::float_type HALF_SQRT2 = 0.5 * SQRT2;
  static const constexpr Grid3DSpatialDef::float_type SQRT3 = sqrt(3.0);

  void SetUp() override {
    const V3 bound_min(0.0, -25.6, -2.0);
    const V3 bound_max(51.2, 25.6, 4.4);
    const C3 voxel_count(256, 256, 32);
    grid_ = Grid3DSpatialDef(bound_min, bound_max, voxel_count);

    ReadPoints(points_raw_, "points.bin");
    ReadPoints(origins_raw_, "ray_origins.bin");
    ReadPoints(frame_counts_, "counts.bin");

    TransformRawData();
  }

  void TransformRawData() {
    // check data consistency
    ASSERT_EQ(points_raw_.size() % 3, 0);
    ASSERT_EQ(origins_raw_.size() % 3, 0);
    const auto num_frames = origins_raw_.size() / 3;
    const auto num_points = points_raw_.size() / 3;
    ASSERT_EQ(num_frames, frame_counts_.size());
    const auto total_from_index =
        std::accumulate(frame_counts_.cbegin(), frame_counts_.cend(), 0);
    ASSERT_EQ(total_from_index, num_points);

    origins_.resize(num_frames);
    points_.resize(num_frames);
    for (std::size_t i = 0, j = 0; j < origins_.size(); ++j, i += 3) {
      origins_[j] = V3f(origins_raw_.data() + i);
    }
    std::size_t batch_count = 0;
    std::size_t current_frame = 0;
    for (const auto frame_count : frame_counts_) {
      auto& pp = points_[current_frame];
      pp.resize(frame_count);
      for (std::size_t i = 0, j = 0; j < frame_count; ++j, i += 3) {
        pp[j] = V3f(points_raw_.data() + batch_count * 3 + i);
      }
      batch_count += frame_count;
      ++current_frame;
    }
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
};

TEST_F(TestOnData, calcIntersections) {
  double t_min, t_max;
  {
    //
    const auto ray = Ray::fromEndpoints({-0.5, -0.5, -0.5}, {2.0, 2.0, 2.0});
    const auto intersect = rayBoxIntersection(ray, grid_, t_min, t_max);
    EXPECT_TRUE(true);
  }
}
