cmake_minimum_required(VERSION 3.21)
project(voxel_traversal)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(voxel_traversal voxel_traversal.cpp)
target_link_libraries(voxel_traversal Eigen3::Eigen)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
        test_voxel_traversal
        test_voxel_traversal.cpp
        voxel_traversal.cpp
)
target_link_libraries(
        test_voxel_traversal
        Eigen3::Eigen
        gtest_main
)

add_executable(
        test_on_data
        test_on_data.cpp
        voxel_traversal.cpp
)
target_link_libraries(
        test_on_data
        Eigen3::Eigen
        gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_voxel_traversal)
gtest_discover_tests(test_on_data PROPERTIES ENVIRONMENT "DATADIR=${CMAKE_CURRENT_SOURCE_DIR}/data/")
