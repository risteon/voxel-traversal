cmake_minimum_required(VERSION 3.21)

# allow for project VERSION option
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()
if (POLICY CMP0003)
  cmake_policy(SET CMP0003 NEW)
endif()

project(VoxelTraversal VERSION 1.0 LANGUAGES CXX)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CXX_CLANG_COMPILE_OPTIONS
  "-march=native"
  "-stdlib=libstdc++"
  "-Weverything"
  "-Wno-c++98-compat"
  "-Wno-c++98-c++11-c++14-compat"
  )
set(CXX_GCC_COMPILE_OPTIONS
  "-march=native"
  "-Wall"
  "-Wno-unknown-pragmas"
  )

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ------------------ INSTALL CONFIGURATION -----------------
# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKEPACKAGE_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/VoxelTraversal/cmake"
  CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where voxelTraversal.cmake is installed"
  )

set(voxel_traversal_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
set(voxel_traversal_ROOT_DIR ${CMAKE_INSTALL_PREFIX})


# ------------- INSTALL CONFIGURATION (RPATH) --------------
# see https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# ------------------- FETCH GTEST ---------------------
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# ------------- VOXEL TRAVERSAL LIBRARY -----------------
add_library(voxel_traversal SHARED src/voxel_traversal.cpp)

set_target_properties(voxel_traversal PROPERTIES EXPORT_NAME VoxelTraversal)
set_target_properties(voxel_traversal PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(voxel_traversal PROPERTIES SOVERSION ${PROJECT_VERSION})
set_target_properties(voxel_traversal PROPERTIES PUBLIC_HEADER "src/grid.h;src/ray.h;src/voxel_traversal.h")

target_include_directories(voxel_traversal PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(voxel_traversal Eigen3::Eigen)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(voxel_traversal PRIVATE ${CXX_CLANG_COMPILE_OPTIONS})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(voxel_traversal PRIVATE ${CXX_GCC_COMPILE_OPTIONS})
endif()

set(voxel_traversal_INCLUDE_INSTALL_DIR
  "${CMAKE_INSTALL_INCLUDEDIR}/VoxelTraversal"
  CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where voxel_traversal header files are installed"
  )

install(TARGETS voxel_traversal EXPORT voxelTraversalTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${voxel_traversal_INCLUDE_INSTALL_DIR}
  )

# The VoxelTraversal target will be located in the VoxelTraversal namespace.
# Other CMake targets can refer to it using VoxelTraversal::<lib>.
export(TARGETS voxel_traversal NAMESPACE VoxelTraversal:: FILE VoxelTraversalTargets.cmake)

# ------------------- TESTS -------------------
enable_testing()
##################################
add_executable(
  test_voxel_traversal
  test/test_voxel_traversal.cpp
)
target_include_directories(test_voxel_traversal PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/VoxelTraversal>
  )
target_link_libraries(
  test_voxel_traversal
  voxel_traversal
  Eigen3::Eigen
  gtest_main
)
##################################
add_executable(
  test_voxel_counter
  test/test_voxel_counter.cpp
)
target_include_directories(test_voxel_counter PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
  $<INSTALL_INTERFACE:include/VoxelTraversal>
  )
target_link_libraries(
  test_voxel_counter
  voxel_traversal
  Eigen3::Eigen
  gtest_main
)
##################################
add_executable(
  test_voxel_intersect
  test/test_voxel_intersect.cpp
)
target_include_directories(test_voxel_intersect PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
  $<INSTALL_INTERFACE:include/VoxelTraversal>
  )
target_link_libraries(
  test_voxel_intersect
  voxel_traversal
  Eigen3::Eigen
  gtest_main
)
##################################
add_executable(
  test_on_data
  test/test_on_data.cpp
)
target_include_directories(test_on_data PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
  $<INSTALL_INTERFACE:include/VoxelTraversal>
  )
target_link_libraries(
  test_on_data
  voxel_traversal
  Eigen3::Eigen
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_voxel_traversal)
gtest_discover_tests(test_voxel_counter)
gtest_discover_tests(test_voxel_intersect)
gtest_discover_tests(test_on_data PROPERTIES ENVIRONMENT "DATADIR=${CMAKE_CURRENT_SOURCE_DIR}/data/")


# ------------------------ INSTALL -------------------------
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VoxelTraversalConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/VoxelTraversalConfig.cmake
  PATH_VARS voxel_traversal_INCLUDE_DIRS voxel_traversal_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
)

write_basic_package_version_file (VoxelTraversalConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

export(PACKAGE ${PROJECT_NAME})

install(EXPORT voxelTraversalTargets NAMESPACE VoxelTraversal:: DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/VoxelTraversalConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/VoxelTraversalConfigVersion.cmake
  DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})
