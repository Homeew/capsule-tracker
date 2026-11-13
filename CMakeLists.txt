cmake_minimum_required(VERSION 3.18)
project(capsule-tracker VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)
FetchContent_Declare(Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    SOURCE_SUBDIR  cmake)   # header-only: skip Eigen's own build
FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.10)
FetchContent_MakeAvailable(Eigen Catch2)

add_library(tracker_lib
    src/trajectory.cpp
    src/lm_baseline.cpp)
target_include_directories(tracker_lib PUBLIC include ${eigen_SOURCE_DIR})
target_compile_options(tracker_lib PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
add_executable(tests tests/test_tracker.cpp)
target_compile_options(tests PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(tests PRIVATE tracker_lib Catch2::Catch2)
add_test(NAME tracker_tests COMMAND tests)
