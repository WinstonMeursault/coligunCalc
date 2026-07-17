cmake_minimum_required(VERSION 3.20)
project(coilgun VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

cmake_policy(SET CMP0135 NEW)
cmake_policy(SET CMP0167 OLD)

include(FetchContent)

# --- Eigen (header-only linear algebra) ---
# Try system package first; download tarball as fallback.
find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found -- downloading v3.4.0 (one-time)")
  FetchContent_Declare(
    Eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
  )
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(Eigen)
endif()

# --- Boost (header-only math functions) ---
# Try system package first; download tarball as fallback.
find_package(Boost QUIET)
if(NOT Boost_FOUND)
  message(STATUS "Boost not found -- downloading v1.86.0 headers (one-time, ~85 MB)")
  FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-cmake.tar.gz
    URL_HASH SHA256=c62ce6e64d34414864fef946363db91cea89c1b90360eabed0515f0eda74c75c
  )
  FetchContent_MakeAvailable(Boost)
endif()

# --- Library target ---
add_subdirectory(src)

# --- Tests ---
option(COILGUN_BUILD_TESTS "Build unit tests" ON)
if(COILGUN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# --- Table generator (optional) ---
option(COILGUN_BUILD_GENERATOR "Build T(q,p) table generator" OFF)
if(COILGUN_BUILD_GENERATOR)
  add_executable(generate_t_table tools/generate_t_table.cpp)
  target_link_libraries(generate_t_table PRIVATE coilgun)
  target_include_directories(generate_t_table
    PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
endif()
