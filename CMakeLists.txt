cmake_minimum_required(VERSION 3.20)
project(coilgunCalc VERSION 4.1.0 LANGUAGES CXX C)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
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

find_package(OpenMP REQUIRED)
target_link_libraries(coilgun PUBLIC OpenMP::OpenMP_CXX)

# --- CUDA acceleration (optional) ---
option(COILGUN_ENABLE_CUDA "Build GPU-accelerated backend" OFF)
option(COILGUN_ENABLE_CPU_PHASE_TIMING
       "Enable low-overhead CPU derivative phase instrumentation" OFF)
if(COILGUN_ENABLE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_ARCHITECTURES "native" CACHE STRING "Target CUDA architectures")
    find_package(CUDAToolkit REQUIRED COMPONENTS cublas cusolver)
    add_subdirectory(src/cuda)
endif()

# --- Tests ---
option(COILGUN_BUILD_TESTS "Build unit tests" ON)
if(COILGUN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# --- Table generator (optional) ---
option(COILGUN_BUILD_GENERATOR "Build T(q,p) table generator" OFF)
option(COILGUN_BUILD_GENERATOR_TESTS "Build T(q,p) generator unit tests" OFF)
if(COILGUN_BUILD_GENERATOR)
  add_executable(generate_t_table tools/generate_t_table.cpp)
  target_link_libraries(generate_t_table PRIVATE coilgun)
  target_include_directories(generate_t_table
    PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
endif()

# --- Installation and package export ---
include(CMakePackageConfigHelpers)

install(TARGETS coilgun
  EXPORT coilgunTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES include/coilgun/coilgun.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coilgun
)
install(FILES include/coilgun/core/types.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coilgun/core
)
install(FILES
  include/coilgun/components/armature.hpp
  include/coilgun/components/driving_coil.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coilgun/components
)
install(FILES
  include/coilgun/physics/cache.hpp
  include/coilgun/physics/constants.hpp
  include/coilgun/physics/elliptic.hpp
  include/coilgun/physics/lookup_tables.hpp
  include/coilgun/physics/mutual_inductance.hpp
  include/coilgun/physics/quadrature.hpp
  include/coilgun/physics/self_inductance.hpp
  include/coilgun/physics/struve.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coilgun/physics
)
install(FILES
  include/coilgun/simulation/cpu_phase_timing.hpp
  include/coilgun/simulation/derivative_workspace.hpp
  include/coilgun/simulation/excitation.hpp
  include/coilgun/simulation/excitation_snapshot.hpp
  include/coilgun/simulation/integration_state.hpp
  include/coilgun/simulation/multi_stage_result.hpp
  include/coilgun/simulation/multi_stage_sim.hpp
  include/coilgun/simulation/sim_result.hpp
  include/coilgun/simulation/single_stage_sim.hpp
  include/coilgun/simulation/termination.hpp
  include/coilgun/simulation/time_stepper.hpp
  include/coilgun/simulation/trigger_config.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coilgun/simulation
)
install(DIRECTORY include/coilgun/optimization/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coilgun/optimization
  FILES_MATCHING PATTERN "*.hpp"
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/coilgunConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/coilgunConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/coilgunConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/coilgun"
)
install(EXPORT coilgunTargets
  FILE coilgunTargets.cmake
  NAMESPACE coilgun::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/coilgun"
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/coilgunConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/coilgunConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/coilgun"
)
