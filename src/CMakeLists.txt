add_library(coilgun STATIC
  coilgun.cpp
  physics/constants.cpp
  physics/elliptic.cpp
  physics/quadrature.cpp
  physics/lookup_tables.cpp
  physics/struve.cpp
  physics/self_inductance.cpp
  physics/mutual_inductance.cpp
  components/driving_coil.cpp
  components/armature.cpp
  simulation/sim_result.cpp
  simulation/termination.cpp
  simulation/cpu_phase_timing.cpp
  simulation/excitation.cpp
  simulation/integration_state.cpp
  simulation/single_stage_sim.cpp
  simulation/multi_stage_result.cpp
  simulation/multi_stage_sim.cpp
)

target_include_directories(coilgun
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_definitions(coilgun
  PUBLIC
    EIGEN_MAX_ALIGN_BYTES=32
    COILGUN_ENABLE_CPU_PHASE_TIMING=$<BOOL:${COILGUN_ENABLE_CPU_PHASE_TIMING}>
)

target_link_libraries(coilgun
  PUBLIC
    Eigen3::Eigen
    Boost::headers
)

# Measurement-only benchmark. It is excluded from the default build because
# execution is machine-dependent and phase timing is opt-in.
add_executable(bench_cpu_sim EXCLUDE_FROM_ALL
  ${CMAKE_SOURCE_DIR}/tests/bench_cpu_sim.cpp)
target_link_libraries(bench_cpu_sim PRIVATE coilgun)
target_include_directories(bench_cpu_sim PRIVATE ${CMAKE_SOURCE_DIR}/include)
