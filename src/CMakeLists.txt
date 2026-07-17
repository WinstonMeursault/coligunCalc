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
  simulation/excitation.cpp
  simulation/single_stage_sim.cpp
  simulation/multi_stage_result.cpp
  simulation/multi_stage_sim.cpp
)

target_include_directories(coilgun
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(coilgun
  PUBLIC
    Eigen3::Eigen
    Boost::headers
)
