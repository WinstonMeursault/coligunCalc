# doctest -- try system package first; download tarball as fallback
find_package(doctest QUIET)
if(NOT doctest_FOUND)
  message(STATUS "doctest not found -- downloading v2.5.3 (one-time)")
  FetchContent_Declare(
    doctest
    URL https://github.com/doctest/doctest/archive/refs/tags/v2.5.3.tar.gz
    URL_HASH SHA256=174ebc4e769928959614789c5b4e9c3d0a0f81a62bb608756b127bfebfb21331
  )
  FetchContent_MakeAvailable(doctest)
endif()

function(add_coilgun_test name)
  add_executable(${name} ${name}.cpp test_main.cpp)
  target_link_libraries(${name} PRIVATE coilgun doctest::doctest)
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_coilgun_test(test_elliptic)
add_coilgun_test(test_lookup)
add_coilgun_test(test_struve)
add_coilgun_test(test_quadrature)
add_coilgun_test(test_self_inductance)
add_coilgun_test(test_mutual_inductance)
add_coilgun_test(test_driving_coil)
add_coilgun_test(test_armature)
add_coilgun_test(test_integration)
add_coilgun_test(test_single_stage_sim)
add_coilgun_test(test_multi_stage_sim)

if(COILGUN_ENABLE_CUDA)
    # Tests with direct kernel launches (<<<>>>) — must compile with nvcc.
    function(add_gpu_test_cuda name)
        add_executable(${name} ${name}.cpp test_main.cpp)
        set_source_files_properties(${name}.cpp PROPERTIES LANGUAGE CUDA)
        set_target_properties(${name} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
        target_compile_definitions(${name} PRIVATE EIGEN_MAX_ALIGN_BYTES=32)
        target_link_libraries(${name} PRIVATE coilgun_cuda coilgun doctest::doctest)
        target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
        add_test(NAME ${name} COMMAND ${name})
    endfunction()

    # Tests using GPU API classes (no kernel launches) — compile with g++.
    function(add_gpu_test name)
        add_executable(${name} ${name}.cpp test_main.cpp)
        target_link_libraries(${name} PRIVATE coilgun_cuda coilgun doctest::doctest
                              CUDA::cudart)
        target_compile_definitions(${name} PRIVATE EIGEN_MAX_ALIGN_BYTES=32)
        target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
        add_test(NAME ${name} COMMAND ${name})
    endfunction()

    add_gpu_test_cuda(test_gpu_elliptic)
    add_gpu_test_cuda(test_gpu_filament)
    add_gpu_test_cuda(test_gpu_coil_pair)
    add_gpu_test(test_gpu_vs_cpu_single)
    add_gpu_test(test_gpu_vs_cpu_multi)
    add_gpu_test(test_gpu_batch)
    add_gpu_test(test_gpu_sim_batch)
endif()
