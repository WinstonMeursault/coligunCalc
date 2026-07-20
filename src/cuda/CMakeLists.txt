add_library(coilgun_cuda STATIC
    gpu_elliptic.cu
    gpu_mutual_inductance.cu
    gpu_adaptor.cu
    gpu_single_stage_sim.cu
    gpu_multi_stage_sim.cu
    sim_batch.cu
    persistent_kernel.cu
    gpu_execution_context.cu
    gpu_solver.cu
    gpu_mutual_pipeline.cu
    gpu_state_kernels.cu
    gpu_thermal.cu
    gpu_engine.cu
    gpu_graph.cu
)

target_include_directories(coilgun_cuda
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

target_compile_definitions(coilgun_cuda PUBLIC COILGUN_CUDA_AVAILABLE=1)

target_link_libraries(coilgun_cuda
    PUBLIC
        coilgun
        Eigen3::Eigen
        Boost::headers
    PUBLIC
        CUDA::cudart
        CUDA::cublas
        CUDA::cusolver
)

target_compile_definitions(coilgun_cuda
    PRIVATE BOOST_MATH_ENABLE_CUDA EIGEN_MAX_ALIGN_BYTES=32
)

target_compile_options(coilgun_cuda
    PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler=-march=native>
)

# Compile/link smoke for the context itself without adding a test-suite file.
# It stays out of the default build because execution requires a CUDA device.
add_executable(coilgun_cuda_context_smoke EXCLUDE_FROM_ALL gpu_execution_context.cu)
target_compile_definitions(coilgun_cuda_context_smoke PRIVATE COILGUN_GPU_CONTEXT_SMOKE)
target_link_libraries(coilgun_cuda_context_smoke PRIVATE CUDA::cudart CUDA::cublas CUDA::cusolver)
target_include_directories(coilgun_cuda_context_smoke PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Measurement-only benchmark. It is excluded from the default build because
# execution requires a CUDA device and its timings are machine-dependent.
add_executable(bench_gpu_engine EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/tests/bench_gpu_engine.cu)
target_compile_definitions(bench_gpu_engine PRIVATE EIGEN_MAX_ALIGN_BYTES=32)
target_link_libraries(bench_gpu_engine PRIVATE coilgun_cuda coilgun CUDA::cudart)
target_include_directories(bench_gpu_engine PRIVATE ${CMAKE_SOURCE_DIR}/include)
