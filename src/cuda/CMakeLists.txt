add_library(coilgun_cuda STATIC
    gpu_elliptic.cu
    gpu_mutual_inductance.cu
    gpu_adaptor.cu
    gpu_single_stage_sim.cu
    gpu_multi_stage_sim.cu
)

target_include_directories(coilgun_cuda
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

target_link_libraries(coilgun_cuda
    PUBLIC
        coilgun
        Eigen3::Eigen
        Boost::headers
    PRIVATE
        CUDA::cudart
)

target_compile_definitions(coilgun_cuda
    PRIVATE BOOST_MATH_ENABLE_CUDA EIGEN_MAX_ALIGN_BYTES=32
)

target_compile_options(coilgun_cuda
    PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler=-march=native>
)
