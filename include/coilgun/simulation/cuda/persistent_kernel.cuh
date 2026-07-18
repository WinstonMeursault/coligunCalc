#pragma once
#include <cstddef>

namespace coilgun::simulation::cuda {
struct PersistentBuffers {
    int *batch_id=nullptr, *active_pairs=nullptr, *shutdown=nullptr, *doorbell=nullptr;
    double *seps=nullptr, *out_M=nullptr, *out_dM=nullptr;
    int *d_batch_id=nullptr, *d_active_pairs=nullptr, *d_shutdown=nullptr, *d_doorbell=nullptr;
    double *d_seps=nullptr, *d_out_M=nullptr, *d_out_dM=nullptr;
};
}

#ifdef __CUDACC__
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include <cuda_runtime.h>
#include <cstring>
#include <iostream>

namespace coilgun::simulation::cuda {

inline bool init_persistent_buffers(PersistentBuffers& p, int N_b, GpuBackend& be) {
    (void)be;
    // WC memory for host→device direction (seps, doorbell)
    auto ah_wc=[&](auto& x,size_t n)->bool{
        cudaError_t e=cudaHostAlloc(&x,n*sizeof(std::remove_pointer_t<decltype(x)>),cudaHostAllocMapped|cudaHostAllocWriteCombined);
        if(e!=cudaSuccess){ std::cerr<<"[coilgun] cudaHostAllocMapped(WC): "<<cudaGetErrorString(e)<<std::endl; return false; }
        return true;
    };
    // Regular mapped memory for device→host direction (out_M, out_dM) — avoid CPU cache coherency issues
    auto ah_reg=[&](auto& x,size_t n)->bool{
        cudaError_t e=cudaHostAlloc(&x,n*sizeof(std::remove_pointer_t<decltype(x)>),cudaHostAllocMapped);
        if(e!=cudaSuccess){ std::cerr<<"[coilgun] cudaHostAllocMapped: "<<cudaGetErrorString(e)<<std::endl; return false; }
        return true;
    };
    auto gdp=[](auto& d,auto h){ cudaHostGetDevicePointer(&d,h,0); };
    if(!ah_wc(p.batch_id,1))return false; gdp(p.d_batch_id,p.batch_id);
    if(!ah_wc(p.active_pairs,1))return false; gdp(p.d_active_pairs,p.active_pairs);
    if(!ah_wc(p.shutdown,1))return false; gdp(p.d_shutdown,p.shutdown);
    if(!ah_reg(p.doorbell,N_b))return false; gdp(p.d_doorbell,p.doorbell);
    if(!ah_wc(p.seps,N_b))return false; gdp(p.d_seps,p.seps);
    if(!ah_reg(p.out_M,N_b))return false; gdp(p.d_out_M,p.out_M);
    if(!ah_reg(p.out_dM,N_b))return false; gdp(p.d_out_dM,p.out_dM);
    memset(p.doorbell,0,N_b*sizeof(int)); memset(p.seps,0,N_b*sizeof(double));
    *p.batch_id=0; *p.active_pairs=0;*p.shutdown=0;
    return true;
}

inline void free_persistent_buffers(PersistentBuffers& p) {
    auto fh=[](void*& x){ if(x){ cudaFreeHost(x); x=nullptr; } };
    fh((void*&)p.batch_id);fh((void*&)p.active_pairs);fh((void*&)p.shutdown);fh((void*&)p.doorbell);
    fh((void*&)p.seps);fh((void*&)p.out_M);fh((void*&)p.out_dM);
}

// launch_persistent_kernel is defined in persistent_kernel.cu
void launch_persistent_kernel(const PersistentBuffers&, const GpuAdaptor&, int, int, int,
                               GpuOptLevel opt_level = GpuOptLevel::Full);

}
#endif
