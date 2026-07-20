#pragma once
#include <cstddef>
#include <cstdint>

namespace coilgun::simulation::cuda {
enum class PersistentStatus : std::uint8_t {
    Uninitialized,
    Ready,
    ShutdownRequested,
    ShutdownComplete,
    Fallback,
};

struct PersistentBuffers {
    // batch_id is retained as the source-compatible name for the generation token.
    int *batch_id=nullptr, *generation=nullptr, *active_pairs=nullptr;
    int *shutdown=nullptr, *shutdown_ack=nullptr, *doorbell=nullptr;
    double *seps=nullptr, *out_M=nullptr, *out_dM=nullptr;
    int *d_batch_id=nullptr, *d_generation=nullptr, *d_active_pairs=nullptr;
    int *d_shutdown=nullptr, *d_shutdown_ack=nullptr, *d_doorbell=nullptr;
    double *d_seps=nullptr, *d_out_M=nullptr, *d_out_dM=nullptr;
    PersistentStatus status=PersistentStatus::Uninitialized;
};
#ifdef __CUDACC__
inline void free_persistent_buffers(PersistentBuffers& p);
#endif
}

#ifdef __CUDACC__
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include <cuda_runtime.h>
#include <cstring>
#include <iostream>
#include <type_traits>

namespace coilgun::simulation::cuda {

inline bool init_persistent_buffers(PersistentBuffers& p, int N_b, GpuBackend& be) {
    (void)be;
    free_persistent_buffers(p);
    p.status = PersistentStatus::Uninitialized;

    auto ah_wc=[&](auto& x,size_t n)->bool {
        cudaError_t e=cudaHostAlloc(&x,n*sizeof(std::remove_pointer_t<decltype(x)>),cudaHostAllocMapped|cudaHostAllocWriteCombined);
        if(e!=cudaSuccess){ return false; }
        return true;
    };
    auto ah_reg=[&](auto& x,size_t n)->bool {
        cudaError_t e=cudaHostAlloc(&x,n*sizeof(std::remove_pointer_t<decltype(x)>),cudaHostAllocMapped);
        if(e!=cudaSuccess){ return false; }
        return true;
    };
    auto gdp=[&](auto& d,auto h)->bool { return cudaHostGetDevicePointer(&d,h,0)==cudaSuccess; };
    if(!ah_wc(p.batch_id,1) || !gdp(p.d_batch_id,p.batch_id)) goto fail;
    p.generation=p.batch_id; p.d_generation=p.d_batch_id;
    if(!ah_wc(p.active_pairs,1) || !gdp(p.d_active_pairs,p.active_pairs)) goto fail;
    if(!ah_wc(p.shutdown,1) || !gdp(p.d_shutdown,p.shutdown)) goto fail;
    if(!ah_reg(p.shutdown_ack,1) || !gdp(p.d_shutdown_ack,p.shutdown_ack)) goto fail;
    if(!ah_reg(p.doorbell,N_b) || !gdp(p.d_doorbell,p.doorbell)) goto fail;
    if(!ah_wc(p.seps,N_b) || !gdp(p.d_seps,p.seps)) goto fail;
    if(!ah_reg(p.out_M,N_b) || !gdp(p.d_out_M,p.out_M)) goto fail;
    if(!ah_reg(p.out_dM,N_b) || !gdp(p.d_out_dM,p.out_dM)) goto fail;
    std::memset(p.doorbell,0,N_b*sizeof(int)); std::memset(p.seps,0,N_b*sizeof(double));
    *p.batch_id=0; *p.active_pairs=0; *p.shutdown=0; *p.shutdown_ack=0;
    p.status=PersistentStatus::Ready;
    return true;

fail:
    p.status=PersistentStatus::Fallback;
    free_persistent_buffers(p);
    p.status=PersistentStatus::Fallback;
    return false;
}

inline void free_persistent_buffers(PersistentBuffers& p) {
    if (p.shutdown) {
        *p.shutdown=1;
        p.status=PersistentStatus::ShutdownRequested;
        cudaDeviceSynchronize();
    }
    auto fh=[](void*& x){ if(x){ cudaFreeHost(x); x=nullptr; } };
    fh((void*&)p.batch_id); p.generation=nullptr; p.d_batch_id=nullptr; p.d_generation=nullptr;
    fh((void*&)p.active_pairs);fh((void*&)p.shutdown);fh((void*&)p.shutdown_ack);
    fh((void*&)p.doorbell);
    fh((void*&)p.seps);fh((void*&)p.out_M);fh((void*&)p.out_dM);
    p.d_active_pairs=nullptr; p.d_shutdown=nullptr; p.d_shutdown_ack=nullptr;
    p.d_doorbell=nullptr; p.d_seps=nullptr; p.d_out_M=nullptr; p.d_out_dM=nullptr;
    if (p.status != PersistentStatus::Fallback)
        p.status=PersistentStatus::ShutdownComplete;
}

// launch_persistent_kernel is defined in persistent_kernel.cu
void launch_persistent_kernel(const PersistentBuffers&, const GpuAdaptor&, int, int, int,
                               GpuOptLevel opt_level = GpuOptLevel::Full);

}
#endif
