/**
 * @file persistent_kernel.cu
 * @brief Persistent kernel — compiled once, linked with other .cu files.
 * @author Winston Meursault
 */
#include "coilgun/simulation/cuda/persistent_kernel.cuh"
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/physics/mutual_inductance.cuh"

namespace coilgun::simulation::cuda {

__global__ void persistent_batch_kernel(
    const CoilGeo* coils, const FilGeo* fils,
    const double* gl_nodes, const double* gl_weights,
    volatile int* generation, volatile int* active_pairs, volatile int* shutdown,
    volatile int* shutdown_ack, volatile double* seps, volatile int* doorbell, volatile double* out_M,
    volatile double* out_dM, int n_s, int n_f, int n_n) {

    (void)active_pairs;
    int idx=blockIdx.x, si=idx/n_f, fi=idx%n_f;
    const auto& c=coils[si]; const auto& f=fils[fi];
    double ra_m=0.5*(c.re+c.ri), ra_h=0.5*(c.re-c.ri);
    double rb_m=0.5*(f.re+f.ri), rb_h=0.5*(f.re-f.ri);
    double la_h=0.5*c.len, lb_h=0.5*f.len, pref=c.turns/16.0;
    int n2=n_n*n_n, n4=n2*n2, tid=threadIdx.x, tot=blockDim.x;

    int last_generation = 0;
    while(*shutdown==0){
        while(*generation == last_generation && *shutdown == 0) {}
        if(*shutdown != 0) break;
        last_generation = *generation;
        while(doorbell[idx]==0 && *shutdown==0){}
        if(*shutdown!=0) break;
        double M=0,dM=0;
        if(idx < n_s * n_f){
            double sep=seps[idx];
            __shared__ double sM[512],sdM[512];
            sM[tid]=0;sdM[tid]=0;
            for(int pt=tid;pt<n4;pt+=tot){
                int i1=pt/(n_n*n_n*n_n);int r1=pt%(n_n*n_n*n_n);
                int j1=r1/(n_n*n_n);int r2=r1%(n_n*n_n);
                int i2=r2/n_n;int j2=r2%n_n;
                double w=gl_weights[i1]*gl_weights[j1]*gl_weights[i2]*gl_weights[j2];
                double ra=ra_m+ra_h*gl_nodes[i1],rb=rb_m+rb_h*gl_nodes[i2];
                double za=la_h*gl_nodes[j1],zb=sep+lb_h*gl_nodes[j2];
                const auto pair=physics::mutual_inductance_filament_pair_device(ra,rb,zb-za);
                sM[tid]+=w*pair.mutual;
                sdM[tid]+=w*pair.gradient;
            }
            __syncthreads();
            for(int s=tot/2;s>0;s/=2){ if(tid<s){ sM[tid]+=sM[tid+s];sdM[tid]+=sdM[tid+s]; } __syncthreads(); }
            if(tid==0){ M=pref*sM[0];dM=pref*sdM[0]; }
        }
        if(tid==0){ out_M[idx]=M;out_dM[idx]=dM;__threadfence_system();doorbell[idx]=0; }
    }
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        *shutdown_ack = 1;
        __threadfence_system();
    }
}

// FP32 variant — uses float integrand with double accumulation
__global__ void persistent_batch_kernel_f32(
    const CoilGeo* coils, const FilGeo* fils,
    const double* gl_nodes_d, const double* gl_weights_d,
    volatile int* generation, volatile int* active_pairs, volatile int* shutdown,
    volatile int* shutdown_ack, volatile double* seps, volatile int* doorbell, volatile double* out_M,
    volatile double* out_dM, int n_s, int n_f, int n_n) {

    int idx=blockIdx.x, si=idx/n_f, fi=idx%n_f;
    const auto& c=coils[si]; const auto& f=fils[fi];
    float ra_m=0.5f*(c.re+c.ri), ra_h=0.5f*(c.re-c.ri);
    float rb_m=0.5f*(f.re+f.ri), rb_h=0.5f*(f.re-f.ri);
    float la_h=0.5f*c.len, lb_h=0.5f*f.len, pref=c.turns/16.0f;
    int n2=n_n*n_n, n4=n2*n2, tid=threadIdx.x, tot=blockDim.x;

    // Convert GL nodes/weights to float for faster access
    float gl_nodes[9], gl_weights[9];
    for(int i=0;i<9;i++){ gl_nodes[i]=(float)gl_nodes_d[i]; gl_weights[i]=(float)gl_weights_d[i]; }

    int last_generation = 0;
    while(*shutdown==0){
        while(*generation == last_generation && *shutdown == 0) {}
        if(*shutdown != 0) break;
        last_generation = *generation;
        while(doorbell[idx]==0 && *shutdown==0){}
        if(*shutdown!=0) break;
        double M=0,dM=0;
        if(idx < n_s * n_f){
            float sep=(float)seps[idx];
            __shared__ double sM[512],sdM[512];
            sM[tid]=0;sdM[tid]=0;
            for(int pt=tid;pt<n4;pt+=tot){
                int i1=pt/(n_n*n_n*n_n);int r1=pt%(n_n*n_n*n_n);
                int j1=r1/(n_n*n_n);int r2=r1%(n_n*n_n);
                int i2=r2/n_n;int j2=r2%n_n;
                float w=gl_weights[i1]*gl_weights[j1]*gl_weights[i2]*gl_weights[j2];
                float ra=ra_m+ra_h*gl_nodes[i1],rb=rb_m+rb_h*gl_nodes[i2];
                float za=la_h*gl_nodes[j1],zb=sep+lb_h*gl_nodes[j2];
                const auto pair=physics::mutual_inductance_filament_pair_f32(ra,rb,zb-za);
                sM[tid]+=w*pair.mutual;
                sdM[tid]+=w*pair.gradient;
            }
            __syncthreads();
            for(int s=tot/2;s>0;s/=2){ if(tid<s){ sM[tid]+=sM[tid+s];sdM[tid]+=sdM[tid+s]; } __syncthreads(); }
            if(tid==0){ M=pref*sM[0];dM=pref*sdM[0]; }
        }
        if(tid==0){ out_M[idx]=M;out_dM[idx]=dM;__threadfence_system();doorbell[idx]=0; }
    }
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        *shutdown_ack = 1;
        __threadfence_system();
    }
}

void launch_persistent_kernel(const PersistentBuffers& p, const GpuAdaptor& a, int N_b, int tpb, int nn,
                               GpuOptLevel opt_level) {
    cudaGetLastError(); if(tpb>512)tpb=512;
    if(opt_level == GpuOptLevel::Aggressive) {
        persistent_batch_kernel_f32<<<N_b,tpb>>>(a.d_coils(),a.d_fils(),a.d_nodes(),a.d_weights(),
            p.d_generation,p.d_active_pairs,p.d_shutdown,p.d_shutdown_ack,p.d_seps,p.d_doorbell,p.d_out_M,p.d_out_dM,
            a.n_stages(),N_b/a.n_stages(),nn);
    } else {
        persistent_batch_kernel<<<N_b,tpb>>>(a.d_coils(),a.d_fils(),a.d_nodes(),a.d_weights(),
            p.d_generation,p.d_active_pairs,p.d_shutdown,p.d_shutdown_ack,p.d_seps,p.d_doorbell,p.d_out_M,p.d_out_dM,
            a.n_stages(),N_b/a.n_stages(),nn);
    }
}

}
