/**
 * @file gpu_mutual_inductance.cu
 * @brief 4D Gauss-Legendre integration kernel for (coil, filament) pair.
 * @author Winston Meursault
 *
 * Each block processes one (driving_coil, armature_filament) pair.
 * 512 threads per block, each thread loops over its share of the
 * n_nodes^4 integration points, with shared-memory tree reduction.
 */

#include "coilgun/physics/mutual_inductance.cuh"
#include "coilgun/physics/quadrature.hpp"
#include <cuda_runtime.h>
#include <cmath>

namespace coilgun::physics {

/**
 * @brief Device-side constant memory for GL nodes.
 *
 * Uploaded once per simulation via upload_gl_nodes().
 */
__constant__ double d_gl_nodes[9];
__constant__ double d_gl_weights[9];

/**
 * @brief Upload Gauss-Legendre nodes and weights to device constant memory.
 * @param n_nodes Number of GL nodes per dimension (max 9).
 */
void upload_gl_nodes(int n_nodes) {
    const auto& gl = gauss_legendre(n_nodes);
    cudaMemcpyToSymbol(d_gl_nodes, gl.nodes.data(), n_nodes * sizeof(double));
    cudaMemcpyToSymbol(d_gl_weights, gl.weights.data(), n_nodes * sizeof(double));
}

/**
 * @brief Map [-1,1] coordinate to physical domain — device version.
 */
__device__ inline double map_coord_device(double centre, double half_width, double xi) {
    return centre + half_width * xi;
}

/**
 * @brief Kernel: compute mutual inductance and gradient for one (coil, filament) pair.
 *
 * Uses 4D Gauss-Legendre tensor-product quadrature. Each thread loops
 * over a subset of the n_nodes^4 integration points. Results are summed
 * via shared-memory tree reduction. Both M and dM/dz are computed
 * simultaneously to avoid redundant elliptic integral evaluations.
 *
 * @param rai Coil inner radius, m.
 * @param rae Coil outer radius, m.
 * @param la  Coil axial length, m.
 * @param na  Coil turns.
 * @param rbi Filament inner radius, m.
 * @param rbe Filament outer radius, m.
 * @param lb  Filament axial length, m.
 * @param nb  Filament turns (1 for armature filament).
 * @param separation Centre-to-centre axial distance between coil and filament, m.
 * @param n_nodes Number of GL nodes per dimension.
 * @param out_M Output pointer for mutual inductance (blockIdx.x).
 * @param out_dM Output pointer for M gradient (blockIdx.x).
 */
__global__ void mutual_inductance_coil_pair_kernel(
        double rai, double rae, double la, int na,
        double rbi, double rbe, double lb, int nb,
        double separation, int n_nodes,
        double* out_M, double* out_dM) {

    const double ra_mid  = 0.5 * (rae + rai);
    const double ra_half = 0.5 * (rae - rai);
    const double rb_mid  = 0.5 * (rbe + rbi);
    const double rb_half = 0.5 * (rbe - rbi);
    const double la_half = 0.5 * la;
    const double lb_half = 0.5 * lb;
    const double prefactor = (na * nb) / 16.0;

    int n2 = n_nodes * n_nodes;
    int n4 = n2 * n2;
    int tid = threadIdx.x;
    int total_threads = blockDim.x;

    __shared__ double sM[512];
    __shared__ double sdM[512];
    sM[tid]  = 0.0;
    sdM[tid] = 0.0;

    for (int idx = tid; idx < n4; idx += total_threads) {
        int i1 = idx / (n_nodes * n_nodes * n_nodes);
        int rem1 = idx % (n_nodes * n_nodes * n_nodes);
        int j1 = rem1 / (n_nodes * n_nodes);
        int rem2 = rem1 % (n_nodes * n_nodes);
        int i2 = rem2 / n_nodes;
        int j2 = rem2 % n_nodes;

        double w = d_gl_weights[i1] * d_gl_weights[j1]
                 * d_gl_weights[i2] * d_gl_weights[j2];

        double ra = map_coord_device(ra_mid, ra_half, d_gl_nodes[i1]);
        double rb = map_coord_device(rb_mid, rb_half, d_gl_nodes[i2]);
        double za = map_coord_device(0.0, la_half, d_gl_nodes[j1]);
        double zb = map_coord_device(separation, lb_half, d_gl_nodes[j2]);

        double abs_sep = fabs(zb - za);
        sM[tid]  += w * mutual_inductance_filament_device(ra, rb, abs_sep);
        sdM[tid] += w * mutual_inductance_gradient_filament_device(ra, rb, zb - za);
    }

    __syncthreads();

    for (int stride = total_threads / 2; stride > 0; stride /= 2) {
        if (tid < stride) {
            sM[tid]  += sM[tid + stride];
            sdM[tid] += sdM[tid + stride];
        }
        __syncthreads();
    }

    if (tid == 0) {
        out_M[blockIdx.x]  = prefactor * sM[0];
        out_dM[blockIdx.x] = prefactor * sdM[0];
    }
}

} // namespace coilgun::physics
