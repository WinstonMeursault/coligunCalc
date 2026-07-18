/**
 * @file gpu_adaptor.hpp
 * @brief Manages device memory for coilgun simulation geometry.
 * @author Winston Meursault
 *
 * Uploads invariant simulation geometry (coil dimensions, filament
 * discretisation, GL quadrature nodes) to device memory once at
 * initialisation. Provides per-step host-to-device upload of armature
 * position and device-to-host download of M and dM results.
 */

#pragma once

#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include <cuda_runtime.h>
#include <cstdint>
#include <vector>

namespace coilgun::simulation::cuda {

/**
 * @brief Packed coil geometry for device transfer.
 */
struct CoilGeo {
    double ri;     ///< Inner winding radius, m.
    double re;     ///< Outer winding radius, m.
    double len;    ///< Axial winding length, m.
    double pos;    ///< Centre position along barrel, m.
    int    turns;  ///< Total winding turns.
};

/**
 * @brief Packed filament geometry for device transfer.
 */
struct FilGeo {
    double ri;   ///< Inner radius, m.
    double re;   ///< Outer radius, m.
    double len;  ///< Axial length, m.
};

/**
 * @brief GPU device buffer manager.
 *
 * Owns all device allocations. Move-only to prevent double-free.
 */
class GpuAdaptor {
public:
    GpuAdaptor() = default;
    ~GpuAdaptor();

    GpuAdaptor(const GpuAdaptor&) = delete;
    GpuAdaptor& operator=(const GpuAdaptor&) = delete;
    GpuAdaptor(GpuAdaptor&&) noexcept;
    GpuAdaptor& operator=(GpuAdaptor&&) noexcept;

    /// Set up device buffers for a single simulation.
    void setup(const std::vector<components::DrivingCoil>& coils,
               const components::Armature& armature,
               int n_nodes);

    /// Set up device buffers for a batch of geometrically identical simulations.
    void setup_batch(const std::vector<components::DrivingCoil>& coils,
                     const components::Armature& armature,
                     int num_sims, int n_nodes);

    /// Upload per-step separation vector for all (stage, filament) pairs.
    void upload_separation(const std::vector<double>& seps);

    /// Download M and dM results.
    void download_results(std::vector<double>& M_out,
                          std::vector<double>& dM_out,
                          int n_pairs);

    /// Upload per-step separation vectors for all simulations in batch mode.
    /// @param seps Flat vector [num_sims][n_stages * n_fil], sim-major.
    void upload_batch_separations(const std::vector<double>& seps);

    /// Download M and dM results for all simulations in batch mode.
    /// @param M_out,dM_out Flat vectors [num_sims][n_pairs], sim-major.
    void download_batch_results(std::vector<double>& M_out,
                                std::vector<double>& dM_out);

    /// @name Device pointers
    /// @{
    const CoilGeo*  d_coils()   const { return d_coils_; }
    const FilGeo*   d_fils()    const { return d_fils_; }
    const double*   d_nodes()    const { return d_nodes_; }
    const double*   d_weights()  const { return d_weights_; }
    double*         d_results_M()  { return d_results_M_; }
    double*         d_results_dM() { return d_results_dM_; }
    double*         d_seps()       { return d_seps_; }
    double*         d_batch_seps()        { return d_batch_seps_; }
    double*         d_batch_results_M()   { return d_batch_results_M_; }
    double*         d_batch_results_dM()  { return d_batch_results_dM_; }
    /// @}

    /// @name Dimensions
    /// @{
    int n_stages() const { return n_stages_; }
    int n_fil()    const { return n_fil_; }
    int n_nodes()  const { return n_nodes_; }
    int batch_size() const { return batch_size_; }
    /// @}

private:
    void free_all();

    CoilGeo* d_coils_    = nullptr;
    FilGeo*  d_fils_     = nullptr;
    double*  d_nodes_    = nullptr;
    double*  d_weights_  = nullptr;
    double*  d_results_M_  = nullptr;
    double*  d_results_dM_ = nullptr;
    double*  d_seps_       = nullptr;

    double*  d_batch_seps_       = nullptr;
    double*  d_batch_results_M_  = nullptr;
    double*  d_batch_results_dM_ = nullptr;

    int n_stages_  = 0;
    int n_fil_     = 0;
    int n_nodes_   = 9;
    int batch_size_ = 1;
};

} // namespace coilgun::simulation::cuda
