/**
 * @file gpu_adaptor.cu
 * @brief GPU device buffer manager — implementation.
 * @author Winston Meursault
 */

#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/physics/quadrature.hpp"
#include <stdexcept>
#include <string>

namespace coilgun::simulation::cuda {

GpuAdaptor::~GpuAdaptor() { free_all(); }

GpuAdaptor::GpuAdaptor(GpuAdaptor&& other) noexcept
    : d_coils_(other.d_coils_), d_fils_(other.d_fils_)
    , d_nodes_(other.d_nodes_), d_weights_(other.d_weights_)
    , d_results_M_(other.d_results_M_), d_results_dM_(other.d_results_dM_)
    , d_seps_(other.d_seps_)
    , d_batch_seps_(other.d_batch_seps_)
    , d_batch_results_M_(other.d_batch_results_M_)
    , d_batch_results_dM_(other.d_batch_results_dM_)
    , n_stages_(other.n_stages_), n_fil_(other.n_fil_)
    , n_nodes_(other.n_nodes_), batch_size_(other.batch_size_) {
    other.d_coils_ = nullptr;
    other.d_fils_ = nullptr;
    other.d_nodes_ = nullptr;
    other.d_weights_ = nullptr;
    other.d_results_M_ = nullptr;
    other.d_results_dM_ = nullptr;
    other.d_seps_ = nullptr;
    other.d_batch_seps_ = nullptr;
    other.d_batch_results_M_ = nullptr;
    other.d_batch_results_dM_ = nullptr;
}

GpuAdaptor& GpuAdaptor::operator=(GpuAdaptor&& other) noexcept {
    if (this != &other) {
        free_all();
        d_coils_ = other.d_coils_; other.d_coils_ = nullptr;
        d_fils_ = other.d_fils_; other.d_fils_ = nullptr;
        d_nodes_ = other.d_nodes_; other.d_nodes_ = nullptr;
        d_weights_ = other.d_weights_; other.d_weights_ = nullptr;
        d_results_M_ = other.d_results_M_; other.d_results_M_ = nullptr;
        d_results_dM_ = other.d_results_dM_; other.d_results_dM_ = nullptr;
        d_seps_ = other.d_seps_; other.d_seps_ = nullptr;
        d_batch_seps_ = other.d_batch_seps_; other.d_batch_seps_ = nullptr;
        d_batch_results_M_ = other.d_batch_results_M_; other.d_batch_results_M_ = nullptr;
        d_batch_results_dM_ = other.d_batch_results_dM_; other.d_batch_results_dM_ = nullptr;
        n_stages_ = other.n_stages_; n_fil_ = other.n_fil_;
        n_nodes_ = other.n_nodes_; batch_size_ = other.batch_size_;
    }
    return *this;
}

void GpuAdaptor::free_all() {
    auto safe_free = [](void*& p) {
        if (p) { cudaFree(p); p = nullptr; }
    };
    safe_free((void*&)d_coils_);
    safe_free((void*&)d_fils_);
    safe_free((void*&)d_nodes_);
    safe_free((void*&)d_weights_);
    safe_free((void*&)d_results_M_);
    safe_free((void*&)d_results_dM_);
    safe_free((void*&)d_seps_);
    safe_free((void*&)d_batch_seps_);
    safe_free((void*&)d_batch_results_M_);
    safe_free((void*&)d_batch_results_dM_);
}

void GpuAdaptor::setup(const std::vector<components::DrivingCoil>& coils,
                       const components::Armature& armature,
                       int n_nodes) {
    n_stages_ = static_cast<int>(coils.size());
    n_fil_    = armature.total_filaments();
    n_nodes_  = n_nodes;

    std::vector<CoilGeo> h_coils(n_stages_);
    for (int i = 0; i < n_stages_; ++i) {
        h_coils[i] = {coils[i].inner_radius(), coils[i].outer_radius(),
                      coils[i].length(), coils[i].position(), coils[i].turns()};
    }

    int nr = armature.radial_filaments();
    int na = armature.axial_filaments();
    std::vector<FilGeo> h_fils(n_fil_);
    for (int k = 0; k < n_fil_; ++k) {
        int j = k % nr + 1;
        h_fils[k] = {armature.filament_inner_radius(j),
                     armature.filament_outer_radius(j),
                     armature.length() / static_cast<double>(na)};
    }

    const auto& gl = physics::gauss_legendre(n_nodes_);
    std::vector<double> h_nodes(gl.nodes.begin(), gl.nodes.end());
    std::vector<double> h_weights(gl.weights.begin(), gl.weights.end());

    size_t pair_count = static_cast<size_t>(n_stages_) * n_fil_;

    auto alloc = [](auto& d, size_t n, const char* name) {
        using ElemType = std::remove_pointer_t<std::decay_t<decltype(d)>>;
        cudaError_t err = cudaMalloc(&d, n * sizeof(ElemType));
        if (err != cudaSuccess)
            throw std::runtime_error(std::string("cudaMalloc failed: ") + name);
    };

    free_all();
    alloc(d_coils_,  n_stages_, "coils");
    alloc(d_fils_,   n_fil_,    "fils");
    alloc(d_nodes_,  n_nodes_,  "nodes");
    alloc(d_weights_, n_nodes_, "weights");
    alloc(d_results_M_,  pair_count, "results_M");
    alloc(d_results_dM_, pair_count, "results_dM");
    alloc(d_seps_,       pair_count, "seps");

    cudaMemcpy(d_coils_,  h_coils.data(),  n_stages_ * sizeof(CoilGeo), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fils_,   h_fils.data(),   n_fil_    * sizeof(FilGeo),  cudaMemcpyHostToDevice);
    cudaMemcpy(d_nodes_,  h_nodes.data(),  n_nodes_  * sizeof(double),  cudaMemcpyHostToDevice);
    cudaMemcpy(d_weights_, h_weights.data(), n_nodes_* sizeof(double),  cudaMemcpyHostToDevice);
}

void GpuAdaptor::setup_batch(const std::vector<components::DrivingCoil>& coils,
                             const components::Armature& armature,
                             int num_sims, int n_nodes) {
    batch_size_ = num_sims;
    n_stages_   = static_cast<int>(coils.size());
    n_fil_      = armature.total_filaments();
    n_nodes_    = n_nodes;

    std::vector<CoilGeo> h_coils(n_stages_);
    for (int i = 0; i < n_stages_; ++i) {
        h_coils[i] = {coils[i].inner_radius(), coils[i].outer_radius(),
                      coils[i].length(), coils[i].position(), coils[i].turns()};
    }

    int nr = armature.radial_filaments();
    int na = armature.axial_filaments();
    std::vector<FilGeo> h_fils(n_fil_);
    for (int k = 0; k < n_fil_; ++k) {
        int j = k % nr + 1;
        h_fils[k] = {armature.filament_inner_radius(j),
                     armature.filament_outer_radius(j),
                     armature.length() / static_cast<double>(na)};
    }

    const auto& gl = physics::gauss_legendre(n_nodes_);
    std::vector<double> h_nodes(gl.nodes.begin(), gl.nodes.end());
    std::vector<double> h_weights(gl.weights.begin(), gl.weights.end());

    size_t pair_count = static_cast<size_t>(n_stages_) * n_fil_;
    size_t batch_count = pair_count * batch_size_;

    auto alloc = [](auto& d, size_t n, const char* name) {
        using ElemType = std::remove_pointer_t<std::decay_t<decltype(d)>>;
        cudaError_t err = cudaMalloc(&d, n * sizeof(ElemType));
        if (err != cudaSuccess)
            throw std::runtime_error(std::string("cudaMalloc failed: ") + name);
    };

    free_all();
    alloc(d_coils_,  n_stages_, "coils");
    alloc(d_fils_,   n_fil_,    "fils");
    alloc(d_nodes_,  n_nodes_,  "nodes");
    alloc(d_weights_, n_nodes_, "weights");
    alloc(d_results_M_,  pair_count, "results_M");
    alloc(d_results_dM_, pair_count, "results_dM");
    alloc(d_seps_,       pair_count, "seps");

    alloc(d_batch_seps_,       batch_count, "batch_seps");
    alloc(d_batch_results_M_,  batch_count, "batch_results_M");
    alloc(d_batch_results_dM_, batch_count, "batch_results_dM");

    cudaMemcpy(d_coils_,  h_coils.data(),  n_stages_ * sizeof(CoilGeo), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fils_,   h_fils.data(),   n_fil_    * sizeof(FilGeo),  cudaMemcpyHostToDevice);
    cudaMemcpy(d_nodes_,  h_nodes.data(),  n_nodes_  * sizeof(double),  cudaMemcpyHostToDevice);
    cudaMemcpy(d_weights_, h_weights.data(), n_nodes_* sizeof(double),  cudaMemcpyHostToDevice);
}

void GpuAdaptor::upload_batch_separations(const std::vector<double>& seps) {
    cudaMemcpy(d_batch_seps_, seps.data(), seps.size() * sizeof(double),
               cudaMemcpyHostToDevice);
}

void GpuAdaptor::download_batch_results(std::vector<double>& M_out,
                                        std::vector<double>& dM_out) {
    size_t n = static_cast<size_t>(n_stages_) * n_fil_ * batch_size_;
    M_out.resize(n);
    dM_out.resize(n);
    cudaMemcpy(M_out.data(),  d_batch_results_M_,  n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dM_out.data(), d_batch_results_dM_, n * sizeof(double), cudaMemcpyDeviceToHost);
}

void GpuAdaptor::upload_separation(const std::vector<double>& seps) {
    cudaMemcpy(d_seps_, seps.data(), seps.size() * sizeof(double), cudaMemcpyHostToDevice);
}

void GpuAdaptor::download_results(std::vector<double>& M_out,
                                  std::vector<double>& dM_out,
                                  int n_pairs) {
    M_out.resize(n_pairs);
    dM_out.resize(n_pairs);
    cudaMemcpy(M_out.data(),  d_results_M_,  n_pairs * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dM_out.data(), d_results_dM_, n_pairs * sizeof(double), cudaMemcpyDeviceToHost);
}

} // namespace coilgun::simulation::cuda
