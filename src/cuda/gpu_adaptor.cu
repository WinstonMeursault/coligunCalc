#include "coilgun/simulation/cuda/gpu_adaptor.hpp"

#include "coilgun/physics/quadrature.hpp"

#include <cuda_runtime.h>

#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace coilgun::simulation::cuda {
namespace {

void check_cuda(cudaError_t error, const char* operation) {
    if (error != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " +
                                 cudaGetErrorString(error));
}

std::size_t checked_product(std::size_t left, std::size_t right,
                            const char* operation) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
        throw std::invalid_argument(std::string(operation) + " size overflow");
    return left * right;
}

template<typename Pointer>
void allocate(Pointer& pointer, std::size_t count, const char* name) {
    if (count == 0) throw std::invalid_argument(std::string(name) + " count must be positive");
    using Element = std::remove_pointer_t<Pointer>;
    check_cuda(cudaMalloc(reinterpret_cast<void**>(&pointer),
                          checked_product(count, sizeof(Element), name)),
               name);
}

void require_device(int expected, const char* operation) {
    int current = -1;
    check_cuda(cudaGetDevice(&current), "cudaGetDevice");
    if (current != expected)
        throw std::invalid_argument(std::string(operation) +
                                    " used on a different CUDA device");
}

} // namespace

GpuAdaptor::~GpuAdaptor() { free_all(); }

GpuAdaptor::GpuAdaptor(GpuAdaptor&& other) noexcept
    : d_coils_(other.d_coils_), d_fils_(other.d_fils_),
      d_nodes_(other.d_nodes_), d_weights_(other.d_weights_),
      d_results_M_(other.d_results_M_), d_results_dM_(other.d_results_dM_),
      d_seps_(other.d_seps_), d_batch_seps_(other.d_batch_seps_),
      d_batch_results_M_(other.d_batch_results_M_),
      d_batch_results_dM_(other.d_batch_results_dM_),
      n_stages_(other.n_stages_), n_fil_(other.n_fil_),
      n_nodes_(other.n_nodes_), batch_size_(other.batch_size_),
      configured_(other.configured_), device_id_(other.device_id_),
      single_pair_capacity_(other.single_pair_capacity_),
      batch_pair_capacity_(other.batch_pair_capacity_) {
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
    other.configured_ = false;
    other.device_id_ = -1;
    other.single_pair_capacity_ = 0;
    other.batch_pair_capacity_ = 0;
}

GpuAdaptor& GpuAdaptor::operator=(GpuAdaptor&& other) noexcept {
    if (this == &other) return *this;
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
    n_stages_ = other.n_stages_;
    n_fil_ = other.n_fil_;
    n_nodes_ = other.n_nodes_;
    batch_size_ = other.batch_size_;
    configured_ = other.configured_;
    device_id_ = other.device_id_;
    single_pair_capacity_ = other.single_pair_capacity_;
    batch_pair_capacity_ = other.batch_pair_capacity_;
    other.configured_ = false;
    other.device_id_ = -1;
    other.single_pair_capacity_ = 0;
    other.batch_pair_capacity_ = 0;
    return *this;
}

void GpuAdaptor::free_all() {
    if (device_id_ >= 0) cudaSetDevice(device_id_);
    auto release = [](auto*& pointer) {
        if (pointer) cudaFree(pointer);
        pointer = nullptr;
    };
    release(d_coils_);
    release(d_fils_);
    release(d_nodes_);
    release(d_weights_);
    release(d_results_M_);
    release(d_results_dM_);
    release(d_seps_);
    release(d_batch_seps_);
    release(d_batch_results_M_);
    release(d_batch_results_dM_);
    configured_ = false;
    device_id_ = -1;
    single_pair_capacity_ = 0;
    batch_pair_capacity_ = 0;
}

void GpuAdaptor::setup(const std::vector<components::DrivingCoil>& coils,
                       const components::Armature& armature,
                       int n_nodes) {
    if (coils.empty()) throw std::invalid_argument("GpuAdaptor coils must not be empty");
    if (armature.total_filaments() <= 0)
        throw std::invalid_argument("GpuAdaptor filament count must be positive");
    if (n_nodes != 9)
        throw std::invalid_argument("CUDA mutual pipeline requires n_nodes == 9");

    int current_device = -1;
    check_cuda(cudaGetDevice(&current_device), "cudaGetDevice(setup)");
    const auto stages = coils.size();
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    const auto pair_count = checked_product(stages, filaments, "single pair capacity");

    std::vector<CoilGeo> host_coils(stages);
    for (std::size_t stage = 0; stage < stages; ++stage) {
        host_coils[stage] = {coils[stage].inner_radius(), coils[stage].outer_radius(),
                             coils[stage].length(), coils[stage].position(),
                             coils[stage].turns()};
    }
    std::vector<FilGeo> host_filaments(filaments);
    for (std::size_t filament = 0; filament < filaments; ++filament) {
        const int radial = static_cast<int>(filament) % armature.radial_filaments() + 1;
        host_filaments[filament] = {
            armature.filament_inner_radius(radial),
            armature.filament_outer_radius(radial),
            armature.length() / static_cast<double>(armature.axial_filaments())};
    }
    const auto& quadrature = physics::gauss_legendre_cached(9);

    free_all();
    device_id_ = current_device;
    n_stages_ = static_cast<int>(stages);
    n_fil_ = static_cast<int>(filaments);
    n_nodes_ = 9;
    batch_size_ = 1;
    single_pair_capacity_ = pair_count;
    try {
        allocate(d_coils_, stages, "cudaMalloc(coils)");
        allocate(d_fils_, filaments, "cudaMalloc(filaments)");
        allocate(d_nodes_, 9, "cudaMalloc(nodes)");
        allocate(d_weights_, 9, "cudaMalloc(weights)");
        allocate(d_results_M_, pair_count, "cudaMalloc(mutual results)");
        allocate(d_results_dM_, pair_count, "cudaMalloc(gradient results)");
        allocate(d_seps_, pair_count, "cudaMalloc(separations)");
        check_cuda(cudaMemcpy(d_coils_, host_coils.data(), stages * sizeof(CoilGeo),
                              cudaMemcpyHostToDevice), "cudaMemcpy(coils)");
        check_cuda(cudaMemcpy(d_fils_, host_filaments.data(), filaments * sizeof(FilGeo),
                              cudaMemcpyHostToDevice), "cudaMemcpy(filaments)");
        check_cuda(cudaMemcpy(d_nodes_, quadrature.nodes.data(), 9 * sizeof(double),
                              cudaMemcpyHostToDevice), "cudaMemcpy(nodes)");
        check_cuda(cudaMemcpy(d_weights_, quadrature.weights.data(), 9 * sizeof(double),
                              cudaMemcpyHostToDevice), "cudaMemcpy(weights)");
        configured_ = true;
    } catch (...) {
        free_all();
        throw;
    }
}

void GpuAdaptor::setup_batch(const std::vector<components::DrivingCoil>& coils,
                             const components::Armature& armature,
                             int num_sims, int n_nodes) {
    if (num_sims <= 0) throw std::invalid_argument("GpuAdaptor batch size must be positive");
    setup(coils, armature, n_nodes);
    batch_size_ = num_sims;
    batch_pair_capacity_ = checked_product(
        single_pair_capacity_, static_cast<std::size_t>(num_sims),
        "batch pair capacity");
    try {
        allocate(d_batch_seps_, batch_pair_capacity_, "cudaMalloc(batch separations)");
        allocate(d_batch_results_M_, batch_pair_capacity_, "cudaMalloc(batch mutual results)");
        allocate(d_batch_results_dM_, batch_pair_capacity_, "cudaMalloc(batch gradient results)");
    } catch (...) {
        free_all();
        throw;
    }
}

void GpuAdaptor::upload_separation(const std::vector<double>& separations) {
    if (!configured_) throw std::invalid_argument("GpuAdaptor is not configured");
    require_device(device_id_, "upload_separation");
    if (separations.size() != single_pair_capacity_)
        throw std::invalid_argument("single separation count does not match configured capacity");
    check_cuda(cudaMemcpy(d_seps_, separations.data(),
                          separations.size() * sizeof(double), cudaMemcpyHostToDevice),
               "cudaMemcpy(single separations)");
}

void GpuAdaptor::download_results(std::vector<double>& mutual,
                                  std::vector<double>& gradient,
                                  int pair_count) {
    if (pair_count < 0) throw std::invalid_argument("pair_count must be non-negative");
    if (!configured_) throw std::invalid_argument("GpuAdaptor is not configured");
    require_device(device_id_, "download_results");
    if (static_cast<std::size_t>(pair_count) != single_pair_capacity_)
        throw std::invalid_argument("result count does not match configured capacity");
    mutual.resize(single_pair_capacity_);
    gradient.resize(single_pair_capacity_);
    check_cuda(cudaMemcpy(mutual.data(), d_results_M_,
                          single_pair_capacity_ * sizeof(double), cudaMemcpyDeviceToHost),
               "cudaMemcpy(single mutual results)");
    check_cuda(cudaMemcpy(gradient.data(), d_results_dM_,
                          single_pair_capacity_ * sizeof(double), cudaMemcpyDeviceToHost),
               "cudaMemcpy(single gradient results)");
}

void GpuAdaptor::upload_batch_separations(const std::vector<double>& separations) {
    if (!configured_ || batch_pair_capacity_ == 0)
        throw std::invalid_argument("GpuAdaptor batch mode is not configured");
    require_device(device_id_, "upload_batch_separations");
    if (separations.size() != batch_pair_capacity_)
        throw std::invalid_argument("batch separation count does not match configured capacity");
    check_cuda(cudaMemcpy(d_batch_seps_, separations.data(),
                          batch_pair_capacity_ * sizeof(double), cudaMemcpyHostToDevice),
               "cudaMemcpy(batch separations)");
}

void GpuAdaptor::download_batch_results(std::vector<double>& mutual,
                                        std::vector<double>& gradient) {
    if (!configured_ || batch_pair_capacity_ == 0)
        throw std::invalid_argument("GpuAdaptor batch mode is not configured");
    require_device(device_id_, "download_batch_results");
    mutual.resize(batch_pair_capacity_);
    gradient.resize(batch_pair_capacity_);
    check_cuda(cudaMemcpy(mutual.data(), d_batch_results_M_,
                          batch_pair_capacity_ * sizeof(double), cudaMemcpyDeviceToHost),
               "cudaMemcpy(batch mutual results)");
    check_cuda(cudaMemcpy(gradient.data(), d_batch_results_dM_,
                          batch_pair_capacity_ * sizeof(double), cudaMemcpyDeviceToHost),
               "cudaMemcpy(batch gradient results)");
}

} // namespace coilgun::simulation::cuda
