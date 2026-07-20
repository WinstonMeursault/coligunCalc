#include "coilgun/simulation/cuda/gpu_thermal.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <memory>
#include <cstdint>
#include <limits>

namespace coilgun::simulation::cuda {
namespace {

constexpr double T_REFERENCE = 293.0;

__host__ __device__ double cp_value(int material, double temperature) {
    return material == static_cast<int>(ThermalMaterial::Copper)
        ? 333.0 * exp(3.917e-4 * temperature)
        : 819.0 * exp(3.07e-4 * temperature);
}

__host__ __device__ double rho_value(int material, double temperature) {
    return material == static_cast<int>(ThermalMaterial::Copper)
        ? -3.5e-9 + 7.2e-11 * temperature
        : -6.57e-9 + 1.2e-10 * temperature;
}

template <typename T>
__device__ T table_interpolate(const double* temperatures, const double* values,
                               std::size_t count, double minimum, double maximum,
                               T temperature, bool log_linear) {
    const T clamped = max(static_cast<T>(minimum), min(static_cast<T>(maximum), temperature));
    const T position = (clamped - static_cast<T>(minimum)) /
                       (static_cast<T>(maximum - minimum)) * static_cast<T>(count - 1);
    const std::size_t lower = min(count - 2, static_cast<std::size_t>(position));
    const T fraction = position - static_cast<T>(lower);
    T left = static_cast<T>(values[lower]);
    T right = static_cast<T>(values[lower + 1]);
    if (log_linear) {
        left = log(left);
        right = log(right);
        return exp(left + fraction * (right - left));
    }
    return left + fraction * (right - left);
}

__global__ void thermal_kernel(const double* table_t, const double* cp_al, const double* cp_cu,
                               const double* rho_al, const double* rho_cu, std::size_t table_count,
                               double minimum, double maximum, int precision, std::size_t count,
                               const double* currents, const double* masses,
                               const double* reference_resistances, const int* materials, double dt,
                           double* temperatures, double* resistivities, double* resistances,
                                double* joule_energy) {
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    const int material = materials[index];
    const double old_temperature = temperatures[index];
    const double old_resistance = resistances[index] == 0.0
        ? reference_resistances[index] : resistances[index];
    const double current = currents[index];
    const double joule = current * current * old_resistance * dt;

    if (precision == static_cast<int>(ThermalPrecision::Aggressive)) {
        const float t = static_cast<float>(old_temperature);
        const float i = static_cast<float>(current);
        const float m = static_cast<float>(masses[index]);
        const float r = static_cast<float>(old_resistance);
        const float cp = table_interpolate(table_t, material ? cp_cu : cp_al, table_count,
                                           minimum, maximum, t, false);
        const float next_t = t + i * i * r * static_cast<float>(dt) / (m * cp);
        temperatures[index] = next_t;
        resistivities[index] = static_cast<float>(table_interpolate(
            table_t, material ? rho_cu : rho_al, table_count, minimum, maximum, next_t, false));
        resistances[index] = static_cast<float>(reference_resistances[index] *
            (1.0f + (material ? 4.1e-3f : 4.2e-3f) * (next_t - static_cast<float>(T_REFERENCE))));
        joule_energy[index] = i * i * r * static_cast<float>(dt);
        return;
    }

    const bool standard = precision == static_cast<int>(ThermalPrecision::Standard);
    const double cp = standard
        ? cp_value(material, old_temperature)
        : table_interpolate(table_t, material ? cp_cu : cp_al, table_count,
                            minimum, maximum, old_temperature, false);
    const double next_temperature = old_temperature +
        current * current * old_resistance * dt / (masses[index] * cp);
    temperatures[index] = next_temperature;
    resistivities[index] = table_interpolate(table_t, material ? rho_cu : rho_al, table_count,
                                             minimum, maximum, next_temperature, false);
    resistances[index] = reference_resistances[index] *
        (1.0 + (material ? 4.1e-3 : 4.2e-3) * (next_temperature - T_REFERENCE));
    joule_energy[index] = joule;
}

void check_cuda(cudaError_t error, const char* operation) {
    if (error != cudaSuccess) throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(error));
}

template <typename T>
void alloc_copy(T*& destination, const T* source, std::size_t count) {
    check_cuda(cudaMalloc(reinterpret_cast<void**>(&destination), count * sizeof(T)), "cudaMalloc");
    check_cuda(cudaMemcpy(destination, source, count * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy");
}

} // namespace

struct ThermalWorkspace::Impl {
    double *t = nullptr, *cp_al = nullptr, *cp_cu = nullptr, *rho_al = nullptr, *rho_cu = nullptr;
    double *i = nullptr, *m = nullptr, *r0 = nullptr, *temp = nullptr, *rho = nullptr, *r = nullptr, *q = nullptr;
    int* material = nullptr;
    std::size_t table_count = 0, value_count = 0;
    ~Impl() {
        cudaFree(t); cudaFree(cp_al); cudaFree(cp_cu); cudaFree(rho_al); cudaFree(rho_cu);
        cudaFree(i); cudaFree(m); cudaFree(r0); cudaFree(temp); cudaFree(rho); cudaFree(r); cudaFree(q);
        cudaFree(material);
    }
};

ThermalWorkspace::ThermalWorkspace(const MaterialTables& tables, std::size_t count) { initialize(tables, count); }
ThermalWorkspace::~ThermalWorkspace() = default;

void ThermalWorkspace::initialize(const MaterialTables& tables, std::size_t count) {
    if (tables.temperatures.size() < 2 || tables.cp_aluminum.size() != tables.temperatures.size() ||
        tables.cp_copper.size() != tables.temperatures.size() ||
        tables.rho_aluminum.size() != tables.temperatures.size() ||
        tables.rho_copper.size() != tables.temperatures.size() ||
        count == 0 || !(tables.minimum_temperature < tables.maximum_temperature))
        throw std::invalid_argument("invalid thermal workspace dimensions or material tables");
    if (!impl_) impl_ = std::make_unique<Impl>();
    if (impl_->table_count == tables.temperatures.size() && impl_->value_count == count) return;
    impl_.reset(new Impl{});
    const auto n = tables.temperatures.size();
    alloc_copy(impl_->t, tables.temperatures.data(), n);
    alloc_copy(impl_->cp_al, tables.cp_aluminum.data(), n);
    alloc_copy(impl_->cp_cu, tables.cp_copper.data(), n);
    alloc_copy(impl_->rho_al, tables.rho_aluminum.data(), n);
    alloc_copy(impl_->rho_cu, tables.rho_copper.data(), n);
    alloc_copy(impl_->i, std::vector<double>(count).data(), count);
    alloc_copy(impl_->m, std::vector<double>(count).data(), count);
    alloc_copy(impl_->r0, std::vector<double>(count).data(), count);
    alloc_copy(impl_->material, std::vector<int>(count).data(), count);
    check_cuda(cudaMalloc(&impl_->temp, count * sizeof(double)), "cudaMalloc");
    check_cuda(cudaMalloc(&impl_->rho, count * sizeof(double)), "cudaMalloc");
    check_cuda(cudaMalloc(&impl_->r, count * sizeof(double)), "cudaMalloc");
    check_cuda(cudaMalloc(&impl_->q, count * sizeof(double)), "cudaMalloc");
    impl_->table_count = n;
    impl_->value_count = count;
    allocation_count_ = n * 5 + count * 8;
}

void ThermalWorkspace::update(const MaterialTables& tables, ThermalPrecision precision,
                              std::size_t batch_count, std::size_t filament_count,
                              const double* currents, const double* masses,
                              const double* reference_resistances, const int* materials,
                              double dt, double* temperatures, double* resistivities,
                              double* resistances, double* joule_energy, cudaStream_t stream) {
    if (!batch_count || !filament_count || !currents || !masses || !reference_resistances ||
        !materials || !temperatures || !resistivities || !resistances || !joule_energy ||
        filament_count > std::numeric_limits<std::size_t>::max() / batch_count ||
        !std::isfinite(dt) || dt <= 0.0)
        throw std::invalid_argument("invalid thermal workspace arguments");
    const auto count = batch_count * filament_count;
    if (count > std::numeric_limits<std::size_t>::max() - 255)
        throw std::invalid_argument("thermal launch size overflow");
    cudaDeviceProp properties{};
    int device = 0;
    if (cudaGetDevice(&device) != cudaSuccess || cudaGetDeviceProperties(&properties, device) != cudaSuccess ||
        (count + 255) / 256 > static_cast<std::size_t>(properties.maxGridSize[0]))
        throw std::invalid_argument("thermal launch exceeds device grid limits");
    for (std::size_t i = 0; i < count; ++i) {
        if (materials[i] != 0 && materials[i] != 1 ||
            !std::isfinite(currents[i]) || !std::isfinite(masses[i]) || masses[i] <= 0.0 ||
            !std::isfinite(reference_resistances[i]) || reference_resistances[i] <= 0.0 ||
            !std::isfinite(temperatures[i]) || !std::isfinite(resistances[i]) || resistances[i] < 0.0)
            throw std::invalid_argument("invalid thermal workspace state");
    }
    initialize(tables, count);
    const auto n = tables.temperatures.size();
    check_cuda(cudaMemcpyAsync(impl_->i, currents, count * sizeof(double), cudaMemcpyHostToDevice, stream), "current copy");
    check_cuda(cudaMemcpyAsync(impl_->m, masses, count * sizeof(double), cudaMemcpyHostToDevice, stream), "mass copy");
    check_cuda(cudaMemcpyAsync(impl_->r0, reference_resistances, count * sizeof(double), cudaMemcpyHostToDevice, stream), "resistance copy");
    check_cuda(cudaMemcpyAsync(impl_->material, materials, count * sizeof(int), cudaMemcpyHostToDevice, stream), "material copy");
    check_cuda(cudaMemcpyAsync(impl_->temp, temperatures, count * sizeof(double), cudaMemcpyHostToDevice, stream), "temperature copy");
    check_cuda(cudaMemcpyAsync(impl_->r, resistances, count * sizeof(double), cudaMemcpyHostToDevice, stream), "state resistance copy");
    thermal_kernel<<<static_cast<unsigned>((count + 255) / 256), 256, 0, stream>>>(impl_->t, impl_->cp_al, impl_->cp_cu,
        impl_->rho_al, impl_->rho_cu, n, tables.minimum_temperature, tables.maximum_temperature,
        static_cast<int>(precision), count, impl_->i, impl_->m, impl_->r0, impl_->material, dt,
        impl_->temp, impl_->rho, impl_->r, impl_->q);
    check_cuda(cudaGetLastError(), "thermal kernel");
    check_cuda(cudaMemcpyAsync(temperatures, impl_->temp, count * sizeof(double), cudaMemcpyDeviceToHost, stream), "temperature result");
    check_cuda(cudaMemcpyAsync(resistivities, impl_->rho, count * sizeof(double), cudaMemcpyDeviceToHost, stream), "rho result");
    check_cuda(cudaMemcpyAsync(resistances, impl_->r, count * sizeof(double), cudaMemcpyDeviceToHost, stream), "resistance result");
    check_cuda(cudaMemcpyAsync(joule_energy, impl_->q, count * sizeof(double), cudaMemcpyDeviceToHost, stream), "joule result");
    check_cuda(cudaStreamSynchronize(stream), "thermal synchronize");
}

std::vector<std::uintptr_t> ThermalWorkspace::device_addresses() const {
    if (!impl_) return {};
    return {reinterpret_cast<std::uintptr_t>(impl_->t), reinterpret_cast<std::uintptr_t>(impl_->cp_al),
            reinterpret_cast<std::uintptr_t>(impl_->cp_cu), reinterpret_cast<std::uintptr_t>(impl_->rho_al),
            reinterpret_cast<std::uintptr_t>(impl_->rho_cu), reinterpret_cast<std::uintptr_t>(impl_->i),
            reinterpret_cast<std::uintptr_t>(impl_->m), reinterpret_cast<std::uintptr_t>(impl_->r0),
            reinterpret_cast<std::uintptr_t>(impl_->temp), reinterpret_cast<std::uintptr_t>(impl_->rho),
            reinterpret_cast<std::uintptr_t>(impl_->r), reinterpret_cast<std::uintptr_t>(impl_->q),
            reinterpret_cast<std::uintptr_t>(impl_->material)};
}

MaterialTables generate_material_tables(std::size_t sample_count, double minimum_temperature,
                                        double maximum_temperature) {
    if (sample_count < 2 || !(minimum_temperature < maximum_temperature))
        throw std::invalid_argument("material table requires at least two ordered temperatures");
    MaterialTables result;
    result.minimum_temperature = minimum_temperature;
    result.maximum_temperature = maximum_temperature;
    result.temperatures.resize(sample_count);
    result.cp_aluminum.resize(sample_count);
    result.cp_copper.resize(sample_count);
    result.rho_aluminum.resize(sample_count);
    result.rho_copper.resize(sample_count);
    for (std::size_t i = 0; i < sample_count; ++i) {
        const double t = minimum_temperature + (maximum_temperature - minimum_temperature) *
                         static_cast<double>(i) / static_cast<double>(sample_count - 1);
        result.temperatures[i] = t;
        result.cp_aluminum[i] = cp_value(0, t);
        result.cp_copper[i] = cp_value(1, t);
        result.rho_aluminum[i] = rho_value(0, t);
        result.rho_copper[i] = rho_value(1, t);
    }
    return result;
}

double interpolate_material_cp(const MaterialTables& tables, ThermalMaterial material,
                               double temperature, ThermalPrecision precision) {
    if (precision == ThermalPrecision::Standard)
        return cp_value(static_cast<int>(material), temperature);
    const auto& values = material == ThermalMaterial::Copper ? tables.cp_copper : tables.cp_aluminum;
    const double t = std::clamp(temperature, tables.minimum_temperature, tables.maximum_temperature);
    const double position = (t - tables.minimum_temperature) /
                            (tables.maximum_temperature - tables.minimum_temperature) * (values.size() - 1);
    const std::size_t lower = std::min(values.size() - 2, static_cast<std::size_t>(position));
    return values[lower] + (position - lower) * (values[lower + 1] - values[lower]);
}

double interpolate_material_resistivity(const MaterialTables& tables, ThermalMaterial material,
                                        double temperature, ThermalPrecision precision) {
    if (precision == ThermalPrecision::Standard)
        return rho_value(static_cast<int>(material), temperature);
    const auto& values = material == ThermalMaterial::Copper ? tables.rho_copper : tables.rho_aluminum;
    const double t = std::clamp(temperature, tables.minimum_temperature, tables.maximum_temperature);
    const double position = (t - tables.minimum_temperature) /
                            (tables.maximum_temperature - tables.minimum_temperature) * (values.size() - 1);
    const std::size_t lower = std::min(values.size() - 2, static_cast<std::size_t>(position));
    return values[lower] + (position - lower) * (values[lower + 1] - values[lower]);
}

void update_thermal_batch(const MaterialTables& tables, ThermalPrecision precision,
                          std::size_t batch_count, std::size_t filament_count, const double* currents,
                          const double* masses, const double* reference_resistances, const int* materials,
                           double dt, double* temperatures, double* resistivities, double* resistances,
                           double* joule_energy, cudaStream_t stream) {
    if (!batch_count || !filament_count || !currents || !masses || !reference_resistances || !materials ||
        !temperatures || !resistivities || !resistances || !joule_energy || tables.temperatures.size() < 2 ||
        !std::isfinite(dt) || dt <= 0.0 || filament_count > std::numeric_limits<std::size_t>::max() / batch_count ||
        tables.cp_aluminum.size() != tables.temperatures.size() ||
        tables.cp_copper.size() != tables.temperatures.size() ||
        tables.rho_aluminum.size() != tables.temperatures.size() ||
        tables.rho_copper.size() != tables.temperatures.size())
        throw std::invalid_argument("invalid thermal batch arguments");
    for (std::size_t i = 0; i < batch_count * filament_count; ++i) {
        if (materials[i] != 0 && materials[i] != 1) throw std::invalid_argument("invalid thermal material");
        if (!std::isfinite(currents[i]) || !std::isfinite(masses[i]) || masses[i] <= 0.0 ||
            !std::isfinite(reference_resistances[i]) || reference_resistances[i] <= 0.0 ||
            !std::isfinite(temperatures[i]) || !std::isfinite(resistances[i]) || resistances[i] < 0.0)
            throw std::invalid_argument("invalid thermal state");
    }
    static thread_local ThermalWorkspace workspace;
    workspace.update(tables, precision, batch_count, filament_count, currents, masses,
                     reference_resistances, materials, dt, temperatures, resistivities,
                     resistances, joule_energy, stream);
    return;
}

void update_thermal_batch_cpu(const MaterialTables& tables, ThermalPrecision precision,
                              std::size_t batch_count, std::size_t filament_count,
                              const double* currents, const double* masses,
                              const double* reference_resistances, const int* materials,
                              double dt, double* temperatures, double* resistivities,
                              double* resistances, double* joule_energy) {
    if (!batch_count || !filament_count || !currents || !masses || !reference_resistances ||
        !materials || !temperatures || !resistivities || !resistances || !joule_energy ||
        !std::isfinite(dt) || dt <= 0.0 ||
        filament_count > std::numeric_limits<std::size_t>::max() / batch_count)
        throw std::invalid_argument("invalid CPU thermal batch arguments");
    for (std::size_t i = 0; i < batch_count * filament_count; ++i) {
        if (materials[i] != 0 && materials[i] != 1) throw std::invalid_argument("invalid thermal material");
        if (!std::isfinite(currents[i]) || !std::isfinite(masses[i]) || masses[i] <= 0.0 ||
            !std::isfinite(reference_resistances[i]) || reference_resistances[i] <= 0.0 ||
            !std::isfinite(temperatures[i]) || !std::isfinite(resistances[i]) || resistances[i] < 0.0)
            throw std::invalid_argument("invalid CPU thermal state");
        const auto material = materials[i] == 1 ? ThermalMaterial::Copper : ThermalMaterial::Aluminum;
        const double r = resistances[i] == 0.0 ? reference_resistances[i] : resistances[i];
        const double q = currents[i] * currents[i] * r * dt;
        const double cp = interpolate_material_cp(tables, material, temperatures[i], precision);
        temperatures[i] += q / (masses[i] * cp);
        resistivities[i] = interpolate_material_resistivity(tables, material, temperatures[i], precision);
        resistances[i] = reference_resistances[i] *
            (1.0 + (material == ThermalMaterial::Copper ? 4.1e-3 : 4.2e-3) *
                       (temperatures[i] - T_REFERENCE));
        joule_energy[i] = q;
        if (!std::isfinite(temperatures[i]) || !std::isfinite(resistivities[i]) ||
            !std::isfinite(resistances[i]) || !std::isfinite(joule_energy[i]))
            throw std::runtime_error("CPU thermal update produced a non-finite result");
    }
}

} // namespace coilgun::simulation::cuda
