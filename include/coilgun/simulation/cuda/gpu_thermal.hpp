#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>
#include <cuda_runtime_api.h>

namespace coilgun::simulation::cuda {

enum class ThermalPrecision { Standard, Full, Aggressive };
enum class ThermalMaterial { Aluminum = 0, Copper = 1 };

struct MaterialTables {
    std::size_t version = 1;
    double minimum_temperature = 293.0;
    double maximum_temperature = 2000.0;
    std::vector<double> temperatures;
    std::vector<double> cp_aluminum;
    std::vector<double> cp_copper;
    std::vector<double> rho_aluminum;
    std::vector<double> rho_copper;
};

struct ThermalWorkspaceKey {
    int device_id = -1;
    std::size_t table_version = 0;
    double min_temperature = 0.0;
    double max_temperature = 0.0;
    std::size_t table_length = 0;
    std::size_t value_count = 0;
    bool operator==(const ThermalWorkspaceKey&) const = default;
};

class ThermalWorkspace {
public:
    ThermalWorkspace() = default;
    ThermalWorkspace(const MaterialTables& tables, std::size_t count);
    ~ThermalWorkspace();
    ThermalWorkspace(const ThermalWorkspace&) = delete;
    ThermalWorkspace& operator=(const ThermalWorkspace&) = delete;

    void initialize(const MaterialTables& tables, std::size_t count);
    void update(const MaterialTables& tables, ThermalPrecision precision,
                std::size_t batch_count, std::size_t filament_count,
                const double* currents, const double* masses,
                const double* reference_resistances, const int* materials,
                double dt, double* temperatures, double* resistivities,
                double* resistances, double* joule_energy,
                cudaStream_t stream = nullptr);
    void initialize_device_state(const MaterialTables& tables,
                                 std::size_t batch_count,
                                 std::size_t filament_count,
                                 const double* masses,
                                 const double* reference_resistances,
                                 const int* materials,
                                 const double* temperatures,
                                 const double* resistances,
                                 cudaStream_t stream = nullptr);
    cudaError_t launch_device(ThermalPrecision precision,
                              std::size_t batch_count,
                              std::size_t stage_count,
                              std::size_t filament_count,
                              const double* currents,
                              const std::uint8_t* active_mask,
                              double dt,
                              cudaStream_t stream = nullptr) noexcept;
    cudaError_t download_device_state(double* temperatures,
                                      double* resistivities,
                                      double* resistances,
                                      double* joule_energy,
                                      cudaStream_t stream = nullptr) noexcept;
    double* device_resistances() const noexcept;
    std::size_t allocation_count() const noexcept { return allocation_count_; }
    const ThermalWorkspaceKey& key() const noexcept { return key_; }
    std::vector<std::uintptr_t> device_addresses() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    std::size_t allocation_count_ = 0;
    ThermalWorkspaceKey key_;
};

MaterialTables generate_material_tables(std::size_t sample_count = 1024,
                                        double minimum_temperature = 293.0,
                                        double maximum_temperature = 2000.0);

double interpolate_material_cp(const MaterialTables&, ThermalMaterial, double temperature,
                               ThermalPrecision = ThermalPrecision::Standard);
double interpolate_material_resistivity(const MaterialTables&, ThermalMaterial, double temperature,
                                        ThermalPrecision = ThermalPrecision::Standard);

void update_thermal_batch(const MaterialTables&, ThermalPrecision, std::size_t batch_count,
                          std::size_t filament_count, const double* currents,
                          const double* masses, const double* reference_resistances,
                           const int* materials, double dt, double* temperatures,
                           double* resistivities, double* resistances, double* joule_energy,
                           cudaStream_t stream = nullptr);

void update_thermal_batch_cpu(const MaterialTables&, ThermalPrecision, std::size_t batch_count,
                              std::size_t filament_count, const double* currents,
                              const double* masses, const double* reference_resistances,
                              const int* materials, double dt, double* temperatures,
                              double* resistivities, double* resistances, double* joule_energy);

} // namespace coilgun::simulation::cuda
