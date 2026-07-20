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
    double minimum_temperature = 293.0;
    double maximum_temperature = 2000.0;
    std::vector<double> temperatures;
    std::vector<double> cp_aluminum;
    std::vector<double> cp_copper;
    std::vector<double> rho_aluminum;
    std::vector<double> rho_copper;
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
    std::size_t allocation_count() const noexcept { return allocation_count_; }
    std::vector<std::uintptr_t> device_addresses() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    std::size_t allocation_count_ = 0;
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
