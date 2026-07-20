/**
 * @file gpu_state_layout.hpp
 * @brief Fixed row-major host index mapping for batched GPU state arrays.
 */

#pragma once

#include <cstddef>
#include <stdexcept>
#include <limits>

namespace coilgun::simulation::cuda {

/**
 * @brief Immutable dimensions and index mapping for GPU state buffers.
 *
 * Each named state buffer is a separate flat allocation. Indices are
 * row-major: the rightmost logical coordinate varies fastest. Active and
 * trigger masks only select work; they never compact or renumber physical
 * state positions.
 */
class GpuStateLayout {
public:
    const std::size_t B;
    const std::size_t S;
    const std::size_t F;
    const std::size_t D;

    GpuStateLayout(std::size_t batch_size,
                   std::size_t stage_count,
                   std::size_t filament_count)
        : B(batch_size), S(stage_count), F(filament_count), D(checked_sum(stage_count, filament_count)) {
        if (B == 0 || S == 0 || F == 0) {
            throw std::invalid_argument("GPU state dimensions must be positive");
        }
        checked_product(B, D, "currents");
        checked_product(B, S, F, "state");
        checked_product(B, F, "temperatures");
        checked_product(B, S, "trigger mask");
        checked_product(B, D, D, "system matrix");
    }

    std::size_t batch_size() const { return B; }
    std::size_t stage_count() const { return S; }
    std::size_t filament_count() const { return F; }
    std::size_t current_dimension() const { return D; }

    std::size_t currents(std::size_t batch, std::size_t dimension) const {
        return currents_offset(batch, dimension);
    }

    std::size_t currents_offset(std::size_t batch, std::size_t dimension) const {
        check(batch, B, "batch");
        check(dimension, D, "current dimension");
        return batch * D + dimension;
    }

    std::size_t m1(std::size_t batch, std::size_t stage, std::size_t filament) const {
        return m1_offset(batch, stage, filament);
    }

    std::size_t m1_offset(std::size_t batch,
                          std::size_t stage,
                          std::size_t filament) const {
        return state_index(batch, stage, filament);
    }

    std::size_t dm1(std::size_t batch, std::size_t stage, std::size_t filament) const {
        return dm1_offset(batch, stage, filament);
    }

    std::size_t dm1_offset(std::size_t batch,
                           std::size_t stage,
                           std::size_t filament) const {
        return state_index(batch, stage, filament);
    }

    std::size_t temperatures(std::size_t batch, std::size_t filament) const {
        return temperatures_offset(batch, filament);
    }

    std::size_t temperatures_offset(std::size_t batch, std::size_t filament) const {
        check(batch, B, "batch");
        check(filament, F, "filament");
        return batch * F + filament;
    }

    std::size_t system_matrix(std::size_t batch, std::size_t row, std::size_t column) const {
        check(batch, B, "batch"); check(row, D, "matrix row"); check(column, D, "matrix column");
        return (batch * D + row) * D + column;
    }

    std::size_t rhs(std::size_t batch, std::size_t row) const {
        check(batch, B, "batch"); check(row, D, "rhs row");
        return batch * D + row;
    }

    std::size_t active_mask(std::size_t batch) const {
        return active_mask_offset(batch);
    }

    std::size_t active_mask_offset(std::size_t batch) const {
        check(batch, B, "batch");
        return batch;
    }

    std::size_t trigger_mask(std::size_t batch, std::size_t stage) const {
        return trigger_mask_offset(batch, stage);
    }

    std::size_t trigger_mask_offset(std::size_t batch, std::size_t stage) const {
        check(batch, B, "batch");
        check(stage, S, "stage");
        return batch * S + stage;
    }

    std::size_t currents_size() const { return B * D; }
    std::size_t state_size() const { return B * S * F; }
    std::size_t temperatures_size() const { return B * F; }
    std::size_t active_mask_size() const { return B; }
    std::size_t trigger_mask_size() const { return B * S; }
    std::size_t system_matrix_size() const { return B * D * D; }
    std::size_t rhs_size() const { return B * D; }

private:
    static std::size_t checked_sum(std::size_t left, std::size_t right) {
        if (right > std::numeric_limits<std::size_t>::max() - left) {
            throw std::overflow_error("GPU layout dimension overflow");
        }
        return left + right;
    }

    static std::size_t checked_product(std::size_t left, std::size_t right, const char* name) {
        if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
            throw std::overflow_error(name);
        }
        return left * right;
    }

    static std::size_t checked_product(std::size_t first, std::size_t second,
                                       std::size_t third, const char* name) {
        return checked_product(checked_product(first, second, name), third, name);
    }

    static void check(std::size_t index, std::size_t extent, const char* name) {
        if (index >= extent) {
            throw std::out_of_range(name);
        }
    }

    std::size_t state_index(std::size_t batch,
                            std::size_t stage,
                            std::size_t filament) const {
        check(batch, B, "batch");
        check(stage, S, "stage");
        check(filament, F, "filament");
        return (batch * S + stage) * F + filament;
    }
};

} // namespace coilgun::simulation::cuda
