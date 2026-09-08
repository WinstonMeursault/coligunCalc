/**
 * @file coilgun.hpp
 * @brief Convenience umbrella header — includes the entire coilgun physics
 *        and components library.
 * @author Winston Meursault
 *
 * @par Modules
 *   - core/types     — common type forward declarations
 *   - physics/constants — μ0, Cu/Al material properties, cp(T), rho(T)
 *   - physics/elliptic  — K(m), E(m) complete elliptic integrals
 *   - physics/struve    — Struve H0(x), H1(x) (SciPy three-regime strategy)
 *   - physics/quadrature — Gauss-Legendre & Gauss-Laguerre quadrature
 *   - physics/cache     — LRU cache template (4096 entries default)
 *   - physics/lookup_tables — T(q,p) shape-factor dense table (2.9M entries)
 *   - physics/self_inductance — self-inductance (exact + fast table lookup)
 *   - physics/mutual_inductance — filament-level & coil-level mutual inductance
 *   - components/driving_coil — multi-turn helical driving coil
 *   - components/armature — solid cylindrical armature (m x n filaments)
 *   - optimization — variable schemas, evaluators, optimizers, results, and selectors
 */

#pragma once

#include "coilgun/core/types.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/elliptic.hpp"
#include "coilgun/physics/struve.hpp"
#include "coilgun/physics/quadrature.hpp"
#include "coilgun/physics/cache.hpp"
#include "coilgun/physics/lookup_tables.hpp"
#include "coilgun/physics/self_inductance.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/excitation_snapshot.hpp"
#include "coilgun/simulation/integration_state.hpp"
#include "coilgun/simulation/cpu_phase_timing.hpp"
#include "coilgun/simulation/derivative_workspace.hpp"
#include "coilgun/simulation/time_stepper.hpp"
#include "coilgun/simulation/sim_result.hpp"
#include "coilgun/simulation/termination.hpp"
#include "coilgun/simulation/single_stage_sim.hpp"
#include "coilgun/simulation/trigger_config.hpp"
#include "coilgun/simulation/multi_stage_result.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"

#include "coilgun/optimization/cache.hpp"
#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/comparator.hpp"
#include "coilgun/optimization/config.hpp"
#include "coilgun/optimization/constraint.hpp"
#include "coilgun/optimization/evaluator.hpp"
#include "coilgun/optimization/genetic_operators.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"
#include "coilgun/optimization/nsga2.hpp"
#include "coilgun/optimization/objective.hpp"
#include "coilgun/optimization/population.hpp"
#include "coilgun/optimization/problem.hpp"
#include "coilgun/optimization/result.hpp"
#include "coilgun/optimization/selectors.hpp"
#include "coilgun/optimization/statistics.hpp"
#include "coilgun/optimization/termination.hpp"
#include "coilgun/optimization/types.hpp"
#include "coilgun/optimization/variables.hpp"
