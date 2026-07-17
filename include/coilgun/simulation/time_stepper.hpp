/**
 * @file time_stepper.hpp
 * @brief Time-integration policy classes for ODE systems.
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::simulation {

/**
 * @brief Forward Euler time stepper (1st-order explicit).
 *
 * @f$ x_{n+1} = x_n + \Delta t \cdot f(x_n) @f$
 *
 * Simple and fast; suitable when the time step is small enough to
 * guarantee stability.
 */
struct EulerStepper {
    /**
     * @brief Advance state by one time step.
     * @tparam State State vector type (must support arithmetic scalar-vector ops).
     * @tparam DerivativeFunc Callable mapping State -> State.
     * @param dt Time step, s.
     * @param state Current state.
     * @param f Derivative function.
     * @return State at the next time level.
     */
    template<typename State, typename DerivativeFunc>
    State advance(double dt, const State& state, DerivativeFunc&& f) const {
        return state + dt * f(state);
    }
};

/**
 * @brief Classical 4th-order Runge-Kutta time stepper.
 *
 * Four derivative evaluations per step; 4th-order accurate for smooth
 * right-hand sides.
 */
struct RK4Stepper {
    /**
     * @brief Advance state by one time step.
     * @tparam State State vector type (must support arithmetic scalar-vector ops).
     * @tparam DerivativeFunc Callable mapping State -> State.
     * @param dt Time step, s.
     * @param state Current state.
     * @param f Derivative function.
     * @return State at the next time level.
     */
    template<typename State, typename DerivativeFunc>
    State advance(double dt, const State& state, DerivativeFunc&& f) const {
        State k1 = f(state);
        State k2 = f(state + 0.5 * dt * k1);
        State k3 = f(state + 0.5 * dt * k2);
        State k4 = f(state + dt * k3);
        return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
};

} // namespace coilgun::simulation
