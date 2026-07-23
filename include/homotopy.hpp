#ifndef GSPICE_HOMOTOPY_HPP
#define GSPICE_HOMOTOPY_HPP

// ---------------------------------------------------------------------------
// Pseudo-Transient Continuation (PTC) — Complete Pseudo-Time ODE Homotopy.
//
// Background:
//   Source stepping and gmin stepping are GSPICE's recovery paths when direct
//   Newton fails for difficult DC operating points. Both can fail on highly
//   nonlinear circuits (feedback amplifiers, oscillators, SMPS) because the
//   homotopy parameter must be tuned to avoid divergence or stall.
//
// Full Pseudo-Transient Formulation:
//   PTC converts the steady-state equation F(x) = 0 into a pseudo-time ODE:
//        C_ptc * dx / dτ + F(x) = 0
//   where τ is pseudo-time. Discretizing via Backward Euler gives:
//        (C_ptc / Δτ) * (x_k - x_{k-1}) + F(x_k) = 0
//
//   Jacobian:  J_aug = J_F + (C_ptc / Δτ) * I
//   Residual:  R_aug = R_F - (C_ptc / Δτ) * (x_k - x_{k-1})
//
//   As pseudo-time τ → ∞ (or as C_ptc / Δτ → 0), x_k converges to the true
//   DC operating point solution F(x) = 0.
// ---------------------------------------------------------------------------

#include "matrix.hpp"
#include "sparse_matrix.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace gspice {

/// Parameters that control pseudo-transient continuation.
struct PtcParams {
    double c_start   = 1e-9;   // Starting virtual capacitance (Farads)
    double c_min     = 1e-20;  // Threshold to consider PTC bias fully removed
    double tau       = 0.1;    // Reduction factor per successful pseudo step
    double dt_init   = 1e-12;  // Initial pseudo timestep (seconds)
    double dt_min    = 1e-18;  // Minimum pseudo timestep
    double dt_growth = 2.0;    // Pseudo timestep growth factor
    double g_max     = 1e6;    // Siemens cap on diagonal conductance
};

/// Stateful controller for full PTC solves.
class PtcController {
public:
    explicit PtcController(const PtcParams& p = PtcParams{})
        : params_(p), c_ptc_(p.c_start), dt_(p.dt_init) {}

    bool active() const { return c_ptc_ > params_.c_min; }

    /// Conductance C_ptc / Δτ added to diagonal of J.
    double diagonalConductance() const {
        return std::min(c_ptc_ / std::max(dt_, 1e-30), params_.g_max);
    }

    /// Update previous pseudo-time state vector x_{k-1}.
    void updatePreviousState(const VectorReal& x) {
        prev_x_ = x;
        has_prev_ = true;
    }

    /// Stamp full PTC diagonal conductance and previous-state residual:
    ///   J_ii += g_ptc
    ///   R_i  += g_ptc * (x_k[i] - x_{k-1}[i])
    void stampPtc(SparseMatrixReal& J, VectorReal& rhs, const VectorReal& x, int num_nodes) const {
        if (!active()) return;
        const double g_ptc = diagonalConductance();
        if (g_ptc <= 0.0) return;

        const int n = std::min(num_nodes, x.getSize());
        for (int i = 0; i < n; ++i) {
            J.add(i, i, g_ptc);
            if (has_prev_ && i < prev_x_.getSize()) {
                const double delta_x = x[i] - prev_x_[i];
                rhs.add(i, -g_ptc * delta_x);
            }
        }
    }

    void notifyConverged() {
        c_ptc_ *= params_.tau;
        dt_ = std::min(dt_ * params_.dt_growth,
                       std::max(dt_, 1.0 / std::max(diagonalConductance(), 1e-30)));
    }

    void notifyDiverged() {
        dt_ = std::max(params_.dt_min, dt_ * 0.5);
    }

    double currentCapacitance() const { return c_ptc_; }
    double currentDt()          const { return dt_; }

private:
    PtcParams  params_;
    double     c_ptc_;
    double     dt_;
    VectorReal prev_x_;
    bool       has_prev_ = false;
};

} // namespace gspice

#endif // GSPICE_HOMOTOPY_HPP
