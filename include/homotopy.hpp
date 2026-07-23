#ifndef GSPICE_HOMOTOPY_HPP
#define GSPICE_HOMOTOPY_HPP

// ---------------------------------------------------------------------------
// Pseudo-Transient Continuation (PTC) — Pillar 4 homotopy method.
//
// Background:
//   Source stepping and gmin stepping are GSPICE's existing recovery paths
//   when direct Newton fails for the DC operating point. Both can fail on
//   highly nonlinear circuits (e.g. feedback amplifiers, oscillators,
//   switched-mode power supplies) because the homotopy parameter must be
//   tuned to avoid both divergence and stall.
//
//   Pseudo-transient continuation adds a virtual capacitor C_ptc from every
//   node to ground. This converts the OP Newton system into a stable ODE
//   that drives the circuit toward steady state. As iterations proceed,
//   C_ptc is reduced toward zero so the DC solution is recovered.
//
//   This technique is used by Spectre, HSPICE, and Eldo as a tertiary
//   recovery strategy. It is particularly effective for:
//   - Oscillators (where source-stepping finds a metastable non-physical OP)
//   - Digital circuits with feedback (inverter chains, SRAM cells)
//   - RF mixers with strong LO signal
//
// Algorithm (simplified):
//   1. Start with C_ptc_max (e.g. 1e-9 F) added to each node row of J.
//   2. Each Newton iteration: J_aug = J + (C_ptc / dt) × I where dt is the
//      pseudo time step (adapts to convergence).
//   3. If Newton converges with C_ptc > C_ptc_min, reduce C_ptc by factor τ.
//   4. Repeat until C_ptc < C_ptc_min, then accept as DC OP.
//
// Integration into main.cpp:
//   PtcController is created in solve_dc_with_recovery() as the third
//   fallback after source stepping and gmin stepping.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace gspice {

/// Parameters that control the pseudo-transient continuation.
struct PtcParams {
    /// Starting capacitance added to every node (Farads).
    double c_start = 1e-9;
    /// Final capacitance; below this threshold the DA bias is removed.
    double c_min = 1e-20;
    /// Reduction factor per successful Newton convergence. 0 < tau < 1.
    double tau = 0.1;
    /// Initial pseudo time step (seconds).
    double dt_init = 1e-12;
    /// Minimum pseudo time step.
    double dt_min = 1e-18;
    /// Growth factor applied to dt after convergence.
    double dt_growth = 2.0;
    /// Maximum capacitance-to-time ratio added to diagonal — acts as a
    /// cap on the conductance injection to avoid ill-conditioning.
    double g_max = 1e6; // siemens
};

/// Stateful controller for one PTC solve.
class PtcController {
public:
    explicit PtcController(const PtcParams& p = PtcParams{})
        : params_(p), c_ptc_(p.c_start), dt_(p.dt_init) {}

    /// True while the PTC bias is still active.
    bool active() const { return c_ptc_ > params_.c_min; }

    /// Conductance to stamp on the diagonal of every voltage node's row.
    /// Equivalent to C_ptc / dt added to g_ii.
    double diagonalConductance() const {
        return std::min(c_ptc_ / std::max(dt_, 1e-30), params_.g_max);
    }

    /// Notify the controller that Newton converged at the current bias.
    /// Reduces C_ptc and increases dt for the next iteration.
    void notifyConverged() {
        c_ptc_ *= params_.tau;
        dt_ = std::min(dt_ * params_.dt_growth,
                       std::max(dt_, 1.0 / std::max(diagonalConductance(), 1e-30)));
    }

    /// Notify the controller that Newton diverged at the current bias.
    /// Increases C_ptc slightly to stabilise, reduces dt.
    void notifyDiverged() {
        dt_ = std::max(params_.dt_min, dt_ * 0.5);
        // Avoid reducing C below where it already is; let the solver try
        // the same bias with a smaller dt (larger diagonal conductance).
    }

    double currentCapacitance() const { return c_ptc_; }
    double currentDt()          const { return dt_; }

private:
    PtcParams params_;
    double    c_ptc_;
    double    dt_;
};

} // namespace gspice

#endif // GSPICE_HOMOTOPY_HPP
