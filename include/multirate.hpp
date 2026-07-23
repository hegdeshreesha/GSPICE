#ifndef GSPICE_MULTIRATE_HPP
#define GSPICE_MULTIRATE_HPP

// ---------------------------------------------------------------------------
// Dynamic Multi-Rate Timestepping Controller (Feature B).
//
// Background:
//   In large mixed-signal or multi-stage analog/digital circuits, different
//   nodes operate at vastly different frequencies (e.g. 2 GHz clock core vs
//   1 kHz bias generator). Uniform global timestepping forces the simulator
//   to evaluate slow nodes thousands of times per nanosecond for no physical
//   benefit.
//
// Multi-Rate Concept:
//   1. Classify nodes into Latency Clusters based on local rate of change
//      |dx_i / dt| and local truncation error (LTE).
//   2. Nodes with small derivative |dx_i / dt| < tol are assigned to a SLOW
//      rate group (timestep multiplier M_i > 1).
//   3. Fast nodes step at Δt_fine; slow nodes step at Δt_coarse = M_i * Δt_fine.
//   4. During intermediate fast steps, slow nodes are held constant or linearly
//      interpolated from their previous and projected state.
//   5. If a slow node experiences an event (excitation |ΔV| > vntol), it is
//      immediately promoted back to the FAST rate group.
// ---------------------------------------------------------------------------

#include "matrix.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace gspice {

enum class NodeActivityClass {
    Fast,   // Evaluated at every dt_fine
    Medium, // Evaluated every 2*dt_fine
    Slow    // Evaluated every 4*dt_fine
};

struct MultiRateParams {
    bool   enabled               = false;
    double dv_dt_fast_threshold  = 1e6;  // V/s — above this is Fast
    double dv_dt_slow_threshold  = 1e3;  // V/s — below this is Slow
    double vntol                 = 1e-6; // Voltage absolute tolerance
    double reltol                = 1e-3; // Voltage relative tolerance
    int    max_ratio             = 4;    // Max ratio between coarse and fine step
};

class MultiRateController {
public:
    explicit MultiRateController(const MultiRateParams& params = MultiRateParams{})
        : params_(params) {}

    void initialize(int num_nodes) {
        num_nodes_ = num_nodes;
        activity_.assign(static_cast<std::size_t>(num_nodes), NodeActivityClass::Fast);
        last_x_.assign(static_cast<std::size_t>(num_nodes), 0.0);
        last_t_ = 0.0;
        step_counter_ = 0;
    }

    bool isEnabled() const { return params_.enabled; }
    void setEnabled(bool enable) { params_.enabled = enable; }

    /// Observe current solution x at time t to re-classify node activity classes.
    void observe(const VectorReal& x, double t) {
        if (!params_.enabled || num_nodes_ <= 0) return;
        const double dt = t - last_t_;
        if (dt <= 1e-18) return;

        const std::size_t size = std::min(static_cast<std::size_t>(x.getSize()),
                                          static_cast<std::size_t>(num_nodes_));
        for (std::size_t i = 0; i < size; ++i) {
            const double dv = std::abs(x[static_cast<int>(i)] - last_x_[i]);
            const double dv_dt = dv / dt;

            if (dv_dt > params_.dv_dt_fast_threshold || dv > params_.vntol * 10.0) {
                activity_[i] = NodeActivityClass::Fast;
            } else if (dv_dt < params_.dv_dt_slow_threshold && dv < params_.vntol) {
                activity_[i] = NodeActivityClass::Slow;
            } else {
                activity_[i] = NodeActivityClass::Medium;
            }
            last_x_[i] = x[static_cast<int>(i)];
        }
        last_t_ = t;
        ++step_counter_;
    }

    /// Check if node `node_idx` should be evaluated at current sub-step.
    bool shouldEvaluateNode(int node_idx) const {
        if (!params_.enabled || node_idx < 0 || node_idx >= num_nodes_) return true;
        const auto cls = activity_[static_cast<std::size_t>(node_idx)];
        if (cls == NodeActivityClass::Fast) return true;
        if (cls == NodeActivityClass::Medium) return (step_counter_ % 2 == 0);
        if (cls == NodeActivityClass::Slow)   return (step_counter_ % 4 == 0);
        return true;
    }

    /// Interpolate value for a slow node at intermediate time `t_interp`
    /// given prev_x (at t_prev) and next_x (at t_next).
    static double interpolateNode(double val_prev, double val_next,
                                   double t_prev, double t_next, double t_interp) {
        const double span = t_next - t_prev;
        if (span <= 1e-18) return val_next;
        const double alpha = std::clamp((t_interp - t_prev) / span, 0.0, 1.0);
        return val_prev + alpha * (val_next - val_prev);
    }

    int activeFastNodes() const {
        int count = 0;
        for (const auto& a : activity_) {
            if (a == NodeActivityClass::Fast) ++count;
        }
        return count;
    }

    int activeSlowNodes() const {
        int count = 0;
        for (const auto& a : activity_) {
            if (a == NodeActivityClass::Slow) ++count;
        }
        return count;
    }

private:
    MultiRateParams params_;
    int num_nodes_ = 0;
    double last_t_ = 0.0;
    long long step_counter_ = 0;
    std::vector<NodeActivityClass> activity_;
    std::vector<double> last_x_;
};

} // namespace gspice

#endif // GSPICE_MULTIRATE_HPP
