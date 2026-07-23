#ifndef GSPICE_FASTSPICE_ENGINE_HPP
#define GSPICE_FASTSPICE_ENGINE_HPP

// ---------------------------------------------------------------------------
// FastSpiceEngine — Event-Driven Selective Matrix Assembly (Feature E).
//
// Concept:
//   In large ICs (memory arrays, logic circuits, multi-stage amplifiers),
//   90% of sub-circuits are in a steady/dormant state at any given instant.
//   Standard SPICE re-evaluates every device's nonlinear equations and re-stamps
//   their Jacobians every Newton iteration, wasting 70-90% of total CPU time.
//
// FastSPICE Event-Driven selective assembly:
//   1. Partition circuit into Device Sub-Blocks.
//   2. Track local terminal voltage delta:
//        ΔV = max |v_k(t) - v_k(t_last)|
//   3. If ΔV < reltol * V + vntol, mark the sub-block LATENT.
//   4. LATENT sub-blocks skip expensive device model evaluation (OSDI eval / MOS
//      equations) and re-use their cached Jacobian/RHS stamp.
//   5. Active sub-blocks are fully evaluated and re-stamped.
//   6. Wake up latent sub-blocks instantly when terminal excitation occurs.
// ---------------------------------------------------------------------------

#include "device.hpp"
#include "matrix.hpp"
#include "sparse_matrix.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

namespace gspice {

struct FastSpiceParams {
    bool   enabled            = false;
    double latency_vntol      = 1e-5; // V — voltage tolerance for latent status
    double latency_reltol     = 1e-3; // relative tolerance for latent status
    int    max_latent_steps   = 20;   // max consecutive steps a device can stay latent before forced refresh
};

struct DeviceLatencyStatus {
    bool is_latent = false;
    int  latent_steps = 0;
    std::vector<double> last_terminal_voltages;
};

class FastSpiceEngine {
public:
    struct Stats {
        long long total_evaluations = 0;
        long long bypassed_evaluations = 0;
        long long active_evaluations = 0;
        double bypass_ratio() const {
            return total_evaluations > 0
                ? static_cast<double>(bypassed_evaluations) / static_cast<double>(total_evaluations)
                : 0.0;
        }
    };

    explicit FastSpiceEngine(const FastSpiceParams& params = FastSpiceParams{})
        : params_(params) {}

    void initialize(std::size_t num_devices) {
        status_.resize(num_devices);
        for (auto& s : status_) {
            s.is_latent = false;
            s.latent_steps = 0;
            s.last_terminal_voltages.clear();
        }
    }

    bool isEnabled() const { return params_.enabled; }
    void setEnabled(bool enable) { params_.enabled = enable; }

    /// Check if device `dev_idx` is latent (can bypass evaluation).
    bool isDeviceLatent(std::size_t dev_idx, const std::vector<double>& current_voltages) {
        if (!params_.enabled || dev_idx >= status_.size()) return false;
        auto& st = status_[dev_idx];
        ++stats_.total_evaluations;

        if (st.last_terminal_voltages.empty() || st.latent_steps >= params_.max_latent_steps) {
            st.is_latent = false;
            st.latent_steps = 0;
            st.last_terminal_voltages = current_voltages;
            ++stats_.active_evaluations;
            return false;
        }

        // Compute max terminal delta
        double max_delta = 0.0;
        const std::size_t n = std::min(current_voltages.size(), st.last_terminal_voltages.size());
        for (std::size_t i = 0; i < n; ++i) {
            const double v = current_voltages[i];
            const double v_last = st.last_terminal_voltages[i];
            const double tol = params_.latency_reltol * std::abs(v) + params_.latency_vntol;
            const double diff = std::abs(v - v_last);
            if (diff > tol) {
                max_delta = std::max(max_delta, diff);
                break;
            }
        }

        if (max_delta == 0.0) { // All within tolerance
            st.is_latent = true;
            ++st.latent_steps;
            ++stats_.bypassed_evaluations;
            return true;
        }

        // Active
        st.is_latent = false;
        st.latent_steps = 0;
        st.last_terminal_voltages = current_voltages;
        ++stats_.active_evaluations;
        return false;
    }

    Stats getStats() const { return stats_; }
    void resetStats() { stats_ = Stats{}; }

private:
    FastSpiceParams params_;
    std::vector<DeviceLatencyStatus> status_;
    Stats stats_;
};

} // namespace gspice

#endif // GSPICE_FASTSPICE_ENGINE_HPP
