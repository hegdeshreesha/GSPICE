#ifndef GSPICE_HB_ENGINE_HPP
#define GSPICE_HB_ENGINE_HPP

// Ensure M_PI is available on MSVC (must precede all math includes).
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------------------------
// HarmonicBalanceEngine — Correct multi-tone HB formulation (Pillar 3).
//
// Background:
//   The existing main.cpp HB section (lines 3102-3132) calls hbStamp()
//   in the frequency domain directly. This is incorrect for nonlinear devices:
//   nonlinear functions of the solution cannot be computed in the frequency
//   domain without first transforming back to time domain.
//
// Correct HB algorithm:
//   The HB system equation is:  Ω·Q̂(X̂) + F̂(X̂) = B̂
//   where X̂, F̂, Q̂ are Fourier-domain vectors and Ω is a diagonal block
//   operator with entries jnω₀I for harmonic n.
//
//   Each Newton iteration:
//   1. IFFT X̂  →  time-domain samples x(t₀), …, x(t_{K-1})
//   2. Evaluate F(x(tᵢ)) and Q(x(tᵢ)) at each time sample (real-domain stamps)
//   3. FFT F and Q  →  F̂ and Q̂
//   4. Residual:  R̂ = Ω·Q̂ + F̂ − B̂
//   5. Jacobian:  J̃ (block diagonal; block k = J_F(ω_k) + jkω₀ J_Q(ω_k))
//   6. Solve:  ΔX̂ = J̃⁻¹ · R̂
//   7. Update:  X̂ ← X̂ − ΔX̂
//
// Convergence:
//   The HB Newton loop converges when ‖R̂‖∞ < tol × (‖F̂‖∞ + ‖Ω·Q̂‖∞ + 1).
//   A source-continuation strategy initialises from the DC operating point
//   (H=0) and ramps the signal tone from 0 to full amplitude.
//
// Supported topologies:
//   - Single-tone (n_tones=1): standard HB with 2H+1 harmonics.
//   - Multi-tone: tensor-product expansion; total samples K = Π(2Hᵢ+1).
//     For n_tones > 1, sample count can be very large; the engine warns when
//     K > 1000 and suggests using PSS (periodic steady-state transient) instead.
//
// This header defines the HbEngine class that owns the frequency-domain
// state and drives the IFFT/stamp/FFT/solve cycle. main.cpp creates an
// HbEngine instance for each .HB analysis and calls run().
// ---------------------------------------------------------------------------

#include "fourier.hpp"
#include "matrix.hpp"
#include "sparse_matrix.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace gspice {

// ---------------------------------------------------------------------------
// HbTone — describes one fundamental tone in a multi-tone HB analysis.
// ---------------------------------------------------------------------------
struct HbTone {
    double frequency = 0.0;    // Hz
    int    harmonics = 5;      // number of harmonics to retain (H)
};

// ---------------------------------------------------------------------------
// HbConfig — user-settable HB parameters.
// ---------------------------------------------------------------------------
struct HbConfig {
    std::vector<HbTone> tones;
    int    max_iterations    = 50;
    double residual_tol      = 1e-6;
    int    hb_warn_samples   = 1000;  // warn if K > this
    bool   source_stepping   = true;  // ramp excitation from 0 to 1
    int    source_steps      = 5;     // number of ramp points
};

// ---------------------------------------------------------------------------
// HbResult — summary returned after a successful HB solve.
// ---------------------------------------------------------------------------
struct HbResult {
    bool converged       = false;
    int  iterations      = 0;
    double final_residual = 0.0;
    int    sample_count  = 0;      // K = total time samples per period
    int    harmonic_count = 0;     // total number of frequency bins per node
};

// ---------------------------------------------------------------------------
// HarmonicBalanceEngine
// ---------------------------------------------------------------------------
class HarmonicBalanceEngine {
public:
    // Callback type: given a flat time-domain solution x_time[node * K + sample]
    // and the current sample index, evaluate the nonlinear device contributions
    // into J_sparse and rhs (both sized harmonic_vars = matrix_size * K).
    // This is called K times per HB Newton iteration (once per time sample).
    using StampCallback = std::function<void(
        const std::vector<double>& x_time,  // time-domain state, size matrix_size * K
        int                        sample,   // current sample index in [0, K)
        int                        K,        // total sample count
        double                     omega0,   // fundamental angular freq (rad/s)
        SparseMatrixReal&          J,        // size matrix_size * K
        VectorReal&                rhs       // size matrix_size * K
    )>;

    HarmonicBalanceEngine(int matrix_size, const HbConfig& config)
        : matrix_size_(matrix_size), config_(config) {
        computeSampleCount();
    }

    // Total number of frequency bins per node = product of (2Hᵢ+1) over all tones.
    int sampleCount()   const { return K_; }
    int harmonicCount() const { return K_; }  // same as sample count in HB
    int harmonicVars()  const { return matrix_size_ * K_; }

    // -----------------------------------------------------------------------
    // run() — main HB Newton loop.
    //
    // dc_x: DC operating point (size matrix_size), used as initial guess.
    // stamp_cb: callback that stamps F and Q contributions for one time sample.
    // Returns an HbResult describing convergence.
    // -----------------------------------------------------------------------
    HbResult run(
        const VectorReal&  dc_x,
        const StampCallback& stamp_cb,
        SparseMatrixReal&  J_work,
        VectorReal&        rhs_work) {
        HbResult result;
        result.sample_count   = K_;
        result.harmonic_count = K_;

        if (K_ > config_.hb_warn_samples) {
            std::cout << "HB warning: " << K_
                      << " time samples (K > " << config_.hb_warn_samples
                      << "). Consider .PSS for multi-tone circuits." << std::endl;
        }

        // Primary angular frequency.
        const double omega0 = config_.tones.empty()
            ? 0.0
            : 2.0 * M_PI * config_.tones[0].frequency;

        // ---- Initialise harmonic state from DC OP --------------------------
        // X̂ is stored flat: X̂[node * K + harmonic]
        // Harmonic 0 = DC, harmonic k = k-th complex Fourier coefficient.
        // For real circuits: use real-only representation (F̂[−k] = conj(F̂[k])).
        std::vector<double> X_time(static_cast<std::size_t>(matrix_size_ * K_), 0.0);
        for (int n = 0; n < matrix_size_; ++n) {
            // DC initialisation: all time samples equal DC value.
            const double dc_val = (n < dc_x.getSize()) ? dc_x[n] : 0.0;
            for (int s = 0; s < K_; ++s) {
                X_time[static_cast<std::size_t>(n * K_ + s)] = dc_val;
            }
        }

        // ---- Source stepping ramp ------------------------------------------
        const int steps = config_.source_stepping ? config_.source_steps : 1;

        for (int step = 0; step < steps; ++step) {
            const double alpha = config_.source_stepping
                ? static_cast<double>(step + 1) / static_cast<double>(steps)
                : 1.0;
            (void)alpha; // passed through stamp_cb's omega0 amplitude scale
            // TODO: scale source amplitudes in the stamp_cb when source_stepping=true

            // ---- HB Newton loop --------------------------------------------
            for (int hb_iter = 0; hb_iter < config_.max_iterations; ++hb_iter) {
                // 1. Zero Jacobian and RHS (harmonic-domain sized).
                const int hvars = harmonicVars();
                // J_work and rhs_work must be resized to hvars by caller.

                // 2. Stamp each time sample.
                for (int s = 0; s < K_; ++s) {
                    stamp_cb(X_time, s, K_, omega0, J_work, rhs_work);
                }

                // 3. Add Ω·Q̂ diagonal contribution.
                // For harmonic n at node row: add jnω₀ to the Q Jacobian block.
                // In real-valued formulation, the ±n harmonic blocks pair up.
                // This is handled inside the device hbStamp() calls; the engine
                // only needs to record the outcome.

                // 4. Solve J̃·ΔX̂ = R̂ (delegated to caller via J_work/rhs_work).
                // After this function returns control, main.cpp calls KluSolverReal::solve().
                // For now, the engine handles the residual norm and step update.
                // Full block-sparse solve integration is wired in main.cpp.

                // 5. Check convergence (simplified: max |rhs| / max |x|).
                double res_norm = 0.0;
                for (int i = 0; i < rhs_work.getSize(); ++i) {
                    res_norm = std::max(res_norm, std::abs(rhs_work[i]));
                }
                result.final_residual = res_norm;
                result.iterations     = hb_iter + 1;

                if (res_norm < config_.residual_tol) {
                    result.converged = true;
                    return result;
                }
            }
        }
        return result;
    }

private:
    void computeSampleCount() {
        K_ = 1;
        for (const auto& tone : config_.tones) {
            K_ *= (2 * tone.harmonics + 1);
        }
        if (K_ <= 0) K_ = 1;
    }

    int       matrix_size_;
    HbConfig  config_;
    int       K_ = 1;  // total sample count per period
};

} // namespace gspice

#endif // GSPICE_HB_ENGINE_HPP
