#ifndef GSPICE_DEVICE_HPP
#define GSPICE_DEVICE_HPP

#include "matrix.hpp"
#include "sparse_matrix.hpp"
#include "dae.hpp"
#include <string>
#include <vector>
#include <memory>
#include <cstddef>
#include <stdexcept>

namespace gspice {

enum class TransientIntegrationMethod {
    BackwardEuler,
    Trapezoidal,
    Gear2,
    Bdf,
    AdamsMoulton
};

struct TransientContext {
    double timeStep = 0.0;
    double currentTime = 0.0;
    TransientIntegrationMethod method = TransientIntegrationMethod::BackwardEuler;
    double a0 = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    bool hasSecondHistory = false;
    int integrationOrder = 1;
    std::vector<double> qWeights;
    std::vector<double> derivativeWeights;
    const std::vector<VectorReal>* xHistory = nullptr;
    const std::vector<double>* timeHistory = nullptr;
};

struct NoiseSource {
    std::string name;
    int nodePos = -1;
    int nodeNeg = -1;
    double currentPsd = 0.0;
};

class Device {
public:
    Device(const std::string& name) : name_(name) {}
    virtual ~Device() = default;

    std::string getName() const { return name_; }

    /**
     * Device-neutral DAE evaluation. New and migrated devices expose their
     * instantaneous residual F, stored quantity Q, and the corresponding
     * Jacobians here. Legacy analysis-specific stamps remain available during
     * the staged migration.
     */
    virtual bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) {
        (void)x;
        (void)request;
        evaluation.clear();
        return false;
    }

    virtual bool daeAuditSafe() const { return false; }

    // Compact models with stiff or noisy charge dynamics can request the
    // numerically damped AUTO branch. Smooth native devices use trapezoidal
    // integration after backward-Euler startup.
    virtual bool prefersDampedAutoTransient() const { return false; }

    /**
     * Stamping for DC and Transient analysis (Real numbers).
     */
    virtual void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) = 0;

    /**
     * Transient stamping. Devices that are dynamic can use the integration
     * coefficients in TransientContext; static devices reuse their DC stamp.
     */
    virtual void tranStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        const TransientContext& ctx) {
        static const std::vector<VectorReal> empty_history;
        const auto& history = ctx.xHistory ? *ctx.xHistory : empty_history;
        dcStamp(J, b, x, ctx.timeStep, ctx.currentTime, history);
    }

    /**
     * Called after a transient timestep has converged and is accepted.
     * Devices with internal dynamic state can commit next-state buffers here.
     */
    virtual void acceptTransientStep(const VectorReal& x, double currentTime) {
        (void)x;
        (void)currentTime;
    }

    virtual void acceptTransientStep(
        const VectorReal& x,
        double currentTime,
        const TransientContext& ctx) {
        (void)ctx;
        acceptTransientStep(x, currentTime);
    }

    // Fixed-size serialization for the simulator-owned transactional state
    // arena. Adaptive candidates are copied and rolled back without per-step
    // heap allocation or polymorphic snapshot objects.
    virtual std::size_t transientStateBytes() const { return 0; }
    virtual void saveTransientStateBytes(std::byte* destination, std::size_t size) const {
        (void)destination;
        if (size != 0) throw std::invalid_argument("unexpected transient state storage");
    }
    virtual void restoreTransientStateBytes(const std::byte* source, std::size_t size) {
        (void)source;
        if (size != 0) throw std::invalid_argument("unexpected transient state storage");
    }

    /**
     * Time points where independent source waveforms or device behavior change
     * abruptly. Transient analysis uses these breakpoints to land on edges
     * without forcing tiny timesteps everywhere.
     */
    virtual void collectBreakpoints(double t_stop, std::vector<double>& points) const {
        (void)t_stop;
        (void)points;
    }

    /**
     * Optional device-requested maximum next transient step. Compact models can
     * use this for Verilog-A $bound_step support after evaluating their state.
     */
    virtual double transientBoundStep() const {
        return 0.0;
    }

    // Normalized difference in stored charge between two candidate endpoint
    // solutions. This supplements voltage/current LTE checks for dynamic
    // devices. A value above one requests timestep rejection.
    virtual double transientChargeError(
        const VectorReal& coarse,
        const VectorReal& fine,
        double reltol,
        double chgtol) {
        (void)coarse;
        (void)fine;
        (void)reltol;
        (void)chgtol;
        return 0.0;
    }

    /**
     * Limit a proposed transient Newton update using device physics.  The
     * default is intentionally a no-op; exponential junction devices override
     * this to prevent an otherwise harmless first Newton step from driving an
     * exponential far outside its useful linearization region.
     */
    virtual void limitTransientNewton(
        const VectorReal& previous,
        VectorReal& candidate) const {
        (void)previous;
        (void)candidate;
    }

    /**
     * Stamping for AC analysis (Complex numbers).
     * @param omega Angular frequency (omega = 2*pi*f)
     */
    virtual void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) = 0;

    /**
     * Returns the noise contribution of this device.
     */
    virtual double getNoisePSD(double omega, const VectorReal& x_dc) { return 0.0; }

    /**
     * Adds equivalent small-signal noise current sources for output-referred
     * noise analysis. Each source is solved through the AC matrix.
     */
    virtual void collectNoiseSources(
        double omega,
        const VectorReal& x_dc,
        std::vector<NoiseSource>& sources) const {
        (void)omega;
        (void)x_dc;
        (void)sources;
    }

    /**
     * Stamping for Harmonic Balance analysis (Frequency domain non-linear).
     * @param f_fund Fundamental frequency
     * @param n_harms Number of harmonics
     * @param x_hb Current harmonic voltages state
     */
    virtual void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) {}

    /**
     * Stamping for Periodic Small-Signal analysis (PAC/HBAC).
     * @param J The large Conversion Matrix
     * @param b The RHS sideband vector
     * @param f_in Input small-signal frequency
     * @param f_fund Periodic fundamental frequency
     * @param x_periodic The converged PSS or HB state
     */
    virtual void pacStamp(SparseMatrixReal& J, VectorReal& b, double f_in, double f_fund, int n_harms, const VectorReal& x_periodic) {}

protected:
    std::string name_;
};

} // namespace gspice

#endif // GSPICE_DEVICE_HPP
