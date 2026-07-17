#ifndef GSPICE_DEVICE_HPP
#define GSPICE_DEVICE_HPP

#include "matrix.hpp"
#include "sparse_matrix.hpp"
#include <string>
#include <vector>

namespace gspice {

enum class TransientIntegrationMethod {
    BackwardEuler,
    Trapezoidal,
    Gear2
};

struct TransientContext {
    double timeStep = 0.0;
    double currentTime = 0.0;
    TransientIntegrationMethod method = TransientIntegrationMethod::BackwardEuler;
    double a0 = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    bool hasSecondHistory = false;
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
