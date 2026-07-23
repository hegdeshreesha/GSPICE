#ifndef GSPICE_DIODE_HPP
#define GSPICE_DIODE_HPP

#include "device.hpp"
#include "fourier.hpp"
#include <algorithm>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdint>

namespace gspice {

class Diode : public Device {
public:
    Diode(
        const std::string& name,
        int nodePos,
        int nodeNeg,
        double Is = 1e-14,
        double emissionCoeff = 1.0,
        double cjo = 0.0)
        : Device(name),
          nodePos_(nodePos),
          nodeNeg_(nodeNeg),
          Is_(Is),
          emissionCoeff_(std::max(emissionCoeff, 1e-6)),
          Vt_(emissionCoeff_ * thermalVoltage_),
          cjo_(std::max(cjo, 0.0)) {}

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        const double rawVoltage = terminalVoltage(x);
        const double limited = std::clamp(rawVoltage, -2.0, 0.8);
        double current = diodeCurrentFromLimitedVoltage(limited) + 1e-12 * limited;
        const double conductance = diodeConductanceFromLimitedVoltage(limited) + 1e-12;
        // Preserve a consistent Newton linearization when the exponential is
        // evaluated at a limited voltage but the global solution still holds
        // the raw proposal.
        current += conductance * (rawVoltage - limited);
        if (request.staticResidual) {
            evaluation.staticResidual.push_back({nodePos_, current});
            evaluation.staticResidual.push_back({nodeNeg_, -current});
        }
        if (request.staticJacobian) {
            evaluation.staticJacobian.push_back({nodePos_, nodePos_, conductance});
            evaluation.staticJacobian.push_back({nodePos_, nodeNeg_, -conductance});
            evaluation.staticJacobian.push_back({nodeNeg_, nodePos_, -conductance});
            evaluation.staticJacobian.push_back({nodeNeg_, nodeNeg_, conductance});
        }
        if (request.dynamicResidual && cjo_ > 0.0) {
            const double charge = cjo_ * rawVoltage;
            evaluation.dynamicResidual.push_back({nodePos_, charge, 0});
            evaluation.dynamicResidual.push_back({nodeNeg_, -charge, 0});
        }
        if (request.dynamicJacobian && cjo_ > 0.0) {
            evaluation.dynamicJacobian.push_back({nodePos_, nodePos_, cjo_, 0});
            evaluation.dynamicJacobian.push_back({nodePos_, nodeNeg_, -cjo_, 0});
            evaluation.dynamicJacobian.push_back({nodeNeg_, nodePos_, -cjo_, 0});
            evaluation.dynamicJacobian.push_back({nodeNeg_, nodeNeg_, cjo_, 0});
        }
        return true;
    }

    bool daeAuditSafe() const override { return true; }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        double Vd = limitedVoltage(x);
        double Id = diodeCurrentFromLimitedVoltage(Vd);
        double gd = diodeConductanceFromLimitedVoltage(Vd);
        double gmin = 1e-12; gd += gmin; Id += gmin * Vd;
        double Ieq = Id - gd * Vd;
        J.add(nodePos_, nodePos_, gd);
        J.add(nodeNeg_, nodeNeg_, gd);
        J.add(nodePos_, nodeNeg_, -gd);
        J.add(nodeNeg_, nodePos_, -gd);
        b.add(nodePos_, -Ieq);
        b.add(nodeNeg_, Ieq);
        if (timeStep > 0.0 && !x_hist.empty() && cjo_ > 0.0) {
            stampBackwardEulerCap(J, b, timeStep, x_hist.back());
        }
    }

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, const TransientContext& ctx) override {
        static const std::vector<VectorReal> empty_history;
        dcStamp(J, b, x, 0.0, ctx.currentTime, empty_history);
        stampIntegratedCap(J, b, ctx);
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        (void)currentTime;
        if (cjo_ <= 0.0 || ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;
        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && prevCapCurrentValid_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }
        const double vNow = terminalVoltage(x);
        const double vPrev = terminalVoltage(ctx.xHistory->back());
        const double vPrev2 = useSecond
            ? terminalVoltage((*ctx.xHistory)[ctx.xHistory->size() - 2]) : vPrev;
        if (useTrap) {
            prevCapCurrent_ = cjo_ * (a0 * vNow + a1 * vPrev) - prevCapCurrent_;
        } else {
            prevCapCurrent_ = cjo_ * (a0 * vNow + a1 * vPrev + (useSecond ? a2 * vPrev2 : 0.0));
        }
        prevCapCurrentValid_ = true;
    }

    std::size_t transientStateBytes() const override {
        return sizeof(prevCapCurrent_) + sizeof(std::uint8_t);
    }

    void saveTransientStateBytes(std::byte* destination, std::size_t size) const override {
        if (size != transientStateBytes()) throw std::invalid_argument("diode transient state size");
        std::memcpy(destination, &prevCapCurrent_, sizeof(prevCapCurrent_));
        const std::uint8_t valid = prevCapCurrentValid_ ? 1u : 0u;
        std::memcpy(destination + sizeof(prevCapCurrent_), &valid, sizeof(valid));
    }

    void restoreTransientStateBytes(const std::byte* source, std::size_t size) override {
        if (size != transientStateBytes()) throw std::invalid_argument("diode transient state size");
        std::memcpy(&prevCapCurrent_, source, sizeof(prevCapCurrent_));
        std::uint8_t valid = 0;
        std::memcpy(&valid, source + sizeof(prevCapCurrent_), sizeof(valid));
        prevCapCurrentValid_ = valid != 0;
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        double Vd = limitedVoltage(x_dc);
        double gd = diodeConductanceFromLimitedVoltage(Vd);
        std::complex<double> y = {gd, omega * cjo_};
        J.add(nodePos_, nodePos_, y);
        J.add(nodeNeg_, nodeNeg_, y);
        J.add(nodePos_, nodeNeg_, -y);
        J.add(nodeNeg_, nodePos_, -y);
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc, std::vector<NoiseSource>& sources) const override {
        (void)omega;
        const double q = 1.602176634e-19;
        const double current = std::abs(diodeCurrentFromLimitedVoltage(limitedVoltage(x_dc)));
        const double psd = 2.0 * q * current;
        if (psd > 0.0) {
            sources.push_back({name_ + ".shot", nodePos_, nodeNeg_, psd});
        }
    }

    double transientChargeError(
        const VectorReal& coarse,
        const VectorReal& fine,
        double reltol,
        double chgtol) override {
        if (cjo_ <= 0.0) return 0.0;
        const double qCoarse = cjo_ * (nodeVoltage(coarse, nodePos_) - nodeVoltage(coarse, nodeNeg_));
        const double qFine = cjo_ * (nodeVoltage(fine, nodePos_) - nodeVoltage(fine, nodeNeg_));
        const double tol = chgtol + reltol * std::max(std::abs(qCoarse), std::abs(qFine));
        return std::abs(qFine - qCoarse) / std::max(tol, 1e-30);
    }

    void limitTransientNewton(const VectorReal& previous, VectorReal& candidate) const override {
        const double oldVoltage = terminalVoltage(previous);
        const double newVoltage = terminalVoltage(candidate);
        applyLimitedTerminalVoltage(candidate, pnjlim(newVoltage, oldVoltage));
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        int K = 2 * n_harms + 1;
        int N_samples = 4 * n_harms; // Over-sampling for FFT accuracy

        // 1. Extract frequency domain node voltages
        std::vector<std::complex<double>> Vd_freq(n_harms + 1);
        Vd_freq[0] = { ((nodePos_ >= 0) ? x_hb[nodePos_ * K] : 0.0) - ((nodeNeg_ >= 0) ? x_hb[nodeNeg_ * K] : 0.0), 0.0 };
        for (int h = 1; h <= n_harms; ++h) {
            double v_cos = ((nodePos_ >= 0) ? x_hb[nodePos_ * K + 2*h - 1] : 0.0) - ((nodeNeg_ >= 0) ? x_hb[nodeNeg_ * K + 2*h - 1] : 0.0);
            double v_sin = ((nodePos_ >= 0) ? x_hb[nodePos_ * K + 2*h] : 0.0) - ((nodeNeg_ >= 0) ? x_hb[nodeNeg_ * K + 2*h] : 0.0);
            Vd_freq[h] = {v_cos, v_sin};
        }

        // 2. IDFT to get time-domain voltages
        std::vector<double> Vd_time = Fourier::inverse(Vd_freq, N_samples);

        // 3. Evaluate non-linear Id and gd at each time sample
        std::vector<double> Id_time(N_samples);
        std::vector<double> gd_time(N_samples);
        for (int i = 0; i < N_samples; ++i) {
            double vd = Vd_time[i];
            if (vd > 0.8) vd = 0.8;
            double expV = std::exp(vd / Vt_);
            Id_time[i] = Is_ * (expV - 1.0);
            gd_time[i] = (Is_ / Vt_) * expV;
        }

        // 4. DFT back to frequency domain
        auto Id_freq = Fourier::forward(Id_time);
        auto gd_freq = Fourier::forward(gd_time);

        // 5. Stamp the HB matrix and RHS
        // DC Component (k=0)
        b.add(nodePos_ * K, -Id_freq[0].real());
        b.add(nodeNeg_ * K, Id_freq[0].real());
        
        // This is a simplified convolution for the Jacobian
        // Real simulators use more complex frequency-shifting here.
        for (int h = 0; h < K; ++h) {
            J.add(nodePos_ * K + h, nodePos_ * K + h, gd_freq[0].real());
        }
    }

private:
    double nodeVoltage(const VectorReal& x, int node) const {
        return node >= 0 ? x[node] : 0.0;
    }

    double limitedVoltage(const VectorReal& x) const {
        double Vd = nodeVoltage(x, nodePos_) - nodeVoltage(x, nodeNeg_);
        if (Vd > 0.8) Vd = 0.8;
        if (Vd < -2.0) Vd = -2.0;
        return Vd;
    }

    double diodeCurrentFromLimitedVoltage(double Vd) const {
        return Is_ * (std::exp(Vd / Vt_) - 1.0);
    }

    double diodeConductanceFromLimitedVoltage(double Vd) const {
        return (Is_ / Vt_) * std::exp(Vd / Vt_);
    }

    double terminalVoltage(const VectorReal& x) const {
        return nodeVoltage(x, nodePos_) - nodeVoltage(x, nodeNeg_);
    }

    double pnjlim(double proposed, double previous) const {
        const double vt = std::max(Vt_, 1e-12);
        const double vcrit = vt * std::log(vt / (std::sqrt(2.0) * std::max(Is_, 1e-30)));
        if (proposed > vcrit && std::abs(proposed - previous) > 2.0 * vt) {
            if (previous > 0.0) {
                const double arg = 1.0 + (proposed - previous) / vt;
                proposed = arg > 0.0 ? previous + vt * std::log(arg) : vcrit;
            } else {
                proposed = vt * std::log(std::max(proposed, vt) / vt);
            }
        }
        return proposed;
    }

    void applyLimitedTerminalVoltage(VectorReal& candidate, double voltage) const {
        const double current = terminalVoltage(candidate);
        const double correction = voltage - current;
        if (nodePos_ >= 0 && nodeNeg_ >= 0) {
            candidate[nodePos_] += 0.5 * correction;
            candidate[nodeNeg_] -= 0.5 * correction;
        } else if (nodePos_ >= 0) {
            candidate[nodePos_] += correction;
        } else if (nodeNeg_ >= 0) {
            candidate[nodeNeg_] -= correction;
        }
    }

    void stampIntegratedCap(SparseMatrixReal& J, VectorReal& b, const TransientContext& ctx) const {
        if (cjo_ <= 0.0 || ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;
        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && prevCapCurrentValid_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }
        const double vPrev = terminalVoltage(ctx.xHistory->back());
        const double vPrev2 = useSecond
            ? terminalVoltage((*ctx.xHistory)[ctx.xHistory->size() - 2]) : vPrev;
        const double geq = cjo_ * a0;
        double ieq = -cjo_ * (a1 * vPrev + (useSecond ? a2 * vPrev2 : 0.0));
        if (useTrap) ieq += prevCapCurrent_;
        J.add(nodePos_, nodePos_, geq);
        J.add(nodeNeg_, nodeNeg_, geq);
        J.add(nodePos_, nodeNeg_, -geq);
        J.add(nodeNeg_, nodePos_, -geq);
        b.add(nodePos_, ieq);
        b.add(nodeNeg_, -ieq);
    }

    void stampBackwardEulerCap(
        SparseMatrixReal& J,
        VectorReal& b,
        double timeStep,
        const VectorReal& x_prev) const {
        const double Geq = cjo_ / std::max(timeStep, 1e-30);
        const double Vprev = nodeVoltage(x_prev, nodePos_) - nodeVoltage(x_prev, nodeNeg_);
        const double Ieq = Geq * Vprev;
        J.add(nodePos_, nodePos_, Geq);
        J.add(nodeNeg_, nodeNeg_, Geq);
        J.add(nodePos_, nodeNeg_, -Geq);
        J.add(nodeNeg_, nodePos_, -Geq);
        b.add(nodePos_, Ieq);
        b.add(nodeNeg_, -Ieq);
    }

    static constexpr double thermalVoltage_ = 0.025852;
    int nodePos_;
    int nodeNeg_;
    double Is_;
    double emissionCoeff_;
    double Vt_;
    double cjo_;
    double prevCapCurrent_ = 0.0;
    bool prevCapCurrentValid_ = false;
};

} // namespace gspice

#endif // GSPICE_DIODE_HPP
