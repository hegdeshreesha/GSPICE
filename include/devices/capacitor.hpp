#ifndef GSPICE_CAPACITOR_HPP
#define GSPICE_CAPACITOR_HPP

#include "device.hpp"
#include <string>
#include <cstring>
#include <cstdint>

namespace gspice {

class Capacitor : public Device {
public:
    Capacitor(const std::string& name, int nodePos, int nodeNeg, double value)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value) {}

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        if (request.dynamicResidual) {
            const double charge = value_ * voltage(x);
            evaluation.dynamicResidual.push_back({nodePos_, charge, 0});
            evaluation.dynamicResidual.push_back({nodeNeg_, -charge, 0});
        }
        if (request.dynamicJacobian) {
            evaluation.dynamicJacobian.push_back({nodePos_, nodePos_, value_, 0});
            evaluation.dynamicJacobian.push_back({nodePos_, nodeNeg_, -value_, 0});
            evaluation.dynamicJacobian.push_back({nodeNeg_, nodePos_, -value_, 0});
            evaluation.dynamicJacobian.push_back({nodeNeg_, nodeNeg_, value_, 0});
        }
        return true;
    }

    bool daeAuditSafe() const override { return true; }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        (void)currentTime;
        if (timeStep <= 0) return;
        if (x_hist.empty()) return;
        DaeRequest request;
        request.analysis = DaeAnalysis::Transient;
        request.staticResidual = false;
        request.staticJacobian = false;
        request.dynamicResidual = true;
        request.dynamicJacobian = true;
        DaeEvaluation current;
        DaeEvaluation previous;
        evaluateDae(x, request, current);
        evaluateDae(x_hist.back(), request, previous);
        DaeHistory history;
        appendScaledDaeResidual(history, previous.dynamicResidual, -1.0 / timeStep);
        stampDaeTransient(current, x, 1.0 / timeStep, history, J, b);
    }

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, const TransientContext& ctx) override {
        (void)x;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;

        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && prevCurrentValid_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }

        DaeRequest request;
        request.analysis = DaeAnalysis::Transient;
        request.time = ctx.currentTime;
        request.staticResidual = false;
        request.staticJacobian = false;
        request.dynamicResidual = true;
        request.dynamicJacobian = true;
        DaeEvaluation current;
        DaeEvaluation previous;
        DaeEvaluation previous2;
        evaluateDae(x, request, current);
        evaluateDae((*ctx.xHistory)[ctx.xHistory->size() - 1], request, previous);
        DaeHistory history;
        appendScaledDaeResidual(history, previous.dynamicResidual, a1);
        if (useSecond) {
            evaluateDae((*ctx.xHistory)[ctx.xHistory->size() - 2], request, previous2);
            appendScaledDaeResidual(history, previous2.dynamicResidual, a2);
        }
        if (useTrap) {
            history.push_back({nodePos_, -prevCurrent_});
            history.push_back({nodeNeg_, prevCurrent_});
        }
        stampDaeTransient(current, x, a0, history, J, b);
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        (void)currentTime;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;

        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && prevCurrentValid_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }

        const double vNow = voltage(x);
        const double vPrev = voltage((*ctx.xHistory)[ctx.xHistory->size() - 1]);
        const double vPrev2 = useSecond ? voltage((*ctx.xHistory)[ctx.xHistory->size() - 2]) : vPrev;
        if (useTrap) {
            prevCurrent_ = value_ * (a0 * vNow + a1 * vPrev) - prevCurrent_;
        } else {
            prevCurrent_ = value_ * (a0 * vNow + a1 * vPrev + (useSecond ? a2 * vPrev2 : 0.0));
        }
        prevCurrentValid_ = true;
    }

    std::size_t transientStateBytes() const override {
        return sizeof(prevCurrent_) + sizeof(std::uint8_t);
    }

    void saveTransientStateBytes(std::byte* destination, std::size_t size) const override {
        if (size != transientStateBytes()) throw std::invalid_argument("capacitor transient state size");
        std::memcpy(destination, &prevCurrent_, sizeof(prevCurrent_));
        const std::uint8_t valid = prevCurrentValid_ ? 1u : 0u;
        std::memcpy(destination + sizeof(prevCurrent_), &valid, sizeof(valid));
    }

    void restoreTransientStateBytes(const std::byte* source, std::size_t size) override {
        if (size != transientStateBytes()) throw std::invalid_argument("capacitor transient state size");
        std::memcpy(&prevCurrent_, source, sizeof(prevCurrent_));
        std::uint8_t valid = 0;
        std::memcpy(&valid, source + sizeof(prevCurrent_), sizeof(valid));
        prevCurrentValid_ = valid != 0;
    }

    double transientChargeError(
        const VectorReal& coarse,
        const VectorReal& fine,
        double reltol,
        double chgtol) override {
        const double qCoarse = value_ * voltage(coarse);
        const double qFine = value_ * voltage(fine);
        const double tol = chgtol + reltol * std::max(std::abs(qCoarse), std::abs(qFine));
        return std::abs(qFine - qCoarse) / std::max(tol, 1e-30);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        DaeRequest request;
        request.analysis = DaeAnalysis::SmallSignal;
        request.staticResidual = false;
        request.staticJacobian = false;
        request.dynamicResidual = false;
        request.dynamicJacobian = true;
        DaeEvaluation evaluation;
        evaluateDae(x_dc, request, evaluation);
        stampDaeSmallSignal(evaluation, omega, J);
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        int K = 2 * n_harms + 1;
        for (int h = 1; h <= n_harms; ++h) {
            double omega_h = 2.0 * 3.14159265358979 * f_fund * h;
            double val = omega_h * value_;
            int cos_idx = 2 * h - 1; int sin_idx = 2 * h;
            J.add(nodePos_ * K + cos_idx, nodePos_ * K + sin_idx, -val);
            J.add(nodePos_ * K + sin_idx, nodePos_ * K + cos_idx, val);
            J.add(nodeNeg_ * K + cos_idx, nodeNeg_ * K + sin_idx, val);
            J.add(nodeNeg_ * K + sin_idx, nodeNeg_ * K + cos_idx, -val);
            J.add(nodePos_ * K + cos_idx, nodeNeg_ * K + sin_idx, val);
            J.add(nodePos_ * K + sin_idx, nodeNeg_ * K + cos_idx, -val);
            J.add(nodeNeg_ * K + cos_idx, nodePos_ * K + sin_idx, -val);
            J.add(nodeNeg_ * K + sin_idx, nodePos_ * K + cos_idx, val);
        }
    }

    void pacStamp(SparseMatrixReal& J, VectorReal& b, double f_in, double f_fund, int n_harms, const VectorReal& x_periodic) override {
        int K = 2 * n_harms + 1;
        // Sideband frequencies: f_k = f_in + k*f_fund
        for (int k = -n_harms; k <= n_harms; ++k) {
            double fk = f_in + k * f_fund;
            double omega_k = 2.0 * 3.14159265358979 * fk;
            double val = omega_k * value_;
            
            // Map k to our 2*N+1 indices (0=DC, 1=H1_Cos, 2=H1_Sin, ...)
            int idx;
            if (k == 0) idx = 0;
            else if (k > 0) idx = 2 * k - 1;
            else idx = 2 * std::abs(k);

            // Stamp complex admittance val into the conversion matrix
            // This is a simplified version (diagonal sidebands only)
            J.add(nodePos_ * K + idx, nodePos_ * K + idx, val); // Imaginary parts would need cross-terms
        }
    }

private:
    double voltage(const VectorReal& x) const {
        return ((nodePos_ >= 0) ? x[nodePos_] : 0.0) - ((nodeNeg_ >= 0) ? x[nodeNeg_] : 0.0);
    }

    int nodePos_;
    int nodeNeg_;
    double value_;
    double prevCurrent_ = 0.0;
    bool prevCurrentValid_ = false;
};

} // namespace gspice

#endif // GSPICE_CAPACITOR_HPP
