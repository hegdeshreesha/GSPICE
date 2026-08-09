#ifndef GSPICE_MOSVAR_CAPACITOR_HPP
#define GSPICE_MOSVAR_CAPACITOR_HPP

#include "device.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace gspice {

class MosvarCapacitor : public Device {
public:
    MosvarCapacitor(
        const std::string& name,
        int nodeGate,
        int nodeWell,
        double cAccum,
        double cMinimum,
        double vfb,
        double slope)
        : Device(name),
          nodeGate_(nodeGate),
          nodeWell_(nodeWell),
          cAccum_(std::max(cAccum, 1e-24)),
          cMinimum_(std::clamp(cMinimum, 0.0, std::max(cAccum, 1e-24))),
          vfb_(vfb),
          slope_(std::max(std::abs(slope), 1e-3)) {}

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        const double v = voltage(x);
        const double q = charge(v);
        const double c = capacitance(v);
        if (request.dynamicResidual) {
            evaluation.dynamicResidual.push_back({nodeGate_, q, 0});
            evaluation.dynamicResidual.push_back({nodeWell_, -q, 0});
        }
        if (request.dynamicJacobian) {
            evaluation.dynamicJacobian.push_back({nodeGate_, nodeGate_, c, 0});
            evaluation.dynamicJacobian.push_back({nodeGate_, nodeWell_, -c, 0});
            evaluation.dynamicJacobian.push_back({nodeWell_, nodeGate_, -c, 0});
            evaluation.dynamicJacobian.push_back({nodeWell_, nodeWell_, c, 0});
        }
        return true;
    }

    bool daeAuditSafe() const override { return true; }

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)J;
        (void)b;
        (void)x;
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        DaeRequest request;
        request.analysis = DaeAnalysis::SmallSignal;
        request.dynamicJacobian = true;
        DaeEvaluation evaluation;
        evaluateDae(x_dc, request, evaluation);
        stampDaeSmallSignal(evaluation, omega, J);
    }

    double transientChargeError(
        const VectorReal& coarse,
        const VectorReal& fine,
        double reltol,
        double chgtol) override {
        const double qCoarse = charge(voltage(coarse));
        const double qFine = charge(voltage(fine));
        const double tol = chgtol + reltol * std::max(std::abs(qCoarse), std::abs(qFine));
        return std::abs(qFine - qCoarse) / std::max(tol, 1e-30);
    }

private:
    static double softplus(double x) {
        if (x > 40.0) return x;
        if (x < -40.0) return std::exp(x);
        return std::log1p(std::exp(x));
    }

    double voltage(const VectorReal& x) const {
        return ((nodeGate_ >= 0) ? x[nodeGate_] : 0.0) -
               ((nodeWell_ >= 0) ? x[nodeWell_] : 0.0);
    }

    double capacitance(double v) const {
        const double s = 1.0 / (1.0 + std::exp(std::clamp((v - vfb_) / slope_, -40.0, 40.0)));
        return cMinimum_ + (cAccum_ - cMinimum_) * s;
    }

    double charge(double v) const {
        const double delta = cAccum_ - cMinimum_;
        return cMinimum_ * v - delta * slope_ * softplus(-(v - vfb_) / slope_);
    }

    int nodeGate_;
    int nodeWell_;
    double cAccum_;
    double cMinimum_;
    double vfb_;
    double slope_;
};

} // namespace gspice

#endif // GSPICE_MOSVAR_CAPACITOR_HPP
