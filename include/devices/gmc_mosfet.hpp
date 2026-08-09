#ifndef GSPICE_GMC_MOSFET_HPP
#define GSPICE_GMC_MOSFET_HPP

#include "devices/mosfet.hpp"
#include "gdi.hpp"
#include "gsdi.hpp"

#include <memory>
#include <vector>

namespace gspice {

// Native GMC bridge for the existing primitive MOS equations. The local
// four-terminal Mosfet keeps the GDI model independent of global node IDs.
class GmcMosfetInstance final : public GdiInstanceBase {
public:
    GmcMosfetInstance(
        int type,
        double width,
        double length,
        double threshold,
        double kp,
        double lambda,
        double gamma,
        double phi)
        : model_("GMC_MOS", 0, 1, 2, 3, type, width, length,
                 threshold, kp, lambda, gamma, phi) {}

    bool evaluate(
        const double* voltages,
        const GdiEvaluationRequest& request,
        GdiEvaluationResult& result) override {
        result.clear();
        VectorReal x(4);
        for (int i = 0; i < 4; ++i) x[i] = voltages[i];

        DaeRequest dae_request;
        dae_request.analysis = request.dynamic_residual || request.dynamic_jacobian
            ? DaeAnalysis::Transient : DaeAnalysis::OperatingPoint;
        dae_request.time = request.time;
        dae_request.staticResidual = request.static_residual;
        dae_request.staticJacobian = request.static_jacobian;
        dae_request.dynamicResidual = request.dynamic_residual;
        dae_request.dynamicJacobian = request.dynamic_jacobian;
        dae_request.allowBypass = request.allow_bypass;

        DaeEvaluation evaluation;
        if (!model_.evaluateDae(x, dae_request, evaluation)) return false;
        for (const auto& term : evaluation.staticResidual) {
            result.static_residual.push_back({term.equation, term.value});
        }
        for (const auto& term : evaluation.staticJacobian) {
            result.static_jacobian.push_back({term.equation, term.unknown, term.value});
        }
        for (const auto& term : evaluation.dynamicResidual) {
            result.dynamic_residual.push_back({term.equation, term.value});
        }
        for (const auto& term : evaluation.dynamicJacobian) {
            result.dynamic_jacobian.push_back({term.equation, term.unknown, term.value});
        }
        result.limiting_applied = evaluation.limitingApplied;
        result.bypassed = evaluation.bypassed;
        return true;
    }

    std::size_t numTerminals() const override { return 4; }
    std::size_t numInternalNodes() const override { return 0; }
    std::size_t numStates() const override { return 0; }

private:
    Mosfet model_;
};

} // namespace gspice

#endif // GSPICE_GMC_MOSFET_HPP
