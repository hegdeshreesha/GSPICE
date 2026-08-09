#ifndef GSPICE_GSDI_DEVICE_ADAPTER_HPP
#define GSPICE_GSDI_DEVICE_ADAPTER_HPP

#include "device.hpp"
#include "gsdi.hpp"

#include <memory>

namespace gspice {

class GsdiDaeDeviceAdapter final : public GsdiInstance {
public:
    GsdiDaeDeviceAdapter(std::unique_ptr<Device> device, std::size_t terminal_count)
        : device_(std::move(device)), terminal_count_(terminal_count) {}

    bool evaluate(const GsdiEvalRequest& request, GsdiEvalResult& result) override {
        if (!device_ || !request.solution) return false;
        VectorReal x(static_cast<int>(request.solution_size));
        for (std::size_t i = 0; i < request.solution_size; ++i) {
            x[static_cast<int>(i)] = request.solution[i];
        }

        DaeRequest dae;
        dae.analysis = request.analysis == GsdiAnalysis::Transient ? DaeAnalysis::Transient :
                       request.analysis == GsdiAnalysis::Ac ? DaeAnalysis::SmallSignal :
                       DaeAnalysis::OperatingPoint;
        dae.time = request.time;
        dae.staticResidual = request.residual;
        dae.staticJacobian = request.jacobian;
        dae.dynamicResidual = request.dynamic_residual;
        dae.dynamicJacobian = request.dynamic_jacobian;
        dae.readOnlyState = request.read_only_state;
        dae.allowBypass = request.allow_bypass;

        DaeEvaluation evaluation;
        if (!device_->evaluateDae(x, dae, evaluation)) return false;
        result.clear();
        for (const auto& term : evaluation.staticResidual) {
            result.static_residual.push_back({term.equation, term.value, term.conservationGroup});
        }
        for (const auto& term : evaluation.staticJacobian) {
            result.static_jacobian.push_back({term.equation, term.unknown, term.value, term.conservationGroup});
        }
        for (const auto& term : evaluation.dynamicResidual) {
            result.dynamic_residual.push_back({term.equation, term.value, term.conservationGroup});
        }
        for (const auto& term : evaluation.dynamicJacobian) {
            result.dynamic_jacobian.push_back({term.equation, term.unknown, term.value, term.conservationGroup});
        }
        if (request.noise) {
            std::vector<NoiseSource> sources;
            device_->collectNoiseSources(request.omega, x, sources);
            for (const auto& source : sources) {
                result.noise.push_back({source.nodePos, source.nodeNeg,
                                        source.currentPsd, source.name});
            }
        }
        result.limiting_applied = evaluation.limitingApplied;
        result.bypassed = evaluation.bypassed;
        return true;
    }

    std::size_t terminalCount() const override { return terminal_count_; }

private:
    std::unique_ptr<Device> device_;
    std::size_t terminal_count_ = 0;
};

} // namespace gspice

#endif // GSPICE_GSDI_DEVICE_ADAPTER_HPP
