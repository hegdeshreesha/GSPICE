#ifndef GSPICE_DAE_AUDIT_HPP
#define GSPICE_DAE_AUDIT_HPP

#include "device.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace gspice {

struct DaeAuditOptions {
    double relativeTolerance = 2e-4;
    double staticAbsoluteTolerance = 1e-12;
    double dynamicAbsoluteTolerance = 1e-18;
    double perturbationScale = 2e-6;
};

struct DaeAuditReport {
    bool supported = false;
    bool finite = true;
    std::size_t staticDerivativesChecked = 0;
    std::size_t dynamicDerivativesChecked = 0;
    double worstStaticDerivative = 0.0;
    double worstDynamicDerivative = 0.0;
    int worstStaticEquation = -1;
    int worstStaticUnknown = -1;
    int worstDynamicEquation = -1;
    int worstDynamicUnknown = -1;
    double worstStaticNumerical = 0.0;
    double worstStaticExpected = 0.0;
    double worstDynamicNumerical = 0.0;
    double worstDynamicExpected = 0.0;
    double worstChargeImbalance = 0.0;
    double worstChargeJacobianImbalance = 0.0;

    bool passed() const {
        return supported && finite && worstStaticDerivative <= 1.0 &&
            worstDynamicDerivative <= 1.0 && worstChargeImbalance <= 1.0 &&
            worstChargeJacobianImbalance <= 1.0;
    }
};

namespace detail {

using DaeMatrixKey = std::pair<int, int>;

class DaeAuditStateGuard {
public:
    explicit DaeAuditStateGuard(Device& device)
        : device_(device), bytes_(device.transientStateBytes()) {
        if (!bytes_.empty()) device_.saveTransientStateBytes(bytes_.data(), bytes_.size());
    }

    ~DaeAuditStateGuard() noexcept {
        try {
            restore();
        } catch (...) {
        }
    }

    void restore() {
        if (!bytes_.empty()) device_.restoreTransientStateBytes(bytes_.data(), bytes_.size());
    }

private:
    Device& device_;
    std::vector<std::byte> bytes_;
};

inline std::map<int, double> aggregateResidual(const std::vector<DaeResidualTerm>& terms) {
    std::map<int, double> result;
    for (const auto& term : terms) result[term.equation] += term.value;
    return result;
}

inline std::map<DaeMatrixKey, double> aggregateJacobian(
    const std::vector<DaeJacobianTerm>& terms) {
    std::map<DaeMatrixKey, double> result;
    for (const auto& term : terms) result[{term.equation, term.unknown}] += term.value;
    return result;
}

inline double normalizedDifference(double actual, double expected, double relative, double absolute) {
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return std::abs(actual - expected) / std::max(absolute + relative * scale, 1e-300);
}

} // namespace detail

inline DaeAuditReport auditDaeDevice(
    Device& device,
    const VectorReal& solution,
    const DaeAuditOptions& options = {}) {
    DaeAuditReport report;
    if (!device.daeAuditSafe()) return report;
    detail::DaeAuditStateGuard stateGuard(device);

    DaeRequest fullRequest;
    fullRequest.analysis = DaeAnalysis::Transient;
    fullRequest.dynamicResidual = true;
    fullRequest.dynamicJacobian = true;
    fullRequest.enableLimiting = false;
    fullRequest.readOnlyState = true;
    fullRequest.highPrecision = true;
    DaeEvaluation baseline;
    stateGuard.restore();
    if (!device.evaluateDae(solution, fullRequest, baseline)) return report;
    report.supported = true;
    report.finite = baseline.finite();
    if (!report.finite) return report;

    const auto staticJacobian = detail::aggregateJacobian(baseline.staticJacobian);
    const auto dynamicJacobian = detail::aggregateJacobian(baseline.dynamicJacobian);
    std::set<int> unknowns;
    std::set<int> equations;
    for (const auto& entry : staticJacobian) {
        equations.insert(entry.first.first);
        unknowns.insert(entry.first.second);
    }
    for (const auto& entry : dynamicJacobian) {
        equations.insert(entry.first.first);
        unknowns.insert(entry.first.second);
    }

    // Keep Jacobian flags enabled while finite-differencing residuals. Some
    // compact-model ABIs share intermediate calculations between these output
    // requests; the numerical derivative remains independent because the
    // returned perturbed Jacobians are ignored.
    DaeRequest residualRequest = fullRequest;
    for (int unknown : unknowns) {
        if (unknown < 0 || unknown >= solution.getSize()) continue;
        VectorReal plus = solution;
        VectorReal minus = solution;
        const double delta = options.perturbationScale * std::max(1.0, std::abs(solution[unknown]));
        plus[unknown] += delta;
        minus[unknown] -= delta;
        DaeEvaluation plusEvaluation;
        DaeEvaluation minusEvaluation;
        stateGuard.restore();
        const bool plusSupported = device.evaluateDae(plus, residualRequest, plusEvaluation);
        stateGuard.restore();
        const bool minusSupported = device.evaluateDae(minus, residualRequest, minusEvaluation);
        if (!plusSupported || !minusSupported) {
            report.supported = false;
            return report;
        }
        const auto plusStatic = detail::aggregateResidual(plusEvaluation.staticResidual);
        const auto minusStatic = detail::aggregateResidual(minusEvaluation.staticResidual);
        const auto plusDynamic = detail::aggregateResidual(plusEvaluation.dynamicResidual);
        const auto minusDynamic = detail::aggregateResidual(minusEvaluation.dynamicResidual);
        for (int equation : equations) {
            const double numericalStatic =
                (plusStatic.count(equation) ? plusStatic.at(equation) : 0.0) -
                (minusStatic.count(equation) ? minusStatic.at(equation) : 0.0);
            const double numericalDynamic =
                (plusDynamic.count(equation) ? plusDynamic.at(equation) : 0.0) -
                (minusDynamic.count(equation) ? minusDynamic.at(equation) : 0.0);
            const double staticDerivative = numericalStatic / (2.0 * delta);
            const double dynamicDerivative = numericalDynamic / (2.0 * delta);
            const auto staticIt = staticJacobian.find({equation, unknown});
            const auto dynamicIt = dynamicJacobian.find({equation, unknown});
            const double expectedStatic = staticIt == staticJacobian.end() ? 0.0 : staticIt->second;
            const double expectedDynamic = dynamicIt == dynamicJacobian.end() ? 0.0 : dynamicIt->second;
            const double staticError = detail::normalizedDifference(
                staticDerivative, expectedStatic,
                options.relativeTolerance, options.staticAbsoluteTolerance);
            const double dynamicError = detail::normalizedDifference(
                dynamicDerivative, expectedDynamic,
                options.relativeTolerance, options.dynamicAbsoluteTolerance);
            if (staticError > report.worstStaticDerivative) {
                report.worstStaticDerivative = staticError;
                report.worstStaticEquation = equation;
                report.worstStaticUnknown = unknown;
                report.worstStaticNumerical = staticDerivative;
                report.worstStaticExpected = expectedStatic;
            }
            if (dynamicError > report.worstDynamicDerivative) {
                report.worstDynamicDerivative = dynamicError;
                report.worstDynamicEquation = equation;
                report.worstDynamicUnknown = unknown;
                report.worstDynamicNumerical = dynamicDerivative;
                report.worstDynamicExpected = expectedDynamic;
            }
            ++report.staticDerivativesChecked;
            ++report.dynamicDerivativesChecked;
        }
    }

    std::map<int, double> groupCharge;
    std::map<std::pair<int, int>, double> groupJacobian;
    std::map<int, double> groupMagnitude;
    std::map<std::pair<int, int>, double> groupJacobianMagnitude;
    for (const auto& term : baseline.dynamicResidual) {
        if (term.conservationGroup < 0) continue;
        groupCharge[term.conservationGroup] += term.value;
        groupMagnitude[term.conservationGroup] += std::abs(term.value);
    }
    for (const auto& term : baseline.dynamicJacobian) {
        if (term.conservationGroup < 0) continue;
        const auto key = std::make_pair(term.conservationGroup, term.unknown);
        groupJacobian[key] += term.value;
        groupJacobianMagnitude[key] += std::abs(term.value);
    }
    for (const auto& group : groupCharge) {
        const double tolerance = options.dynamicAbsoluteTolerance +
            options.relativeTolerance * groupMagnitude[group.first];
        report.worstChargeImbalance = std::max(
            report.worstChargeImbalance, std::abs(group.second) / std::max(tolerance, 1e-300));
    }
    for (const auto& group : groupJacobian) {
        const double tolerance = options.dynamicAbsoluteTolerance +
            options.relativeTolerance * groupJacobianMagnitude[group.first];
        report.worstChargeJacobianImbalance = std::max(
            report.worstChargeJacobianImbalance,
            std::abs(group.second) / std::max(tolerance, 1e-300));
    }
    return report;
}

} // namespace gspice

#endif // GSPICE_DAE_AUDIT_HPP
