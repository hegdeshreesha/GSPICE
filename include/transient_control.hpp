#ifndef GSPICE_TRANSIENT_CONTROL_HPP
#define GSPICE_TRANSIENT_CONTROL_HPP

#include "integration_formula.hpp"
#include "matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace gspice {

struct PolynomialPrediction {
    VectorReal value;
    std::vector<double> weights;
    int order = 0;
    bool valid = false;
};

inline std::vector<double> polynomialExtrapolationWeights(
    const std::vector<double>& chronologicalTimes,
    double targetTime,
    int order) {
    if (order < 1 || chronologicalTimes.size() < static_cast<std::size_t>(order + 1)) {
        return {};
    }
    const std::size_t first = chronologicalTimes.size() - static_cast<std::size_t>(order + 1);
    std::vector<double> weights(static_cast<std::size_t>(order + 1), 1.0);
    const double scale = std::max(
        std::abs(targetTime - chronologicalTimes.back()),
        std::numeric_limits<double>::min());
    for (int i = 0; i <= order; ++i) {
        const double ti = (chronologicalTimes[first + static_cast<std::size_t>(i)] - targetTime) / scale;
        for (int j = 0; j <= order; ++j) {
            if (i == j) continue;
            const double tj = (chronologicalTimes[first + static_cast<std::size_t>(j)] - targetTime) / scale;
            const double denominator = ti - tj;
            if (!std::isfinite(denominator) || std::abs(denominator) < 1e-14) return {};
            weights[static_cast<std::size_t>(i)] *= -tj / denominator;
        }
    }
    return weights;
}

inline PolynomialPrediction polynomialPredict(
    const std::vector<VectorReal>& chronologicalValues,
    const std::vector<double>& chronologicalTimes,
    double targetTime,
    int order) {
    PolynomialPrediction result;
    if (chronologicalValues.size() != chronologicalTimes.size() || chronologicalValues.empty()) {
        return result;
    }
    result.weights = polynomialExtrapolationWeights(chronologicalTimes, targetTime, order);
    if (result.weights.empty()) return result;
    const std::size_t first = chronologicalValues.size() - result.weights.size();
    const int vectorSize = chronologicalValues.back().getSize();
    result.value = VectorReal(vectorSize);
    for (int component = 0; component < vectorSize; ++component) {
        double predicted = 0.0;
        for (std::size_t i = 0; i < result.weights.size(); ++i) {
            if (chronologicalValues[first + i].getSize() != vectorSize) return PolynomialPrediction{};
            predicted += result.weights[i] * chronologicalValues[first + i][component];
        }
        if (!std::isfinite(predicted) || std::abs(predicted) > 1e100) return PolynomialPrediction{};
        result.value[component] = predicted;
    }
    result.order = order;
    result.valid = true;
    return result;
}

// For a method of order p, both the implicit corrector and a degree-p
// polynomial extrapolator first fail on a polynomial of degree p+1.  Applying
// both algorithms to that normalized polynomial yields the exact multiplier
// that maps (corrector - predictor) to the corrector's local error.  This is
// valid for nonuniform history and avoids hard-coded uniform-step constants.
inline double predictorCorrectorErrorFactor(
    const IntegrationFormula& corrector,
    const std::vector<double>& newestFirstFormulaTimes,
    const std::vector<double>& chronologicalPredictionTimes,
    double targetTime) {
    const int order = corrector.order;
    if (order < 1 || corrector.qWeights.empty() ||
        newestFirstFormulaTimes.size() < corrector.qWeights.size() ||
        chronologicalPredictionTimes.size() < static_cast<std::size_t>(order + 1)) {
        return 0.0;
    }
    const double h = targetTime - chronologicalPredictionTimes.back();
    if (!(h > 0.0) || !std::isfinite(h)) return 0.0;
    const auto monomial = [order](double u) {
        return std::pow(u, order + 1);
    };
    const auto derivative = [order, h](double u) {
        return static_cast<double>(order + 1) * std::pow(u, order) / h;
    };

    double residual = 0.0;
    for (std::size_t i = 0; i < corrector.qWeights.size(); ++i) {
        const double u = (newestFirstFormulaTimes[i] - targetTime) / h;
        residual += corrector.qWeights[i] * monomial(u);
    }
    for (std::size_t i = 0; i < corrector.derivativeWeights.size(); ++i) {
        const std::size_t timeIndex = i + 1;
        if (timeIndex >= newestFirstFormulaTimes.size()) return 0.0;
        const double u = (newestFirstFormulaTimes[timeIndex] - targetTime) / h;
        residual += corrector.derivativeWeights[i] * derivative(u);
    }
    residual -= derivative(0.0);
    const double correctorError = -residual / corrector.qWeights.front();

    const auto predictorWeights = polynomialExtrapolationWeights(
        chronologicalPredictionTimes, targetTime, order);
    if (predictorWeights.empty()) return 0.0;
    const std::size_t first = chronologicalPredictionTimes.size() - predictorWeights.size();
    double predictorError = 0.0;
    for (std::size_t i = 0; i < predictorWeights.size(); ++i) {
        const double u = (chronologicalPredictionTimes[first + i] - targetTime) / h;
        predictorError += predictorWeights[i] * monomial(u);
    }
    const double denominator = correctorError - predictorError;
    if (!std::isfinite(denominator) || std::abs(denominator) < 1e-30) return 0.0;
    return std::clamp(std::abs(correctorError / denominator), 1e-6, 1.0);
}

struct AdaptiveOrderDecision {
    int order = 1;
    double projectedFactor = 1.0;
};

inline double projectedStepFactor(double normalizedError, int order) {
    if (!std::isfinite(normalizedError) || normalizedError <= 1e-30) return 2.5;
    return std::clamp(
        0.9 * std::pow(normalizedError, -1.0 / static_cast<double>(order + 1)),
        0.2,
        2.5);
}

inline AdaptiveOrderDecision chooseAdaptiveOrder(
    int currentOrder,
    int maximumOrder,
    double currentError,
    const std::optional<double>& lowerError,
    const std::optional<double>& higherError) {
    AdaptiveOrderDecision best{currentOrder, projectedStepFactor(currentError, currentOrder)};
    if (currentOrder > 1 && lowerError && std::isfinite(*lowerError)) {
        const double factor = projectedStepFactor(*lowerError, currentOrder - 1);
        if (currentError > 0.8 || factor > best.projectedFactor * 1.20) {
            best = {currentOrder - 1, factor};
        }
    }
    if (currentOrder < maximumOrder && higherError && std::isfinite(*higherError)) {
        const double factor = projectedStepFactor(*higherError, currentOrder + 1);
        // Hysteresis prevents order chatter when two formulas predict nearly
        // identical work.  Raising order also requires comfortable LTE margin.
        const bool asymptoticallySmooth = currentError < 0.05 && *higherError < 0.05 &&
            factor >= best.projectedFactor * 0.95;
        if ((*higherError < 0.8 && factor > best.projectedFactor * 1.12) ||
            asymptoticallySmooth) {
            best = {currentOrder + 1, factor};
        }
    }
    return best;
}

inline bool detectTrapezoidalRinging(
    const std::vector<VectorReal>& chronologicalHistory,
    int voltageUnknowns,
    double absoluteTolerance,
    double relativeTolerance) {
    if (chronologicalHistory.size() < 4) return false;
    const VectorReal& x0 = chronologicalHistory[chronologicalHistory.size() - 4];
    const VectorReal& x1 = chronologicalHistory[chronologicalHistory.size() - 3];
    const VectorReal& x2 = chronologicalHistory[chronologicalHistory.size() - 2];
    const VectorReal& x3 = chronologicalHistory[chronologicalHistory.size() - 1];
    const int n = std::min({voltageUnknowns, x0.getSize(), x1.getSize(), x2.getSize(), x3.getSize()});
    for (int i = 0; i < n; ++i) {
        const double d0 = x1[i] - x0[i];
        const double d1 = x2[i] - x1[i];
        const double d2 = x3[i] - x2[i];
        const double scale = std::max({std::abs(x0[i]), std::abs(x1[i]), std::abs(x2[i]), std::abs(x3[i])});
        const double floor = 4.0 * (absoluteTolerance + relativeTolerance * scale);
        if (std::abs(d1) <= floor || std::abs(d2) <= floor) continue;
        const bool alternating = d0 * d1 < 0.0 && d1 * d2 < 0.0;
        const bool weaklyDamped = std::abs(d2) >= 0.70 * std::abs(d1);
        if (alternating && weaklyDamped) return true;
    }
    return false;
}

class AutomaticTransientMethodController {
public:
    explicit AutomaticTransientMethodController(bool permanentlyDamped = false)
        : permanentlyDamped_(permanentlyDamped) {}

    bool useTrapezoidal() const {
        return !permanentlyDamped_ && dampedStepsRemaining_ == 0;
    }

    bool observe(
        const std::vector<VectorReal>& history,
        int voltageUnknowns,
        double absoluteTolerance,
        double relativeTolerance) {
        if (permanentlyDamped_) return false;
        const bool ringing = detectTrapezoidalRinging(
            history, voltageUnknowns, absoluteTolerance, relativeTolerance);
        if (ringing) {
            const bool switched = dampedStepsRemaining_ == 0;
            dampedStepsRemaining_ = 6;
            stableSteps_ = 0;
            return switched;
        }
        if (dampedStepsRemaining_ > 0) {
            ++stableSteps_;
            if (stableSteps_ >= 4) {
                dampedStepsRemaining_ = 0;
                stableSteps_ = 0;
                return true;
            }
            --dampedStepsRemaining_;
        }
        return false;
    }

    void restart() {
        dampedStepsRemaining_ = 0;
        stableSteps_ = 0;
    }

private:
    bool permanentlyDamped_ = false;
    int dampedStepsRemaining_ = 0;
    int stableSteps_ = 0;
};

} // namespace gspice

#endif // GSPICE_TRANSIENT_CONTROL_HPP
