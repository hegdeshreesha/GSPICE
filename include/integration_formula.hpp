#ifndef GSPICE_INTEGRATION_FORMULA_HPP
#define GSPICE_INTEGRATION_FORMULA_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace gspice {

enum class IntegrationFamily {
    Bdf,
    AdamsMoulton
};

// Represents
//   derivative(now) = sum(qWeights[i] * q[i])
//                   + sum(derivativeWeights[i] * derivativeHistory[i])
// where q[0] is the new unknown Q and q[1...] are accepted Q values.
struct IntegrationFormula {
    IntegrationFamily family = IntegrationFamily::Bdf;
    int order = 1;
    std::vector<double> qWeights;
    std::vector<double> derivativeWeights;

    double leading() const {
        return qWeights.empty() ? 0.0 : qWeights.front();
    }

    double differentiate(
        double current,
        const std::vector<double>& qHistory,
        const std::vector<double>& derivativeHistory = {}) const {
        if (qWeights.empty()) throw std::logic_error("empty integration formula");
        if (qHistory.size() + 1 < qWeights.size() ||
            derivativeHistory.size() < derivativeWeights.size()) {
            throw std::invalid_argument("insufficient integration history");
        }
        double result = qWeights[0] * current;
        for (std::size_t i = 1; i < qWeights.size(); ++i) {
            result += qWeights[i] * qHistory[i - 1];
        }
        for (std::size_t i = 0; i < derivativeWeights.size(); ++i) {
            result += derivativeWeights[i] * derivativeHistory[i];
        }
        return result;
    }
};

namespace integration_detail {

inline std::vector<double> solveDense(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs) {
    const std::size_t size = rhs.size();
    if (matrix.size() != size) throw std::invalid_argument("integration matrix size mismatch");
    for (std::size_t column = 0; column < size; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < size; ++row) {
            if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
        }
        if (std::abs(matrix[pivot][column]) < 1e-14) {
            throw std::invalid_argument("singular integration time history");
        }
        std::swap(matrix[pivot], matrix[column]);
        std::swap(rhs[pivot], rhs[column]);
        const double diagonal = matrix[column][column];
        for (std::size_t entry = column; entry < size; ++entry) matrix[column][entry] /= diagonal;
        rhs[column] /= diagonal;
        for (std::size_t row = 0; row < size; ++row) {
            if (row == column) continue;
            const double factor = matrix[row][column];
            if (factor == 0.0) continue;
            for (std::size_t entry = column; entry < size; ++entry) {
                matrix[row][entry] -= factor * matrix[column][entry];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    return rhs;
}

inline std::vector<double> normalizedTimes(const std::vector<double>& times) {
    if (times.size() < 2) throw std::invalid_argument("integration requires two time points");
    const double step = times[0] - times[1];
    if (!(step > 0.0) || !std::isfinite(step)) {
        throw std::invalid_argument("integration time points must be newest first and distinct");
    }
    std::vector<double> normalized(times.size(), 0.0);
    for (std::size_t i = 0; i < times.size(); ++i) {
        normalized[i] = (times[i] - times[0]) / step;
        if (!std::isfinite(normalized[i])) throw std::invalid_argument("non-finite integration time");
    }
    return normalized;
}

inline std::vector<double> multiplyPolynomial(
    const std::vector<double>& polynomial,
    double root) {
    std::vector<double> result(polynomial.size() + 1, 0.0);
    for (std::size_t degree = 0; degree < polynomial.size(); ++degree) {
        result[degree] -= root * polynomial[degree];
        result[degree + 1] += polynomial[degree];
    }
    return result;
}

inline double integratePolynomial(
    const std::vector<double>& polynomial,
    double lower,
    double upper) {
    double integral = 0.0;
    for (std::size_t degree = 0; degree < polynomial.size(); ++degree) {
        const double power = static_cast<double>(degree + 1);
        integral += polynomial[degree] / power *
            (std::pow(upper, power) - std::pow(lower, power));
    }
    return integral;
}

} // namespace integration_detail

inline IntegrationFormula makeBdfFormula(const std::vector<double>& newestFirstTimes) {
    const auto normalized = integration_detail::normalizedTimes(newestFirstTimes);
    const std::size_t count = normalized.size();
    std::vector<std::vector<double>> consistency(count, std::vector<double>(count, 0.0));
    std::vector<double> target(count, 0.0);
    if (count > 1) target[1] = 1.0;
    for (std::size_t power = 0; power < count; ++power) {
        for (std::size_t point = 0; point < count; ++point) {
            consistency[power][point] = std::pow(normalized[point], static_cast<int>(power));
        }
    }
    auto weights = integration_detail::solveDense(consistency, target);
    const double step = newestFirstTimes[0] - newestFirstTimes[1];
    for (double& weight : weights) weight /= step;
    return {IntegrationFamily::Bdf, static_cast<int>(count - 1), std::move(weights), {}};
}

inline IntegrationFormula makeAdamsMoultonFormula(
    const std::vector<double>& newestFirstDerivativeTimes) {
    const auto normalized = integration_detail::normalizedTimes(newestFirstDerivativeTimes);
    const std::size_t count = normalized.size();
    const double step = newestFirstDerivativeTimes[0] - newestFirstDerivativeTimes[1];
    std::vector<double> integrationWeights(count, 0.0);
    for (std::size_t basis = 0; basis < count; ++basis) {
        std::vector<double> polynomial{1.0};
        double denominator = 1.0;
        for (std::size_t other = 0; other < count; ++other) {
            if (other == basis) continue;
            polynomial = integration_detail::multiplyPolynomial(polynomial, normalized[other]);
            denominator *= normalized[basis] - normalized[other];
        }
        for (double& coefficient : polynomial) coefficient /= denominator;
        integrationWeights[basis] = step *
            integration_detail::integratePolynomial(polynomial, normalized[1], normalized[0]);
    }
    if (std::abs(integrationWeights[0]) < 1e-30) {
        throw std::invalid_argument("invalid implicit Adams coefficient");
    }
    IntegrationFormula formula;
    formula.family = IntegrationFamily::AdamsMoulton;
    formula.order = static_cast<int>(count);
    const double inverseLeading = 1.0 / integrationWeights[0];
    formula.qWeights = {inverseLeading, -inverseLeading};
    for (std::size_t i = 1; i < integrationWeights.size(); ++i) {
        formula.derivativeWeights.push_back(-integrationWeights[i] * inverseLeading);
    }
    return formula;
}

} // namespace gspice

#endif // GSPICE_INTEGRATION_FORMULA_HPP
