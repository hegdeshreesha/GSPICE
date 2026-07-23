#ifndef GSPICE_DAE_HPP
#define GSPICE_DAE_HPP

#include "matrix.hpp"
#include "sparse_matrix.hpp"

#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace gspice {

// GSPICE's device-neutral modified-nodal DAE contract:
//
//     F(x, t) + d Q(x, t) / dt = 0
//
// F contains instantaneous flow terms and Q contains conserved quantities
// such as terminal charge and branch flux.  Analyses decide how Q is
// differentiated; devices never need to manufacture an analysis-specific
// companion model in order to expose their physical storage model.
enum class DaeAnalysis {
    OperatingPoint,
    Transient,
    SmallSignal,
    HarmonicBalance
};

struct DaeRequest {
    DaeAnalysis analysis = DaeAnalysis::OperatingPoint;
    double time = 0.0;
    bool staticResidual = true;
    bool dynamicResidual = false;
    bool staticJacobian = true;
    bool dynamicJacobian = false;
    bool highPrecision = false;
    bool readOnlyState = false;
    bool enableLimiting = true;
    bool nodeset = false;
    bool allowBypass = false;
    double bypassRelativeTolerance = 0.0;
    double bypassAbsoluteTolerance = 0.0;
    std::uint64_t evaluationEpoch = 0;
};

struct DaeResidualTerm {
    int equation = -1;
    double value = 0.0;
    int conservationGroup = -1;
};

struct DaeJacobianTerm {
    int equation = -1;
    int unknown = -1;
    double value = 0.0;
    int conservationGroup = -1;
};

struct DaeEvaluation {
    std::vector<DaeResidualTerm> staticResidual;
    std::vector<DaeResidualTerm> dynamicResidual;
    std::vector<DaeJacobianTerm> staticJacobian;
    std::vector<DaeJacobianTerm> dynamicJacobian;
    bool limitingApplied = false;
    bool bypassed = false;
    double maximumTimeStep = 0.0;

    void clear() {
        staticResidual.clear();
        dynamicResidual.clear();
        staticJacobian.clear();
        dynamicJacobian.clear();
        limitingApplied = false;
        bypassed = false;
        maximumTimeStep = 0.0;
    }

    bool finite() const {
        const auto residualsFinite = [](const std::vector<DaeResidualTerm>& terms) {
            for (const auto& term : terms) {
                if (!std::isfinite(term.value)) return false;
            }
            return true;
        };
        const auto jacobianFinite = [](const std::vector<DaeJacobianTerm>& terms) {
            for (const auto& term : terms) {
                if (!std::isfinite(term.value)) return false;
            }
            return true;
        };
        return residualsFinite(staticResidual) && residualsFinite(dynamicResidual) &&
            jacobianFinite(staticJacobian) && jacobianFinite(dynamicJacobian);
    }
};

using DaeHistory = std::vector<DaeResidualTerm>;

inline void appendScaledDaeResidual(
    DaeHistory& destination,
    const std::vector<DaeResidualTerm>& source,
    double scale) {
    if (scale == 0.0) return;
    destination.reserve(destination.size() + source.size());
    for (const auto& term : source) {
        destination.push_back({term.equation, scale * term.value, term.conservationGroup});
    }
}

inline double daeUnknownValue(const VectorReal& x, int unknown) {
    return unknown >= 0 && unknown < x.getSize() ? x[unknown] : 0.0;
}

inline void stampDaeStatic(
    const DaeEvaluation& evaluation,
    const VectorReal& x,
    SparseMatrixReal& jacobian,
    VectorReal& rhs) {
    if (!evaluation.finite()) {
        throw std::runtime_error("non-finite value in DAE device evaluation");
    }
    for (const auto& term : evaluation.staticResidual) {
        rhs.add(term.equation, -term.value);
    }
    for (const auto& term : evaluation.staticJacobian) {
        jacobian.add(term.equation, term.unknown, term.value);
        rhs.add(term.equation, term.value * daeUnknownValue(x, term.unknown));
    }
}

// history is the already weighted, known part of dQ/dt.  For a BDF formula,
// for example, dQ/dt = alpha * Q(x_now) + history.
inline void stampDaeTransient(
    const DaeEvaluation& evaluation,
    const VectorReal& x,
    double alpha,
    const DaeHistory& history,
    SparseMatrixReal& jacobian,
    VectorReal& rhs) {
    if (!evaluation.finite() || !std::isfinite(alpha)) {
        throw std::runtime_error("non-finite value in transient DAE assembly");
    }
    for (const auto& term : evaluation.staticResidual) {
        rhs.add(term.equation, -term.value);
    }
    for (const auto& term : evaluation.dynamicResidual) {
        rhs.add(term.equation, -alpha * term.value);
    }
    for (const auto& term : history) {
        if (!std::isfinite(term.value)) {
            throw std::runtime_error("non-finite DAE history contribution");
        }
        rhs.add(term.equation, -term.value);
    }
    for (const auto& term : evaluation.staticJacobian) {
        jacobian.add(term.equation, term.unknown, term.value);
        rhs.add(term.equation, term.value * daeUnknownValue(x, term.unknown));
    }
    for (const auto& term : evaluation.dynamicJacobian) {
        const double value = alpha * term.value;
        jacobian.add(term.equation, term.unknown, value);
        rhs.add(term.equation, value * daeUnknownValue(x, term.unknown));
    }
}

inline void stampDaeSmallSignal(
    const DaeEvaluation& evaluation,
    double omega,
    SparseMatrixComplex& jacobian) {
    if (!evaluation.finite() || !std::isfinite(omega)) {
        throw std::runtime_error("non-finite value in small-signal DAE assembly");
    }
    for (const auto& term : evaluation.staticJacobian) {
        jacobian.add(term.equation, term.unknown, {term.value, 0.0});
    }
    for (const auto& term : evaluation.dynamicJacobian) {
        jacobian.add(term.equation, term.unknown, {0.0, omega * term.value});
    }
}

inline double daeResidualAt(
    const std::vector<DaeResidualTerm>& residual,
    int equation) {
    double value = 0.0;
    for (const auto& term : residual) {
        if (term.equation == equation) value += term.value;
    }
    return value;
}

} // namespace gspice

#endif // GSPICE_DAE_HPP
