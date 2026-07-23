#ifndef GSPICE_INDUCTOR_HPP
#define GSPICE_INDUCTOR_HPP

#include "device.hpp"
#include <string>
#include <cstring>
#include <cstdint>

namespace gspice {

class Inductor : public Device {
public:
    Inductor(const std::string& name, int nodePos, int nodeNeg, double value, int branchIndex = -1)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value), branchIndex_(branchIndex) {}

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        if (branchIndex_ < 0) return true;
        const double current = branchCurrent(x);
        const double voltage = nodeVoltage(x, nodePos_) - nodeVoltage(x, nodeNeg_);
        if (request.staticResidual) {
            evaluation.staticResidual.push_back({nodePos_, current});
            evaluation.staticResidual.push_back({nodeNeg_, -current});
            evaluation.staticResidual.push_back({branchIndex_, voltage});
        }
        if (request.staticJacobian) {
            evaluation.staticJacobian.push_back({nodePos_, branchIndex_, 1.0});
            evaluation.staticJacobian.push_back({nodeNeg_, branchIndex_, -1.0});
            evaluation.staticJacobian.push_back({branchIndex_, nodePos_, 1.0});
            evaluation.staticJacobian.push_back({branchIndex_, nodeNeg_, -1.0});
        }
        // The branch equation is Vp - Vn - d(L*I)/dt = 0, therefore its
        // conserved DAE quantity is -flux = -L*I.
        if (request.dynamicResidual) {
            evaluation.dynamicResidual.push_back({branchIndex_, -value_ * current});
        }
        if (request.dynamicJacobian) {
            evaluation.dynamicJacobian.push_back({branchIndex_, branchIndex_, -value_});
        }
        return true;
    }

    bool daeAuditSafe() const override { return true; }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        if (branchIndex_ < 0) return;
        if (timeStep <= 0) {
            DaeRequest request;
            DaeEvaluation evaluation;
            evaluateDae(x, request, evaluation);
            stampDaeStatic(evaluation, x, J, b);
            return;
        }
        if (x_hist.empty()) return;
        DaeRequest request;
        request.analysis = DaeAnalysis::Transient;
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
        if (branchIndex_ < 0) return;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) {
            dcStamp(J, b, x, 0.0, ctx.currentTime, {});
            return;
        }

        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && prevVoltageValid_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }

        DaeRequest request;
        request.analysis = DaeAnalysis::Transient;
        request.time = ctx.currentTime;
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
            history.push_back({branchIndex_, prevVoltage_});
        }
        stampDaeTransient(current, x, a0, history, J, b);
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        (void)currentTime;
        (void)ctx;
        prevVoltage_ = nodeVoltage(x, nodePos_) - nodeVoltage(x, nodeNeg_);
        prevVoltageValid_ = true;
    }

    std::size_t transientStateBytes() const override {
        return sizeof(prevVoltage_) + sizeof(std::uint8_t);
    }

    void saveTransientStateBytes(std::byte* destination, std::size_t size) const override {
        if (size != transientStateBytes()) throw std::invalid_argument("inductor transient state size");
        std::memcpy(destination, &prevVoltage_, sizeof(prevVoltage_));
        const std::uint8_t valid = prevVoltageValid_ ? 1u : 0u;
        std::memcpy(destination + sizeof(prevVoltage_), &valid, sizeof(valid));
    }

    void restoreTransientStateBytes(const std::byte* source, std::size_t size) override {
        if (size != transientStateBytes()) throw std::invalid_argument("inductor transient state size");
        std::memcpy(&prevVoltage_, source, sizeof(prevVoltage_));
        std::uint8_t valid = 0;
        std::memcpy(&valid, source + sizeof(prevVoltage_), sizeof(valid));
        prevVoltageValid_ = valid != 0;
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        if (branchIndex_ < 0) return;
        (void)b;
        DaeRequest request;
        request.analysis = DaeAnalysis::SmallSignal;
        request.staticResidual = false;
        request.dynamicResidual = false;
        request.dynamicJacobian = true;
        DaeEvaluation evaluation;
        evaluateDae(x_dc, request, evaluation);
        stampDaeSmallSignal(evaluation, omega, J);
    }

    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }

private:
    static double nodeVoltage(const VectorReal& x, int node) {
        return node >= 0 ? x[node] : 0.0;
    }

    double branchCurrent(const VectorReal& x) const {
        return branchIndex_ >= 0 ? x[branchIndex_] : 0.0;
    }

    int nodePos_;
    int nodeNeg_;
    double value_;
    int branchIndex_;
    double prevVoltage_ = 0.0;
    bool prevVoltageValid_ = false;
};

} // namespace gspice

#endif // GSPICE_INDUCTOR_HPP
