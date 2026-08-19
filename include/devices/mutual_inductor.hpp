#ifndef GSPICE_MUTUAL_INDUCTOR_HPP
#define GSPICE_MUTUAL_INDUCTOR_HPP

#include "device.hpp"
#include <cmath>
#include <string>

namespace gspice {

class MutualInductor : public Device {
public:
    MutualInductor(
        const std::string& name,
        std::string primaryName,
        std::string secondaryName,
        double coupling)
        : Device(name),
          primaryName_(std::move(primaryName)),
          secondaryName_(std::move(secondaryName)),
          coupling_(coupling) {}

    const std::string& primaryName() const { return primaryName_; }
    const std::string& secondaryName() const { return secondaryName_; }
    double coupling() const { return coupling_; }

    void setResolvedBranches(int primaryBranch, int secondaryBranch, double mutualInductance) {
        primaryBranch_ = primaryBranch;
        secondaryBranch_ = secondaryBranch;
        mutualInductance_ = mutualInductance;
    }

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        if (primaryBranch_ < 0 || secondaryBranch_ < 0 || mutualInductance_ == 0.0) return true;
        if (request.dynamicResidual) {
            evaluation.dynamicResidual.push_back(
                {primaryBranch_, -mutualInductance_ * branchCurrent(x, secondaryBranch_)});
            evaluation.dynamicResidual.push_back(
                {secondaryBranch_, -mutualInductance_ * branchCurrent(x, primaryBranch_)});
        }
        if (request.dynamicJacobian) {
            evaluation.dynamicJacobian.push_back({primaryBranch_, secondaryBranch_, -mutualInductance_});
            evaluation.dynamicJacobian.push_back({secondaryBranch_, primaryBranch_, -mutualInductance_});
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
        (void)currentTime;
        if (timeStep <= 0.0 || x_hist.empty()) return;
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
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;
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
        appendScaledDaeResidual(history, previous.dynamicResidual, ctx.a1);
        if (ctx.hasSecondHistory && ctx.xHistory->size() >= 2) {
            evaluateDae((*ctx.xHistory)[ctx.xHistory->size() - 2], request, previous2);
            appendScaledDaeResidual(history, previous2.dynamicResidual, ctx.a2);
        }
        stampDaeTransient(current, x, ctx.a0, history, J, b);
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

private:
    static double branchCurrent(const VectorReal& x, int branch) {
        return branch >= 0 && branch < x.getSize() ? x[branch] : 0.0;
    }

    std::string primaryName_;
    std::string secondaryName_;
    double coupling_ = 0.0;
    int primaryBranch_ = -1;
    int secondaryBranch_ = -1;
    double mutualInductance_ = 0.0;
};

} // namespace gspice

#endif // GSPICE_MUTUAL_INDUCTOR_HPP
