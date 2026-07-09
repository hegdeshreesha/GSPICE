#ifndef GSPICE_INDUCTOR_HPP
#define GSPICE_INDUCTOR_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class Inductor : public Device {
public:
    Inductor(const std::string& name, int nodePos, int nodeNeg, double value, int branchIndex = -1)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value), branchIndex_(branchIndex) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        if (branchIndex_ < 0) return;
        if (timeStep <= 0) {
            J.add(nodePos_, branchIndex_, 1.0);
            J.add(nodeNeg_, branchIndex_, -1.0);
            J.add(branchIndex_, nodePos_, 1.0);
            J.add(branchIndex_, nodeNeg_, -1.0);
            return;
        }
        double Req = value_ / timeStep;
        const VectorReal& x_prev = x_hist.back();
        double Ibr_prev = x_prev[branchIndex_];
        J.add(nodePos_, branchIndex_, 1.0);
        J.add(nodeNeg_, branchIndex_, -1.0);
        J.add(branchIndex_, nodePos_, 1.0);
        J.add(branchIndex_, nodeNeg_, -1.0);
        J.add(branchIndex_, branchIndex_, -Req);
        b.add(branchIndex_, -Req * Ibr_prev);
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

        const double iPrev = branchCurrent((*ctx.xHistory)[ctx.xHistory->size() - 1]);
        const double iPrev2 = useSecond ? branchCurrent((*ctx.xHistory)[ctx.xHistory->size() - 2]) : iPrev;
        const double Req = value_ * a0;
        double rhs = value_ * (a1 * iPrev + (useSecond ? a2 * iPrev2 : 0.0));
        if (useTrap) {
            rhs -= prevVoltage_;
        }

        J.add(nodePos_, branchIndex_, 1.0);
        J.add(nodeNeg_, branchIndex_, -1.0);
        J.add(branchIndex_, nodePos_, 1.0);
        J.add(branchIndex_, nodeNeg_, -1.0);
        J.add(branchIndex_, branchIndex_, -Req);
        b.add(branchIndex_, rhs);
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        (void)currentTime;
        (void)ctx;
        prevVoltage_ = nodeVoltage(x, nodePos_) - nodeVoltage(x, nodeNeg_);
        prevVoltageValid_ = true;
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        if (branchIndex_ < 0) return;
        std::complex<double> Z = {0.0, omega * value_};
        J.add(nodePos_, branchIndex_, {1.0, 0.0});
        J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
        J.add(branchIndex_, nodePos_, {1.0, 0.0});
        J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
        J.add(branchIndex_, branchIndex_, -Z);
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
