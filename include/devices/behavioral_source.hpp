#ifndef GSPICE_BEHAVIORAL_SOURCE_HPP
#define GSPICE_BEHAVIORAL_SOURCE_HPP

#include "device.hpp"
#include "expression.hpp"
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace gspice {

class BehavioralSource : public Device {
public:
    enum class Mode {
        Current,
        Voltage
    };

    BehavioralSource(
        const std::string& name,
        int nodePos,
        int nodeNeg,
        Mode mode,
        BehavioralExpression expression,
        int branchIndex = -1)
        : Device(name),
          nodePos_(nodePos),
          nodeNeg_(nodeNeg),
          mode_(mode),
          expression_(std::move(expression)),
          branchIndex_(branchIndex) {}

    bool isVoltageMode() const { return mode_ == Mode::Voltage; }
    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }

    const std::vector<BehavioralExpression::CurrentRef>& currentRefs() const {
        return expression_.currentRefs();
    }

    void bindBranch(const std::string& name, int branchIndex) {
        expression_.bindBranch(name, branchIndex);
    }

    bool allBranchesBound(std::string& missing) const {
        return expression_.allBranchesBound(missing);
    }

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)timeStep;
        (void)x_hist;
        if (mode_ == Mode::Voltage && branchIndex_ < 0) return;

        const double value = expression_.evaluate(x, currentTime);
        std::vector<std::pair<int, double>> derivs = numericalDerivatives(x, currentTime);

        if (mode_ == Mode::Current) {
            double ieq = value;
            for (const auto& [idx, deriv] : derivs) {
                ieq -= deriv * x[idx];
                J.add(nodePos_, idx, deriv);
                J.add(nodeNeg_, idx, -deriv);
            }
            b.add(nodePos_, -ieq);
            b.add(nodeNeg_, ieq);
        } else {
            J.add(nodePos_, branchIndex_, 1.0);
            J.add(nodeNeg_, branchIndex_, -1.0);
            J.add(branchIndex_, nodePos_, 1.0);
            J.add(branchIndex_, nodeNeg_, -1.0);

            double rhs = value;
            for (const auto& [idx, deriv] : derivs) {
                rhs -= deriv * x[idx];
                J.add(branchIndex_, idx, -deriv);
            }
            b.add(branchIndex_, rhs);
        }
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        if (mode_ == Mode::Voltage && branchIndex_ < 0) return;
        std::vector<std::pair<int, double>> derivs = numericalDerivatives(x_dc, 0.0);

        if (mode_ == Mode::Current) {
            for (const auto& [idx, deriv] : derivs) {
                J.add(nodePos_, idx, {deriv, 0.0});
                J.add(nodeNeg_, idx, {-deriv, 0.0});
            }
        } else {
            J.add(nodePos_, branchIndex_, {1.0, 0.0});
            J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
            J.add(branchIndex_, nodePos_, {1.0, 0.0});
            J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
            for (const auto& [idx, deriv] : derivs) {
                J.add(branchIndex_, idx, {-deriv, 0.0});
            }
        }
    }

private:
    int nodePos_;
    int nodeNeg_;
    Mode mode_;
    BehavioralExpression expression_;
    int branchIndex_;

    std::vector<std::pair<int, double>> numericalDerivatives(const VectorReal& x, double time) const {
        std::vector<std::pair<int, double>> derivs;
        const std::vector<int> deps = expression_.dependencyIndices();
        for (int idx : deps) {
            if (idx < 0 || idx >= x.getSize()) continue;
            const double base = x[idx];
            const double h = std::max(1e-9, std::abs(base) * 1e-6);
            VectorReal xp = x;
            VectorReal xm = x;
            xp[idx] = base + h;
            xm[idx] = base - h;
            const double fp = expression_.evaluate(xp, time);
            const double fm = expression_.evaluate(xm, time);
            const double deriv = (fp - fm) / (2.0 * h);
            if (std::isfinite(deriv) && std::abs(deriv) > 1e-18) {
                derivs.push_back({idx, deriv});
            }
        }
        return derivs;
    }
};

} // namespace gspice

#endif // GSPICE_BEHAVIORAL_SOURCE_HPP
