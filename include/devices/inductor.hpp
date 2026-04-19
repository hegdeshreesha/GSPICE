#ifndef GSPICE_INDUCTOR_HPP
#define GSPICE_INDUCTOR_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class Inductor : public Device {
public:
    Inductor(const std::string& name, int nodePos, int nodeNeg, double value, int branchIndex = -1)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value), branchIndex_(branchIndex) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
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
    int nodePos_;
    int nodeNeg_;
    double value_;
    int branchIndex_;
};

} // namespace gspice

#endif // GSPICE_INDUCTOR_HPP
