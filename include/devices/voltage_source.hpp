#ifndef GSPICE_VOLTAGE_SOURCE_HPP
#define GSPICE_VOLTAGE_SOURCE_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class VoltageSource : public Device {
public:
    VoltageSource(const std::string& name, int nodePos, int nodeNeg, double voltage, int branchIndex = -1)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), voltage_(voltage), branchIndex_(branchIndex) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        if (branchIndex_ < 0) return;
        J.add(nodePos_, branchIndex_, 1.0);
        J.add(nodeNeg_, branchIndex_, -1.0);
        J.add(branchIndex_, nodePos_, 1.0);
        J.add(branchIndex_, nodeNeg_, -1.0);
        b.add(branchIndex_, voltage_);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        if (branchIndex_ < 0) return;
        J.add(nodePos_, branchIndex_, {1.0, 0.0});
        J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
        J.add(branchIndex_, nodePos_, {1.0, 0.0});
        J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
        b.add(branchIndex_, {1.0, 0.0});
    }

    void pacStamp(SparseMatrixReal& J, VectorReal& b, double f_in, double f_fund, int n_harms, const VectorReal& x_periodic) override {
        hbStamp(J, b, f_fund, n_harms, x_periodic); 
        int K = 2 * n_harms + 1;
        b.add(branchIndex_ * K, 1.0);
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        if (branchIndex_ < 0) return;
        int K = 2 * n_harms + 1;
        for (int k = 0; k < K; ++k) {
            J.add(nodePos_ * K + k, branchIndex_ * K + k, 1.0);
            J.add(nodeNeg_ * K + k, branchIndex_ * K + k, -1.0);
            J.add(branchIndex_ * K + k, nodePos_ * K + k, 1.0);
            J.add(branchIndex_ * K + k, nodeNeg_ * K + k, -1.0);
        }
        b.add(branchIndex_ * K, voltage_);
    }

    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }

private:
    int nodePos_, nodeNeg_;
    double voltage_;
    int branchIndex_;
};

} // namespace gspice

#endif // GSPICE_VOLTAGE_SOURCE_HPP
