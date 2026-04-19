#ifndef GSPICE_CURRENT_SOURCE_HPP
#define GSPICE_CURRENT_SOURCE_HPP

#include "device.hpp"
#include <string>
#include <cmath>

namespace gspice {

class CurrentSource : public Device {
public:
    /**
     * @param name Name of source
     * @param nodePos Positive node (current flows FROM this node)
     * @param nodeNeg Negative node (current flows TO this node)
     * @param dcValue DC current value
     */
    CurrentSource(const std::string& name, int nodePos, int nodeNeg, double dcValue)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), dcValue_(dcValue) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        // Current Source RHS Stamp:
        // Ax = b. A current source is a constant 'b'.
        // Flowing out of nodePos: RHS[nodePos] -= dcValue
        // Flowing into nodeNeg:  RHS[nodeNeg] += dcValue
        b.add(nodePos_, -dcValue_);
        b.add(nodeNeg_, dcValue_);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        // AC Magnitude is 1A at phase 0 by default
        b.add(nodePos_, {-1.0, 0.0});
        b.add(nodeNeg_, {1.0, 0.0});
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        int K = 2 * n_harms + 1;
        // RHS: DC component at index 0
        b.add(nodePos_ * K, -dcValue_);
        b.add(nodeNeg_ * K, dcValue_);
    }

private:
    int nodePos_;
    int nodeNeg_;
    double dcValue_;
};

} // namespace gspice

#endif // GSPICE_CURRENT_SOURCE_HPP
