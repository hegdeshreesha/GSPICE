#ifndef GSPICE_PROBE_HPP
#define GSPICE_PROBE_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class StabilityProbe : public Device {
public:
    /**
     * Tian Probe for STB Analysis
     * @param name Device name
     * @param nodePos Input of loop
     * @param nodeNeg Output of loop
     */
    StabilityProbe(const std::string& name, int nodePos, int nodeNeg, int branchIndex = -1)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), branchIndex_(branchIndex) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        if (branchIndex_ < 0) return;
        // In DC, probe is a short circuit (0V source)
        J.add(nodePos_, branchIndex_, 1.0);
        J.add(nodeNeg_, branchIndex_, -1.0);
        J.add(branchIndex_, nodePos_, 1.0);
        J.add(branchIndex_, nodeNeg_, -1.0);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        if (branchIndex_ < 0) return;
        // In standard AC, probe is a short circuit
        J.add(nodePos_, branchIndex_, {1.0, 0.0});
        J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
        J.add(branchIndex_, nodePos_, {1.0, 0.0});
        J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
    }

    /**
     * Custom stamp for the two STB passes
     * @param pass 1 for Voltage excitation, 2 for Current excitation
     */
    void stbStamp(SparseMatrixComplex& J, VectorComplex& b, int pass) {
        if (branchIndex_ < 0) return;
        if (pass == 1) {
            // Voltage Pass: Probe is a 1V AC source
            J.add(nodePos_, branchIndex_, {1.0, 0.0});
            J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
            J.add(branchIndex_, nodePos_, {1.0, 0.0});
            J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
            b.add(branchIndex_, {1.0, 0.0});
        } else {
            // Current Pass: Probe is a 1A AC current source in parallel with the short
            // (Tian math actually uses a more complex nullor setup, 
            // but for GSPICE we will implement the source injection here).
            J.add(nodePos_, branchIndex_, {1.0, 0.0});
            J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
            J.add(branchIndex_, nodePos_, {1.0, 0.0});
            J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
            b.add(nodePos_, {-1.0, 0.0});
            b.add(nodeNeg_, {1.0, 0.0});
        }
    }

    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }
    int getNodePos() const { return nodePos_; }
    int getNodeNeg() const { return nodeNeg_; }

private:
    int nodePos_, nodeNeg_;
    int branchIndex_;
};

} // namespace gspice

#endif // GSPICE_PROBE_HPP
