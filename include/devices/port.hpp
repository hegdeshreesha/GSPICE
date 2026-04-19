#ifndef GSPICE_PORT_HPP
#define GSPICE_PORT_HPP

#include "device.hpp"
#include <string>
#include <complex>

namespace gspice {

class Port : public Device {
public:
    Port(const std::string& name, int nodePos, int nodeNeg, int portNum, double z0 = 50.0)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), portNum_(portNum), z0_(z0) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        double G = 1.0 / z0_;
        J.add(nodePos_, nodePos_, G);
        J.add(nodeNeg_, nodeNeg_, G);
        J.add(nodePos_, nodeNeg_, -G);
        J.add(nodeNeg_, nodePos_, -G);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        if (branchIndex_ < 0) return;
        J.add(nodePos_, branchIndex_, {1.0, 0.0});
        J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
        J.add(branchIndex_, nodePos_, {1.0, 0.0});
        J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
        J.add(branchIndex_, branchIndex_, {-z0_, 0.0});
    }

    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }
    double getZ0() const { return z0_; }
    int getPortNum() const { return portNum_; }
    int getNodePos() const { return nodePos_; }
    int getNodeNeg() const { return nodeNeg_; }

private:
    int nodePos_;
    int nodeNeg_;
    int portNum_;
    double z0_;
    int branchIndex_ = -1;
};

} // namespace gspice

#endif // GSPICE_PORT_HPP
