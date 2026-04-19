#ifndef GSPICE_RESISTOR_HPP
#define GSPICE_RESISTOR_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class Resistor : public Device {
public:
    Resistor(const std::string& name, int nodePos, int nodeNeg, double value)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value) {}

    double getConductance() const { return 1.0 / value_; }
    int getNodePos() const { return nodePos_; }
    int getNodeNeg() const { return nodeNeg_; }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        double G = getConductance();
        J.add(nodePos_, nodePos_, G);
        J.add(nodeNeg_, nodeNeg_, G);
        J.add(nodePos_, nodeNeg_, -G);
        J.add(nodeNeg_, nodePos_, -G);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        double G = getConductance();
        J.add(nodePos_, nodePos_, {G, 0.0});
        J.add(nodeNeg_, nodeNeg_, {G, 0.0});
        J.add(nodePos_, nodeNeg_, {-G, 0.0});
        J.add(nodeNeg_, nodePos_, {-G, 0.0});
    }

    double getNoisePSD(double omega, const VectorReal& x_dc) override {
        double k = 1.380649e-23;
        double T = 298.15;
        double G = getConductance();
        return 4.0 * k * T * G;
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        double G = getConductance();
        int K = 2 * n_harms + 1; // Size of harmonic block per node

        // Stamp G for DC and every Cos/Sin harmonic component
        for (int k = 0; k < K; ++k) {
            J.add(nodePos_ * K + k, nodePos_ * K + k, G);
            J.add(nodeNeg_ * K + k, nodeNeg_ * K + k, G);
            J.add(nodePos_ * K + k, nodeNeg_ * K + k, -G);
            J.add(nodeNeg_ * K + k, nodePos_ * K + k, -G);
        }
    }
    int nodePos_;
    int nodeNeg_;
    double value_;
};

} // namespace gspice

#endif // GSPICE_RESISTOR_HPP
