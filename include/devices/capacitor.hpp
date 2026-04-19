#ifndef GSPICE_CAPACITOR_HPP
#define GSPICE_CAPACITOR_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class Capacitor : public Device {
public:
    Capacitor(const std::string& name, int nodePos, int nodeNeg, double value)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        if (timeStep <= 0) return; 
        double Geq = value_ / timeStep;
        if (x_hist.empty()) return;
        const VectorReal& x_prev = x_hist.back();
        double Vd_prev = ((nodePos_ >= 0) ? x_prev[nodePos_] : 0.0) - ((nodeNeg_ >= 0) ? x_prev[nodeNeg_] : 0.0);
        double Ieq = Geq * Vd_prev;
        J.add(nodePos_, nodePos_, Geq);
        J.add(nodeNeg_, nodeNeg_, Geq);
        J.add(nodePos_, nodeNeg_, -Geq);
        J.add(nodeNeg_, nodePos_, -Geq);
        b.add(nodePos_, Ieq);
        b.add(nodeNeg_, -Ieq);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        std::complex<double> Y = {0.0, omega * value_};
        J.add(nodePos_, nodePos_, Y);
        J.add(nodeNeg_, nodeNeg_, Y);
        J.add(nodePos_, nodeNeg_, -Y);
        J.add(nodeNeg_, nodePos_, -Y);
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        int K = 2 * n_harms + 1;
        // Capacitor is open circuit for DC (k=0), so skip
        for (int h = 1; h <= n_harms; ++h) {
            double omega_h = 2.0 * 3.14159265358979 * f_fund * h;
            double val = omega_h * value_;
            
            int cos_idx = 2 * h - 1;
            int sin_idx = 2 * h;

            // Stamping dQ/dt (imaginary part)
            // [ 0   -w ] [ Vcos ] = [ Icos ]
            // [ w    0 ] [ Vsin ] = [ Isin ]
            J.add(nodePos_ * K + cos_idx, nodePos_ * K + sin_idx, -val);
            J.add(nodePos_ * K + sin_idx, nodePos_ * K + cos_idx, val);
            
            J.add(nodeNeg_ * K + cos_idx, nodeNeg_ * K + sin_idx, val);
            J.add(nodeNeg_ * K + sin_idx, nodeNeg_ * K + cos_idx, -val);

            // Mutual terms
            J.add(nodePos_ * K + cos_idx, nodeNeg_ * K + sin_idx, val);
            J.add(nodePos_ * K + sin_idx, nodeNeg_ * K + cos_idx, -val);
            J.add(nodeNeg_ * K + cos_idx, nodePos_ * K + sin_idx, -val);
            J.add(nodeNeg_ * K + sin_idx, nodePos_ * K + cos_idx, val);
        }
    }

private:
    int nodePos_;
    int nodeNeg_;
    double value_;
};

} // namespace gspice

#endif // GSPICE_CAPACITOR_HPP
