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
        for (int h = 1; h <= n_harms; ++h) {
            double omega_h = 2.0 * 3.14159265358979 * f_fund * h;
            double val = omega_h * value_;
            int cos_idx = 2 * h - 1; int sin_idx = 2 * h;
            J.add(nodePos_ * K + cos_idx, nodePos_ * K + sin_idx, -val);
            J.add(nodePos_ * K + sin_idx, nodePos_ * K + cos_idx, val);
            J.add(nodeNeg_ * K + cos_idx, nodeNeg_ * K + sin_idx, val);
            J.add(nodeNeg_ * K + sin_idx, nodeNeg_ * K + cos_idx, -val);
            J.add(nodePos_ * K + cos_idx, nodeNeg_ * K + sin_idx, val);
            J.add(nodePos_ * K + sin_idx, nodeNeg_ * K + cos_idx, -val);
            J.add(nodeNeg_ * K + cos_idx, nodePos_ * K + sin_idx, -val);
            J.add(nodeNeg_ * K + sin_idx, nodePos_ * K + cos_idx, val);
        }
    }

    void pacStamp(SparseMatrixReal& J, VectorReal& b, double f_in, double f_fund, int n_harms, const VectorReal& x_periodic) override {
        int K = 2 * n_harms + 1;
        // Sideband frequencies: f_k = f_in + k*f_fund
        for (int k = -n_harms; k <= n_harms; ++k) {
            double fk = f_in + k * f_fund;
            double omega_k = 2.0 * 3.14159265358979 * fk;
            double val = omega_k * value_;
            
            // Map k to our 2*N+1 indices (0=DC, 1=H1_Cos, 2=H1_Sin, ...)
            int idx;
            if (k == 0) idx = 0;
            else if (k > 0) idx = 2 * k - 1;
            else idx = 2 * std::abs(k);

            // Stamp complex admittance val into the conversion matrix
            // This is a simplified version (diagonal sidebands only)
            J.add(nodePos_ * K + idx, nodePos_ * K + idx, val); // Imaginary parts would need cross-terms
        }
    }

private:
    int nodePos_;
    int nodeNeg_;
    double value_;
};

} // namespace gspice

#endif // GSPICE_CAPACITOR_HPP
