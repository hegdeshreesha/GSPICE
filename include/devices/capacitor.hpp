#ifndef GSPICE_CAPACITOR_HPP
#define GSPICE_CAPACITOR_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class Capacitor : public Device {
public:
    Capacitor(const std::string& name, int nodePos, int nodeNeg, double value)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
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

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, const TransientContext& ctx) override {
        (void)x;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;

        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && prevCurrentValid_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }

        const double vPrev = voltage((*ctx.xHistory)[ctx.xHistory->size() - 1]);
        const double vPrev2 = useSecond ? voltage((*ctx.xHistory)[ctx.xHistory->size() - 2]) : vPrev;
        const double Geq = value_ * a0;
        double Ieq = -value_ * (a1 * vPrev + (useSecond ? a2 * vPrev2 : 0.0));
        if (useTrap) {
            Ieq += prevCurrent_;
        }

        J.add(nodePos_, nodePos_, Geq);
        J.add(nodeNeg_, nodeNeg_, Geq);
        J.add(nodePos_, nodeNeg_, -Geq);
        J.add(nodeNeg_, nodePos_, -Geq);
        b.add(nodePos_, Ieq);
        b.add(nodeNeg_, -Ieq);
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        (void)currentTime;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;

        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && prevCurrentValid_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }

        const double vNow = voltage(x);
        const double vPrev = voltage((*ctx.xHistory)[ctx.xHistory->size() - 1]);
        const double vPrev2 = useSecond ? voltage((*ctx.xHistory)[ctx.xHistory->size() - 2]) : vPrev;
        if (useTrap) {
            prevCurrent_ = value_ * (a0 * vNow + a1 * vPrev) - prevCurrent_;
        } else {
            prevCurrent_ = value_ * (a0 * vNow + a1 * vPrev + (useSecond ? a2 * vPrev2 : 0.0));
        }
        prevCurrentValid_ = true;
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
    double voltage(const VectorReal& x) const {
        return ((nodePos_ >= 0) ? x[nodePos_] : 0.0) - ((nodeNeg_ >= 0) ? x[nodeNeg_] : 0.0);
    }

    int nodePos_;
    int nodeNeg_;
    double value_;
    double prevCurrent_ = 0.0;
    bool prevCurrentValid_ = false;
};

} // namespace gspice

#endif // GSPICE_CAPACITOR_HPP
