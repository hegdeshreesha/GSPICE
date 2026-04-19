#ifndef GSPICE_DIODE_HPP
#define GSPICE_DIODE_HPP

#include "device.hpp"
#include "fourier.hpp"
#include <cmath>
#include <string>

namespace gspice {

class Diode : public Device {
public:
    Diode(const std::string& name, int nodePos, int nodeNeg, double Is = 1e-14, double Vt = 0.026)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), Is_(Is), Vt_(Vt) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        double Vd = ((nodePos_ >= 0) ? x[nodePos_] : 0.0) - ((nodeNeg_ >= 0) ? x[nodeNeg_] : 0.0);
        if (Vd > 0.8) Vd = 0.8; if (Vd < -2.0) Vd = -2.0;
        double expV = std::exp(Vd / Vt_);
        double Id = Is_ * (expV - 1.0);
        double gd = (Is_ / Vt_) * expV;
        double gmin = 1e-12; gd += gmin; Id += gmin * Vd;
        double Ieq = Id - gd * Vd;
        J.add(nodePos_, nodePos_, gd);
        J.add(nodeNeg_, nodeNeg_, gd);
        J.add(nodePos_, nodeNeg_, -gd);
        J.add(nodeNeg_, nodePos_, -gd);
        b.add(nodePos_, -Ieq);
        b.add(nodeNeg_, Ieq);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        double Vd = ((nodePos_ >= 0) ? x_dc[nodePos_] : 0.0) - ((nodeNeg_ >= 0) ? x_dc[nodeNeg_] : 0.0);
        double gd = (Is_ / Vt_) * std::exp(Vd / Vt_);
        J.add(nodePos_, nodePos_, {gd, 0.0});
        J.add(nodeNeg_, nodeNeg_, {gd, 0.0});
        J.add(nodePos_, nodeNeg_, {-gd, 0.0});
        J.add(nodeNeg_, nodePos_, {-gd, 0.0});
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        int K = 2 * n_harms + 1;
        int N_samples = 4 * n_harms; // Over-sampling for FFT accuracy

        // 1. Extract frequency domain node voltages
        std::vector<std::complex<double>> Vd_freq(n_harms + 1);
        Vd_freq[0] = { ((nodePos_ >= 0) ? x_hb[nodePos_ * K] : 0.0) - ((nodeNeg_ >= 0) ? x_hb[nodeNeg_ * K] : 0.0), 0.0 };
        for (int h = 1; h <= n_harms; ++h) {
            double v_cos = ((nodePos_ >= 0) ? x_hb[nodePos_ * K + 2*h - 1] : 0.0) - ((nodeNeg_ >= 0) ? x_hb[nodeNeg_ * K + 2*h - 1] : 0.0);
            double v_sin = ((nodePos_ >= 0) ? x_hb[nodePos_ * K + 2*h] : 0.0) - ((nodeNeg_ >= 0) ? x_hb[nodeNeg_ * K + 2*h] : 0.0);
            Vd_freq[h] = {v_cos, v_sin};
        }

        // 2. IDFT to get time-domain voltages
        std::vector<double> Vd_time = Fourier::inverse(Vd_freq, N_samples);

        // 3. Evaluate non-linear Id and gd at each time sample
        std::vector<double> Id_time(N_samples);
        std::vector<double> gd_time(N_samples);
        for (int i = 0; i < N_samples; ++i) {
            double vd = Vd_time[i];
            if (vd > 0.8) vd = 0.8;
            double expV = std::exp(vd / Vt_);
            Id_time[i] = Is_ * (expV - 1.0);
            gd_time[i] = (Is_ / Vt_) * expV;
        }

        // 4. DFT back to frequency domain
        auto Id_freq = Fourier::forward(Id_time);
        auto gd_freq = Fourier::forward(gd_time);

        // 5. Stamp the HB matrix and RHS
        // DC Component (k=0)
        b.add(nodePos_ * K, -Id_freq[0].real());
        b.add(nodeNeg_ * K, Id_freq[0].real());
        
        // This is a simplified convolution for the Jacobian
        // Real simulators use more complex frequency-shifting here.
        for (int h = 0; h < K; ++h) {
            J.add(nodePos_ * K + h, nodePos_ * K + h, gd_freq[0].real());
        }
    }

private:
    int nodePos_;
    int nodeNeg_;
    double Is_;
    double Vt_;
};

} // namespace gspice

#endif // GSPICE_DIODE_HPP
