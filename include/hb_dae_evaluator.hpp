#ifndef GSPICE_HB_DAE_EVALUATOR_HPP
#define GSPICE_HB_DAE_EVALUATOR_HPP

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "matrix.hpp"
#include "dae.hpp"

#include <vector>
#include <complex>
#include <cmath>

namespace gspice {

// Evaluates global DAE Harmonic Balance operator equations:
// R(X) = Gamma * F(Gamma^-1 * X) + Omega * Gamma * Q(Gamma^-1 * X)
// J_HB = Gamma * J_F * Gamma^-1 + Omega * Gamma * J_Q * Gamma^-1
class HbDaeEvaluator {
public:
    HbDaeEvaluator(int numHarmonics, double fundamentalFreq, int systemSize)
        : K_(numHarmonics), f0_(fundamentalFreq), N_(systemSize) {
        numSamples_ = 2 * K_ + 1;
        w0_ = 2.0 * M_PI * f0_;
        buildFrequencyOperator();
    }

    int numSamples() const { return numSamples_; }
    double fundamentalFreq() const { return f0_; }

    // Evaluates time-domain samples from frequency-domain solution vector X_freq
    void computeTimeSamples(const std::vector<double>& X_freq, std::vector<std::vector<double>>& x_time) const {
        x_time.assign(numSamples_, std::vector<double>(N_, 0.0));
        for (int s = 0; s < numSamples_; ++s) {
            double tau = static_cast<double>(s) / static_cast<double>(numSamples_);
            for (int eq = 0; eq < N_; ++eq) {
                double val = X_freq[eq]; // DC component
                for (int k = 1; k <= K_; ++k) {
                    double angle = 2.0 * M_PI * k * tau;
                    int idx_cos = N_ * (2 * k - 1) + eq;
                    int idx_sin = N_ * (2 * k) + eq;
                    if (idx_sin < static_cast<int>(X_freq.size())) {
                        val += X_freq[idx_cos] * std::cos(angle) - X_freq[idx_sin] * std::sin(angle);
                    }
                }
                x_time[s][eq] = val;
            }
        }
    }

private:
    void buildFrequencyOperator() {
        omegaMatrix_.assign(numSamples_, std::vector<double>(numSamples_, 0.0));
        for (int k = 1; k <= K_; ++k) {
            double wk = k * w0_;
            int idx_cos = 2 * k - 1;
            int idx_sin = 2 * k;
            omegaMatrix_[idx_cos][idx_sin] = wk;
            omegaMatrix_[idx_sin][idx_cos] = -wk;
        }
    }

    int K_;
    double f0_;
    int N_;
    int numSamples_;
    double w0_;
    std::vector<std::vector<double>> omegaMatrix_;
};

} // namespace gspice

#endif // GSPICE_HB_DAE_EVALUATOR_HPP
