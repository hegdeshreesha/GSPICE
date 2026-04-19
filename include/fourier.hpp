#ifndef GSPICE_FOURIER_HPP
#define GSPICE_FOURIER_HPP

#include <vector>
#include <complex>
#include <cmath>

namespace gspice {

class Fourier {
public:
    /**
     * Discrete Fourier Transform (DFT)
     * Converts N time-domain samples to N/2 harmonics.
     */
    static std::vector<std::complex<double>> forward(const std::vector<double>& time_samples) {
        int N = static_cast<int>(time_samples.size());
        int K = N / 2;
        std::vector<std::complex<double>> harmonics(K + 1, {0.0, 0.0});

        for (int k = 0; k <= K; ++k) {
            for (int n = 0; n < N; ++n) {
                double angle = 2.0 * 3.14159265358979 * k * n / N;
                harmonics[k] += std::complex<double>(time_samples[n] * std::cos(angle), 
                                                    -time_samples[n] * std::sin(angle));
            }
            harmonics[k] /= (double)N;
            if (k > 0) harmonics[k] *= 2.0; // Account for negative frequencies
        }
        return harmonics;
    }

    /**
     * Inverse Discrete Fourier Transform (IDFT)
     * Converts K harmonics back to N time-domain samples.
     */
    static std::vector<double> inverse(const std::vector<std::complex<double>>& harmonics, int N) {
        std::vector<double> time_samples(N, 0.0);
        int K = static_cast<int>(harmonics.size()) - 1;

        for (int n = 0; n < N; ++n) {
            // DC Component
            time_samples[n] = harmonics[0].real();
            // Harmonics
            for (int k = 1; k <= K; ++k) {
                double angle = 2.0 * 3.14159265358979 * k * n / N;
                time_samples[n] += harmonics[k].real() * std::cos(angle) - 
                                   harmonics[k].imag() * std::sin(angle);
            }
        }
        return time_samples;
    }
};

} // namespace gspice

#endif // GSPICE_FOURIER_HPP
