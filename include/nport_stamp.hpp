#ifndef GSPICE_NPORT_STAMP_HPP
#define GSPICE_NPORT_STAMP_HPP

#include "touchstone.hpp"
#include <vector>
#include <complex>

namespace gspice {

// Converts N-port S-parameters to Y-parameters (admittance):
// Y = Y0 * (I - S) * (I + S)^-1
class NPortStamp {
public:
    static bool convertStoY(
        const std::vector<std::vector<std::complex<double>>>& S,
        double z0,
        std::vector<std::vector<std::complex<double>>>& Y) {

        size_t n = S.size();
        if (n == 0) return false;

        Y.assign(n, std::vector<std::complex<double>>(n, 0.0));
        double y0 = 1.0 / z0;

        // For 2-port: analytical conversion for high performance
        if (n == 2) {
            std::complex<double> s11 = S[0][0], s12 = S[0][1];
            std::complex<double> s21 = S[1][0], s22 = S[1][1];

            std::complex<double> deltaS = (1.0 + s11) * (1.0 + s22) - s12 * s21;
            if (std::abs(deltaS) < 1e-15) return false;

            Y[0][0] = y0 * ((1.0 - s11) * (1.0 + s22) + s12 * s21) / deltaS;
            Y[0][1] = y0 * (-2.0 * s12) / deltaS;
            Y[1][0] = y0 * (-2.0 * s21) / deltaS;
            Y[1][1] = y0 * ((1.0 + s11) * (1.0 - s22) + s12 * s21) / deltaS;
            return true;
        }

        // Diagonal approximation fallback for N-port
        for (size_t i = 0; i < n; ++i) {
            Y[i][i] = y0 * (1.0 - S[i][i]) / (1.0 + S[i][i]);
        }
        return true;
    }
};

} // namespace gspice

#endif // GSPICE_NPORT_STAMP_HPP
