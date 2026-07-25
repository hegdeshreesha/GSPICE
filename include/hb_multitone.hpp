#ifndef GSPICE_HB_MULTITONE_HPP
#define GSPICE_HB_MULTITONE_HPP

#include <vector>
#include <cmath>
#include <utility>

namespace gspice {

struct FrequencyMix {
    int k1 = 0;
    int k2 = 0;
    double freq = 0.0;
};

// Generates 2D frequency lattice for multi-tone Harmonic Balance: f = k1*f1 + k2*f2
class HbMultiToneLattice {
public:
    static std::vector<FrequencyMix> generateLattice(
        double f1, int maxK1,
        double f2, int maxK2) {

        std::vector<FrequencyMix> lattice;
        lattice.push_back({0, 0, 0.0}); // DC

        for (int k1 = -maxK1; k1 <= maxK1; ++k1) {
            for (int k2 = -maxK2; k2 <= maxK2; ++k2) {
                if (k1 == 0 && k2 == 0) continue;
                double freq = k1 * f1 + k2 * f2;
                if (freq > 0.0) { // Positive frequencies only
                    lattice.push_back({k1, k2, freq});
                }
            }
        }
        return lattice;
    }
};

} // namespace gspice

#endif // GSPICE_HB_MULTITONE_HPP
