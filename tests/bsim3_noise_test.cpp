#include "devices/bsim3_model.hpp"

#include <cassert>
#include <cmath>
#include <vector>

int main() {
    gspice::Bsim3Mosfet device("M1", 0, 1, 2, 3, 1, 1e-6, 1e-6, 0.4,
                               120e-6, 0.05, 1e5, 0.5, 1.0, 1e-8,
                               1e-20, 1.0, 1.0);
    gspice::VectorReal x(4);
    x[0] = 0.8;
    x[1] = 1.0;

    const double low_frequency_psd = device.getNoisePSD(2.0 * 3.141592653589793 * 1.0e3, x);
    const double psd = device.getNoisePSD(2.0 * 3.141592653589793 * 1.0e9, x);
    assert(low_frequency_psd > psd);
    assert(std::isfinite(psd));
    assert(psd >= 0.0);

    std::vector<gspice::NoiseSource> sources;
    device.collectNoiseSources(0.0, x, sources);
    assert(sources.size() == (psd > 0.0 ? 1u : 0u));
    if (!sources.empty()) {
        assert(sources[0].nodePos == 0);
        assert(sources[0].nodeNeg == 2);
        assert(std::isfinite(sources[0].currentPsd));
        assert(sources[0].currentPsd > 0.0);
    }
}
