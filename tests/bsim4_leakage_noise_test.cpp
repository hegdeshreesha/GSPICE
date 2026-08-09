#include "devices/bsim4_model.hpp"

#include <cassert>
#include <cmath>
#include <vector>

int main() {
    const auto model = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"U0", 500.0}, {"IS", 1.0e-14},
        {"JS", 1.0e-2}, {"JSW", 1.0e-3},
        {"NJ", 1.0}, {"AIGC", 1.0e-6}, {"BIGC", 1.0},
        {"CIGC", 0.0}, {"AGIDL", 1.0e-6}, {"BGIDL", 1.0e-8},
        {"CGIDL", 0.5}, {"EGIDL", 0.1}, {"KF", 1.0e-20}
    }).prepare(1.0e-6, 1.0e-6);
    assert(model.js == 1.0e-2 && model.jsw == 1.0e-3);
    assert(model.agidl == 1.0e-6 && model.bgidl == 1.0e-8);
    gspice::Bsim4Mosfet device("M1", 0, 1, 2, 3, model);
    gspice::VectorReal x(4);
    x[0] = 0.8;
    x[1] = 0.0;
    x[3] = 0.1;
    gspice::DaeRequest request;
    request.staticResidual = true;
    gspice::DaeEvaluation evaluation;
    assert(device.evaluateDae(x, request, evaluation));
    double sum = 0.0;
    for (const auto& term : evaluation.staticResidual) sum += term.value;
    assert(std::abs(sum) < 1.0e-18);
    const double psd = device.getNoisePSD(2.0 * 3.141592653589793 * 1.0e9, x);
    assert(std::isfinite(psd) && psd >= 0.0);
    std::vector<gspice::NoiseSource> sources;
    device.collectNoiseSources(0.0, x, sources);
    assert(sources.size() == (psd > 0.0 ? 1u : 0u));
}
