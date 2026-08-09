#include "devices/bsim3_dc_core.hpp"
#include <cassert>
#include <cmath>

int main() {
    auto model = gspice::Bsim3ParameterSet::from({
        {"LEVEL", 49.0}, {"VTH0", 0.4}, {"KP", 120e-6}, {"TOXE", 1e-8},
        {"U0", 0.05}, {"VSAT", 1e5}, {"NFACTOR", 1.4}
    }).prepare(1e-6, 1e-6);
    const auto off = gspice::bsim3EvaluateDc(model, {0.8, 0.0, 0.0, 0.0});
    const auto on = gspice::bsim3EvaluateDc(model, {0.8, 0.8, 0.0, 0.0});
    assert(off.valid && on.valid);
    assert(std::abs(off.current[0] + off.current[2]) < 1e-24);
    assert(std::abs(on.current[0] + on.current[2]) < 1e-24);
    assert(on.current[0] > off.current[0]);
    assert(std::isfinite(on.jacobian[0][1]));

    auto leakageModel = gspice::Bsim3ParameterSet::from({
        {"LEVEL", 49.0}, {"VTH0", 0.4}, {"KP", 120e-6}, {"TOXE", 1e-8},
        {"U0", 0.05}, {"VSAT", 1e5}, {"IS", 1e-14}, {"JS", 1e-2},
        {"JSW", 1e-3}, {"NJ", 1.0}, {"AIGC", 1e-6}, {"BIGC", 1.0},
        {"CIGC", 0.0}, {"AGIDL", 1e-6}, {"BGIDL", 1e-8},
        {"CGIDL", 0.5}, {"EGIDL", 0.1}
    }).prepare(1e-6, 1e-6);
    const auto leakage = gspice::bsim3EvaluateDc(leakageModel, {0.8, 0.0, 0.0, 0.1});
    assert(leakage.valid);
    double sum = 0.0;
    for (double current : leakage.current) sum += current;
    assert(std::abs(sum) < 1e-18);
    for (const auto& row : leakage.jacobian)
        for (double entry : row)
            assert(std::isfinite(entry));
}
