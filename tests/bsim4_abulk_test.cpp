#include "devices/bsim4_dc_core.hpp"

#include <cassert>
#include <cmath>

int main() {
    const auto base = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"U0", 500.0}, {"TOXE", 1.0e-8}
    }).prepare(1.0e-6, 1.0e-6);
    const auto altered = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"U0", 500.0}, {"TOXE", 1.0e-8},
        {"A0", 2.0}, {"B0", 0.2e-6}, {"B1", 1.0e-6}, {"AGS", 0.2}
    }).prepare(1.0e-6, 1.0e-6);
    const auto a = gspice::bsim4EvaluateDc(base, {0.8, 0.8, 0.0, 0.0});
    const auto b = gspice::bsim4EvaluateDc(altered, {0.8, 0.8, 0.0, 0.0});
    assert(a.valid && b.valid);
    assert(std::isfinite(a.current[0]) && std::isfinite(b.current[0]));
    assert(std::abs(a.current[0] - b.current[0]) > 1.0e-18);
    assert(std::abs(b.current[0] + b.current[2]) < 1.0e-18);
    for (const auto& row : b.jacobian)
        for (double value : row) assert(std::isfinite(value));
}
