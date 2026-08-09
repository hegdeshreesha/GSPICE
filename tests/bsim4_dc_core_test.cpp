#include "devices/bsim4_dc_core.hpp"

#include <cassert>
#include <cmath>

int main() {
    const auto model = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"VTH0", 0.4}, {"KP", 120.0e-6},
        {"U0", 500.0}, {"VSAT", 1.0e5}, {"TOXE", 1.0e-8}
    }).prepare(1.0e-6, 1.0e-6);
    const auto off = gspice::bsim4EvaluateDc(model, {0.8, 0.0, 0.0, 0.0});
    const auto on = gspice::bsim4EvaluateDc(model, {0.8, 1.0, 0.0, 0.0});
    const auto noSeriesResistance = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"VTH0", 0.4}, {"KP", 120.0e-6},
        {"U0", 500.0}, {"VSAT", 1.0e5}, {"TOXE", 1.0e-8},
        {"RDSMOD", 1.0}
    }).prepare(1.0e-6, 1.0e-6);
    const auto onNoSeries = gspice::bsim4EvaluateDc(
        noSeriesResistance, {0.8, 1.0, 0.0, 0.0});
    assert(off.valid && on.valid);
    assert(onNoSeries.valid && onNoSeries.current[0] > on.current[0]);
    assert(on.current[0] > off.current[0]);
    assert(std::abs(on.current[0] + on.current[2]) < 1.0e-18);
    for (const auto& row : on.jacobian)
        for (double value : row) assert(std::isfinite(value));
}
