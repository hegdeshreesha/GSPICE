#include "devices/psp103_temperature.hpp"

#include <cassert>
#include <cmath>

int main() {
    const gspice::Psp103SelfHeatingInputs dc{2.0, 0.0, 0.0, 10.0, 1.0};
    const auto dc_result = gspice::psp103SelfHeatingResidual(dc, 20.0);
    assert(std::abs(dc_result.residual_w) < 1e-12);
    assert(std::abs(dc_result.dresidual_drise_w_per_k - 0.1) < 1e-12);

    const gspice::Psp103SelfHeatingInputs transient{0.0, 5.0, 2.0, 10.0, 4.0};
    const auto transient_result = gspice::psp103SelfHeatingResidual(transient, 9.0);
    assert(std::abs(transient_result.residual_w - 8.9) < 1e-12);
    assert(std::abs(transient_result.dresidual_drise_w_per_k - 2.1) < 1e-12);
    return 0;
}
