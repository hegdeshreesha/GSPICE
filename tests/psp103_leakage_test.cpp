#include "devices/psp103_leakage.hpp"

#include <cassert>
#include <cmath>

int main() {
    const gspice::Psp103GidlInputs gidl{-1.0, 1.0, 0.0, 2.0, 10.0};
    assert(gspice::psp103GidlCurrent(gidl) > 0.0);
    assert(gspice::psp103GidlCurrent({1.0, 1.0, 0.0, 2.0, 10.0}) == 0.0);

    const gspice::Psp103GateOverlapInputs overlap{1.0, 0.1, -0.5, 1e-12, 1.0, 0.2, 3.1};
    assert(std::isfinite(gspice::psp103GateOverlapCurrent(overlap)));

    const gspice::Psp103GateChannelInputs channel{1.0, 0.8, 0.1, 0.026, 1.0, 0.2,
        0.4, 1.0, 1.0, 0.2, 1.0, 1.0, 0.1, 0.2, 3.1, 1e-12};
    const auto gate = gspice::psp103GateChannelCurrent(channel);
    assert(std::isfinite(gate.igc) && std::isfinite(gate.igd));
    assert(std::abs(gate.igs + gate.igd - gate.igc) < 1e-18);

    gspice::GsdiEvalResult result;
    gspice::psp103AppendCurrentBranch(result, 1, 2, 1, 2, 3.0, 4.0, 7);
    assert(result.static_residual.size() == 2);
    assert(result.static_jacobian.size() == 4);
    assert(result.static_residual[0].value + result.static_residual[1].value == 0.0);
    assert(result.finite());
    return 0;
}
