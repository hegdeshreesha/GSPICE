#include "devices/bsim3_parameters.hpp"
#include <cassert>
#include <cmath>

int main() {
    const auto params = gspice::Bsim3ParameterSet::from({
        {"LEVEL", 49.0}, {"VTH0", 0.42}, {"KP", 120.0e-6},
        {"TOXE", 1.0e-8}, {"DL", 1.0e-8}, {"DW", 2.0e-8}
    });
    const auto prepared = params.prepare(1.0e-6, 1.0e-6, 27.0);
    assert(prepared.validation && prepared.leff > 0.0 && prepared.weff > 0.0);
    assert(prepared.beta > 0.0 && prepared.cox > 0.0);

    const auto cardUnits = gspice::Bsim3ParameterSet::from({{"LEVEL", 49.0}, {"U0", 500.0}})
                               .prepare(1.0e-6, 1.0e-6);
    assert(cardUnits.validation && cardUnits.u0 == 0.05);

    const auto binned = gspice::Bsim3ParameterSet::from({
        {"LEVEL", 49.0}, {"VTH0", 0.4}, {"LVTH0", 1.0e-7},
        {"NFACTOR", 1.0}, {"LNFACTOR", 2.0e-7}, {"JS", 1.0e-14},
        {"LJS", 1.0e-20}, {"AGIDL", 1.0e-6}, {"LAGIDL", 1.0e-12}
    }).prepare(1.0e-6, 1.0e-6);
    assert(binned.validation);
    assert(std::abs(binned.vth0 - 0.5) < 1.0e-12);
    assert(std::abs(binned.nfactor - 1.2) < 1.0e-12);
    assert(std::abs(binned.js - 2.0e-14) < 1.0e-20);
    assert(std::abs(binned.agidl - 2.0e-6) < 1.0e-18);

    const auto hotJunction = gspice::Bsim3ParameterSet::from({
        {"LEVEL", 49.0}, {"JS", 1.0e-14}, {"JSW", 1.0e-12},
        {"XTI", 3.0}, {"EG", 1.11}
    }).prepare(1.0e-6, 1.0e-6, 100.0);
    assert(hotJunction.validation);
    assert(hotJunction.js > 1.0e-14);
    assert(hotJunction.jsw > 1.0e-12);

    const auto invalid = gspice::Bsim3ParameterSet::from({{"LEVEL", 54.0}}).prepare(1.0e-6, 1.0e-6);
    assert(!invalid.validation && invalid.validation.reason == "only BSIM3 level 49 is supported");
}
