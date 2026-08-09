#include "devices/bsim_translator.hpp"
#include <cassert>

int main() {
    const auto b3 = gspice::BsimModelTranslator::translate({
        {"level", 49.0}, {"vto", 0.42}, {"beta", 100.0e-6}, {"tox", 1.0e-8}
    });
    assert(b3 && b3.family == gspice::BsimFamily::Bsim3);
    assert(b3.parameters.at("VTH0") == 0.42);
    assert(b3.parameters.at("KP") == 100.0e-6);
    assert(b3.parameters.at("TOXE") == 1.0e-8);

    const auto b4 = gspice::BsimModelTranslator::translate({
        {"LEVEL", 54.0}, {"VTH0", 0.45}, {"KP", 200.0e-6}
    });
    assert(b4 && b4.family == gspice::BsimFamily::Bsim4);

    const auto invalid = gspice::BsimModelTranslator::translate({{"LEVEL", 1.0}});
    assert(!invalid && !invalid.reason.empty());
}
