#include "devices/bjt.hpp"
#include "devices/diode.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gspice;

namespace {

double diodeCurrent(const Diode& diode, double voltage) {
    VectorReal x(2);
    x[0] = voltage;
    double current = 0.0;
    assert(diode.probeCurrent(x, current));
    return current;
}

std::array<double, 3> bjtCurrents(const Bjt& bjt, double vc, double vb, double ve) {
    VectorReal x(3);
    x[0] = vc;
    x[1] = vb;
    x[2] = ve;
    DaeEvaluation eval;
    DaeRequest request;
    request.staticResidual = true;
    assert(const_cast<Bjt&>(bjt).evaluateDae(x, request, eval));
    std::array<double, 3> currents{};
    for (const auto& stamp : eval.staticResidual) {
        currents[static_cast<std::size_t>(stamp.equation)] = stamp.value;
    }
    return currents;
}

} // namespace

int main() {
    Diode ideal("Dideal", 0, 1, 1e-14, 1.0, 0.0, 0.0);
    Diode series("Drs", 0, 1, 1e-14, 1.0, 0.0, 100.0);
    assert(diodeCurrent(ideal, 0.8) > 100.0 * diodeCurrent(series, 0.8));

    Diode noBreakdown("Dnb", 0, 1, 1e-14, 1.0, 0.0, 0.0, 0.0);
    Diode breakdown("Dbv", 0, 1, 1e-14, 1.0, 0.0, 0.0, 0.5, 1e-6, 1.0);
    assert(std::abs(diodeCurrent(breakdown, -1.0)) > 1e6 * std::abs(diodeCurrent(noBreakdown, -1.0)));

    Bjt nominal("Qnom", 0, 1, 2, 1, 1e-16, 100.0, 1.0, 1.0, 1.0);
    Bjt highInjection("Qikf", 0, 1, 2, 1, 1e-16, 100.0, 1.0, 1.0, 1.0,
                      1.0, 0.0, 0.0, 0.0, 1e-6);
    const auto nominalHighBias = bjtCurrents(nominal, 1.0, 0.9, 0.0);
    const auto ikfHighBias = bjtCurrents(highInjection, 1.0, 0.9, 0.0);
    assert(std::abs(nominalHighBias[0]) > 10.0 * std::abs(ikfHighBias[0]));

    Bjt early("Qvaf", 0, 1, 2, 1, 1e-16, 100.0, 1.0, 1.0, 1.0,
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0);
    const auto lowVce = bjtCurrents(early, 0.5, 0.72, 0.0);
    const auto highVce = bjtCurrents(early, 2.0, 0.72, 0.0);
    assert(std::abs(highVce[0]) > 1.05 * std::abs(lowVce[0]));

    Bjt noLeak("Qnol", 0, 1, 2, 1, 1e-16, 100.0, 1.0, 1.0, 1.0);
    Bjt leak("Qleak", 0, 1, 2, 1, 1e-16, 100.0, 1.0, 1.0, 1.0,
             1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-13, 0.0, 1.5);
    const auto noLeakCurrents = bjtCurrents(noLeak, 1.0, 0.65, 0.0);
    const auto leakCurrents = bjtCurrents(leak, 1.0, 0.65, 0.0);
    assert(std::abs(leakCurrents[1]) > 10.0 * std::abs(noLeakCurrents[1]));

    std::cout << "primitive model parameter probe passed\n";
    return 0;
}
