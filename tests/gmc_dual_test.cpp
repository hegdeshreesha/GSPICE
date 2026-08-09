#include "gmc_dual.hpp"

#include <cassert>
#include <cmath>

namespace {

double fdExpSqrtSum(double vg, double vd, double vb) {
    return std::exp(vg - vd) * std::sqrt(vd) + std::tanh(vb) * vb;
}

}  // namespace

int main() {
    using gspice::GmcDual4;
    const auto vgs = GmcDual4::variable(0.7, 0);
    const auto vds = GmcDual4::variable(0.4, 1);
    const auto expression = gspice::gmcExp(vgs) * gspice::gmcSqrt(vds) + 2.0;

    assert(std::abs(expression.value - (std::exp(0.7) * std::sqrt(0.4) + 2.0)) < 1e-12);
    assert(std::abs(expression.derivative[0] - std::exp(0.7) * std::sqrt(0.4)) < 1e-12);
    assert(std::abs(expression.derivative[1] - std::exp(0.7) / (2.0 * std::sqrt(0.4))) < 1e-12);
    assert(expression.derivative[2] == 0.0);
    assert(expression.derivative[3] == 0.0);

    // BSIM-class sizing: eight local unknowns (4 terminals + internal nodes).
    using gspice::GmcDual8;
    const auto vg = GmcDual8::variable(0.8, 0);
    const auto vd = GmcDual8::variable(0.4, 1);
    const auto vb = GmcDual8::variable(-0.3, 2);
    const auto model = gspice::gmcExp(vg - vd) * gspice::gmcSqrt(vd) + gspice::gmcTanh(vb) * vb;

    assert(model.value == fdExpSqrtSum(0.8, 0.4, -0.3));

    const double h = 1e-6;
    assert(std::abs(model.derivative[0] - (fdExpSqrtSum(0.8 + h, 0.4, -0.3) -
                                           fdExpSqrtSum(0.8 - h, 0.4, -0.3)) / (2.0 * h)) < 1e-6);
    assert(std::abs(model.derivative[1] - (fdExpSqrtSum(0.8, 0.4 + h, -0.3) -
                                           fdExpSqrtSum(0.8, 0.4 - h, -0.3)) / (2.0 * h)) < 1e-6);
    assert(std::abs(model.derivative[2] - (fdExpSqrtSum(0.8, 0.4, -0.3 + h) -
                                           fdExpSqrtSum(0.8, 0.4, -0.3 - h)) / (2.0 * h)) < 1e-6);

    // Limiting and branch-selection helpers used by compact models.
    const auto clamped = gspice::gmcLim(vg, 0.5, 0.6);
    assert(clamped.value == 0.6);
    assert(clamped.derivative[0] == 0.0);
    const auto in_range = gspice::gmcLim(GmcDual8::variable(0.55, 0), 0.5, 0.6);
    assert(in_range.value == 0.55);
    assert(in_range.derivative[0] == 1.0);
    assert(gspice::gmcAbs(-vb).value == 0.3);
    assert(gspice::gmcMax(GmcDual8::constant(0.2), vd).value == 0.4);
    assert(gspice::gmcMin(GmcDual8::constant(0.2), vd).value == 0.2);

    // Extended transcendental set used by the deepened GMC Verilog-A subset.
    const auto x = GmcDual4::variable(0.7, 1);
    const double hdg = 2e-6;
    const auto fd_deriv = [&](auto f, double v) {
        return (f(v + hdg) - f(v - hdg)) / (2.0 * hdg);
    };
    auto sin_x = gspice::gmcSin(x);
    assert(std::abs(sin_x.value - std::sin(0.7)) < 1e-12);
    assert(std::abs(sin_x.derivative[1] - std::cos(0.7)) < 1e-9);
    auto cos_x = gspice::gmcCos(x);
    assert(std::abs(cos_x.value - std::cos(0.7)) < 1e-12);
    assert(std::abs(cos_x.derivative[1] - fd_deriv([](double v) { return std::cos(v); }, 0.7)) < 1e-9);
    auto tan_x = gspice::gmcTan(GmcDual4::variable(0.3, 0));
    assert(std::abs(tan_x.value - std::tan(0.3)) < 1e-12);
    assert(std::abs(tan_x.derivative[0] - (1.0 + std::tan(0.3) * std::tan(0.3))) < 1e-9);
    auto sinh_x = gspice::gmcSinh(x);
    assert(std::abs(sinh_x.value - std::sinh(0.7)) < 1e-12);
    assert(std::abs(sinh_x.derivative[1] - std::cosh(0.7)) < 1e-9);
    auto cosh_x = gspice::gmcCosh(x);
    assert(std::abs(cosh_x.value - std::cosh(0.7)) < 1e-12);
    assert(std::abs(cosh_x.derivative[1] - std::sinh(0.7)) < 1e-9);
    return 0;
}
