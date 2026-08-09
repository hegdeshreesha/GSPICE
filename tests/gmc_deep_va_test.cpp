#include "devices/gdi_device.hpp"
#include "gmc_deep.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace {

// Mirror of the deepened Verilog-A surface in tests/gmc_deep.va. Keeping the
// expected model in plain doubles makes the gate a genuine cross-check of the
// generated C++ and its GmcDual derivatives.
constexpr double kG = 2e-3;
constexpr double kCj = 3e-6;
constexpr double kNwhite = 4e-12;
constexpr double kNflicker = 9e-12;
constexpr double kTemperature = 300.0;
constexpr double kPi = 3.14159265358979323846;

// Mirror of the analog functions in tests/gmc_deep.va.
double gclamp(double x, double lo, double hi) {
    return std::min(std::max(x, lo), hi);
}

double gsoft(double x) {
    return 1.0 + 0.5 * gclamp(x, 0.0, 1.0);
}

double expectedW(double v) {
    return std::exp(v) + std::log(v + 1.0) + std::sqrt(std::abs(v) + 1.0)
        + std::tanh(v) + std::sin(v) + std::cosh(v) + std::sinh(v) / 4.0
        + (v > 0.5 ? v * v : 0.5)
        + std::max(v - 0.5, 0.0) + std::min(v, 1.5)
        + std::log10(v * 10.0) + std::log(v + 2.0)
        + std::exp(std::clamp(v - 2.0, -40.0, 40.0))
        + kTemperature / 300.0 - 1.0
        + std::pow(v, 2) + gsoft(v - 0.25);
}

double expectedI(double v) {
    // sel = 2 selects the "2, 3:" case item -> 1e-9 * v (default item unused).
    return kG * expectedW(v) + 1e-9 * v + (v > 0.0 ? 1e-12 * v : 0.0);
}

double expectedDIdv(double v) {
    const double h = 2e-6;
    return (expectedI(v + h) - expectedI(v - h)) / (2.0 * h);
}

}  // namespace

int main() {
    gspice::GmcDeepProbeModel model;

    const auto& descriptor = model.descriptor();
    assert(descriptor.model_type == "deep_probe");
    assert(descriptor.terminal_count == 2);
    assert(descriptor.nodes.size() == 2);
    assert(descriptor.parameters.size() == 4);
    assert(descriptor.parameters[0].name == "G");
    assert(descriptor.parameters[0].default_value == 2e-3);
    assert(descriptor.parameters[1].name == "Cj");
    assert(descriptor.parameters[1].default_value == 3e-6);
    assert(descriptor.supports_noise);
    assert(descriptor.jacobian_pattern.size() == 4);
    assert(descriptor.supports_transient);

    gspice::GsdiModelCard card;
    card.name = "DEEP1";
    card.type = "deep_probe";
    card.parameters["G"] = kG;
    card.parameters["Cj"] = kCj;
    card.parameters["Nwhite"] = kNwhite;
    card.parameters["Nflicker"] = kNflicker;
    auto instance = model.createInstance(card, {3, 5});
    assert(instance != nullptr);
    assert(instance->terminalCount() == 2);
    assert(instance->internalNodeCount() == 0);

    const double v = 1.0;
    const double voltages[] = {1.5, 0.5};
    gspice::GsdiEvalRequest request;
    request.solution = voltages;
    request.solution_size = 2;
    request.residual = true;
    request.jacobian = true;
    request.dynamic_residual = true;
    request.dynamic_jacobian = true;
    gspice::GsdiEvalResult result;
    assert(instance->evaluate(request, result));
    assert(result.finite());

    assert(result.static_residual.size() == 2);
    assert(result.static_jacobian.size() == 4);
    assert(result.dynamic_residual.size() == 2);
    assert(result.dynamic_jacobian.size() == 4);

    const double expected_i = expectedI(v);
    assert(std::abs(result.static_residual[0].value - expected_i) < 1e-12);
    assert(std::abs(result.static_residual[1].value + expected_i) < 1e-12);
    assert(std::abs(result.dynamic_residual[0].value - kCj * v) < 1e-21);
    assert(std::abs(result.dynamic_residual[1].value + kCj * v) < 1e-21);

    // Analytical GmcDual derivative vs central finite differences of the
    // plain-double mirror model.
    double j00 = 0.0, q00 = 0.0;
    for (const auto& term : result.static_jacobian) {
        if (term.equation == 0 && term.unknown == 0) j00 += term.value;
    }
    for (const auto& term : result.dynamic_jacobian) {
        if (term.equation == 0 && term.unknown == 0) q00 += term.value;
    }
    assert(std::abs(j00 - expectedDIdv(v)) < 1e-9);
    assert(std::abs(q00 - kCj) < 1e-18);

    {
        gspice::GsdiEvalRequest noise_req;
        noise_req.solution = voltages;
        noise_req.solution_size = 2;
        noise_req.residual = false;
        noise_req.jacobian = false;
        noise_req.noise = true;
        noise_req.omega = 2.0 * kPi * 3.0;
        gspice::GsdiEvalResult noise_res;
        assert(instance->evaluate(noise_req, noise_res));
        assert(noise_res.noise.size() == 2);
        assert(noise_res.noise[0].node_pos == 0);
        assert(noise_res.noise[0].node_neg == 1);
        assert(noise_res.noise[0].name == "thermal");
        assert(std::abs(noise_res.noise[0].spectral_density - kNwhite) < 1e-24);
        assert(noise_res.noise[1].name == "flicker");
        assert(std::abs(noise_res.noise[1].spectral_density - kNflicker / 3.0) < 1e-24);
    }

    // Re-evaluate in the region where gsoft(v - 0.25) sits against its lower
    // clamp bound, so the analog-function branch is exercised off-corner.
    // v = 0.1 -> gsoft(-0.15) = 1.0 (clamped), constant derivative (0); keeps
    // log10(v * 10.0) = log10(1) in domain (v must stay > 0).
    {
        const double v2 = 0.1;
        const double voltages2[] = {0.1, 0.0};
        gspice::GsdiEvalRequest req2;
        req2.solution = voltages2;
        req2.solution_size = 2;
        req2.residual = true;
        req2.jacobian = true;
        gspice::GsdiEvalResult res2;
        assert(instance->evaluate(req2, res2));
        assert(res2.finite());
        const double i2 = expectedI(v2);
        assert(std::abs(res2.static_residual[0].value - i2) < 1e-12);
        double j00_2 = 0.0;
        for (const auto& term : res2.static_jacobian) {
            if (term.equation == 0 && term.unknown == 0) j00_2 += term.value;
        }
        assert(std::abs(j00_2 - expectedDIdv(v2)) < 1e-9);
    }

    // The generated GsdiInstance stamps correctly through the GdiDevice
    // adapter into global MNA columns.
    gspice::GdiDevice device("D1", model.createInstance(card, {3, 5}), {3, 5});
    gspice::VectorReal x(6);
    x[3] = 1.5;
    x[5] = 0.5;
    gspice::DaeRequest dae;
    dae.staticResidual = true;
    dae.staticJacobian = true;
    dae.dynamicResidual = true;
    dae.dynamicJacobian = true;
    gspice::DaeEvaluation evaluation;
    assert(device.evaluateDae(x, dae, evaluation));
    double g33 = 0.0, q33 = 0.0;
    for (const auto& term : evaluation.staticJacobian) {
        if (term.equation == 3 && term.unknown == 3) g33 += term.value;
    }
    for (const auto& term : evaluation.dynamicJacobian) {
        if (term.equation == 3 && term.unknown == 3) q33 += term.value;
    }
    assert(std::abs(g33 - expectedDIdv(v)) < 1e-9);
    assert(std::abs(q33 - kCj) < 1e-18);
    return 0;
}
