#include "devices/gdi_device.hpp"
#include "gmc_va_probe.hpp"

#include <cassert>
#include <cmath>

int main() {
    gspice::GmcVaProbeModel model;

    // GSDI-native descriptor metadata emitted by GMC from the .va module.
    const auto& descriptor = model.descriptor();
    assert(descriptor.model_type == "va_probe");
    assert(descriptor.version == "1.0");
    assert(descriptor.terminal_count == 2);
    assert(descriptor.nodes.size() == 2);
    assert(descriptor.nodes[0].name == "p");
    assert(descriptor.nodes[0].role == gspice::GsdiNodeRole::Terminal);
    assert(descriptor.nodes[0].collapse_partner == gspice::GsdiNoCollapse);
    assert(descriptor.nodes[1].name == "n");
    assert(descriptor.parameters.size() == 2);
    assert(descriptor.parameters[0].name == "G");
    assert(descriptor.parameters[0].default_value == 2e-3);
    assert(descriptor.parameters[0].is_model_param);
    assert(descriptor.parameters[1].name == "Cj");
    assert(descriptor.parameters[1].default_value == 3e-12);
    assert(descriptor.jacobian_pattern.size() == 4);
    assert(descriptor.supports_op);
    assert(descriptor.supports_transient);  // ddt(Cj * vd) is present
    assert(descriptor.supports_ac);
    assert(!descriptor.supports_noise);

    // Route the generated model through GsdiModel / GsdiInstance directly.
    gspice::GsdiModelCard card;
    card.name = "VA1";
    card.type = "va_probe";
    card.parameters["G"] = 2e-3;
    card.parameters["Cj"] = 3e-12;
    auto instance = model.createInstance(card, {3, 5});
    assert(instance != nullptr);
    assert(instance->terminalCount() == 2);
    assert(instance->internalNodeCount() == 0);

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

    assert(std::abs(result.static_residual[0].value - 2e-3) < 1e-15);
    assert(std::abs(result.static_residual[1].value + 2e-3) < 1e-15);
    assert(std::abs(result.dynamic_residual[0].value - 3e-12) < 1e-24);
    assert(std::abs(result.dynamic_residual[1].value + 3e-12) < 1e-24);
    assert(std::abs(result.static_jacobian[0].value - 2e-3) < 1e-12);
    assert(std::abs(result.static_jacobian[0].unknown - 0) < 1e-12);
    assert(std::abs(result.dynamic_jacobian[0].value - 3e-12) < 1e-21);

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
    assert(evaluation.staticJacobian.size() == 4);
    assert(evaluation.dynamicJacobian.size() == 4);

    double j33 = 0.0, j35 = 0.0, j53 = 0.0, j55 = 0.0;
    double q33 = 0.0;
    for (const auto& term : evaluation.staticJacobian) {
        if (term.equation == 3 && term.unknown == 3) j33 += term.value;
        if (term.equation == 3 && term.unknown == 5) j35 += term.value;
        if (term.equation == 5 && term.unknown == 3) j53 += term.value;
        if (term.equation == 5 && term.unknown == 5) j55 += term.value;
    }
    for (const auto& term : evaluation.dynamicJacobian) {
        if (term.equation == 3 && term.unknown == 3) q33 += term.value;
    }
    assert(std::abs(j33 - 2e-3) < 1e-12);
    assert(std::abs(j35 + 2e-3) < 1e-12);
    assert(std::abs(j53 + 2e-3) < 1e-12);
    assert(std::abs(j55 - 2e-3) < 1e-12);
    assert(std::abs(q33 - 3e-12) < 1e-21);
    return 0;
}
