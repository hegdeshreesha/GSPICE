#include "devices/gmc_juncap_express.hpp"

#include <cassert>
#include <cmath>

int main() {
    gspice::GmcJuncapExpressModule module;
    const auto instance = module.createInstance(
        {{"PHITD", 0.026}, {"ISATFOR1", 1e-14}, {"MFOR1", 1.0}},
        {{"VJ", 0.2}, {"DVJ", 1.0}});
    double voltages[] = {0.6, 0.0};
    gspice::GdiEvaluationRequest request;
    request.dynamic_residual = true;
    request.dynamic_jacobian = true;
    request.noise = true;
    gspice::GdiEvaluationResult result;
    assert(instance->evaluate(voltages, request, result));
    assert(result.static_residual.size() == 2);
    assert(result.static_jacobian.size() == 4);
    assert(result.dynamic_residual.size() == 2);
    assert(result.dynamic_jacobian.size() == 4);
    assert(result.noise.size() == 1);
    assert(std::abs(result.static_residual[0].value + result.static_residual[1].value) < 1e-30);
    assert(std::abs(result.dynamic_residual[0].value + result.dynamic_residual[1].value) < 1e-30);
    return 0;
}
