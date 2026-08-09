#include "devices/gdi_device.hpp"
#include "devices/psp103_gsdi.hpp"
#include "devices/psp103_model.hpp"
#include "gsdi_device_adapter.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace {

bool near(double a, double b) { return std::abs(a - b) < 1e-10; }

}  // namespace

int main() {
    // --- 1. Descriptor metadata ---
    const auto& descriptor = gspice::psp103GsdiDescriptor();
    assert(descriptor.model_type == "PSP103VA");
    assert(descriptor.terminal_count == 4);
    assert(descriptor.nodes.size() == 4);
    for (std::size_t i = 0; i < descriptor.nodes.size(); ++i) {
        assert(descriptor.nodes[i].role == gspice::GsdiNodeRole::Terminal);
        assert(descriptor.nodes[i].index == static_cast<int>(i));
        assert(descriptor.nodes[i].collapse_partner == gspice::GsdiNoCollapse);
    }
    assert(descriptor.nodes[0].name == "d");
    assert(descriptor.nodes[1].name == "g");
    assert(descriptor.nodes[2].name == "s");
    assert(descriptor.nodes[3].name == "b");
    assert(descriptor.jacobian_pattern.size() == 16);
    bool dense = true;
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            const bool found = std::any_of(
                descriptor.jacobian_pattern.begin(), descriptor.jacobian_pattern.end(),
                [row, column](const std::pair<int, int>& entry) {
                    return entry.first == row && entry.second == column;
                });
            dense = dense && found;
        }
    }
    assert(dense);
    assert(descriptor.supports_op);
    assert(descriptor.supports_transient);
    assert(descriptor.supports_ac);
    assert(descriptor.supports_noise);

    // --- 2. Standard collapse plan on a pure 4-terminal descriptor is the
    //        identity: every local node keeps its own column, nothing merges.
    const gspice::GsdiCollapseMap plan =
        gspice::GsdiCollapseMap::standard(gspice::psp103GsdiDescriptor());
    assert(plan.localCount() == 4);
    assert(plan.unknownCount() == 4);
    for (int local = 0; local < 4; ++local) {
        assert(plan.collapsedIndex(local) == local);
        assert(plan.localNodeOf(local) == local);
        assert(plan.valueSource(local) == local);
    }

    // --- 3. Numeric identity: direct native evaluator on local nodes 0-3
    //        versus the same evaluator routed through GsdiDaeDeviceAdapter +
    //        GdiDevice (the Phase E recast path). Both must produce identical
    //        DAE residuals and Jacobians for the same global solution. The
    //        routed device is wired to global columns {0,1,2,3} with the
    //        identity collapse plan, so local index k == global column k.
    const auto params = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VFB", -0.85}, {"TOX", 2e-9}, {"PHIB", 0.7}, {"BET", 1e-3}
    });
    const gspice::GsdiParamMap instance = {{"W", 1e-6}, {"L", 1e-6}};

    auto direct = std::make_unique<gspice::Psp103Mosfet>(
        "M1", 0, 1, 2, 3, params, instance, 27.0);
    gspice::GdiDevice routed(
        "M1",
        std::make_unique<gspice::GsdiDaeDeviceAdapter>(
            std::make_unique<gspice::Psp103Mosfet>("M1", 0, 1, 2, 3, params, instance, 27.0),
            4),
        std::vector<int>{0, 1, 2, 3},
        gspice::GsdiCollapseMap::standard(gspice::psp103GsdiDescriptor()));

    // A mixed-bias operating point: drain up, gate well above threshold,
    // source and bulk pinned.
    gspice::VectorReal x(4);
    x[0] = 0.5;
    x[1] = 1.2;
    x[2] = 0.0;
    x[3] = 0.0;

    gspice::DaeRequest request;
    request.analysis = gspice::DaeAnalysis::OperatingPoint;
    request.staticResidual = true;
    request.staticJacobian = true;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;

    gspice::DaeEvaluation directEval;
    gspice::DaeEvaluation routedEval;
    assert(direct->evaluateDae(x, request, directEval));
    assert(routed.evaluateDae(x, request, routedEval));

    // Static/dynamic residuals and Jacobians must match term-for-term.
    assert(directEval.staticResidual.size() == routedEval.staticResidual.size());
    for (std::size_t i = 0; i < directEval.staticResidual.size(); ++i) {
        if (!near(directEval.staticResidual[i].value, routedEval.staticResidual[i].value)) {
            assert(false);
            return 1;
        }
    }
    assert(directEval.staticJacobian.size() == routedEval.staticJacobian.size());
    for (std::size_t i = 0; i < directEval.staticJacobian.size(); ++i) {
        if (!near(directEval.staticJacobian[i].value, routedEval.staticJacobian[i].value)) {
            assert(false);
            return 1;
        }
    }
    assert(directEval.dynamicResidual.size() == routedEval.dynamicResidual.size());
    for (std::size_t i = 0; i < directEval.dynamicResidual.size(); ++i) {
        if (!near(directEval.dynamicResidual[i].value, routedEval.dynamicResidual[i].value)) {
            assert(false);
            return 1;
        }
    }
    assert(directEval.dynamicJacobian.size() == routedEval.dynamicJacobian.size());
    for (std::size_t i = 0; i < directEval.dynamicJacobian.size(); ++i) {
        if (!near(directEval.dynamicJacobian[i].value, routedEval.dynamicJacobian[i].value)) {
            assert(false);
            return 1;
        }
    }

    return 0;
}