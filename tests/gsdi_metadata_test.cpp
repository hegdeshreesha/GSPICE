#include "gsdi.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

namespace {

// Fake PSP-like device: 4 terminals (d, g, s, b), one internal node (series
// resistance junction Td), one collapsible node (source-internal Tint).
// Channel current I = KP * exp(V(g) - V(s)) flows from the intrinsic drain Td
// to the source; the drain series resistor carries Ird = (V(Td) - V(d)) / Rd.
class FakePspModel final : public gspice::GsdiModel {
public:
    static const gspice::GsdiModelDescriptor& descriptorRef() {
        static const gspice::GsdiModelDescriptor descriptor = [] {
            gspice::GsdiModelDescriptor d;
            d.model_type = "FAKEPSP";
            d.version = "1.0";
            d.terminal_count = 4;
            d.nodes = {
                {"D", 0, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
                {"G", 1, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
                {"S", 2, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
                {"B", 3, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
                {"Td", 4, gspice::GsdiNodeRole::Internal, gspice::GsdiNoCollapse},
                {"Tint", 5, gspice::GsdiNodeRole::Collapsible, 3},
            };
            d.parameters = {
                {"VTH0", 0.4, "V", "threshold voltage", true, true, 0.0, false, 0.0},
                {"KP", 120e-6, "A/V^2", "transconductance", true, true, 0.0, false, 0.0},
                {"TOXE", 1e-9, "m", "oxide thickness", true, true, 0.0, false, 0.0},
                {"L", 1e-6, "m", "channel length", false, true, 0.0, false, 0.0},
                {"W", 1e-6, "m", "channel width", false, true, 0.0, false, 0.0},
            };
            d.opvars = {
                {"gm", "transconductance"},
                {"gds", "drain conductance"},
                {"vdsat", "saturation voltage"},
            };
            d.jacobian_pattern = {{0, 0}, {0, 4}, {2, 1}, {2, 2},
                                  {4, 0}, {4, 1}, {4, 2}, {4, 4}};
            d.supports_op = true;
            d.supports_transient = true;
            d.supports_ac = true;
            d.supports_noise = false;
            return d;
        }();
        return descriptor;
    }

    const gspice::GsdiModelDescriptor& descriptor() const override { return descriptorRef(); }

    std::unique_ptr<gspice::GsdiInstance> createInstance(
        const gspice::GsdiModelCard& card,
        const std::vector<int>& terminal_nodes) const override {
        (void)terminal_nodes;
        return std::make_unique<Instance>(card);
    }

private:
    class Instance final : public gspice::GsdiInstance {
    public:
        explicit Instance(const gspice::GsdiModelCard& card) {
            const auto it = card.parameters.find("KP");
            kp_ = it != card.parameters.end() ? it->second : 120e-6;
            rd_ = 1000.0;
        }

        bool evaluate(const gspice::GsdiEvalRequest& request,
                      gspice::GsdiEvalResult& result) override {
            if (!request.solution || request.solution_size < 6) return false;
            const double vd = request.solution[0];
            const double vg = request.solution[1];
            const double vs = request.solution[2];
            const double vtint = request.solution[4];

            const double i_channel = kp_ * std::exp(vg - vs);
            const double i_series = (vtint - vd) / rd_;

            result.clear();
            if (request.residual) {
                // KCL: current entering the device at each local node.
                result.static_residual.push_back({0, i_series, 0});
                result.static_residual.push_back({4, -(i_channel + i_series), 0});
                result.static_residual.push_back({2, i_channel, 0});
                result.static_residual.push_back({1, 0.0, 0});
                result.static_residual.push_back({3, 0.0, 0});
                result.static_residual.push_back({5, 0.0, 0});
            }
            if (request.jacobian) {
                const double gmi = i_channel;
                const double grd = 1.0 / rd_;
                result.static_jacobian.push_back({0, 0, -grd, 0});
                result.static_jacobian.push_back({0, 4, grd, 0});
                result.static_jacobian.push_back({4, 0, grd, 0});
                result.static_jacobian.push_back({4, 4, -grd, 0});
                result.static_jacobian.push_back({4, 1, -gmi, 0});
                result.static_jacobian.push_back({4, 2, gmi, 0});
                result.static_jacobian.push_back({2, 1, gmi, 0});
                result.static_jacobian.push_back({2, 2, -gmi, 0});
            }
            if (request.dynamic_residual) {
                result.dynamic_residual.push_back({1, 1e-15 * vg, 0});
            }
            if (request.dynamic_jacobian) {
                result.dynamic_jacobian.push_back({1, 1, 1e-15, 0});
            }
            result.opvars = {i_channel, 1e-9, vg - vs};
            return true;
        }

        std::size_t terminalCount() const override { return 4; }
        std::size_t internalNodeCount() const override { return 2; }

    private:
        double kp_ = 120e-6;
        double rd_ = 1000.0;
    };
};

}  // namespace

int main() {
    const auto& descriptor = FakePspModel::descriptorRef();
    assert(descriptor.model_type == "FAKEPSP");
    assert(descriptor.terminal_count == 4);
    assert(descriptor.nodes.size() == 6);
    assert(descriptor.nodes[0].role == gspice::GsdiNodeRole::Terminal);
    assert(descriptor.nodes[4].role == gspice::GsdiNodeRole::Internal);
    assert(descriptor.nodes[5].role == gspice::GsdiNodeRole::Collapsible);
    assert(descriptor.nodes[5].collapse_partner == 3);
    assert(descriptor.parameters.size() == 5);
    assert(descriptor.parameters[0].name == "VTH0");
    assert(descriptor.parameters[0].is_model_param);
    assert(!descriptor.parameters[3].is_model_param);
    assert(descriptor.parameters[0].has_min && descriptor.parameters[0].min_value == 0.0);
    assert(descriptor.opvars.size() == 3);
    assert(descriptor.opvars[0].name == "gm");
    assert(descriptor.jacobian_pattern.size() == 8);
    assert(descriptor.supports_transient && !descriptor.supports_noise);

    // Default collapse decisions: every node keeps its own column.
    const gspice::GsdiCollapseMap keep_all(descriptor, {});
    assert(keep_all.unknownCount() == 6);
    for (int node = 0; node < 6; ++node) {
        assert(keep_all.collapsedIndex(node) == node);
        assert(keep_all.localNodeOf(node) == node);
    }

    // Collapse the source-internal pair (node 5 into node 3).
    const gspice::GsdiCollapseMap collapse_pair(
        descriptor, {gspice::GsdiNoCollapse, gspice::GsdiNoCollapse,
                     gspice::GsdiNoCollapse, gspice::GsdiNoCollapse,
                     gspice::GsdiNoCollapse, 3});
    assert(collapse_pair.unknownCount() == 5);
    assert(collapse_pair.collapsedIndex(5) == gspice::GsdiNoCollapse);
    assert(collapse_pair.collapsedIndex(3) == 3);
    assert(collapse_pair.localNodeOf(3) == 3);
    assert(collapse_pair.localNodeOf(4) == 4);

    // Collapse the internal node into the reference (ground).
    const gspice::GsdiCollapseMap collapse_ground(
        descriptor, {gspice::GsdiNoCollapse, gspice::GsdiNoCollapse,
                     gspice::GsdiNoCollapse, gspice::GsdiNoCollapse,
                     gspice::GsdiCollapseToGround, gspice::GsdiNoCollapse});
    assert(collapse_ground.unknownCount() == 5);
    assert(collapse_ground.collapsedIndex(4) == gspice::GsdiNoCollapse);
    assert(collapse_ground.collapsedIndex(5) == 4);
    assert(collapse_ground.localNodeOf(4) == 5);

    // Full evaluation with all six local unknowns present.
    FakePspModel model;
    gspice::GsdiModelCard card;
    card.name = "FAKE1";
    card.type = "FAKEPSP";
    card.parameters["KP"] = 120e-6;
    auto instance = model.createInstance(card, {0, 1, 2, 3});

    const double voltages[] = {0.3, 0.8, 0.0, 0.0, 0.45, 0.0};
    gspice::GsdiEvalRequest request;
    request.analysis = gspice::GsdiAnalysis::Transient;
    request.solution = voltages;
    request.solution_size = 6;
    request.residual = true;
    request.jacobian = true;
    request.dynamic_residual = true;
    request.dynamic_jacobian = true;

    gspice::GsdiEvalResult result;
    assert(instance->evaluate(request, result));
    assert(result.finite());

    const double expected_i = 120e-6 * std::exp(0.8);
    const double expected_iseries = (0.45 - 0.3) / 1000.0;
    assert(std::abs(result.static_residual[0].value - expected_iseries) < 1e-20);
    assert(std::abs(result.static_residual[2].value - expected_i) < 1e-20);
    assert(std::abs(result.static_residual[1].value +
                    (expected_i + expected_iseries)) < 1e-20);

    // KCL conservation across all local nodes.
    double total = 0.0;
    for (const auto& term : result.static_residual) total += term.value;
    assert(std::abs(total) < 1e-24);

    // Jacobian row sums must be zero when the device has no terminal currents
    // (all local unknowns shifted together).
    for (std::size_t row = 0; row < 6; ++row) {
        double row_sum = 0.0;
        for (const auto& term : result.static_jacobian) {
            if (term.equation == static_cast<int>(row)) row_sum += term.value;
        }
        assert(std::abs(row_sum) < 1e-20);
    }

    assert(result.dynamic_residual.size() == 1);
    assert(result.dynamic_residual[0].equation == 1);
    assert(std::abs(result.dynamic_residual[0].value - 1e-15 * 0.8) < 1e-27);
    assert(result.dynamic_jacobian.size() == 1);

    assert(result.opvars.size() == 3);
    assert(std::abs(result.opvars[0] - expected_i) < 1e-24);
    assert(result.opvars[2] == 0.8);

    // Bypass and limiting flags stay off for a plain evaluation.
    assert(!result.bypassed);
    assert(!result.limiting_applied);
    return 0;
}
