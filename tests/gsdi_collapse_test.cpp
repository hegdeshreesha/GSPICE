#include "devices/gdi_device.hpp"
#include "gsdi.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

namespace {

// Linear 3-local-node model: two terminals (A, B) plus one collapsible node
// (Ci). Each node pair is bridged by a conductance, so the device is a
// weighted Laplacian (KCL-conservative). A capacitance C shunts Ci-A and is
// charge-conservative. With Ci collapsed into A (V(Ci,A) <+ 0) the branches
// A-Ci carry nothing and the supernode's conductance to B is g01+g12.
struct FakeCoupledCoupled {
    double g01, g02, g12, c;
};

class CoupledFake final : public gspice::GsdiInstance {
public:
    CoupledFake(double g01, double g02, double g12, double c)
        : g01_(g01), g02_(g02), g12_(g12), c_(c) {}

    bool evaluate(const gspice::GsdiEvalRequest& request,
                  gspice::GsdiEvalResult& result) override {
        last_ = std::vector<double>(request.solution,
                                    request.solution + request.solution_size);
        const double v0 = request.solution[0];
        const double v1 = request.solution[1];
        const double v2 = request.solution[2];

        result.clear();
        if (request.residual) {
            result.static_residual.push_back({0, g01_ * (v0 - v1) + g02_ * (v0 - v2), 0});
            result.static_residual.push_back({1, g01_ * (v1 - v0) + g12_ * (v1 - v2), 0});
            result.static_residual.push_back({2, g02_ * (v2 - v0) + g12_ * (v2 - v1), 0});
        }
        if (request.jacobian) {
            result.static_jacobian.push_back({0, 0, g01_ + g02_, 0});
            result.static_jacobian.push_back({0, 1, -g01_, 0});
            result.static_jacobian.push_back({0, 2, -g02_, 0});
            result.static_jacobian.push_back({1, 0, -g01_, 0});
            result.static_jacobian.push_back({1, 1, g01_ + g12_, 0});
            result.static_jacobian.push_back({1, 2, -g12_, 0});
            result.static_jacobian.push_back({2, 0, -g02_, 0});
            result.static_jacobian.push_back({2, 1, -g12_, 0});
            result.static_jacobian.push_back({2, 2, g02_ + g12_, 0});
        }
        if (request.dynamic_residual) {
            result.dynamic_residual.push_back({0, -c_ * (v2 - v0), 0});
            result.dynamic_residual.push_back({2, c_ * (v2 - v0), 0});
        }
        if (request.dynamic_jacobian) {
            result.dynamic_jacobian.push_back({0, 0, c_, 0});
            result.dynamic_jacobian.push_back({0, 2, -c_, 0});
            result.dynamic_jacobian.push_back({2, 0, -c_, 0});
            result.dynamic_jacobian.push_back({2, 2, c_, 0});
        }
        return true;
    }

    std::size_t terminalCount() const override { return 2; }
    std::size_t internalNodeCount() const override { return 1; }

    const std::vector<double>& lastSolution() const { return last_; }

private:
    double g01_, g02_, g12_, c_;
    std::vector<double> last_;
};

gspice::GsdiModelDescriptor makeDescriptor(int partner) {
    gspice::GsdiModelDescriptor d;
    d.model_type = "COUPLED";
    d.terminal_count = 2;
    d.nodes = {
        {"A", 0, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
        {"B", 1, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
        {"Ci", 2, gspice::GsdiNodeRole::Collapsible, partner},
    };
    return d;
}

gspice::GdiDevice makeDevice(std::unique_ptr<CoupledFake> instance,
                             const std::vector<int>& nodes,
                             const gspice::GsdiCollapseMap* collapse) {
    if (collapse) {
        return gspice::GdiDevice("D1", std::move(instance), nodes, *collapse, nullptr);
    }
    return gspice::GdiDevice("D1", std::move(instance), nodes, nullptr);
}

gspice::DaeEvaluation run(gspice::GdiDevice& device, const std::vector<double>& x_values) {
    gspice::VectorReal x(static_cast<int>(x_values.size()));
    for (std::size_t i = 0; i < x_values.size(); ++i) x[static_cast<int>(i)] = x_values[i];
    gspice::DaeRequest request;
    request.analysis = gspice::DaeAnalysis::OperatingPoint;
    request.staticResidual = true;
    request.staticJacobian = true;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    gspice::DaeEvaluation evaluation;
    const bool ok = device.evaluateDae(x, request, evaluation);
    assert(ok);
    return evaluation;
}

// Assemble the stamped static Jacobian into row-major form indexed by position
// of each node inside `nodes`. A reference-node (-1) term is folded into a
// dedicated ground slot at index `n`. Returns false only if a stamp referenced
// a node that is not in `nodes` and not the reference (a leak out of the
// allocated MNA columns).
bool assembleJ(const gspice::DaeEvaluation& evaluation, const std::vector<int>& nodes,
               std::vector<std::vector<double>>& J) {
    const int n = static_cast<int>(nodes.size());
    const std::size_t ground = static_cast<std::size_t>(n);
    J.assign(static_cast<std::size_t>(n) + 1,
             std::vector<double>(static_cast<std::size_t>(n) + 1, 0.0));
    const auto slot = [&](int node) {
        if (node < 0) return static_cast<int>(ground);
        for (int i = 0; i < n; ++i) {
            if (nodes[static_cast<std::size_t>(i)] == node) return i;
        }
        return -1;
    };
    for (const auto& term : evaluation.staticJacobian) {
        const int r = slot(term.equation);
        const int c = slot(term.unknown);
        if (r < 0 || c < 0) return false;
        J[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] += term.value;
    }
    return true;
}

bool assembleB(const gspice::DaeEvaluation& evaluation, const std::vector<int>& nodes,
               std::vector<double>& b) {
    const int n = static_cast<int>(nodes.size());
    const std::size_t ground = static_cast<std::size_t>(n);
    b.assign(static_cast<std::size_t>(n) + 1, 0.0);
    const auto slot = [&](int node) {
        if (node < 0) return static_cast<int>(ground);
        for (int i = 0; i < n; ++i) {
            if (nodes[static_cast<std::size_t>(i)] == node) return i;
        }
        return -1;
    };
    for (const auto& term : evaluation.staticResidual) {
        const int r = slot(term.equation);
        if (r < 0) return false;
        b[static_cast<std::size_t>(r)] += term.value;
    }
    return true;
}

bool near(double a, double b) { return std::abs(a - b) < 1e-12; }

}  // namespace

int main() {
    constexpr double g01 = 0.5;
    constexpr double g02 = 0.25;
    constexpr double g12 = 0.75;
    constexpr double cap = 1e-12;
    const std::vector<int> nodes = {5, 7};

    // --- GsdiCollapseMap API: standard() from the descriptor ---
    const auto descriptor_pair = makeDescriptor(0);  // Ci collapses into A
    const gspice::GsdiCollapseMap standard_pair = gspice::GsdiCollapseMap::standard(descriptor_pair);
    assert(standard_pair.localCount() == 3);
    assert(standard_pair.unknownCount() == 2);
    assert(standard_pair.collapsedIndex(0) == 0);
    assert(standard_pair.collapsedIndex(1) == 1);
    assert(standard_pair.collapsedIndex(2) == gspice::GsdiNoCollapse);
    assert(standard_pair.localNodeOf(0) == 0);
    assert(standard_pair.localNodeOf(1) == 1);
    assert(standard_pair.valueSource(0) == 0);
    assert(standard_pair.valueSource(1) == 1);
    assert(standard_pair.valueSource(2) == 0);  // value read from terminal A's column

    // Chain resolution: node 2 -> 1, node 3 -> 2 -> 1.
    {
        gspice::GsdiModelDescriptor d;
        d.nodes = {
            {"N0", 0, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
            {"N1", 1, gspice::GsdiNodeRole::Terminal, gspice::GsdiNoCollapse},
            {"N2", 2, gspice::GsdiNodeRole::Collapsible, 1},
            {"N3", 3, gspice::GsdiNodeRole::Collapsible, 2},
        };
        const gspice::GsdiCollapseMap chain(
            d, {gspice::GsdiNoCollapse, gspice::GsdiNoCollapse, 1, 2});
        assert(chain.unknownCount() == 2);
        assert(chain.valueSource(2) == 1);
        assert(chain.valueSource(3) == 1);

        // Chain terminating at the reference node.
        const gspice::GsdiCollapseMap grounded(
            d, {gspice::GsdiNoCollapse, gspice::GsdiNoCollapse,
                gspice::GsdiCollapseToGround, 2});
        assert(grounded.unknownCount() == 2);
        assert(grounded.valueSource(2) == gspice::GsdiNoCollapse);
        assert(grounded.valueSource(3) == gspice::GsdiNoCollapse);

        // Self-referential decision (invalid input) must not loop forever.
        const gspice::GsdiCollapseMap self_loop(
            d, {gspice::GsdiNoCollapse, gspice::GsdiNoCollapse, 2, 3});
        assert(self_loop.valueSource(2) == gspice::GsdiNoCollapse);
        assert(self_loop.valueSource(3) == gspice::GsdiNoCollapse);
    }

    // --- 1. Pair collapse through GdiDevice ---
    auto instance = std::make_unique<CoupledFake>(g01, g02, g12, cap);
    auto device = makeDevice(std::move(instance), nodes, &standard_pair);
    const std::vector<double> x_vals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.05};
    const gspice::DaeEvaluation eval = run(device, x_vals);

    // The instance must have seen V(Ci) == V(A) (collapse resolves the value to
    // the kept terminal A's column).
    {
        auto probe = std::make_unique<CoupledFake>(g01, g02, g12, cap);
        CoupledFake* raw = probe.get();
        auto probe_device = makeDevice(std::move(probe), nodes, &standard_pair);
        const gspice::DaeEvaluation probe_eval = run(probe_device, x_vals);
        (void)probe_eval;
        const auto& s = raw->lastSolution();
        assert(s.size() == 3);
        assert(near(s[0], 0.3));
        assert(near(s[1], 0.05));
        assert(near(s[2], 0.3));  // collapsed V(Ci) reads terminal A's column
    }

    std::vector<std::vector<double>> J;
    std::vector<double> b;
    assert(assembleJ(eval, nodes, J));
    assert(assembleB(eval, nodes, b));
    // Shorted supernode A+Ci vs B: single conductance g01+g12.
    assert(near(J[0][0], g01 + g12));
    assert(near(J[0][1], -(g01 + g12)));
    assert(near(J[1][0], -(g01 + g12)));
    assert(near(J[1][1], g01 + g12));
    // Residual: I(A+Ci) = (g01+g12)*(vA-vB), I(B) = -I(A+Ci).
    const double expected_i = (g01 + g12) * (0.3 - 0.05);
    assert(near(b[0], expected_i));
    assert(near(b[1], -expected_i));

    // Dynamic stamps collapse too: the cap shunts the shorted pair (net zero).
    double net_dyn_J = 0.0;
    for (const auto& term : eval.dynamicJacobian) {
        const int r = (term.equation == 5) ? 0 : (term.equation == 7) ? 1 : -1;
        const int c = (term.unknown == 5) ? 0 : (term.unknown == 7) ? 1 : -1;
        assert(r >= 0 && c >= 0);
        net_dyn_J += term.value;
    }
    assert(near(net_dyn_J, 0.0));
    for (const auto& term : eval.dynamicResidual) {
        const int r = (term.equation == 5) ? 0 : (term.equation == 7) ? 1 : -1;
        assert(r >= 0);
        assert(near(term.value, 0.0));
    }

    // --- 2. Without a collapse plan: local index == column, all three unknown ---
    const std::vector<int> full_nodes = {5, 7, 9};
    auto full_instance = std::make_unique<CoupledFake>(g01, g02, g12, cap);
    auto full_device = makeDevice(std::move(full_instance), full_nodes, nullptr);
    const std::vector<double> full_x = {0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.05, 0.0, 0.1};
    const gspice::DaeEvaluation full_eval = run(full_device, full_x);
    std::vector<std::vector<double>> Jf;
    std::vector<double> bf;
    assert(assembleJ(full_eval, full_nodes, Jf));
    assert(assembleB(full_eval, full_nodes, bf));
    assert(near(Jf[0][0], g01 + g02) && near(Jf[0][1], -g01) && near(Jf[0][2], -g02));
    assert(near(Jf[1][0], -g01) && near(Jf[1][1], g01 + g12) && near(Jf[1][2], -g12));
    assert(near(Jf[2][0], -g02) && near(Jf[2][1], -g12) && near(Jf[2][2], g02 + g12));

    // --- 3. Projection invariant: collapsed Jacobian == shorted projection of ---
    //        the full Jacobian. src(local) maps every local node to the kept
    //        column that carries it (A->A, B->B, Ci->A).
    {
        const int src[3] = {0, 1, 0};
        double proj[2][2] = {};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) proj[src[i]][src[j]] += Jf[i][j];
        }
        assert(near(proj[0][0], J[0][0]) && near(proj[0][1], J[0][1]));
        assert(near(proj[1][0], J[1][0]) && near(proj[1][1], J[1][1]));
    }

    // --- 4. Collapse to ground: Ci's value and KCL equation vanish ---
    const auto descriptor_ground = makeDescriptor(gspice::GsdiCollapseToGround);
    const gspice::GsdiCollapseMap ground_map =
        gspice::GsdiCollapseMap::standard(descriptor_ground);
    assert(ground_map.unknownCount() == 2);
    assert(ground_map.valueSource(2) == gspice::GsdiNoCollapse);
    auto ground_instance = std::make_unique<CoupledFake>(g01, g02, g12, cap);
    auto ground_device = makeDevice(std::move(ground_instance), nodes, &ground_map);
    const gspice::DaeEvaluation ground_eval = run(ground_device, x_vals);
    std::vector<std::vector<double>> Jg;
    std::vector<double> bg;
    assert(assembleJ(ground_eval, nodes, Jg));
    assert(assembleB(ground_eval, nodes, bg));
    // No stamp may reference a node outside {5,7} or ground.
    assert(near(Jg[0][0], g01 + g02) && near(Jg[0][1], -g01));
    assert(near(Jg[1][0], -g01) && near(Jg[1][1], g01 + g12));
    // Ci is grounded, so vCi enters the A branch as 0.
    assert(near(bg[0], g01 * (0.3 - 0.05) + g02 * (0.3 - 0.0)));
    assert(near(bg[1], g01 * (0.05 - 0.3) + g12 * (0.05 - 0.0)));
    // The reference slot carries Ci's row/column so the conservation group
    // closes (identical to how Capacitor emits its grounded terminal).
    assert(near(Jg[0][2], -g02));
    assert(near(Jg[2][0], -g02));
    assert(near(Jg[1][2], -g12));
    assert(near(Jg[2][1], -g12));
    assert(near(Jg[2][2], g02 + g12));
    for (int row = 0; row < 3; ++row) {
        double row_sum = 0.0;
        for (int col = 0; col < 3; ++col) {
            row_sum += Jg[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
        }
        assert(std::abs(row_sum) < 1e-12);
    }
    // Ground-row residual balances the two kept rows.
    assert(near(bg[2], -g02 * 0.3 - g12 * 0.05));
    {
        auto probe = std::make_unique<CoupledFake>(g01, g02, g12, cap);
        CoupledFake* raw = probe.get();
        auto probe_device = makeDevice(std::move(probe), nodes, &ground_map);
        const gspice::DaeEvaluation probe_eval = run(probe_device, x_vals);
        (void)probe_eval;
        const auto& s = raw->lastSolution();
        assert(s.size() == 3);
        assert(near(s[0], 0.3));
        assert(near(s[1], 0.05));
        assert(near(s[2], 0.0));  // value of the reference node
    }

    // --- 5. Dynamic ground collapse: cap between Ci(ground) and A ---
    {
        int count = 0;
        double aa = 0.0, ag = 0.0, ga = 0.0, gg = 0.0;
        int max_row = 0, max_col = 0;
        for (const auto& term : ground_eval.dynamicJacobian) {
            const int r = (term.equation == 5) ? 0 : (term.equation == 7) ? 1
                          : (term.equation == -1) ? 2 : -1;
            const int c = (term.unknown == 5) ? 0 : (term.unknown == 7) ? 1
                          : (term.unknown == -1) ? 2 : -1;
            assert(r >= 0 && c >= 0);
            if (r == 0 && c == 0) aa += term.value;
            if (r == 0 && c == 2) ag += term.value;
            if (r == 2 && c == 0) ga += term.value;
            if (r == 2 && c == 2) gg += term.value;
            ++count;
            max_row = (r > max_row) ? r : max_row;
            max_col = (c > max_col) ? c : max_col;
        }
        // Ci's self-capacitance stamps at (A,A) and the negative pair moves to
        // the reference row — mirroring Capacitor, ground-referenced terms are
        // emitted so the charge conservation group closes.
        assert(count == 4);
        assert(near(aa, cap));
        assert(near(ag, -cap));
        assert(near(ga, -cap));
        assert(near(gg, cap));
        assert(max_row == 2 && max_col == 2);
        // Per-column charge sums vanish within the single conservation group.
        assert(near(aa + ga, 0.0));
        assert(near(ag + gg, 0.0));
    }

    return 0;
}