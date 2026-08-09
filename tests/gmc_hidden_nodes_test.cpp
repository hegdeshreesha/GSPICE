#include "devices/gdi_device.hpp"
#include "gmc_hidden.hpp"

#include <cassert>
#include <cmath>

int main() {
    gspice::GmcPspLikeModel model;

    // GSDI-native descriptor: first 4 ports are external terminals, Td is an
    // internal node, Tint is collapsible into s (local index 2).
    const auto& descriptor = model.descriptor();
    assert(descriptor.model_type == "psp_like");
    assert(descriptor.terminal_count == 4);
    assert(descriptor.nodes.size() == 6);
    assert(descriptor.nodes[0].name == "d" &&
           descriptor.nodes[0].role == gspice::GsdiNodeRole::Terminal);
    assert(descriptor.nodes[1].name == "g" &&
           descriptor.nodes[1].role == gspice::GsdiNodeRole::Terminal);
    assert(descriptor.nodes[2].name == "s" &&
           descriptor.nodes[2].role == gspice::GsdiNodeRole::Terminal);
    assert(descriptor.nodes[3].name == "b" &&
           descriptor.nodes[3].role == gspice::GsdiNodeRole::Terminal);
    assert(descriptor.nodes[4].name == "Td" &&
           descriptor.nodes[4].role == gspice::GsdiNodeRole::Internal);
    assert(descriptor.nodes[4].collapse_partner == gspice::GsdiNoCollapse);
    assert(descriptor.nodes[5].name == "Tint" &&
           descriptor.nodes[5].role == gspice::GsdiNodeRole::Collapsible);
    assert(descriptor.nodes[5].collapse_partner == 2);
    assert(descriptor.jacobian_pattern.size() == 36);  // dense 6x6
    assert(descriptor.supports_transient);  // ddt(Cj * (V(Tint) - V(s)))
    assert(descriptor.parameters.size() == 3);
    assert(descriptor.parameters[0].name == "Rd" &&
           descriptor.parameters[0].default_value == 100.0);
    assert(descriptor.parameters[1].name == "KP" &&
           descriptor.parameters[1].default_value == 2e-3);
    assert(descriptor.parameters[2].name == "Cj" &&
           descriptor.parameters[2].default_value == 3e-12);

    // Standard collapse plan derived from the descriptor.
    auto collapse = gspice::GsdiCollapseMap::standard(descriptor);
    assert(collapse.localCount() == 6);
    assert(collapse.unknownCount() == 5);
    assert(collapse.collapsedIndex(0) == 0);
    assert(collapse.collapsedIndex(4) == 4);
    assert(collapse.collapsedIndex(5) == gspice::GsdiNoCollapse);
    assert(collapse.valueSource(0) == 0);
    assert(collapse.valueSource(4) == 4);
    assert(collapse.valueSource(5) == 2);

    // Direct GsdiInstance evaluation over the full local node space.
    gspice::GsdiModelCard card;
    card.name = "PSP1";
    card.type = "psp_like";
    auto instance = model.createInstance(card, {10, 11, 12, 13});
    assert(instance != nullptr);
    assert(instance->terminalCount() == 4);
    assert(instance->internalNodeCount() == 1);

    // Local voltages: d=1.5, g=0.8, s=0.5, b=0.1, Td=1.0, Tint=V(s)=0.5.
    const double voltages[] = {1.5, 0.8, 0.5, 0.1, 1.0, 0.5};
    gspice::GsdiEvalRequest request;
    request.solution = voltages;
    request.solution_size = 6;
    request.residual = true;
    request.jacobian = true;
    request.dynamic_residual = true;
    request.dynamic_jacobian = true;
    gspice::GsdiEvalResult result;
    assert(instance->evaluate(request, result));
    assert(result.finite());

    // Static residual: row0 = I(d,Td) = 5m, row2 = -vch = -0.6m,
    // row4 = -5m + vch = -4.4m. Rows 1, 3, 5 are zero.
    assert(result.static_residual.size() == 6);
    double r0 = 0.0, r2 = 0.0, r4 = 0.0;
    for (const auto& res : result.static_residual) {
        if (res.equation == 0) r0 += res.value;
        if (res.equation == 2) r2 += res.value;
        if (res.equation == 4) r4 += res.value;
    }
    assert(std::abs(r0 - 5e-3) < 1e-15);
    assert(std::abs(r2 + 6e-4) < 1e-15);
    assert(std::abs(r4 + 4.4e-3) < 1e-15);

    // Static jacobian (independent dual per local slot).
    assert(result.static_jacobian.size() == 36);
    double j00 = 0.0, j04 = 0.0, j41 = 0.0, j42 = 0.0, j22 = 0.0;
    for (const auto& elem : result.static_jacobian) {
        if (elem.equation == 0 && elem.unknown == 0) j00 += elem.value;
        if (elem.equation == 0 && elem.unknown == 4) j04 += elem.value;
        if (elem.equation == 4 && elem.unknown == 1) j41 += elem.value;
        if (elem.equation == 4 && elem.unknown == 2) j42 += elem.value;
        if (elem.equation == 2 && elem.unknown == 2) j22 += elem.value;
    }
    assert(std::abs(j00 - 0.01) < 1e-12);
    assert(std::abs(j04 + 0.01) < 1e-12);
    assert(std::abs(j41 - 2e-3) < 1e-12);
    assert(std::abs(j42 + 2e-3) < 1e-12);
    assert(std::abs(j22 - 2e-3) < 1e-12);

    // Dynamic residual: Q(Tint) = Q(s) = 0 at the collapse point V(Tint)=V(s).
    assert(result.dynamic_residual.size() == 6);
    for (const auto& res : result.dynamic_residual) {
        assert(std::abs(res.value) < 1e-30);
    }

    // Dynamic jacobian (dual of the charge) carries Cj entries on the
    // collapsible pair before the global aliasing is applied.
    assert(result.dynamic_jacobian.size() == 36);
    double q55 = 0.0, q52 = 0.0, q22dyn = 0.0;
    for (const auto& elem : result.dynamic_jacobian) {
        if (elem.equation == 5 && elem.unknown == 5) q55 += elem.value;
        if (elem.equation == 5 && elem.unknown == 2) q52 += elem.value;
        if (elem.equation == 2 && elem.unknown == 2) q22dyn += elem.value;
    }
    assert(std::abs(q55 - 3e-12) < 1e-21);
    assert(std::abs(q52 + 3e-12) < 1e-21);
    assert(std::abs(q22dyn - 3e-12) < 1e-21);

    // Stamp through GdiDevice with the collapse plan. nodes_ = kept columns:
    // d=10, g=11, s=12, b=13, Td=14. Tint's row/col alias into s's column.
    gspice::GdiDevice device("D1", model.createInstance(card, {10, 11, 12, 13}),
                             {10, 11, 12, 13, 14}, collapse);
    gspice::VectorReal x(15);
    x[10] = 1.5;
    x[11] = 0.8;
    x[12] = 0.5;
    x[13] = 0.1;
    x[14] = 1.0;
    gspice::DaeRequest dae;
    dae.staticResidual = true;
    dae.staticJacobian = true;
    dae.dynamicResidual = true;
    dae.dynamicJacobian = true;
    gspice::DaeEvaluation evaluation;
    assert(device.evaluateDae(x, dae, evaluation));
    assert(evaluation.staticJacobian.size() == 36);
    assert(evaluation.dynamicJacobian.size() == 36);

    double sr10 = 0.0, sr12 = 0.0, sr14 = 0.0;
    for (const auto& res : evaluation.staticResidual) {
        if (res.equation == 10) sr10 += res.value;
        if (res.equation == 12) sr12 += res.value;
        if (res.equation == 14) sr14 += res.value;
    }
    assert(std::abs(sr10 - 5e-3) < 1e-15);
    assert(std::abs(sr12 + 6e-4) < 1e-15);
    assert(std::abs(sr14 + 4.4e-3) < 1e-15);

    // Collapsed static jacobian: Tint's columns alias into s (column 12), so
    // (12,12) = +2m from vch, (14,12) = -2m, (12,11) = -2m, (14,11) = +2m.
    double g1010 = 0.0, g1014 = 0.0, g1410 = 0.0, g1414 = 0.0;
    double g1411 = 0.0, g1412 = 0.0, g1211 = 0.0, g1212 = 0.0;
    for (const auto& elem : evaluation.staticJacobian) {
        if (elem.equation == 10 && elem.unknown == 10) g1010 += elem.value;
        if (elem.equation == 10 && elem.unknown == 14) g1014 += elem.value;
        if (elem.equation == 14 && elem.unknown == 10) g1410 += elem.value;
        if (elem.equation == 14 && elem.unknown == 14) g1414 += elem.value;
        if (elem.equation == 14 && elem.unknown == 11) g1411 += elem.value;
        if (elem.equation == 14 && elem.unknown == 12) g1412 += elem.value;
        if (elem.equation == 12 && elem.unknown == 11) g1211 += elem.value;
        if (elem.equation == 12 && elem.unknown == 12) g1212 += elem.value;
    }
    assert(std::abs(g1010 - 0.01) < 1e-12);
    assert(std::abs(g1014 + 0.01) < 1e-12);
    assert(std::abs(g1410 + 0.01) < 1e-12);
    assert(std::abs(g1414 - 0.01) < 1e-12);
    assert(std::abs(g1411 - 2e-3) < 1e-12);
    assert(std::abs(g1412 + 2e-3) < 1e-12);
    assert(std::abs(g1211 + 2e-3) < 1e-12);
    assert(std::abs(g1212 - 2e-3) < 1e-12);

    // Chain-rule consistency of the collapse: the four Cj entries of the
    // Tint/s pair alias into (12,12) and must cancel exactly, since the
    // collapsed DAE has V(Tint) == V(s) and the charge vanishes.
    double q1212 = 0.0;
    for (const auto& elem : evaluation.dynamicJacobian) {
        if (elem.equation == 12 && elem.unknown == 12) q1212 += elem.value;
    }
    assert(std::abs(q1212) < 1e-21);

    double dr12 = 0.0;
    for (const auto& res : evaluation.dynamicResidual) {
        if (res.equation == 12) dr12 += res.value;
    }
    assert(std::abs(dr12) < 1e-30);

    // No stamp may target a merged-away column: all equations/unknowns are in
    // nodes_ = {10..14}.
    for (const auto& elem : evaluation.staticJacobian) {
        assert(elem.equation >= 10 && elem.equation <= 14);
        assert(elem.unknown >= 10 && elem.unknown <= 14);
    }
    for (const auto& elem : evaluation.dynamicJacobian) {
        assert(elem.equation >= 10 && elem.equation <= 14);
        assert(elem.unknown >= 10 && elem.unknown <= 14);
    }
    return 0;
}
