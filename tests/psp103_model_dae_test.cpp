#include "devices/psp103_model.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

int main() {
    gspice::Psp103ParameterSet params = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VTO", 0.45}, {"U0", 500.0}, {"TOX", 2e-9}, {"EPSROX", 3.9}
    });
    gspice::Psp103Mosfet mos("M1", 0, 1, 2, 3, params, {{"W", 1e-6}, {"L", 1e-6}}, 27.0);
    assert(mos.prefersDampedAutoTransient());

    gspice::VectorReal x(4);
    x[0] = 0.05;
    x[1] = 1.2;
    x[2] = 0.0;
    x[3] = 0.0;

    gspice::DaeRequest request;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    gspice::DaeEvaluation evaluation;
    assert(mos.evaluateDae(x, request, evaluation));

    double chargeSum = 0.0;
    for (const auto& term : evaluation.dynamicResidual) chargeSum += term.value;
    assert(std::abs(chargeSum) < 1e-24);
    assert(evaluation.dynamicJacobian.size() == 16);

    double gateCap = 0.0;
    for (const auto& term : evaluation.dynamicJacobian) {
        if (term.equation == 1 && term.unknown == 1) gateCap += term.value;
    }
    // Physical bound: gate capacitance cannot exceed Cox*W*L for the gate stack
    // (porous: TOX=2nm, EPSROX=3.9, W=L=1um -> Cox*WL ~ 1.726e-14 F).
    const double coxWL = 8.8541878128e-12 * 3.9 / 2e-9 * 1e-12;
    assert(gateCap > 0.0 && gateCap < coxWL);

    gspice::Psp103ParameterSet leakageParams = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VTO", 0.45}, {"BET", 1e-3}, {"IGINV", 1e-6}, {"IGOVD", 2e-6}
    });
    gspice::Psp103Mosfet leakage("M2", 0, 1, 2, 3, leakageParams, {{"W", 2e-6}, {"L", 1e-6}}, 27.0);
    gspice::DaeRequest staticRequest;
    staticRequest.staticResidual = true;
    staticRequest.staticJacobian = true;
    gspice::DaeEvaluation staticEval;
    assert(leakage.evaluateDae(x, staticRequest, staticEval));
    double gateCurrent = 0.0;
    double gateConductance = 0.0;
    for (const auto& term : staticEval.staticResidual) {
        if (term.equation == 1) gateCurrent += term.value;
    }
    for (const auto& term : staticEval.staticJacobian) {
        if (term.equation == 1 && term.unknown == 1) gateConductance += term.value;
    }
    assert(gateCurrent > 0.0);
    assert(gateConductance > 0.0);

    gspice::VectorReal limited(4);
    for (int i = 0; i < 4; ++i) limited[i] = x[i] + 10.0;
    mos.limitTransientNewton(x, limited);
    for (int i = 0; i < 4; ++i) {
        assert(std::abs(limited[i] - x[i]) <= 0.55);
    }

    gspice::Psp103ParameterSet gidlParams = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VTO", 0.45}, {"BET", 1e-3}, {"AGIDLD", 1e-3}, {"BGIDLD", 0.1}
    });
    gspice::Psp103Mosfet gidl("M3", 0, 1, 2, 3, gidlParams, {{"W", 2e-6}, {"L", 1e-6}}, 27.0);
    gspice::VectorReal xGidl(4);
    xGidl[0] = 1.8;
    xGidl[1] = 0.0;
    xGidl[2] = 0.0;
    xGidl[3] = 0.0;
    gspice::DaeEvaluation gidlEval;
    assert(gidl.evaluateDae(xGidl, staticRequest, gidlEval));
    double drainCurrent = 0.0;
    double bulkCurrent = 0.0;
    for (const auto& term : gidlEval.staticResidual) {
        if (term.equation == 0) drainCurrent += term.value;
        if (term.equation == 3) bulkCurrent += term.value;
    }
    assert(drainCurrent > 0.0);
    assert(bulkCurrent < 0.0);
    assert(std::abs(drainCurrent + bulkCurrent) > 0.0);

    gspice::Psp103ParameterSet junctionParams = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VTO", 0.45}, {"BET", 1e-3},
        {"CJORBOT", 1e-15}, {"CJORBOTD", 2e-15}, {"PBRBOT", 0.8}, {"PBRBOTD", 0.8},
        {"IDSATRBOT", 1e-18}, {"IDSATRBOTD", 2e-18}
    });
    gspice::Psp103Mosfet junction("M4", 0, 1, 2, 3, junctionParams, {{"W", 2e-6}, {"L", 1e-6}}, 27.0);
    gspice::VectorReal xJunc(4);
    xJunc[0] = 0.45;
    xJunc[1] = 0.0;
    xJunc[2] = 0.0;
    xJunc[3] = 0.0;
    gspice::DaeEvaluation juncDynEval;
    assert(junction.evaluateDae(xJunc, request, juncDynEval));
    double drainBulkCap = 0.0;
    for (const auto& term : juncDynEval.dynamicJacobian) {
        if (term.equation == 0 && term.unknown == 0) drainBulkCap += term.value;
    }
    assert(drainBulkCap > 0.0);

    gspice::DaeEvaluation juncStaticEval;
    assert(junction.evaluateDae(xJunc, staticRequest, juncStaticEval));
    double junctionDrainCurrent = 0.0;
    double junctionBulkCurrent = 0.0;
    for (const auto& term : juncStaticEval.staticResidual) {
        if (term.equation == 0) junctionDrainCurrent += term.value;
        if (term.equation == 3) junctionBulkCurrent += term.value;
    }
    assert(junctionDrainCurrent > 0.0);
    assert(junctionBulkCurrent < 0.0);
    return 0;
}
