#include "dae.hpp"
#include "dae_audit.hpp"
#include "devices/bjt.hpp"
#include "devices/capacitor.hpp"
#include "devices/diode.hpp"
#include "devices/inductor.hpp"
#include "devices/mosfet.hpp"
#include "devices/osdi_device.hpp"
#include "devices/resistor.hpp"
#include "integration_formula.hpp"
#include "transient_control.hpp"
#include "transient_state_store.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace {

using namespace gspice;

struct TestOsdiInstance {
    uint32_t mapping[1]{};
    double* resistPointers[1]{};
    double residual = 0.0;
};

int testOsdiEvaluations = 0;

void testOsdiSetupModel(void*, void*, OsdiSimParas*, OsdiInitInfo*) {}
void testOsdiSetupInstance(void*, void*, void*, double, uint32_t, OsdiSimParas*, OsdiInitInfo*) {}
uint32_t testOsdiEval(void*, void* rawInstance, void*, OsdiSimInfo* info) {
    ++testOsdiEvaluations;
    auto* instance = reinterpret_cast<TestOsdiInstance*>(rawInstance);
    instance->residual = info->prev_solve[instance->mapping[0]];
    return 0;
}
void testOsdiLoadResidual(void* rawInstance, void*, double* destination) {
    auto* instance = reinterpret_cast<TestOsdiInstance*>(rawInstance);
    destination[instance->mapping[0]] = instance->residual;
}
void testOsdiWriteJacobian(void*, void*, double* destination) {
    destination[0] = 1.0;
}

OsdiDescriptor testOsdiDescriptor() {
    static OsdiNode node{};
    static OsdiJacobianEntry jacobian{};
    node.name = const_cast<char*>("p");
    node.resist_residual_off = static_cast<uint32_t>(offsetof(TestOsdiInstance, residual));
    node.react_residual_off = UINT32_MAX;
    jacobian.nodes = {0, 0};
    jacobian.react_ptr_off = UINT32_MAX;
    jacobian.flags = JACOBIAN_ENTRY_RESIST;
    OsdiDescriptor descriptor{};
    descriptor.name = const_cast<char*>("test_bypass");
    descriptor.num_nodes = 1;
    descriptor.num_terminals = 1;
    descriptor.nodes = &node;
    descriptor.num_jacobian_entries = 1;
    descriptor.jacobian_entries = &jacobian;
    descriptor.node_mapping_offset = static_cast<uint32_t>(offsetof(TestOsdiInstance, mapping));
    descriptor.jacobian_ptr_resist_offset = static_cast<uint32_t>(offsetof(TestOsdiInstance, resistPointers));
    descriptor.state_idx_off = UINT32_MAX;
    descriptor.bound_step_offset = UINT32_MAX;
    descriptor.instance_size = sizeof(TestOsdiInstance);
    descriptor.model_size = 1;
    descriptor.setup_model = testOsdiSetupModel;
    descriptor.setup_instance = testOsdiSetupInstance;
    descriptor.eval = testOsdiEval;
    descriptor.load_residual_resist = testOsdiLoadResidual;
    descriptor.num_resistive_jacobian_entries = 1;
    descriptor.write_jacobian_array_resist = testOsdiWriteJacobian;
    return descriptor;
}

bool close(double lhs, double rhs, double scale = 1.0) {
    return std::abs(lhs - rhs) <= 1e-12 * std::max(scale, std::max(std::abs(lhs), std::abs(rhs)));
}

void testCapacitorDae() {
    using namespace gspice;
    Capacitor capacitor("C1", 0, 1, 2.0e-12);
    VectorReal x(2);
    x[0] = 1.5;
    x[1] = 0.5;

    DaeRequest request;
    request.analysis = DaeAnalysis::Transient;
    request.staticResidual = false;
    request.staticJacobian = false;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    DaeEvaluation evaluation;
    assert(capacitor.evaluateDae(x, request, evaluation));
    assert(close(daeResidualAt(evaluation.dynamicResidual, 0), 2.0e-12, 2.0e-12));
    assert(close(daeResidualAt(evaluation.dynamicResidual, 1), -2.0e-12, 2.0e-12));
    assert(close(
        daeResidualAt(evaluation.dynamicResidual, 0) +
        daeResidualAt(evaluation.dynamicResidual, 1),
        0.0,
        2.0e-12));

    SparseMatrixReal transientJacobian(2);
    VectorReal transientRhs(2);
    DaeHistory history{{0, -1.0e-3}, {1, 1.0e-3}};
    stampDaeTransient(evaluation, x, 1.0e9, history, transientJacobian, transientRhs);
    const auto dense = transientJacobian.toDense();
    assert(close(dense(0, 0), 2.0e-3));
    assert(close(dense(0, 1), -2.0e-3));
    assert(close(dense(1, 0), -2.0e-3));
    assert(close(dense(1, 1), 2.0e-3));
    assert(close(transientRhs[0], 1.0e-3));
    assert(close(transientRhs[1], -1.0e-3));

    SparseMatrixComplex acJacobian(2);
    stampDaeSmallSignal(evaluation, 5.0e6, acJacobian);
    const auto acDense = acJacobian.toDense();
    assert(close(acDense(0, 0).imag(), 1.0e-5));
    assert(close(acDense(0, 1).imag(), -1.0e-5));
}

void testResistorDae() {
    using namespace gspice;
    Resistor resistor("R1", 0, 1, 1000.0);
    VectorReal x(2);
    x[0] = 1.0;
    x[1] = 0.25;
    DaeEvaluation evaluation;
    DaeRequest request;
    assert(resistor.evaluateDae(x, request, evaluation));
    assert(close(daeResidualAt(evaluation.staticResidual, 0), 7.5e-4));
    assert(close(daeResidualAt(evaluation.staticResidual, 1), -7.5e-4));

    SparseMatrixReal jacobian(2);
    VectorReal rhs(2);
    stampDaeStatic(evaluation, x, jacobian, rhs);
    assert(close(rhs[0], 0.0));
    assert(close(rhs[1], 0.0));
}

void testInductorDae() {
    using namespace gspice;
    Inductor inductor("L1", 0, 1, 4.0e-9, 2);
    VectorReal x(3);
    x[0] = 1.2;
    x[1] = 0.2;
    x[2] = 3.0e-3;

    DaeRequest request;
    request.analysis = DaeAnalysis::Transient;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    DaeEvaluation evaluation;
    assert(inductor.evaluateDae(x, request, evaluation));
    assert(close(daeResidualAt(evaluation.staticResidual, 0), 3.0e-3));
    assert(close(daeResidualAt(evaluation.staticResidual, 1), -3.0e-3));
    assert(close(daeResidualAt(evaluation.staticResidual, 2), 1.0));
    assert(close(daeResidualAt(evaluation.dynamicResidual, 2), -12.0e-12, 12.0e-12));

    SparseMatrixComplex acJacobian(3);
    stampDaeSmallSignal(evaluation, 2.5e8, acJacobian);
    const auto dense = acJacobian.toDense();
    assert(close(dense(2, 2).imag(), -1.0));
    assert(close(dense(2, 0).real(), 1.0));
    assert(close(dense(2, 1).real(), -1.0));

    std::vector<VectorReal> history;
    VectorReal previous(3);
    previous[2] = 2.0e-3;
    history.push_back(previous);
    TransientContext context;
    context.timeStep = 1.0e-9;
    context.currentTime = 1.0e-9;
    context.a0 = 1.0e9;
    context.a1 = -1.0e9;
    context.xHistory = &history;
    SparseMatrixReal transientJacobian(3);
    VectorReal transientRhs(3);
    inductor.tranStamp(transientJacobian, transientRhs, x, context);
    const auto transientDense = transientJacobian.toDense();
    assert(close(transientDense(2, 2), -4.0));
    assert(close(transientRhs[2], -8.0e-3));
}

void testTransactionalStateStore() {
    using namespace gspice;
    TransientStateLayout layout;
    const auto deviceA = layout.allocate(2);
    const auto deviceB = layout.allocate(1);
    layout.seal();

    TransientStateStore store(layout, 3, 2);
    auto initialA = store.initial(deviceA);
    initialA[0] = 1.0;
    initialA[1] = 2.0;
    store.initial(deviceB)[0] = 9.0;

    const auto baseline = store.checkpoint();
    store.prepareCandidate();
    auto firstA = store.candidate(deviceA);
    firstA[0] = 3.0;
    firstA[1] = 4.0;
    store.acceptCandidate();
    assert(close(store.accepted(deviceA)[0], 3.0));
    assert(close(store.accepted(deviceA, 1)[0], 1.0));
    assert(close(store.accepted(deviceB)[0], 9.0));

    store.prepareCandidate();
    store.candidate(deviceA)[0] = 5.0;
    store.acceptCandidate();
    assert(close(store.accepted(deviceA)[0], 5.0));
    store.rollback(baseline);
    assert(close(store.accepted(deviceA)[0], 1.0));
    assert(close(store.accepted(deviceA)[1], 2.0));
    assert(store.acceptedCount() == 0);

    bool rejectedUnpreparedCandidate = false;
    try {
        (void)store.candidate(deviceA);
    } catch (const std::logic_error&) {
        rejectedUnpreparedCandidate = true;
    }
    assert(rejectedUnpreparedCandidate);
}

void testOpaqueTransactionalStateStore() {
    using namespace gspice;
    OpaqueTransientStateLayout layout;
    const auto first = layout.allocate(3);
    const auto second = layout.allocate(1);
    layout.seal();
    OpaqueTransientStateStore store(layout, 3, 2);
    store.initial(first)[0] = std::byte{0x11};
    store.initial(first)[1] = std::byte{0x22};
    store.initial(second)[0] = std::byte{0x7f};

    const auto baseline = store.checkpoint();
    store.prepareCandidate();
    store.candidate(first)[0] = std::byte{0x33};
    store.acceptCandidate();
    assert(store.accepted(first)[0] == std::byte{0x33});
    assert(store.accepted(first, 1)[0] == std::byte{0x11});
    assert(store.accepted(second)[0] == std::byte{0x7f});

    store.rollback(baseline);
    assert(store.accepted(first)[1] == std::byte{0x22});
    assert(store.accepted(second)[0] == std::byte{0x7f});
}

void testNonlinearDaeAudit() {
    using namespace gspice;
    DaeAuditOptions options;
    options.relativeTolerance = 1e-3;

    Diode diode("D1", 0, 1, 1e-14, 1.0, 2e-12);
    VectorReal diodePoint(2);
    diodePoint[0] = 0.35;
    const auto diodeReport = auditDaeDevice(diode, diodePoint, options);
    assert(diodeReport.passed());

    Bjt bjt("Q1", 0, 1, 2, 1, 1e-16, 100.0, 1.0, 1.0, 1.0,
            1.0, 2e-12, 1e-12, 1e-10);
    VectorReal bjtPoint(3);
    bjtPoint[0] = 1.0;
    bjtPoint[1] = 0.65;
    const auto bjtReport = auditDaeDevice(bjt, bjtPoint, options);
    assert(bjtReport.passed());

    Mosfet mos("M1", 0, 1, 2, 3, 1, 2e-6, 1e-6, 0.45, 120e-6);
    VectorReal mosPoint(4);
    mosPoint[0] = 1.0;
    mosPoint[1] = 1.2;
    const auto mosReport = auditDaeDevice(mos, mosPoint, options);
    assert(mosReport.passed());
}

void testVariableStepIntegration() {
    using namespace gspice;
    const auto be = makeBdfFormula({1.0, 0.75});
    assert(be.order == 1);
    assert(close(be.qWeights[0], 4.0));
    assert(close(be.qWeights[1], -4.0));

    const auto bdf2 = makeBdfFormula({1.0, 0.75, 0.5});
    assert(bdf2.order == 2);
    assert(close(bdf2.qWeights[0], 6.0));
    assert(close(bdf2.qWeights[1], -8.0));
    assert(close(bdf2.qWeights[2], 2.0));
    // A second-order formula differentiates a quadratic exactly.
    assert(close(bdf2.differentiate(1.0, {0.75 * 0.75, 0.25}), 2.0));

    const auto variableBdf2 = makeBdfFormula({1.0, 0.8, 0.5});
    assert(close(variableBdf2.differentiate(1.0, {0.64, 0.25}), 2.0));

    const auto trap = makeAdamsMoultonFormula({1.0, 0.75});
    assert(trap.order == 2);
    assert(close(trap.qWeights[0], 8.0));
    assert(close(trap.qWeights[1], -8.0));
    assert(close(trap.derivativeWeights[0], -1.0));
    assert(close(trap.differentiate(1.0, {0.75 * 0.75}, {1.5}), 2.0));

    const auto adams3 = makeAdamsMoultonFormula({1.0, 0.8, 0.5});
    assert(adams3.order == 3);
    // Integrated quadratic derivative history recovers q'=3*t^2 exactly.
    assert(close(adams3.differentiate(1.0, {0.8 * 0.8 * 0.8},
        {3.0 * 0.8 * 0.8, 3.0 * 0.5 * 0.5}), 3.0, 1e-9));
}

void testPredictorCorrectorControl() {
    using namespace gspice;
    std::vector<double> times{0.0, 1.0, 2.0};
    std::vector<VectorReal> values;
    for (double time : times) {
        VectorReal value(1);
        value[0] = time * time;
        values.push_back(value);
    }
    const auto prediction = polynomialPredict(values, times, 3.0, 2);
    assert(prediction.valid);
    assert(close(prediction.value[0], 9.0));

    const auto be = makeBdfFormula({3.0, 2.0});
    assert(close(predictorCorrectorErrorFactor(be, {3.0, 2.0}, times, 3.0), 1.0 / 3.0));

    const auto trap = makeAdamsMoultonFormula({3.0, 2.0});
    assert(close(predictorCorrectorErrorFactor(trap, {3.0, 2.0}, times, 3.0), 1.0 / 13.0));

    const auto bdf2 = makeBdfFormula({3.0, 2.0, 1.0});
    assert(close(predictorCorrectorErrorFactor(bdf2, {3.0, 2.0, 1.0}, times, 3.0), 2.0 / 11.0));

    const auto raise = chooseAdaptiveOrder(2, 5, 0.2, 0.6, 0.01);
    assert(raise.order == 3);
    const auto lower = chooseAdaptiveOrder(3, 5, 0.95, 0.1, std::nullopt);
    assert(lower.order == 2);
}

void testTrapRingingDetection() {
    using namespace gspice;
    std::vector<VectorReal> ringing;
    for (double value : {0.0, 1.0, 0.1, 0.9}) {
        VectorReal point(1);
        point[0] = value;
        ringing.push_back(point);
    }
    assert(detectTrapezoidalRinging(ringing, 1, 1e-9, 1e-6));

    std::vector<VectorReal> smooth;
    for (double value : {0.0, 0.5, 0.75, 0.875}) {
        VectorReal point(1);
        point[0] = value;
        smooth.push_back(point);
    }
    assert(!detectTrapezoidalRinging(smooth, 1, 1e-9, 1e-6));

    AutomaticTransientMethodController controller;
    assert(controller.useTrapezoidal());
    assert(controller.observe(ringing, 1, 1e-9, 1e-6));
    assert(!controller.useTrapezoidal());
}

void testOsdiDeviceBypass() {
    using namespace gspice;
    testOsdiEvaluations = 0;
    OSDIDevice device("NTEST", testOsdiDescriptor(), {0});
    VectorReal point(1);
    point[0] = 1.0;
    DaeRequest request;
    request.allowBypass = true;
    request.bypassRelativeTolerance = 1e-3;
    request.bypassAbsoluteTolerance = 1e-9;
    request.evaluationEpoch = 7;
    DaeEvaluation first;
    assert(device.evaluateDae(point, request, first));
    assert(testOsdiEvaluations == 1);
    point[0] += 1e-5;
    DaeEvaluation second;
    assert(device.evaluateDae(point, request, second));
    assert(second.bypassed);
    assert(testOsdiEvaluations == 1);
    point[0] += 1e-2;
    DaeEvaluation third;
    assert(device.evaluateDae(point, request, third));
    assert(!third.bypassed);
    assert(testOsdiEvaluations == 2);
}

} // namespace

int main() {
    testCapacitorDae();
    testResistorDae();
    testInductorDae();
    testTransactionalStateStore();
    testOpaqueTransactionalStateStore();
    testNonlinearDaeAudit();
    testVariableStepIntegration();
    testPredictorCorrectorControl();
    testTrapRingingDetection();
    testOsdiDeviceBypass();
    return 0;
}
