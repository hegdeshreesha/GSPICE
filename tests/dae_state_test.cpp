#include "dae.hpp"
#include "dae_audit.hpp"
#include "devices/bjt.hpp"
#include "devices/capacitor.hpp"
#include "devices/diode.hpp"
#include "devices/inductor.hpp"
#include "devices/mosfet.hpp"
#include "devices/osdi_device.hpp"
#include "circuit_topology.hpp"
#include "solvers/klu_solver.hpp"
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

struct TestOsdiModel {
    double marker = 0.0;
};

int testOsdiEvaluations = 0;
int testOsdiModelSetups = 0;

void testOsdiSetupModel(void*, void* rawModel, OsdiSimParas*, OsdiInitInfo*) {
    ++testOsdiModelSetups;
    reinterpret_cast<TestOsdiModel*>(rawModel)->marker = 41.0;
}
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
    descriptor.model_size = sizeof(TestOsdiModel);
    descriptor.setup_model = testOsdiSetupModel;
    descriptor.setup_instance = testOsdiSetupInstance;
    descriptor.eval = testOsdiEval;
    descriptor.load_residual_resist = testOsdiLoadResidual;
    descriptor.num_resistive_jacobian_entries = 1;
    descriptor.write_jacobian_array_resist = testOsdiWriteJacobian;
    return descriptor;
}

struct ConformanceOsdiInstance {
    uint32_t mapping[3]{};
    bool collapsed[1]{};
    double resist[2]{};
    double react[2]{};
    double limitResist[2]{};
    double limitReact[2]{};
    double scale = 0.0;
};

int conformanceInstanceSetups = 0;
bool conformanceCollapseInternal = true;
uint32_t conformanceLastFlags = 0;
uint32_t conformanceReturnFlags = 0;
double conformanceSetupGmin = 0.0;
double conformanceSetupMinr = 0.0;
double conformanceEvalGmin = 0.0;
double conformanceEvalMinr = 0.0;
double conformanceSetupTnom = 0.0;
double conformanceSetupScale = 0.0;
double conformanceSetupReltol = 0.0;
double conformanceInstanceTemperature = 0.0;

double conformanceSimparam(const OsdiSimParas& paras, const char* wanted) {
    if (!paras.names || !paras.vals) return 0.0;
    for (std::size_t i = 0; paras.names[i] != nullptr; ++i) {
        if (std::strcmp(paras.names[i], wanted) == 0) return paras.vals[i];
    }
    return 0.0;
}

void conformanceSetupModel(
    void*,
    void* rawModel,
    OsdiSimParas* simParams,
    OsdiInitInfo*) {
    reinterpret_cast<TestOsdiModel*>(rawModel)->marker = 17.0;
    conformanceSetupGmin = conformanceSimparam(*simParams, "gmin");
    conformanceSetupMinr = conformanceSimparam(*simParams, "minr");
    conformanceSetupTnom = conformanceSimparam(*simParams, "tnom");
    conformanceSetupScale = conformanceSimparam(*simParams, "scale");
    conformanceSetupReltol = conformanceSimparam(*simParams, "reltol");
}

void conformanceSetupInstance(
    void*,
    void* rawInstance,
    void*,
    double temperature,
    uint32_t,
    OsdiSimParas*,
    OsdiInitInfo*) {
    auto* instance = reinterpret_cast<ConformanceOsdiInstance*>(rawInstance);
    conformanceInstanceTemperature = temperature;
    instance->collapsed[0] = conformanceCollapseInternal;
    instance->scale = 2.0 + static_cast<double>(conformanceInstanceSetups++);
}

uint32_t conformanceEval(void*, void* rawInstance, void*, OsdiSimInfo* info) {
    auto* instance = reinterpret_cast<ConformanceOsdiInstance*>(rawInstance);
    conformanceLastFlags = info->flags;
    conformanceEvalGmin = conformanceSimparam(info->paras, "gmin");
    conformanceEvalMinr = conformanceSimparam(info->paras, "minr");
    const double input = info->prev_solve[instance->mapping[0]];
    double limited = input;
    uint32_t result = conformanceReturnFlags;
    if ((info->flags & ENABLE_LIM) != 0 && input > 0.5) {
        limited = 0.5;
        result |= EVAL_RET_FLAG_LIM;
    }
    const double delta = limited - input;
    for (int i = 0; i < 2; ++i) {
        instance->resist[i] = 0.5 * instance->scale * limited;
        instance->react[i] = 0.05 * instance->scale * limited;
        instance->limitResist[i] = 0.5 * instance->scale * delta;
        instance->limitReact[i] = 0.05 * instance->scale * delta;
    }
    return result;
}

void conformanceLoadResist(void* rawInstance, void*, double* destination) {
    auto* instance = reinterpret_cast<ConformanceOsdiInstance*>(rawInstance);
    destination[instance->mapping[0]] += instance->resist[0];
    destination[instance->mapping[2]] += instance->resist[1];
}

void conformanceLoadReact(void* rawInstance, void*, double* destination) {
    auto* instance = reinterpret_cast<ConformanceOsdiInstance*>(rawInstance);
    destination[instance->mapping[0]] += instance->react[0];
    destination[instance->mapping[2]] += instance->react[1];
}

void conformanceLoadLimitResist(void* rawInstance, void*, double* destination) {
    auto* instance = reinterpret_cast<ConformanceOsdiInstance*>(rawInstance);
    destination[instance->mapping[0]] += instance->limitResist[0];
    destination[instance->mapping[2]] += instance->limitResist[1];
}

void conformanceLoadLimitReact(void* rawInstance, void*, double* destination) {
    auto* instance = reinterpret_cast<ConformanceOsdiInstance*>(rawInstance);
    destination[instance->mapping[0]] += instance->limitReact[0];
    destination[instance->mapping[2]] += instance->limitReact[1];
}

void conformanceWriteResistJacobian(void* rawInstance, void*, double* destination) {
    const auto* instance = reinterpret_cast<const ConformanceOsdiInstance*>(rawInstance);
    destination[0] = 0.5 * instance->scale;
    destination[1] = 0.5 * instance->scale;
}

void conformanceWriteReactJacobian(void* rawInstance, void*, double* destination) {
    const auto* instance = reinterpret_cast<const ConformanceOsdiInstance*>(rawInstance);
    destination[0] = 0.05 * instance->scale;
    destination[1] = 0.05 * instance->scale;
}

void* conformanceAccess(
    void* rawInstance,
    void*,
    uint32_t id,
    uint32_t flags) {
    if (!rawInstance || id != 0 ||
        (flags & ACCESS_FLAG_INSTANCE) == 0) {
        return nullptr;
    }
    return &reinterpret_cast<ConformanceOsdiInstance*>(rawInstance)->scale;
}

OsdiDescriptor conformanceOsdiDescriptor() {
    static OsdiNode nodes[3]{};
    static OsdiNodePair collapsible[1]{{2, 0}};
    static OsdiJacobianEntry jacobian[2]{};
    static char* opvarNames[1]{const_cast<char*>("gm")};
    static OsdiParamOpvar opvar{};
    for (auto& node : nodes) {
        node.resist_residual_off = UINT32_MAX;
        node.react_residual_off = UINT32_MAX;
    }
    nodes[0].name = const_cast<char*>("p");
    nodes[1].name = const_cast<char*>("n");
    nodes[2].name = const_cast<char*>("internal");
    jacobian[0].nodes = {0, 0};
    jacobian[0].flags = JACOBIAN_ENTRY_RESIST | JACOBIAN_ENTRY_REACT;
    jacobian[0].react_ptr_off = UINT32_MAX;
    jacobian[1].nodes = {2, 0};
    jacobian[1].flags = JACOBIAN_ENTRY_RESIST | JACOBIAN_ENTRY_REACT;
    jacobian[1].react_ptr_off = UINT32_MAX;
    opvar.name = opvarNames;
    opvar.num_alias = 0;
    opvar.description = const_cast<char*>("synthetic transconductance");
    opvar.units = const_cast<char*>("S");
    opvar.flags = PARA_KIND_OPVAR | PARA_TY_REAL;
    opvar.len = 1;

    OsdiDescriptor descriptor{};
    descriptor.name = const_cast<char*>("test_conformance");
    descriptor.num_nodes = 3;
    descriptor.num_terminals = 2;
    descriptor.nodes = nodes;
    descriptor.num_jacobian_entries = 2;
    descriptor.jacobian_entries = jacobian;
    descriptor.num_collapsible = 1;
    descriptor.collapsible = collapsible;
    descriptor.collapsed_offset =
        static_cast<uint32_t>(offsetof(ConformanceOsdiInstance, collapsed));
    descriptor.node_mapping_offset =
        static_cast<uint32_t>(offsetof(ConformanceOsdiInstance, mapping));
    descriptor.jacobian_ptr_resist_offset = UINT32_MAX;
    descriptor.state_idx_off = UINT32_MAX;
    descriptor.bound_step_offset = UINT32_MAX;
    descriptor.instance_size = sizeof(ConformanceOsdiInstance);
    descriptor.model_size = sizeof(TestOsdiModel);
    descriptor.num_opvars = 1;
    descriptor.param_opvar = &opvar;
    descriptor.access = conformanceAccess;
    descriptor.setup_model = conformanceSetupModel;
    descriptor.setup_instance = conformanceSetupInstance;
    descriptor.eval = conformanceEval;
    descriptor.load_residual_resist = conformanceLoadResist;
    descriptor.load_residual_react = conformanceLoadReact;
    descriptor.load_limit_rhs_resist = conformanceLoadLimitResist;
    descriptor.load_limit_rhs_react = conformanceLoadLimitReact;
    descriptor.num_resistive_jacobian_entries = 2;
    descriptor.num_reactive_jacobian_entries = 2;
    descriptor.write_jacobian_array_resist = conformanceWriteResistJacobian;
    descriptor.write_jacobian_array_react = conformanceWriteReactJacobian;
    return descriptor;
}

OsdiDescriptor terminalCollapseOsdiDescriptor() {
    static OsdiNodePair terminalPair[1]{{0, 1}};
    OsdiDescriptor descriptor = conformanceOsdiDescriptor();
    descriptor.name = const_cast<char*>("test_terminal_collapse");
    descriptor.collapsible = terminalPair;
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

void testDaeNewtonLimitingGate() {
    assert(daeNewtonIterationConverged(true, false, false, 0.0, false));
    assert(!daeNewtonIterationConverged(true, false, false, 0.0, true));
    assert(!daeNewtonIterationConverged(true, true, false, 0.0, false));
    assert(!daeNewtonIterationConverged(true, false, true, 1.01, false));
    assert(daeNewtonIterationConverged(true, false, true, 1.0, false));
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

void testOsdiSharedModelLifecycle() {
    using namespace gspice;
    testOsdiModelSetups = 0;
    OSDIDevice first("N1", testOsdiDescriptor(), {0});
    const auto sharedModel = first.sharedModelState();
    assert(sharedModel);
    assert(sharedModel->initialized);
    assert(testOsdiModelSetups == 1);
    assert(close(reinterpret_cast<const TestOsdiModel*>(sharedModel->data())->marker, 41.0));

    OSDIDevice second("N2", testOsdiDescriptor(), {0}, {}, {}, 27.0,
        false, false, false, false, nullptr, sharedModel);
    assert(second.sharedModelState() == sharedModel);
    assert(testOsdiModelSetups == 1);

    std::vector<std::byte> snapshot(first.transientStateBytes());
    first.saveTransientStateBytes(snapshot.data(), snapshot.size());
    reinterpret_cast<TestOsdiModel*>(sharedModel->data())->marker = 99.0;
    first.restoreTransientStateBytes(snapshot.data(), snapshot.size());
    assert(close(reinterpret_cast<const TestOsdiModel*>(sharedModel->data())->marker, 99.0));
}

void testOsdiLimitingCollapseAndFqConformance() {
    using namespace gspice;
    conformanceInstanceSetups = 0;
    conformanceCollapseInternal = true;
    conformanceReturnFlags = 0;
    const OsdiDescriptor descriptor = conformanceOsdiDescriptor();
    OSDIDevice first(
        "NC1", descriptor, {0, -1}, {}, {}, 27.0, true,
        false, false, false, nullptr, nullptr, 7e-9, 3e-8,
        31.0, 2.5, 2e-4, 7e-7, 4e-13, 8e-15, 6e-12);
    assert(close(conformanceSetupGmin, 7e-9));
    assert(close(conformanceSetupMinr, 3e-8));
    assert(close(conformanceSetupTnom, 31.0));
    assert(close(conformanceSetupScale, 2.5));
    assert(close(conformanceSetupReltol, 2e-4));
    assert(close(conformanceInstanceTemperature, 300.15));
    assert(first.metadata().opvar_ids.size() == 1);
    assert(first.metadata().parameters.size() == 1);
    int allocations = 0;
    assert(first.bindInternalUnknowns([&]() { return 1 + allocations++; }) == 1);
    assert(allocations == 1);

    VectorReal firstPoint(1);
    firstPoint[0] = 1.0;
    DaeRequest request;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    request.enableLimiting = false;
    request.simulationGmin = 9e-7;
    DaeEvaluation unlimited;
    assert(first.evaluateDae(firstPoint, request, unlimited));
    assert(close(conformanceEvalGmin, 9e-7));
    assert(close(conformanceEvalMinr, 3e-8));
    assert((conformanceLastFlags & ENABLE_LIM) == 0);
    assert((conformanceLastFlags & INIT_LIM) == 0);

    conformanceReturnFlags = EVAL_RET_FLAG_FINISH | EVAL_RET_FLAG_STOP;
    DaeEvaluation controlRequest;
    assert(first.evaluateDae(firstPoint, request, controlRequest));
    assert(controlRequest.finishRequested);
    assert(controlRequest.stopRequested);
    conformanceReturnFlags = 0;

    request.enableLimiting = true;
    DaeEvaluation limited;
    assert(first.evaluateDae(firstPoint, request, limited));
    assert(limited.limitingApplied);
    assert((conformanceLastFlags & ENABLE_LIM) != 0);
    assert((conformanceLastFlags & INIT_LIM) != 0);
    assert((conformanceLastFlags & CALC_RESIST_LIM_RHS) != 0);
    assert((conformanceLastFlags & CALC_REACT_LIM_RHS) != 0);

    double staticResidual = 0.0;
    double dynamicResidual = 0.0;
    double staticJacobian = 0.0;
    double dynamicJacobian = 0.0;
    for (const auto& term : limited.staticResidual) {
        if (term.equation == 0) staticResidual += term.value;
    }
    for (const auto& term : limited.dynamicResidual) {
        if (term.equation == 0) dynamicResidual += term.value;
    }
    for (const auto& term : limited.staticJacobian) {
        if (term.equation == 0 && term.unknown == 0) staticJacobian += term.value;
    }
    for (const auto& term : limited.dynamicJacobian) {
        if (term.equation == 0 && term.unknown == 0) dynamicJacobian += term.value;
    }
    assert(close(staticResidual, 2.0));
    assert(close(dynamicResidual, 0.2));
    assert(close(staticJacobian, 2.0));
    assert(close(dynamicJacobian, 0.2));

    std::vector<OperatingPointVariable> opvars;
    first.collectOperatingPointVariables(firstPoint, opvars);
    assert(opvars.size() == 1);
    assert(opvars[0].name == "NC1.gm");
    assert(opvars[0].units == "S");
    assert(opvars[0].numericValues.size() == 1);
    assert(close(opvars[0].numericValues[0], 2.0));

    conformanceCollapseInternal = false;
    assert(first.reconfigureInstanceTemperature(57.0));
    assert(close(conformanceInstanceTemperature, 330.15));
    assert(first.topologyRevision() == 1);
    VectorReal expandedPoint(2);
    expandedPoint[0] = 0.25;
    expandedPoint[1] = 0.10;
    DaeEvaluation expandedEvaluation;
    assert(first.evaluateDae(expandedPoint, request, expandedEvaluation));
    bool foundInternalEquation = false;
    for (const auto& term : expandedEvaluation.staticResidual) {
        if (term.equation == 1) foundInternalEquation = true;
    }
    assert(foundInternalEquation);

    DaeEvaluation secondLimitedCall;
    assert(first.evaluateDae(firstPoint, request, secondLimitedCall));
    assert((conformanceLastFlags & INIT_LIM) == 0);

    OSDIDevice second(
        "NC2", descriptor, {1, -1}, {}, {}, 27.0, true,
        false, false, false, nullptr, first.sharedModelState(), 7e-9, 3e-8,
        31.0, 2.5, 2e-4, 7e-7, 4e-13, 8e-15, 6e-12);
    allocations = 0;
    assert(second.bindInternalUnknowns([&]() { return 2 + allocations++; }) == 1);
    assert(allocations == 1);
    VectorReal secondPoint(2);
    secondPoint[1] = 0.25;
    DaeEvaluation secondEvaluation;
    assert(second.evaluateDae(secondPoint, request, secondEvaluation));
    assert(!secondEvaluation.limitingApplied);
    double secondResidual = 0.0;
    for (const auto& term : secondEvaluation.staticResidual) {
        if (term.equation == 1) secondResidual += term.value;
    }
    assert(close(secondResidual, 1.0));

    conformanceCollapseInternal = true;
    OSDIDevice terminalCollapse(
        "NCMERGE", terminalCollapseOsdiDescriptor(), {0, 1}, {}, {}, 27.0, true);
    assert(terminalCollapse.bindInternalUnknowns([]() { return 2; }) == 1);
    std::vector<NodeCollapse> collapses;
    terminalCollapse.collectNodeCollapses(collapses);
    assert(collapses.size() == 1);
    CircuitTopology topology;
    assert(topology.rebuild(2, collapses));
    assert(topology.hasAliases());
    assert(topology.aliasCount() == 1);

    // Two independently connected external nodes become one KCL unknown.
    // The retained alias row enforces equality without changing matrix size.
    SparseMatrixReal mergedMatrix(3);
    VectorReal mergedRhs(3);
    mergedMatrix.add(0, 0, 1.0);
    mergedRhs.add(0, 1.0);
    mergedMatrix.add(1, 1, 1.0);
    mergedRhs.add(1, 3.0);
    mergedMatrix.add(2, 2, 1.0);
    topology.apply(mergedMatrix, mergedRhs);
    const VectorReal mergedSolution = KluSolverReal::solve(mergedMatrix, mergedRhs);
    assert(close(mergedSolution[0], 2.0));
    assert(close(mergedSolution[1], 2.0));
    assert(close(mergedSolution[2], 0.0));

    OSDIDevice terminalGroundCollapse(
        "NCGROUND", terminalCollapseOsdiDescriptor(), {0, -1}, {}, {}, 27.0, true);
    assert(terminalGroundCollapse.bindInternalUnknowns([]() { return 2; }) == 1);
    std::vector<NodeCollapse> groundCollapses;
    terminalGroundCollapse.collectNodeCollapses(groundCollapses);
    assert(groundCollapses.size() == 1);
    CircuitTopology groundTopology;
    assert(groundTopology.rebuild(1, groundCollapses));
    SparseMatrixReal groundedMatrix(2);
    VectorReal groundedRhs(2);
    groundedMatrix.add(0, 0, 2.0);
    groundedRhs.add(0, 5.0);
    groundedMatrix.add(1, 1, 1.0);
    groundTopology.apply(groundedMatrix, groundedRhs);
    const VectorReal groundedSolution =
        KluSolverReal::solve(groundedMatrix, groundedRhs);
    assert(close(groundedSolution[0], 0.0));

    // Re-running setup can split the terminals again. Rebuilding from the
    // original device coordinates restores the uncollapsed system safely.
    conformanceCollapseInternal = false;
    assert(terminalCollapse.reconfigureInstanceTemperature(47.0));
    collapses.clear();
    terminalCollapse.collectNodeCollapses(collapses);
    assert(collapses.empty());
    assert(topology.rebuild(2, collapses));
    assert(!topology.hasAliases());
}

void testSparseMatrixStructureCache() {
    SparseMatrixReal matrix(3);
    matrix.setStructureCacheEnabled(true);
    matrix.add(0, 0, 2.0);
    matrix.add(0, 1, -1.0);
    matrix.add(1, 1, 4.0);
    auto first = matrix.getEntries();
    assert(matrix.structureCacheReady());
    assert(first.size() == 3);
    assert(close(first[0].value, 2.0));

    matrix.clear();
    matrix.add(0, 0, 3.0);
    matrix.add(0, 1, -2.0);
    matrix.add(1, 1, 5.0);
    auto second = matrix.getEntries();
    assert(matrix.structureCacheReady());
    assert(second.size() == 3);
    assert(close(second[0].value, 3.0));
    assert(close(second[1].value, -2.0));
    assert(close(second[2].value, 5.0));

    matrix.clear();
    matrix.add(0, 0, 7.0);
    matrix.add(2, 2, 9.0);
    auto rebuilt = matrix.getEntries();
    assert(matrix.structureCacheReady());
    assert(rebuilt.size() == 4);
    bool found_new_slot = false;
    for (const auto& entry : rebuilt) {
        if (entry.row == 2 && entry.col == 2 && close(entry.value, 9.0)) {
            found_new_slot = true;
        }
    }
    assert(found_new_slot);
}

} // namespace

int main() {
    testSparseMatrixStructureCache();
    testCapacitorDae();
    testDaeNewtonLimitingGate();
    testResistorDae();
    testInductorDae();
    testTransactionalStateStore();
    testOpaqueTransactionalStateStore();
    testNonlinearDaeAudit();
    testVariableStepIntegration();
    testPredictorCorrectorControl();
    testTrapRingingDetection();
    testOsdiDeviceBypass();
    testOsdiSharedModelLifecycle();
    testOsdiLimitingCollapseAndFqConformance();
    return 0;
}
