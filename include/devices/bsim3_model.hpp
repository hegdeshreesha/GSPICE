#pragma once

#include "device.hpp"
#include "dae.hpp"
#include "bsim3_dc_core.hpp"
#include "../gmc_dual.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace gspice {

/**
 * Native C++ BSIM3 (Berkeley Short-Channel IGFET Model) Compact MOSFET Class.
 * Evaluates through the dual-number DC core (bsim3DcEvaluateDual) so OP/DC,
 * transient, AC, and noise share one analytic-derivative path.
 */
class Bsim3Mosfet : public Device {
public:
    // Legacy constructor kept for the standing unit gates; it routes through
    // the same model set used by the GMC factory.
    Bsim3Mosfet(const std::string& name, int nodeD, int nodeG, int nodeS, int nodeB,
                int type = 1,
                double W = 1e-6, double L = 1e-6, double Vth0 = 0.4,
                double Kp = 120e-6, double u0 = 0.05, double vsat = 1e5,
                double k1 = 0.5, double nfactor = 1.0, double tox = 1e-8,
                double kf = 0.0, double af = 1.0, double ef = 1.0)
        : Bsim3Mosfet(name, nodeD, nodeG, nodeS, nodeB,
                      legacyModel(W, L, Vth0, Kp, u0, vsat, k1, nfactor, tox, kf, af, ef),
                      type) {}

    Bsim3Mosfet(const std::string& name, int nodeD, int nodeG, int nodeS, int nodeB,
                Bsim3PreparedModel model, int type = 1)
        : Device(name), nodeD_(nodeD), nodeG_(nodeG), nodeS_(nodeS), nodeB_(nodeB),
          model_(std::move(model)), type_(type) {}

    bool daeAuditSafe() const override { return true; }

    bool evaluateDae(const VectorReal& x, const DaeRequest& request,
                     DaeEvaluation& evaluation) override {
        evaluation.clear();
        const std::array<int, 4> nodes = {nodeD_, nodeG_, nodeS_, nodeB_};
        const auto voltage = [&](int node) { return node >= 0 ? x[node] : 0.0; };
        const std::array<double, 4> v = {voltage(nodeD_), voltage(nodeG_), voltage(nodeS_), voltage(nodeB_)};
        if (request.staticResidual || request.staticJacobian) {
            const auto dc = bsim3EvaluateDc(model_, v, type_);
            if (!dc.valid) return false;
            if (request.staticResidual)
                for (int row = 0; row < 4; ++row) evaluation.staticResidual.push_back({nodes[row], dc.current[row]});
            if (request.staticJacobian)
                for (int row = 0; row < 4; ++row)
                    for (int column = 0; column < 4; ++column)
                        evaluation.staticJacobian.push_back({nodes[row], nodes[column], dc.jacobian[row][column]});
        }
        if (request.dynamicResidual || request.dynamicJacobian) {
            const auto charges = terminalChargesDual(v);
            if (request.dynamicResidual)
                for (int row = 0; row < 4; ++row) evaluation.dynamicResidual.push_back({nodes[row], charges[row].value});
            if (request.dynamicJacobian)
                for (int row = 0; row < 4; ++row)
                    for (int column = 0; column < 4; ++column)
                        evaluation.dynamicJacobian.push_back({nodes[row], nodes[column], charges[row].derivative[column]});
        }
        return model_.validation.valid;
    }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x,
                 double, double, const std::vector<VectorReal>&) override {
        DaeRequest request;
        request.staticResidual = true;
        request.staticJacobian = true;
        DaeEvaluation evaluation;
        evaluateDae(x, request, evaluation);
        stampDaeStatic(evaluation, x, J, b);
    }

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x,
                   const TransientContext& ctx) override {
        DaeRequest request;
        request.analysis = DaeAnalysis::Transient;
        request.time = ctx.currentTime;
        request.staticResidual = true;
        request.staticJacobian = true;
        request.dynamicResidual = ctx.timeStep > 0.0;
        request.dynamicJacobian = ctx.timeStep > 0.0;
        DaeEvaluation current;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) {
            request.dynamicResidual = false;
            request.dynamicJacobian = false;
            evaluateDae(x, request, current);
            stampDaeStatic(current, x, J, b);
            return;
        }
        evaluateDae(x, request, current);
        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }
        DaeEvaluation previous;
        evaluateDae((*ctx.xHistory)[ctx.xHistory->size() - 1], request, previous);
        DaeHistory history;
        appendScaledDaeResidual(history, previous.dynamicResidual, a1);
        if (useSecond) {
            DaeEvaluation previous2;
            evaluateDae((*ctx.xHistory)[ctx.xHistory->size() - 2], request, previous2);
            appendScaledDaeResidual(history, previous2.dynamicResidual, a2);
        }
        stampDaeTransient(current, x, a0, history, J, b);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega,
                 const VectorReal& x_dc) override {
        (void)b;
        DaeRequest request;
        request.analysis = DaeAnalysis::SmallSignal;
        request.staticJacobian = true;
        request.dynamicJacobian = true;
        DaeEvaluation evaluation;
        evaluateDae(x_dc, request, evaluation);
        stampDaeSmallSignal(evaluation, omega, J);
    }

    double getNoisePSD(double omega, const VectorReal& x_dc) override {
        const std::array<double, 4> v = {
            nodeD_ >= 0 ? x_dc[nodeD_] : 0.0,
            nodeG_ >= 0 ? x_dc[nodeG_] : 0.0,
            nodeS_ >= 0 ? x_dc[nodeS_] : 0.0,
            nodeB_ >= 0 ? x_dc[nodeB_] : 0.0};
        const auto result = bsim3EvaluateDc(model_, v, type_);
        const double gm = std::abs(result.jacobian[0][1]);
        const double gds = std::abs(result.jacobian[0][0]);
        // Drain current noise: thermal + flicker.  MOSIS-style channel
        // conductance approximation keeps the OFF state at zero.
        const double k = 1.380649e-23;
        const double temperature = model_.temperature_k;
        const double channelConductance = std::max(gds + (2.0 / 3.0) * gm, 0.0);
        const double thermal = 4.0 * k * temperature * channelConductance;
        const double frequency = std::max(std::abs(omega) / (2.0 * 3.141592653589793), 1.0e-30);
        const double current = std::abs(result.current[0]);
        const double flicker = model_.kf > 0.0
            ? model_.kf * std::pow(std::max(current, 1.0e-30), model_.af) /
                (std::max(3.453133e-11 * model_.weff * model_.leff /
                             std::max(model_.toxe, 1.0e-12),
                          1.0e-30) *
                 std::pow(frequency, model_.ef))
            : 0.0;
        return thermal + flicker;
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc,
                             std::vector<NoiseSource>& sources) const override {
        const std::array<double, 4> v = {
            nodeD_ >= 0 ? x_dc[nodeD_] : 0.0,
            nodeG_ >= 0 ? x_dc[nodeG_] : 0.0,
            nodeS_ >= 0 ? x_dc[nodeS_] : 0.0,
            nodeB_ >= 0 ? x_dc[nodeB_] : 0.0};
        const auto result = bsim3EvaluateDc(model_, v, type_);
        const double k = 1.380649e-23;
        const double temperatureK = model_.temperature_k;
        const double thermal = 4.0 * k * temperatureK *
            std::max(std::abs(result.jacobian[0][0]) +
                         (2.0 / 3.0) * std::abs(result.jacobian[0][1]), 0.0);
        const double frequency = std::max(std::abs(omega) / (2.0 * 3.141592653589793), 1.0e-30);
        const double current = std::abs(result.current[0]);
        const double flicker = model_.kf > 0.0
            ? model_.kf * std::pow(std::max(current, 1.0e-30), model_.af) /
                (std::max(3.453133e-11 * model_.weff * model_.leff /
                             std::max(model_.toxe, 1.0e-12), 1.0e-30) *
                 std::pow(frequency, model_.ef))
            : 0.0;
        const double psd = thermal + flicker;
        if (psd > 0.0)
            sources.push_back({name_ + ".channel", nodeD_, nodeS_, psd});
    }

private:
    static Bsim3PreparedModel legacyModel(double W, double L, double Vth0,
                                          double Kp, double u0, double vsat,
                                          double k1, double nfactor, double tox,
                                          double kf, double af, double ef) {
        std::unordered_map<std::string, double> raw = {
            {"LEVEL", 49.0}, {"VTH0", Vth0}, {"KP", Kp}, {"U0", u0},
            {"VSAT", vsat}, {"K1", k1}, {"NFACTOR", nfactor}, {"TOXE", tox},
            {"KF", kf}, {"AF", af}, {"EF", ef}};
        return Bsim3ParameterSet::from(raw).prepare(W, L);
    }

    std::array<GmcDual4, 4> terminalChargesDual(const std::array<double, 4>& v) const {
        const auto vd = GmcDual4::variable(v[0], 0);
        const auto vg = GmcDual4::variable(v[1], 1);
        const auto vs = GmcDual4::variable(v[2], 2);
        const auto vb = GmcDual4::variable(v[3], 3);
        const auto vgs = type_ * (vg - vs);
        const auto vds = type_ * (vd - vs);
        const auto vsb = bsim3Positive(type_ * (vs - vb));
const double sqrtPhi = std::sqrt(std::max(model_.phi, 1.0e-12));
        const auto vth = model_.vth0 +
            model_.k1 * (gmcSqrt(bsim3Positive(model_.phi + vsb)) - sqrtPhi) -
            model_.k2 * vds;
        const auto vov = gmcMax(GmcDual4::constant(0.0), vgs - vth);
        const double cox = 3.453133e-11 / std::max(model_.toxe, 1.0e-12);
        const double coxwl = cox * model_.weff * model_.leff;
        const auto qinv = -coxwl * vov;  // inversion charge, C
        // Velocity-saturation boundary for the source/drain partition.
        const double esat = 2.0 * model_.vsat / std::max(model_.u0, 1.0e-6);
        const double esatL = esat * std::max(model_.leff, 1.0e-15);
        const auto vdsat = (vov * esatL) / (vov + esatL + 1.0e-30);
        const auto partition = unitInterval(vds / (vdsat + 1.0e-30));
        const auto qg = -qinv;                                  // total gate charge
        const auto qd = qinv * (0.5 + partition / 6.0);         // drain share
        const auto qb = GmcDual4::constant(0.0);                // bulk charge
        const auto qs = -(qg + qd + qb);                        // source share closes the group
        return {qd, qg, qs, qb};
    }

    static GmcDual4 nonnegative(const GmcDual4& value) {
        return value.value > 0.0 ? value : GmcDual4::constant(0.0);
    }
    static GmcDual4 unitInterval(const GmcDual4& value) {
        if (value.value <= 0.0) return GmcDual4::constant(0.0);
        if (value.value >= 1.0) return GmcDual4::constant(1.0);
        return value;
    }

    int nodeD_, nodeG_, nodeS_, nodeB_;
    int type_;
    Bsim3PreparedModel model_;
};

} // namespace gspice
