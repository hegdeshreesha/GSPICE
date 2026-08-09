#pragma once

#include "../device.hpp"
#include "../dae.hpp"
#include "bsim4_dc_core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace gspice {

class Bsim4Mosfet final : public Device {
public:
    Bsim4Mosfet(const std::string& name, int nodeD, int nodeG, int nodeS, int nodeB,
                Bsim4PreparedModel model, int type = 1)
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
            const auto currents = bsim4DcEvaluateDual(model_, v, type_);
            if (request.staticResidual)
                for (int row = 0; row < 4; ++row) evaluation.staticResidual.push_back({nodes[row], currents[row].value});
            if (request.staticJacobian)
                for (int row = 0; row < 4; ++row)
                    for (int column = 0; column < 4; ++column)
                        evaluation.staticJacobian.push_back({nodes[row], nodes[column], currents[row].derivative[column]});
        }
        if (request.dynamicResidual || request.dynamicJacobian) {
            const auto charges = terminalCharges(v);
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
        const auto result = bsim4EvaluateDc(model_, v, type_);
        const double gm = std::abs(result.jacobian[0][1]);
        const double gds = std::max(result.jacobian[0][0], 0.0);
        const double thermal = 4.0 * 1.380649e-23 * model_.temperature_k *
            (gds + (2.0 / 3.0) * gm);
        const double frequency = std::max(std::abs(omega) / (2.0 * 3.141592653589793), 1.0e-30);
        const double current = std::abs(result.current[0]);
        const double flicker = model_.kf > 0.0 ? model_.kf *
            std::pow(std::max(current, 1.0e-30), model_.af) /
            (std::max(model_.cox * model_.weff * model_.leff, 1.0e-30) *
             std::pow(frequency, model_.ef)) : 0.0;
        return thermal + flicker;
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc,
                             std::vector<NoiseSource>& sources) const override {
        const double psd = const_cast<Bsim4Mosfet*>(this)->getNoisePSD(omega, x_dc);
        if (psd > 0.0) sources.push_back({name_ + ".channel", nodeD_, nodeS_, psd});
    }

private:
    std::array<GmcDual4, 4> terminalCharges(const std::array<double, 4>& v) const {
        const auto vd = GmcDual4::variable(v[0], 0);
        const auto vg = GmcDual4::variable(v[1], 1);
        const auto vs = GmcDual4::variable(v[2], 2);
        const auto vb = GmcDual4::variable(v[3], 3);
        const auto vgs = type_ * (vg - vs);
        const auto vds = type_ * (vd - vs);
        const auto vsb = bsim4Nonnegative(type_ * (vs - vb));
        const double vt = 8.617333262e-5 * model_.temperature_k;
        const auto vth = model_.vth0 + model_.k1 *
            (gmcSqrt(bsim4Nonnegative(model_.phi + vsb)) - std::sqrt(std::max(model_.phi, 1.0e-6))) -
            model_.k2 * vsb;
        const auto vov = vgs - vth;
        const double n = std::max(model_.nfactor, 1.0);
        const double mstar = 0.5 + std::atan(model_.minv) / 3.141592653589793;
        const auto softplus = n * vt * gmcLog(1.0 + bsim4LimitedExp(mstar * vov / (n * vt)));
        const auto denom = mstar + n * (model_.cox / std::max(model_.cdep0, 1.0e-30)) *
            bsim4LimitedExp((model_.voff - (1.0 - mstar) * vov) / (n * vt));
        const auto vgtEff = softplus / denom;
        const double coxwl = model_.cox * model_.weff * model_.leff;
        const auto qinv = -coxwl * vgtEff;
        const auto vdsat = (vgtEff + 2.0 * vt) * std::max(2.0 * model_.vsat * model_.leff / std::max(model_.u0, 1.0e-30), 1.0e-15) /
            (vgtEff + 2.0 * vt + std::max(2.0 * model_.vsat * model_.leff / std::max(model_.u0, 1.0e-30), 1.0e-15));
        const auto partition = bsim4Nonnegative(vds) /
            (vdsat + 1.0e-30);
        const auto intrinsicQd = qinv * (0.5 + partition / 6.0);
        const auto intrinsicQb = GmcDual4::constant(0.0);
        const auto intrinsicQg = -qinv;
        const auto qgsOverlap = model_.cgso * model_.weff * vgs;
        const auto qgdOverlap = model_.cgdo * model_.weff * (vg - vd);
        const auto qgbOverlap = model_.cgbo * model_.leff * (vg - vb);
        const auto junctionCharge = [&](const GmcDual4& junctionVoltage) {
            const auto reverseVoltage = junctionVoltage.value < 0.0
                ? -junctionVoltage : GmcDual4::constant(0.0);
            const auto depletion = [&](double capacitance, double builtIn,
                                       double grading) {
                if (capacitance <= 0.0) return GmcDual4::constant(0.0);
                const auto argument = 1.0 + reverseVoltage / builtIn;
                return -capacitance * builtIn / (1.0 - grading) *
                    (gmcPow(argument, 1.0 - grading) - 1.0);
            };
            return model_.weff * model_.leff *
                       depletion(model_.cj, model_.pb, model_.mj) +
                   2.0 * (model_.weff + model_.leff) *
                       depletion(model_.cjsw, model_.pbsw, model_.mjsw);
        };
        const auto qbd = junctionCharge(vb - vd);
        const auto qbs = junctionCharge(vb - vs);
        const auto qg = intrinsicQg + qgsOverlap + qgdOverlap + qgbOverlap;
        const auto qd = intrinsicQd - qgdOverlap - qbd;
        const auto qb = intrinsicQb - qgbOverlap + qbd + qbs;
        const auto qs = -(qg + qd + qb);
        return {qd, qg, qs, qb};
    }

    int nodeD_, nodeG_, nodeS_, nodeB_;
    Bsim4PreparedModel model_;
    int type_;
};

} // namespace gspice
