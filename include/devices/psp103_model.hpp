#ifndef GSPICE_PSP103_MODEL_HPP
#define GSPICE_PSP103_MODEL_HPP

#include "device.hpp"
#include "psp103_core.hpp"
#include "psp103_parameters.hpp"
#include "psp103_temperature.hpp"
#include "../gmc_dual.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace gspice {

class Psp103Mosfet final : public Device {
public:
    // PSP103 charge decomposition (module loadDynamic): intrinsic gate/bulk/
    // drain/source charges plus the gate-source / gate-drain fringe+overlap
    // and the gate-bulk overlap charges. All quantities are after the
    // S/D-interchange undo, in Coulombs per unit CHNL_TYPE*MULT.
    struct ChargeComponents {
        GmcDual4 qg;   // Qg  (intrinsic gate)
        GmcDual4 qb;   // Qb  (intrinsic bulk)
        GmcDual4 qd;   // Qd  (intrinsic drain; undoes S/D-interchange)
        GmcDual4 qs;   // Qs  = -(Qg+Qb+Qd)
        GmcDual4 qfgs; // Qfgs + Qgs_ov (fringe + source overlap)
        GmcDual4 qfgd; // Qfgd + Qgd_ov (fringe + drain overlap)
        GmcDual4 qgb;  // Qgb_ov
        GmcDual4 cox_qm;
    };

    ChargeComponents intrinsicAndOverlapCharges(const VectorReal& x) const {
        const auto vd = GmcDual4::variable(nodeVoltage(x, nodeD_), 0);
        const auto vg = GmcDual4::variable(nodeVoltage(x, nodeG_), 1);
        const auto vs = GmcDual4::variable(nodeVoltage(x, nodeS_), 2);
        const auto vb = GmcDual4::variable(nodeVoltage(x, nodeB_), 3);
        const double type = static_cast<double>(prepared_.polarity);

        auto vgs = type * (vg - vs);
        auto vds = type * (vd - vs);
        auto vsb = type * (vs - vb);
        const auto vgs_prime = vgs;
        const auto vgd_prime = vgs - vds;
        int sigVds = 1;
        if (vds.value < 0.0) {
            sigVds = -1;
            vgs = vgs - vds;
            vsb = vsb + vds;
            vds = -vds;
        }

        Psp103DcCoreInputs in;
        in.phib = setup_.phib_dc;
        in.g0 = setup_.g0_dc;
        in.phit0 = setup_.phit0;
        in.vfb_t = setup_.vfb_t;
        in.kp = setup_.kp;
        in.cf = setup_.cf;
        in.cfd = setup_.cfd;
        in.cfb = setup_.cfb;
        in.psce = setup_.psce;
        in.psce_d = setup_.psce_d;
        in.psce_b = setup_.psce_b;
        in.dnsub = setup_.dnsub;
        in.vnsub = setup_.vnsub;
        in.nslp = setup_.nslp;
        in.rsb = setup_.rsb;
        in.rsg = setup_.rsg;
        in.ther = setup_.ther;
        in.e_eff0 = setup_.e_eff0;
        in.eta_mu = setup_.eta_mu;
        in.eta_mu1 = setup_.eta_mu1;
        in.mue_t = setup_.mue_t;
        in.themu_t = setup_.themu_t;
        in.cs_t = setup_.cs_t;
        in.xcor_t = setup_.xcor_t;
        in.thesat_b = setup_.thesat_b;
        in.thesat_g = setup_.thesat_g;
        in.theta_sat_t = setup_.theta_sat_t;
        in.ax = setup_.ax;
        in.alp = setup_.alp;
        in.alp1 = setup_.alp1;
        in.alp2 = setup_.alp2;
        in.vp = setup_.vp;
        in.bet = setup_.bet_i;
        in.phix = setup_.phix_dc;
        in.aphi = setup_.aphi_dc;
        in.phix1 = setup_.phix1_dc;
        in.vsbnud = setup_.vsbnud;
        in.dvsbnud = setup_.dvsbnud;
        in.gfac_nud = setup_.gfac_nud;
        in.us1 = setup_.us1;
        in.us21 = setup_.us21;
        in.sw_nud = setup_.sw_nud != 0;
        in.pmos = setup_.channel_type == -1;

        const auto core = psp103DcCore(in, vgs, vds, vsb);

        GmcDual4 cox_qm = GmcDual4::constant(setup_.cox_i);
        if (setup_.qq > 0.0) {
            cox_qm = setup_.cox_i /
                (1.0 + setup_.qq * gmcPow(
                    core.qeff1 * core.qeff1 + setup_.qlim2,
                    -psp103_const::oneSixth));
        }

        GmcDual4 qg, qi, qd, qb;
        const double xg_val = core.xg.value;
        const double w_acc = std::clamp((0.001 - xg_val) / 0.002, 0.0, 1.0);
        const double s_acc = w_acc * w_acc * (3.0 - 2.0 * w_acc);

        const auto qg_acc = core.voxm;
        const auto qb_acc = qg_acc;

        const auto fj = 0.5 * (core.dps / core.h);
        const auto fj2 = fj * fj;
        const auto qclm = (1.0 - core.g_delta_l) *
            (core.qim - 0.5 * (core.alpha * core.dps));
        const auto qg_inv = core.voxm + 0.5 * (core.eta_p * core.dps *
            (fj * core.g_delta_l * psp103_const::oneThird - 1.0 + core.g_delta_l));
        const auto ftemp = core.alpha * core.dps * psp103_const::oneSixth;
        const auto qi_inv = core.g_delta_l * (core.qim + ftemp * fj) + qclm;
        const auto qd_inv = 0.5 * (core.g_delta_l * core.g_delta_l *
            (core.qim - ftemp * (1.0 - fj - 0.2 * fj2)) +
            qclm * (1.0 + core.g_delta_l));
        const auto qb_inv = qg_inv - qi_inv;

        qg = (1.0 - s_acc) * qg_inv + s_acc * qg_acc;
        qi = (1.0 - s_acc) * qi_inv;
        qd = (1.0 - s_acc) * qd_inv;
        qb = (1.0 - s_acc) * qb_inv + s_acc * qb_acc;
        ChargeComponents comp;
        comp.cox_qm = cox_qm;
        comp.qg = qg * cox_qm;
        comp.qd = -qd * cox_qm;
        comp.qb = -qb * cox_qm;
        comp.qs = -(comp.qg + comp.qb + comp.qd);
        if (sigVds < 0) {
            const auto temp = comp.qd;
            comp.qd = comp.qs;
            comp.qs = temp;
        }

        GmcDual4 vovs = GmcDual4::constant(0.0);
        GmcDual4 vovd = GmcDual4::constant(0.0);
        if (setup_.cgov_i > 0.0 || setup_.cgovd_i > 0.0) {
            const auto xgs_ov = -vgs_prime * setup_.inv_phit;
            const auto xgd_ov = -vgd_prime * setup_.inv_phit;
            GmcDual4 sp_ov_xg = 0.5 * (xgs_ov + gmcSqrt(xgs_ov * xgs_ov + setup_.spov_eps2_s));
            const auto xs_ov = -sp_ov_xg - setup_.gov2_s * 0.5 +
                setup_.gov_s * gmcSqrt(sp_ov_xg + setup_.gov2_s * 0.25 + setup_.spov_a_s) +
                GmcDual4::constant(setup_.spov_delta1_s);
            sp_ov_xg = 0.5 * (xgd_ov + gmcSqrt(xgd_ov * xgd_ov + setup_.spov_eps2_d));
            const auto xd_ov = -sp_ov_xg - setup_.gov2_d * 0.5 +
                setup_.gov_d * gmcSqrt(sp_ov_xg + setup_.gov2_d * 0.25 + setup_.spov_a_d) +
                GmcDual4::constant(setup_.spov_delta1_d);
            vovs = -setup_.phit * (xgs_ov + xs_ov);
            vovd = -setup_.phit * (xgd_ov + xd_ov);
        }
        const auto qgs_ov = setup_.cgov_i * vovs;
        const auto qgd_ov = setup_.cgovd_i * vovd;
        comp.qgb = setup_.cgbov_i * (vgs + vsb);
        comp.qfgs = setup_.cfr_i * vgs_prime + qgs_ov;
        comp.qfgd = setup_.cfrd_i * vgd_prime + qgd_ov;
        return comp;
    }

    Psp103Mosfet(
        const std::string& name,
        int nodeD,
        int nodeG,
        int nodeS,
        int nodeB,
        const Psp103ParameterSet& model,
        const GsdiParamMap& instance,
        double temperature_c)
        : Device(name),
          nodeD_(nodeD),
          nodeG_(nodeG),
          nodeS_(nodeS),
          nodeB_(nodeB),
          model_(model),
          prepared_(model.prepare(instance, temperature_c)),
          scaled_(psp103ScaleTemperature(model_, prepared_, temperature_c)),
          setup_(psp103PrepareDevice(model_, prepared_, temperature_c)),
          cox_(coxDensity(model)),
          threshold_(model.alias(std::abs(scaled_.vfb) + 0.5 * model.alias(0.7, {"PHIB", "PHIBO"}),
                                 {"VTO", "VT0", "VTH0", "VTH", "VFB", "VFB0"})),
          beta_(std::max(betaDensity(model, scaled_), 1.0e-9) *
                prepared_.width_length_ratio),
          lambda_(std::max(model.alias(0.02, {"LAMBDA", "LAM", "LAMDA", "CLM"}), 0.0)),
          nfactor_(std::max(model.alias(1.4, {"NFACTOR", "N", "N0"}), 1.01)),
          gateInvLeakage_(std::max(scaled_.ig_inv, 0.0) * prepared_.effective_width_m),
          gateSrcOverlapLeakage_(std::max(scaled_.ig_ov, 0.0) * prepared_.effective_width_m),
          gateDrnOverlapLeakage_(std::max(scaled_.ig_ovd, 0.0) * prepared_.effective_width_m),
          sourceGidlScale_(std::max(scaled_.agidl, 0.0) * prepared_.effective_width_m),
          sourceGidlBarrier_(std::max(scaled_.bgidl, 0.0)),
          drainGidlScale_(std::max(scaled_.agidld, 0.0) * prepared_.effective_width_m),
          drainGidlBarrier_(std::max(scaled_.bgidld, 0.0)),
          sourceJunctionCap_(junctionCapacitance(model_, prepared_, false)),
          drainJunctionCap_(junctionCapacitance(model_, prepared_, true)),
          sourceJunctionLeakage_(junctionLeakage(model_, prepared_, false)),
          drainJunctionLeakage_(junctionLeakage(model_, prepared_, true)),
          sourceJunctionVbi_(junctionBuiltIn(model_, false)),
          drainJunctionVbi_(junctionBuiltIn(model_, true)) {}

    bool daeAuditSafe() const override { return true; }

    bool prefersDampedAutoTransient() const override { return true; }

    void limitTransientNewton(
        const VectorReal& previous,
        VectorReal& candidate) const override {
        limitNodeUpdate(previous, candidate, nodeD_);
        limitNodeUpdate(previous, candidate, nodeG_);
        limitNodeUpdate(previous, candidate, nodeS_);
        limitNodeUpdate(previous, candidate, nodeB_);
    }

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        if (!prepared_.valid()) return false;
        const std::array<int, 4> nodes = {nodeD_, nodeG_, nodeS_, nodeB_};
        const auto currents = terminalCurrents(x);
        if (!finiteDualVector(currents)) return false;
        if (request.staticResidual) {
            for (int row = 0; row < 4; ++row) {
                evaluation.staticResidual.push_back({nodes[row], currents[row].value});
            }
        }
        if (request.staticJacobian) {
            for (int row = 0; row < 4; ++row) {
                for (int column = 0; column < 4; ++column) {
                    evaluation.staticJacobian.push_back(
                        {nodes[row], nodes[column], currents[row].derivative[column]});
                }
            }
        }
        if (request.dynamicResidual || request.dynamicJacobian) {
            const auto charges = terminalCharges(x);
            if (!finiteDualVector(charges)) return false;
            if (request.dynamicResidual) {
                for (int row = 0; row < 4; ++row) {
                    evaluation.dynamicResidual.push_back({nodes[row], charges[row].value, 0});
                }
            }
            if (request.dynamicJacobian) {
                for (int row = 0; row < 4; ++row) {
                    for (int column = 0; column < 4; ++column) {
                        evaluation.dynamicJacobian.push_back(
                            {nodes[row], nodes[column], charges[row].derivative[column], 0});
                    }
                }
            }
        }
        return true;
    }

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
        DaeRequest request;
        request.analysis = DaeAnalysis::OperatingPoint;
        DaeEvaluation evaluation;
        if (evaluateDae(x, request, evaluation)) stampDaeStatic(evaluation, x, J, b);
    }

    void tranStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
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
            if (evaluateDae(x, request, current)) stampDaeStatic(current, x, J, b);
            return;
        }
        if (!evaluateDae(x, request, current)) return;

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
        (void)omega;
        DaeRequest request;
        request.analysis = DaeAnalysis::SmallSignal;
        request.staticJacobian = true;
        request.dynamicJacobian = true;
        DaeEvaluation evaluation;
        if (evaluateDae(x_dc, request, evaluation)) stampDaeSmallSignal(evaluation, omega, J);
    }

    // PSP103 noise sources used by .NOISE/.PNOISE: channel thermal, flicker,
    // and shot noise on the gate leakage legs.
    double getNoisePSD(double omega, const VectorReal& x_dc) override {
        return psp103NoisePsd(omega, x_dc);
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc,
                             std::vector<NoiseSource>& sources) const override {
        const auto noise = psp103NoiseComponents(omega, x_dc);
        if (noise.channel_psd > 0.0) {
            sources.push_back({name_ + ".drain", nodeD_, nodeS_, noise.channel_psd});
        }
        if (noise.gate_shot_psd > 0.0 && nodeG_ >= 0 && nodeS_ >= 0) {
            sources.push_back({name_ + ".gate", nodeG_, nodeS_, noise.gate_shot_psd});
        }
    }

    bool probeCurrent(const VectorReal& x, double& current, double time = 0.0) const override {
        (void)time;
        if (!prepared_.valid()) return false;
        const auto currents = terminalCurrents(x);
        if (!finiteDualVector(currents)) return false;
        current = currents[0].value;
        return true;
    }

private:
    struct Psp103NoiseComponents {
        double channel_psd = 0.0;
        double gate_shot_psd = 0.0;
    };

    Psp103NoiseComponents psp103NoiseComponents(double omega, const VectorReal& x_dc) const {
        Psp103NoiseComponents out;
        if (!prepared_.valid() || setup_.cox_i <= 0.0) return out;
        const auto currents = terminalCurrents(x_dc);
        if (!finiteDualVector(currents)) return out;
        const double gm = std::abs(currents[0].derivative[1]);
        const double gds = std::abs(currents[0].derivative[0]);
        const double id = std::abs(currents[0].value);
        const double ig = std::abs(currents[1].value);
        const double k = 1.380649e-23;
        const double temperature = prepared_.temperature.kelvin;
        // Channel thermal: 4kT * (gds + (2/3)gm) covers triode through
        // saturation continuously.
        out.channel_psd = 4.0 * k * temperature * (gds + (2.0 / 3.0) * gm);
        // Flicker: PSP-style SFL/ALPNOI * Id^2 / f^EF, or the MOSIS
        // KF * Id^AF / (f^EF) form when PSP flicker is not given.
        const double frequency = std::max(std::abs(omega) / (2.0 * 3.141592653589793), 1.0e-30);
        const double sfl = std::max(model_.alias(0.0, {"SFL", "SFLN", "ALPNOI"}), 0.0);
        const double kf = std::max(model_.alias(0.0, {"KF", "KF0"}), 0.0);
        const double ef = std::max(model_.alias(1.0, {"EF", "EF0", "EFO", "EFEDGEO"}), 1.0e-30);
        if (sfl > 0.0) {
            out.channel_psd += sfl * id * id / std::pow(frequency, ef);
        } else if (kf > 0.0) {
            const double af = std::max(model_.alias(1.0, {"AF", "AF0"}), 1.0e-30);
            const double cox_wl = 3.453133e-11 * prepared_.effective_width_m *
                                  prepared_.effective_length_m /
                                  std::max(setup_.tox_i, 1.0e-12);
            out.channel_psd += kf * std::pow(std::max(id, 1.0e-30), af) /
                (std::max(cox_wl, 1.0e-30) * std::pow(frequency, ef));
        }
        const double fnt = std::max(model_.sum({"FNT", "FNTO", "FNTEXCL", "FNTEDGEO"}), 0.0);
        if (fnt > 0.0) {
            out.channel_psd += fnt * 4.0 * k * temperature * (gds + gm);
        }
        // Gate leakage shot noise: 2q * |Ig|.
        out.gate_shot_psd = 2.0 * 1.602176634e-19 * ig;
        return out;
    }

    double psp103NoisePsd(double omega, const VectorReal& x_dc) const {
        const auto noise = psp103NoiseComponents(omega, x_dc);
        return noise.channel_psd + noise.gate_shot_psd;
    }

    static void limitNodeUpdate(
        const VectorReal& previous,
        VectorReal& candidate,
        int node) {
        if (node < 0 || node >= previous.getSize() || node >= candidate.getSize()) return;
        const double old_value = previous[node];
        const double delta = candidate[node] - old_value;
        const double bound = std::min(0.25 + 0.25 * std::abs(old_value), 2.0);
        candidate[node] = old_value + std::clamp(delta, -bound, bound);
    }

    static double coxDensity(const Psp103ParameterSet& model) {
        const double tox = model.alias(0.0, {"TOX", "TOXO"});
        const double epsr = model.alias(3.9, {"EPSROX", "EPSROXO"});
        if (tox <= 0.0 || epsr <= 0.0) return 0.0;
        return 8.8541878128e-12 * epsr / tox;
    }

    static double betaDensity(const Psp103ParameterSet& model, const Psp103TemperatureScaled& scaled) {
        const double explicit_beta = model.alias(
            std::numeric_limits<double>::quiet_NaN(), {"BET", "BETN", "BETO", "BETA", "KP"});
        if (std::isfinite(explicit_beta)) return explicit_beta;
        if (scaled.beta > 0.0) return scaled.beta;
        double mobility = model.alias(0.0, {"U0", "UO", "MUE", "MOBILITY"});
        if (mobility > 1.0) mobility *= 1.0e-4;
        const double tox = model.alias(0.0, {"TOX", "TOXO"});
        const double epsr = model.alias(3.9, {"EPSROX", "EPSROXO"});
        if (mobility <= 0.0 || tox <= 0.0 || epsr <= 0.0) return 0.0;
        return mobility * 8.8541878128e-12 * epsr / tox;
    }

    std::array<GmcDual4, 4> terminalCharges(const VectorReal& x) const {
        std::array<GmcDual4, 4> charge{
            GmcDual4::constant(0.0), GmcDual4::constant(0.0),
            GmcDual4::constant(0.0), GmcDual4::constant(0.0)};
        const auto vd = GmcDual4::variable(nodeVoltage(x, nodeD_), 0);
        const auto vs = GmcDual4::variable(nodeVoltage(x, nodeS_), 2);
        const auto vb = GmcDual4::variable(nodeVoltage(x, nodeB_), 3);
        const double type = static_cast<double>(prepared_.polarity);
        if (setup_.cox_i > 0.0) {
            const double mult = type * setup_.mult_i;
            const auto comp = intrinsicAndOverlapCharges(x);
            // I(GP,SI) <+ ddt(CHNL_TYPE*MULT*Qg)
            charge[1] = charge[1] + mult * comp.qg;
            charge[2] = charge[2] - mult * comp.qg;
            // I(BP,SI) <+ ddt(CHNL_TYPE*MULT*Qb)
            charge[3] = charge[3] + mult * comp.qb;
            charge[2] = charge[2] - mult * comp.qb;
            // I(DI,SI) <+ ddt(CHNL_TYPE*MULT*Qd)
            charge[0] = charge[0] + mult * comp.qd;
            charge[2] = charge[2] - mult * comp.qd;
            // I(GP,SI) <+ ddt(CHNL_TYPE*MULT*Qfgs)
            charge[1] = charge[1] + mult * comp.qfgs;
            charge[2] = charge[2] - mult * comp.qfgs;
            // I(GP,DI) <+ ddt(CHNL_TYPE*MULT*Qfgd)
            charge[1] = charge[1] + mult * comp.qfgd;
            charge[0] = charge[0] - mult * comp.qfgd;
            // I(GP,BP) <+ ddt(CHNL_TYPE*MULT*Qgb_ov)
            charge[1] = charge[1] + mult * comp.qgb;
            charge[3] = charge[3] - mult * comp.qgb;
        }
        addJunctionCharge(charge, type * (vs - vb), 2, 3, sourceJunctionCap_, sourceJunctionVbi_);
        addJunctionCharge(charge, type * (vd - vb), 0, 3, drainJunctionCap_, drainJunctionVbi_);
        return charge;
    }

    double nodeVoltage(const VectorReal& x, int node) const {
        return node >= 0 && node < x.getSize() ? x[node] : 0.0;
    }

    static bool finiteDualVector(const std::array<GmcDual4, 4>& values) {
        for (const auto& value : values) {
            if (!std::isfinite(value.value)) return false;
            for (double derivative : value.derivative) {
                if (!std::isfinite(derivative)) return false;
            }
        }
        return true;
    }

    std::array<GmcDual4, 4> terminalCurrents(const VectorReal& x) const {
        return terminalCurrentsUnchecked(x);
    }

    std::array<GmcDual4, 4> terminalCurrentsUnchecked(const VectorReal& x) const {
        const auto vd = GmcDual4::variable(nodeVoltage(x, nodeD_), 0);
        const auto vg = GmcDual4::variable(nodeVoltage(x, nodeG_), 1);
        const auto vs = GmcDual4::variable(nodeVoltage(x, nodeS_), 2);
        const auto vb = GmcDual4::variable(nodeVoltage(x, nodeB_), 3);
        const double type = static_cast<double>(prepared_.polarity);

        // PSP103.4 module lines 2045-2057 (voltage affectations) and
        // 2067-2074 (source-drain interchange)
        auto vgs = type * (vg - vs);
        auto vds = type * (vd - vs);
        auto vsb = type * (vs - vb);
        int sigVds = 1;
        if (vds.value < 0.0) {
            sigVds = -1;
            vgs = vgs - vds;
            vsb = vsb + vds;
            vds = -vds;
        }

        // Map Psp103DeviceSetup -> Psp103DcCoreInputs
        Psp103DcCoreInputs in;
        in.phib = setup_.phib_dc;
        in.g0 = setup_.g0_dc;
        in.phit0 = setup_.phit0;
        in.vfb_t = setup_.vfb_t;
        in.kp = setup_.kp;
        in.cf = setup_.cf;
        in.cfd = setup_.cfd;
        in.cfb = setup_.cfb;
        in.psce = setup_.psce;
        in.psce_d = setup_.psce_d;
        in.psce_b = setup_.psce_b;
        in.dnsub = setup_.dnsub;
        in.vnsub = setup_.vnsub;
        in.nslp = setup_.nslp;
        in.rsb = setup_.rsb;
        in.rsg = setup_.rsg;
        in.ther = setup_.ther;
        in.e_eff0 = setup_.e_eff0;
        in.eta_mu = setup_.eta_mu;
        in.eta_mu1 = setup_.eta_mu1;
        in.mue_t = setup_.mue_t;
        in.themu_t = setup_.themu_t;
        in.cs_t = setup_.cs_t;
        in.xcor_t = setup_.xcor_t;
        in.thesat_b = setup_.thesat_b;
        in.thesat_g = setup_.thesat_g;
        in.theta_sat_t = setup_.theta_sat_t;
        in.ax = setup_.ax;
        in.alp = setup_.alp;
        in.alp1 = setup_.alp1;
        in.alp2 = setup_.alp2;
        in.vp = setup_.vp;
        in.bet = setup_.bet_i;
        in.phix = setup_.phix_dc;
        in.aphi = setup_.aphi_dc;
        in.phix1 = setup_.phix1_dc;
        in.vsbnud = setup_.vsbnud;
        in.dvsbnud = setup_.dvsbnud;
        in.gfac_nud = setup_.gfac_nud;
        in.us1 = setup_.us1;
        in.us21 = setup_.us21;
        in.sw_nud = setup_.sw_nud != 0;
        in.pmos = setup_.channel_type == -1;

        const auto core = psp103DcCore(in, vgs, vds, vsb);
        const auto ids = type * setup_.mult_i * core.ids;

        std::array<GmcDual4, 4> currents{
            GmcDual4::constant(0.0), GmcDual4::constant(0.0),
            GmcDual4::constant(0.0), GmcDual4::constant(0.0)};
        // loadStatic: I(DI,SI) <+ CHNL_TYPE * MULT_i * Ids (module lines 2544-2554)
        if (sigVds > 0) {
            currents[0] = currents[0] + ids;
            currents[2] = currents[2] - ids;
        } else {
            currents[2] = currents[2] + ids;
            currents[0] = currents[0] - ids;
        }
        addGateLeakage(currents, vg - vs, 1, 2, gateInvLeakage_ + gateSrcOverlapLeakage_);
        addGateLeakage(currents, vg - vd, 1, 0, gateDrnOverlapLeakage_);
        addGidlLeakage(currents, type * (vd - vg), type * (vd - vb), 0, 3,
                       type, drainGidlScale_, drainGidlBarrier_);
        addGidlLeakage(currents, type * (vs - vg), type * (vs - vb), 2, 3,
                       type, sourceGidlScale_, sourceGidlBarrier_);
        addJunctionLeakage(currents, type * (vd - vb), 0, 3, type, drainJunctionLeakage_);
        addJunctionLeakage(currents, type * (vs - vb), 2, 3, type, sourceJunctionLeakage_);
        return currents;
    }

    static double junctionCapacitance(
        const Psp103ParameterSet& model,
        const Psp103PreparedModel& prepared,
        bool drain) {
        const double area = std::max(prepared.effective_width_m * prepared.effective_length_m, 0.0);
        const double bottom = model.sum(drain
            ? std::initializer_list<const char*>{"CJORBOTD"}
            : std::initializer_list<const char*>{"CJORBOT"}) * area;
        const double scale = model.alias(1.0, {"SWJUNCAP"});
        return std::max(bottom * std::max(scale, 0.0), 0.0);
    }

    static double junctionLeakage(
        const Psp103ParameterSet& model,
        const Psp103PreparedModel& prepared,
        bool drain) {
        const double area = std::max(prepared.effective_width_m * prepared.effective_length_m, 0.0);
        const double perimeter = std::max(prepared.effective_width_m + prepared.effective_length_m, 0.0);
        const double bottom = model.sum(drain
            ? std::initializer_list<const char*>{"IDSATRBOTD"}
            : std::initializer_list<const char*>{"IDSATRBOT"}) * area;
        const double side = model.sum(drain
            ? std::initializer_list<const char*>{"IDSATRGATD", "IDSATRSTID"}
            : std::initializer_list<const char*>{"IDSATRGAT", "IDSATRSTI"}) * perimeter;
        return std::max(bottom + side, 0.0);
    }

    static double junctionBuiltIn(const Psp103ParameterSet& model, bool drain) {
        const double vbi = model.sum(drain
            ? std::initializer_list<const char*>{"VBIRBOTD", "VBIRGATD", "VBIRSTID"}
            : std::initializer_list<const char*>{"VBIRBOT", "VBIRGAT", "VBIRSTI"});
        const double fallback = model.sum(drain
            ? std::initializer_list<const char*>{"PBRBOTD", "PBRGATD", "PBRSTID"}
            : std::initializer_list<const char*>{"PBRBOT", "PBRGAT", "PBRSTI"});
        const double averaged = vbi > 0.0 ? vbi / 3.0 : fallback / 3.0;
        return std::clamp(averaged, 0.1, 2.0);
    }

    void addJunctionCharge(
        std::array<GmcDual4, 4>& charge,
        const GmcDual4& vdb,
        int diffusion,
        int bulk,
        double cap,
        double vbi) const {
        if (cap <= 0.0) return;
        const auto forward = gmcLim(vdb / std::max(vbi, 0.1), -20.0, 0.95);
        const auto q = cap * std::max(vbi, 0.1) * (1.0 - gmcSqrt(1.0 - forward));
        charge[diffusion] = charge[diffusion] + q;
        charge[bulk] = charge[bulk] - q;
    }

    void addJunctionLeakage(
        std::array<GmcDual4, 4>& currents,
        const GmcDual4& vdb,
        int diffusion,
        int bulk,
        double type,
        double scale) const {
        if (scale <= 0.0) return;
        const auto limited = gmcLim(vdb, -2.0, 0.8);
        const auto v = limited / std::max(2.0 * prepared_.temperature.thermal_voltage, 1.0e-6);
        const auto current = type * scale * (gmcExp(v) - 1.0);
        currents[diffusion] = currents[diffusion] + current;
        currents[bulk] = currents[bulk] - current;
    }

    void addGateLeakage(
        std::array<GmcDual4, 4>& currents,
        const GmcDual4& voltage,
        int gate,
        int other,
        double scale) const {
        if (scale <= 0.0) return;
        const auto limited = gmcLim(voltage, -4.0, 4.0);
        const auto current = scale * (gmcExp(limited / 0.35) - gmcExp(-limited / 4.0));
        currents[gate] = currents[gate] + current;
        currents[other] = currents[other] - current;
    }

    void addGidlLeakage(
        std::array<GmcDual4, 4>& currents,
        const GmcDual4& vdg,
        const GmcDual4& vdb,
        int diffusion,
        int bulk,
        double type,
        double scale,
        double barrier) const {
        if (scale <= 0.0) return;
        const auto field = gmcMax(vdg, 0.0);
        if (field.value <= 0.0) return;
        const auto junction = gmcMax(vdb, 0.0);
        const double b = barrier > 0.0 ? barrier : 1.0;
        const auto current = type * scale * junction * gmcExp(gmcLim(-b / (field + 1.0e-3), -80.0, 0.0));
        currents[diffusion] = currents[diffusion] + current;
        currents[bulk] = currents[bulk] - current;
    }

    int nodeD_ = -1;
    int nodeG_ = -1;
    int nodeS_ = -1;
    int nodeB_ = -1;
    Psp103ParameterSet model_;
    Psp103PreparedModel prepared_;
    Psp103TemperatureScaled scaled_;
    Psp103DeviceSetup setup_;
    double cox_ = 0.0;
    double threshold_ = 0.5;
    double beta_ = 1.0e-3;
    double lambda_ = 0.02;
    double nfactor_ = 1.4;
    double gateInvLeakage_ = 0.0;
    double gateSrcOverlapLeakage_ = 0.0;
    double gateDrnOverlapLeakage_ = 0.0;
    double sourceGidlScale_ = 0.0;
    double sourceGidlBarrier_ = 0.0;
    double drainGidlScale_ = 0.0;
    double drainGidlBarrier_ = 0.0;
    double sourceJunctionCap_ = 0.0;
    double drainJunctionCap_ = 0.0;
    double sourceJunctionLeakage_ = 0.0;
    double drainJunctionLeakage_ = 0.0;
    double sourceJunctionVbi_ = 0.7;
    double drainJunctionVbi_ = 0.7;
};

} // namespace gspice

#endif // GSPICE_PSP103_MODEL_HPP
