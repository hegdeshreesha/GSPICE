#pragma once

#include "../gmc_dual.hpp"
#include "bsim4_parameters.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace gspice {

struct Bsim4DcResult {
    std::array<double, 4> current{};
    std::array<std::array<double, 4>, 4> jacobian{};
    struct State {
        double vth = 0.0;
        double vgtEff = 0.0;
        double vdsat = 0.0;
        double abulk = 0.0;
        double vasat = 0.0;
        double mobility = 0.0;
        double effectiveBeta = 0.0;
        double esatLength = 0.0;
        double lambdaFactor = 0.0;
        double idl = 0.0;
        double clmBoost = 1.0;
        double vadibl = 0.0;
    } state;
    bool valid = false;
    std::string reason;
};

using Bsim4DcState = Bsim4DcResult::State;

inline GmcDual4 bsim4Nonnegative(const GmcDual4& value) {
    return value.value > 0.0 ? value : GmcDual4::constant(0.0);
}

inline GmcDual4 bsim4LimitedExp(const GmcDual4& value) {
    // BSIM4 uses a very small positive MIN_EXP in its weak-inversion
    // asymptotes; -80 would create an artificial off-state floor.
    const double limited = std::clamp(value.value, -700.0, 80.0);
    return gmcExp(value - value.value + limited);
}

inline std::array<GmcDual4, 4> bsim4DcEvaluateDual(
    const Bsim4PreparedModel& model,
    const std::array<double, 4>& voltage,
    int type = 1,
    Bsim4DcState* state = nullptr) {
    const auto vd = GmcDual4::variable(voltage[0], 0);
    const auto vg = GmcDual4::variable(voltage[1], 1);
    const auto vs = GmcDual4::variable(voltage[2], 2);
    const auto vb = GmcDual4::variable(voltage[3], 3);
    const auto vgs = type * (vg - vs);
    const auto vds = type * (vd - vs);
    const auto vsb = bsim4Nonnegative(type * (vs - vb));
    const double phi = std::max(model.phi, 1.0e-6);
    const double sqrtPhi = std::sqrt(phi);
    const double epssub = 1.03594e-10;
    const double epsRatioToxe = (epssub / (3.9 * 8.8541878128e-12)) *
                                std::max(model.toxe, 1.0e-12);
    const double factor1 = std::sqrt(epsRatioToxe);
    // Vbseff smoothing (b4ld.c:1002-1019).  vsb is the nonnegative body bias.
    const double vbsc = -30.0;
    const auto vbse0 = vsb - vbsc - 0.001;
    const auto vbse1 = gmcSqrt(vbse0 * vbse0 - 0.004 * vbsc);
    const auto vbseStage = vbse0.value >= 0.0
        ? vbsc + 0.5 * (vbse0 + vbse1)
        : vbsc * (1.0 - 0.002 / (vbse1.value - vbse0.value));
    const double vbsfwd = 0.95 * phi;
    const auto vbsf0 = vbsfwd - vbseStage - 0.001;
    const auto vbsf1 = gmcSqrt(vbsf0 * vbsf0 + 0.004 * vbsfwd);
    const auto vbseff = vbsfwd - 0.5 * (vbsf0 + vbsf1);
    const auto phis = phi - vbseff;
    const auto sqrtPhis = gmcSqrt(phis);
    const double xdep = model.xdep0 * sqrtPhis.value / sqrtPhi;
    // DVT decay factors (b4ld.c:1037-1094).
    const auto dvtTheta = [&](double dv1, double dv2) -> double {
        const double t0 = dv2 * vbseff.value;
        double t1v;
        if (t0 >= -0.5) { t1v = 1.0 + t0; }
        else { const double t4 = 1.0 / (3.0 + 8.0 * t0); t1v = (1.0 + 3.0 * t0) * t4; }
        const double lt = std::max(factor1 * std::sqrt(xdep) * t1v, 1.0e-30);
        const double t = dv1 * model.leff / lt;
        if (t < 80.0) {
            const double e = std::exp(t);
            return e / ((e - 1.0) * (e - 1.0) + 2.0 * e * 1.0e-300);
        }
        return 1.0 / (1.0e36 - 2.0);
    };
    const double th0 = dvtTheta(model.dvt1, model.dvt2);
    const double th0w = dvtTheta(model.dvt1w, model.dvt2w);
    // DIBL coupling (b4temp.c:1519-1529).
    const double thetaExponent = model.dsub * model.leff /
        std::sqrt(std::max(epsRatioToxe * model.xdep0, 1.0e-30));
    const double thetaExp = std::exp(std::clamp(thetaExponent, -700.0, 80.0));
    const double theta0vb0 =
        thetaExp / ((thetaExp - 1.0) * (thetaExp - 1.0) + 2.0e-280);
    // Official 4.8.3 Vth assembly (b4ld.c:1034-1129).
    const double v0vbi = model.vbi - phi;
    const double deltVth = model.dvt0 * th0 * v0vbi;
    const double delt2w = model.dvt0w * th0w * v0vbi;
    const double tempRatio = model.temperature_k / model.nominal_temperature_k - 1.0;
    const double leffI = std::max(model.leff, 1.0e-15);
    const double lpeVb = std::sqrt(1.0 + model.lpeb / leffI);
    const double vthNarrowW = model.toxe * phi / (model.weff + model.w0);
    const double lpeTerm = model.k1ox *
        (std::sqrt(1.0 + model.lpe0 / leffI) - 1.0) * sqrtPhi +
        (model.kt1 + model.kt1l / leffI + model.kt2 * vbseff.value) * tempRatio;
    double etaD = model.eta0 + model.etab * vbseff.value;
    if (etaD < 1.0e-4) {
        const double t1 = 1.0 / (3.0 - 2.0e4 * etaD);
        etaD = (2.0e-4 - etaD) * t1;
    }
    const auto diblShift = GmcDual4::constant(theta0vb0 * etaD) * vds;
    const auto vth = GmcDual4::constant(
        model.vth0 - deltVth - delt2w + lpeTerm - model.k1 * sqrtPhi * lpeVb +
        model.k3 * vthNarrowW)
        + sqrtPhis * GmcDual4::constant(model.k1ox * lpeVb)
        + vbseff * GmcDual4::constant(-model.k2ox + model.k3b * vthNarrowW)
        - diblShift;
    const auto vov = vgs - vth;
    const double vt = 8.617333262e-5 * model.temperature_k;
    // BSIM4.8.3 subthreshold slope (b4ld.c:1133-1154).
    const double cdscEff = model.cdsc + model.cdscb * vbseff.value +
                           model.cdscd * vds.value;
    const double nFactorTerm =
        (model.nfactor * epssub / std::max(xdep, 1.0e-30) +
         cdscEff * th0 + model.cit) / model.cox;
    const double n = nFactorTerm >= -0.5 ? 1.0 + nFactorTerm
                    : (1.0 + 3.0 * nFactorTerm) / (3.0 + 8.0 * nFactorTerm);
    const double mstar = 0.5 + std::atan(model.minv) / 3.141592653589793;
    // BSIM4.8.3 Vgsteff smoothing (b4ld.c:1236-1296): above the exp
    // threshold the softplus is replaced by the exact linear law
    // (T10 = mstar*Vgst), which the limited-exp clamp alone cannot reach.
    const auto mstarVov = mstar * vov;
    const double vgtExponent = mstarVov.value / (n * vt);
    const auto softplus = [&]() -> GmcDual4 {
        if (vgtExponent > 80.0) return mstarVov;
        if (vgtExponent < -80.0) return GmcDual4::constant(n * vt * 1.0e-300);
        return n * vt * gmcLog(1.0 + gmcExp(mstarVov / (n * vt)));
    }();
    const double voffcbn = model.voff;
    const double denExponent =
        (voffcbn - (1.0 - mstar) * vov.value) / (n * vt);
    const double coxCdep0 = model.cox / std::max(model.cdep0, 1.0e-30);
    const auto vgtDenominator = [&]() -> GmcDual4 {
        if (denExponent > 80.0) return GmcDual4::constant(mstar + n * coxCdep0 * 1.0e36);
        if (denExponent < -80.0) return GmcDual4::constant(mstar + n * coxCdep0 * 1.0e-300);
        return mstar + n * coxCdep0 *
            gmcExp((voffcbn - (1.0 - mstar) * vov) / (n * vt));
    }();
    const auto vgtEff = softplus / vgtDenominator;
    const double length = std::max(model.leff, 1.0e-15);
    const double beta = model.kp * model.weff / length;
    // BSIM4 mobMod=0 uses the effective vertical field
    // (Vgsteff + 2*Vth) / TOXE.  Keep it in SI units: UA/UB/UC are already
    // normalized by the parameter loader.
    const auto field = (vgtEff + 2.0 * vth) /
                       std::max(model.toxe, 1.0e-12);
    const auto mobilityDenominator = 1.0 + field *
        (model.ua + model.uc * vsb + model.ub * field);
    const auto mobility = model.u0 / mobilityDenominator;
    // Quantum-limited effective oxide capacitance (b4ld.c:1790-1805).
    const auto qFactor0 = (vgtEff + model.vtfbphi2) /
                          (2.0e8 * std::max(model.toxp, 1.0e-15));
    const double qExp = std::exp(std::clamp(
        model.bdos * 0.7 * std::log(std::max(qFactor0.value, 1.0e-300)),
        -709.0, 700.0));
    const double qCen = model.ados * 1.9e-9 / (1.0 + qExp);
    const double coxEff = model.coxp * epssub /
        std::max(epssub + model.coxp * qCen, 1.0e-30);
    const auto effectiveBeta = beta * mobility / std::max(model.u0, 1.0e-30) *
                               (coxEff / std::max(model.cox, 1.0e-30));
    const auto esatLength = 2.0 * model.vsat * length /
                            std::max(mobility.value, 1.0e-30);
    const double xj = std::max(model.xj, 1.0e-12);
    const double sqrtXjXdep = std::sqrt(xj * xdep);
    const auto t5 = model.leff / (model.leff + 2.0 * sqrtXjXdep);
    // BSIM4.8.3 Abulk0 (b4ld.c:1338-1376).
    const auto bulkFactor = 0.5 * model.k1ox * lpeVb / std::max(sqrtPhis.value, 1.0e-12)
        + model.k2ox - model.k3b * vthNarrowW;
    const auto bulkGeometry = model.a0 * t5;
    const auto bulkWidth = model.b0 /
        (model.b1 + std::max(model.weff, 1.0e-18));
    const auto abulk0 = 1.0 + bulkFactor * (bulkGeometry + bulkWidth);
    const auto abulkBeforeKeta = abulk0 -
        model.ags * model.a0 * t5 * t5 * t5 * vgtEff;
    const auto ketaTerm = model.keta * vbseff;
    const auto abulk = abulkBeforeKeta /
        (ketaTerm.value >= -0.9 ? 1.0 + ketaTerm : 0.8 + ketaTerm);
    const auto rdsBias = 1.0 + model.prwg * vgtEff;
    const auto rdsBody = model.prwb * sqrtPhis;
    const auto rdsShape = 1.0 / rdsBias + rdsBody;
    const auto rds = model.rdsmod == 1 ? GmcDual4::constant(0.0) :
        model.rdswmin + 0.5 * model.rds0 *
        (rdsShape + gmcSqrt(rdsShape * rdsShape + 0.01));
    const auto vgst2Vt = vgtEff + 2.0 * vt;
    const auto lambdaFactor = [&] {
        if (model.a1 == 0.0) return GmcDual4::constant(model.a2);
        if (model.a1 > 0.0) {
            const double t0 = 1.0 - model.a2;
            const auto t1 = t0 - model.a1 * vgtEff - 0.0001;
            const auto t2 = gmcSqrt(t1 * t1 + 0.0004 * t0);
            return GmcDual4::constant(model.a2 + t0) - 0.5 * (t1 + t2);
        }
        const auto t1 = model.a2 + model.a1 * vgtEff - 0.0001;
        const auto t2 = gmcSqrt(t1 * t1 + 0.0004 * model.a2);
        return 0.5 * (t1 + t2);
    }();
    const auto wvcoxRds = model.weff * model.vsat * model.cox * rds;
    const auto vdsat = [&] {
        if (rds.value <= 0.0 && std::abs(lambdaFactor.value - 1.0) < 1.0e-12) {
            return (vgst2Vt * esatLength) /
                   (abulk * esatLength + vgst2Vt);
        }
        const auto invLambda = 1.0 / bsim4Nonnegative(lambdaFactor);
        const auto t9 = abulk * wvcoxRds;
        const auto t0 = 2.0 * abulk * (t9 - 1.0 + invLambda);
        const auto t1 = vgst2Vt * (2.0 * invLambda - 1.0) +
                        abulk * esatLength + 3.0 * vgst2Vt * t9;
        const auto t2 = vgst2Vt *
                        (esatLength + 2.0 * vgst2Vt * wvcoxRds);
        const auto discriminant = bsim4Nonnegative(t1 * t1 - 2.0 * t0 * t2);
        if (std::abs(t0.value) < 1.0e-30)
            return t2 / bsim4Nonnegative(t1);
        return (t1 - gmcSqrt(discriminant)) / t0;
    }();
    if (state) {
        state->vth = vth.value;
        state->vgtEff = vgtEff.value;
        state->vdsat = vdsat.value;
        state->abulk = abulk.value;
        const auto vasatFactor = 1.0 - 0.5 * abulk * vdsat /
                                 bsim4Nonnegative(vgst2Vt);
        const auto vasat = (esatLength + vdsat +
                            2.0 * wvcoxRds * vgtEff * vasatFactor) /
                           (2.0 / bsim4Nonnegative(lambdaFactor) - 1.0 +
                            abulk * wvcoxRds);
        state->vasat = vasat.value;
    }
    const auto delta = std::max(model.delta, 1.0e-6);
    const auto transition = vdsat - vds - delta;
    const auto smoothing = gmcSqrt(transition * transition + 4.0 * delta * vdsat);
    const auto vdsEffRaw = transition.value >= 0.0
        ? vdsat - 0.5 * (transition + smoothing)
        : vdsat * (1.0 - 2.0 * delta / (smoothing - transition));
    const auto vdsEff = bsim4Nonnegative(vdsEffRaw);
    GmcDual4 ids = GmcDual4::constant(0.0);
    GmcDual4 idl = GmcDual4::constant(0.0);
    GmcDual4 clmBoost = GmcDual4::constant(1.0);
    GmcDual4 vadibl = GmcDual4::constant(0.0);
    if (vds.value > 0.0) {
        const auto channelFactor = vgtEff *
            (1.0 - 0.5 * abulk * vdsEff / vgst2Vt);
        const auto channelConductance = effectiveBeta * channelFactor /
            (1.0 + vdsEff / esatLength);
        idl = channelConductance /
              (1.0 + rds * channelConductance);
        ids = idl;
        const auto diffVds = vds - vdsEff;
        // Official 4.8.3 channel-length-modulation / DIBL-in-Va factors
        // (b4ld.c:1853-2045).  FP and PvagTerm are 1 at their defaults.
        const double clmAllowed = model.pclm > 1.0e-300 && diffVds.value > 1.0e-10;
        if (clmAllowed) {
            const double cclmFactor = (1.0 + rds.value * idl.value) *
                (model.leff + vdsat.value * model.leff / std::max(esatLength, 1.0e-30)) /
                (model.pclm * std::max(model.litl, 1.0e-30));
            const auto vaclm = GmcDual4::constant(cclmFactor) * diffVds;
            const auto vasatLocal = (esatLength + vdsat +
                2.0 * wvcoxRds * vgtEff *
                (1.0 - 0.5 * abulk * vdsat / bsim4Nonnegative(vgst2Vt))) /
                (2.0 / bsim4Nonnegative(lambdaFactor) - 1.0 + abulk * wvcoxRds);
            const auto va = vasatLocal + vaclm;
            clmBoost = 1.0 + gmcLog(va / bsim4Nonnegative(vasatLocal)) /
                std::max(cclmFactor, 1.0e-30);
            if (model.thetaRout > 1.0e-300) {
                const auto t8 = abulk * vdsat;
                vadibl = (vgst2Vt - vgst2Vt * t8 /
                    (vgst2Vt + t8)) / model.thetaRout;
                ids = idl * (1.0 + diffVds / bsim4Nonnegative(vadibl)) * clmBoost;
            } else {
                ids = idl * clmBoost;
            }
        }
    }
    // b4ld.c: cdrain = Ids * Vdseff (Vdseff applied at the very end).
    if (state) {
        state->mobility = mobility.value;
        state->effectiveBeta = effectiveBeta.value;
        state->esatLength = esatLength;
        state->lambdaFactor = lambdaFactor.value;
        state->idl = idl.value;
        state->clmBoost = clmBoost.value;
        state->vadibl = vadibl.value;
    }
    ids = vdsEff * ids;
    ids = type * ids;
    const auto diode = [&](const GmcDual4& junctionVoltage) {
        const double saturationCurrent = model.is +
            model.js * model.weff * model.leff +
            model.jsw * 2.0 * (model.weff + model.leff);
        if (saturationCurrent <= 0.0) return GmcDual4::constant(0.0);
        return saturationCurrent * (bsim4LimitedExp(junctionVoltage /
            (std::max(model.nj, 1.0) * vt)) - 1.0);
    };
    const auto gate = [&](const GmcDual4& gateVoltage) {
        if (model.aigc <= 0.0) return GmcDual4::constant(0.0);
        return model.aigc * model.cox * model.weff * model.leff *
            bsim4LimitedExp(model.bigc * gateVoltage + model.cigc);
    };
    const auto ibd = diode(type * (vb - vd));
    const auto ibs = diode(type * (vb - vs));
    const auto igd = gate(type * (vg - vd));
    const auto igs = gate(type * (vg - vs));
    const auto gidl = [&](const GmcDual4& drainVoltage,
                          const GmcDual4& gateVoltage,
                          double agidl, double bgidl, double cgidl,
                          double egidl) {
        if (agidl <= 0.0 || bgidl <= 0.0 || cgidl <= 0.0) return GmcDual4::constant(0.0);
        const auto field = (drainVoltage - gateVoltage - egidl) /
                           std::max(3.0 * model.toxe, 1.0e-30);
        if (field.value <= 0.0) return GmcDual4::constant(0.0);
        const auto exponent = bgidl / field;
        return agidl * model.weff * field *
            bsim4LimitedExp(-exponent) /
            (1.0 + cgidl * bsim4Nonnegative(type * (vb - vd)));
    };
    const auto igidl = gidl(vds, vgs, model.agidl, model.bgidl,
                            model.cgidl, model.egidl);
    const auto igisl = gidl(-vds, type * (vg - vd), model.agisl, model.bgisl,
                            model.cgisl, model.egisl);
    return {ids - ibd - igd + igidl, igd + igs,
            -ids - ibs - igs + igisl, ibd + ibs - igidl - igisl};
}

inline Bsim4DcResult bsim4EvaluateDc(
    const Bsim4PreparedModel& model,
    const std::array<double, 4>& voltage,
    int type = 1) {
    Bsim4DcResult out;
    if (!model.validation) {
        out.reason = model.validation.reason;
        return out;
    }
    const auto currents = bsim4DcEvaluateDual(model, voltage, type, &out.state);
    for (std::size_t row = 0; row < currents.size(); ++row) {
        out.current[row] = currents[row].value;
        for (std::size_t column = 0; column < currents[row].derivative.size(); ++column)
            out.jacobian[row][column] = currents[row].derivative[column];
    }
    out.valid = true;
    return out;
}

} // namespace gspice
