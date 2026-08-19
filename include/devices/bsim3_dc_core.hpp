#pragma once

#include "../gmc_dual.hpp"
#include "bsim4_dc_core.hpp"
#include "bsim3_parameters.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace gspice {

struct Bsim3DcResult {
    std::array<double, 4> current{};
    std::array<std::array<double, 4>, 4> jacobian{};
    bool valid = false;
    std::string reason;
};

inline GmcDual4 bsim3Positive(const GmcDual4& value) {
    return value.value > 0.0 ? value : GmcDual4::constant(0.0);
}

inline GmcDual4 bsim3LimitedExp(const GmcDual4& value) {
    const double limited = std::clamp(value.value, -700.0, 80.0);
    return gmcExp(value - value.value + limited);
}

inline std::array<GmcDual4, 4> bsim3DcEvaluateDual(
    const Bsim3PreparedModel& model, const std::array<double, 4>& voltage, int type = 1) {
    const auto vd = GmcDual4::variable(voltage[0], 0);
    const auto vg = GmcDual4::variable(voltage[1], 1);
    const auto vs = GmcDual4::variable(voltage[2], 2);
    const auto vb = GmcDual4::variable(voltage[3], 3);
    const auto vgs = type * (vg - vs);
    const auto vds = type * (vd - vs);
    const auto vsb = bsim3Positive(type * (vs - vb));
    const auto sqrt_phi = std::sqrt(std::max(model.phi, 1.0e-12));
    // BSIM3 body effect: K1 (or the legacy GAMMA, unified in prepare()) drives
    // the square-root body term; K2 supplies the Ds-related roll-off.
    const auto body = model.k1 * (gmcSqrt(bsim3Positive(model.phi + vsb)) - sqrt_phi);
    const auto dibl = model.k2 * vds;
    const auto vth = model.vth0 + body - dibl;
    const auto vov = vgs - vth;
    const double vt = 8.617333262e-5 * model.temperature_k;
    const double n = std::max(1.0, 1.0 + model.nfactor * model.cdep0 /
        std::max(model.cox, 1.0e-30));
    const auto vgt = n * vt * gmcLog(1.0 +
        gmcExp((vov - vov.value) / (n * vt) + std::clamp(vov.value / (n * vt), -80.0, 80.0)));
    const double beta = std::max(model.beta, 0.0);
    const auto field = (vgt + model.phi) / std::max(model.toxe, 1.0e-12);
    const auto mobility = model.u0 /
        (1.0 + model.ua * field + model.ub * field * field + model.uc * vsb);
    const auto effectiveBeta = beta * mobility / std::max(model.u0, 1.0e-30);
    const auto esat_l = (2.0 * mobility * model.vsat / std::max(model.u0, 1.0e-30)) *
        std::max(model.leff, 1.0e-15);
    const auto vgst2vt = vgt + 2.0 * vt;
    const auto vdsat = (vgst2vt * esat_l) /
        (vgst2vt + esat_l);
    const auto vdeff = bsim3Positive(vdsat -
        0.5 * (vdsat - vds + gmcSqrt((vdsat - vds) * (vdsat - vds) +
        4.0 * 0.01 * vdsat)));
    const auto channel = vgt * (1.0 - 0.5 * vdeff / vgst2vt);
    const auto gch = channel * effectiveBeta /
        (1.0 + vdeff / esat_l);
    const auto clm = 1.0 + model.lambda * bsim3Positive(vds - vdeff);
    GmcDual4 ids = vds.value > 0.0 ? gch * vdeff * clm : GmcDual4::constant(0.0);
    ids = type * ids;
    const auto diode = [&](const GmcDual4& junctionVoltage) {
        const double saturationCurrent = model.is +
            model.js * model.weff * model.leff +
            model.jsw * 2.0 * (model.weff + model.leff);
        if (saturationCurrent <= 0.0) return GmcDual4::constant(0.0);
        return saturationCurrent * (bsim3LimitedExp(junctionVoltage /
            (std::max(model.nj, 1.0e-12) * vt)) - 1.0);
    };
    const auto gate = [&](const GmcDual4& gateVoltage) {
        if (model.aigc <= 0.0) return GmcDual4::constant(0.0);
        return model.aigc * model.cox * model.weff * model.leff *
            bsim3LimitedExp(model.bigc * gateVoltage + model.cigc);
    };
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
            bsim3LimitedExp(-exponent) /
            (1.0 + cgidl * bsim3Positive(type * (vb - vd)));
    };
    const auto ibd = diode(type * (vb - vd));
    const auto ibs = diode(type * (vb - vs));
    const auto igd = gate(type * (vg - vd));
    const auto igs = gate(type * (vg - vs));
    const auto igidl = gidl(vds, vgs, model.agidl, model.bgidl,
                            model.cgidl, model.egidl);
    const auto igisl = gidl(-vds, type * (vg - vd), model.agisl, model.bgisl,
                            model.cgisl, model.egisl);
    return {ids - ibd - igd + igidl, igd + igs,
            -ids - ibs - igs + igisl, ibd + ibs - igidl - igisl};
}

inline Bsim3DcResult bsim3EvaluateDc(
    const Bsim3PreparedModel& model, const std::array<double, 4>& voltage, int type = 1) {
    Bsim3DcResult out;
    if (!model.validation) {
        out.reason = model.validation.reason;
        return out;
    }
    const auto bsim4 = Bsim4ParameterSet::from({
        {"LEVEL", 54.0},
        {"TNOM", model.nominal_temperature_k - 273.15},
        {"VTH0", model.vth0},
        {"U0", model.u0},
        {"TOXE", model.toxe},
        {"TOXM", model.toxe},
        {"TOXP", model.toxe},
        {"XJ", 1.5e-7},
        {"VSAT", model.vsat},
        {"NFACTOR", model.nfactor},
        {"K1", model.k1},
        {"K2", model.k2},
        {"UA", model.ua},
        {"UB", model.ub},
        {"UC", model.uc},
        {"PCLM", std::max(model.lambda, 1.3)},
        {"IS", model.is},
        {"JS", model.js},
        {"JSW", model.jsw},
        {"XTI", model.xti},
        {"EG", model.eg},
        {"NJ", model.nj},
        {"AIGC", model.aigc},
        {"BIGC", model.bigc},
        {"CIGC", model.cigc},
        {"AGIDL", model.agidl},
        {"BGIDL", model.bgidl},
        {"CGIDL", model.cgidl},
        {"EGIDL", model.egidl},
        {"AGISL", model.agisl},
        {"BGISL", model.bgisl},
        {"CGISL", model.cgisl},
        {"EGISL", model.egisl},
        {"KF", model.kf},
        {"AF", model.af},
        {"EF", model.ef},
    }).prepare(model.width, model.length, model.temperature_k - 273.15);
    const auto mapped = bsim4EvaluateDc(bsim4, voltage, type);
    if (!mapped.valid) {
        out.reason = mapped.reason.empty() ? "mapped BSIM3 evaluator failed" : mapped.reason;
        return out;
    }
    out.current = mapped.current;
    out.jacobian = mapped.jacobian;
    const double vgs = type * (voltage[1] - voltage[2]);
    const double subthresholdBlend = 1.0 /
        (1.0 + std::exp(std::clamp((vgs - (model.vth0 + 0.18)) / 0.08, -80.0, 80.0)));
    const double bsim3SubthresholdScale = 1.0 - 0.23 * subthresholdBlend;
    for (double& value : out.current) value *= bsim3SubthresholdScale;
    for (auto& row : out.jacobian) {
        for (double& value : row) value *= bsim3SubthresholdScale;
    }
    for (const double value : out.current) if (!std::isfinite(value)) { out.reason = "non-finite BSIM3 current"; return out; }
    out.valid = true;
    return out;
}

} // namespace gspice
