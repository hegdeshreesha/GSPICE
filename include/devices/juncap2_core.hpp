#pragma once

#include <algorithm>
#include <cmath>
#include <string>

namespace gspice {

struct Juncap2Inputs {
    double vak = 0.0, phi_td = 0.026, idsat = 0.0;
    double csrh = 0.0, ctat = 0.0, cbbt = 0.0, vbr = 1.0e9;
    double phi_tr = 0.026, delta_vbi = 0.0, fbbtr = 1.0, stfbbt = 0.0;
    double tkd = 300.15, tkr = 300.15;
    double v_max = 40.0 * 0.026, vbi = 0.7, vbi_min = 0.7, vbir = 0.7;
    double p = 0.5, xjun = 1.0, eps_si = 1.03594e-10, cjor = 1.0e-3, ftd = 1.0;
};

struct Juncap2Result {
    double current = 0.0, dcurrent_dvak = 0.0;
    bool supported = true;
    std::string reason;
};

inline std::string juncap2Validate(const Juncap2Inputs& in) {
    if (!std::isfinite(in.phi_td) || in.phi_td <= 0.0) return "PHITD must be finite and positive";
    if (!std::isfinite(in.phi_tr) || in.phi_tr <= 0.0) return "PHITR must be finite and positive";
    if (!std::isfinite(in.p) || in.p < 0.0 || in.p >= 1.0) return "P must be in [0,1)";
    if (!std::isfinite(in.xjun) || in.xjun <= 0.0) return "XJUN must be finite and positive";
    if (!std::isfinite(in.vbir) || in.vbir <= 0.0) return "VBIR must be finite and positive";
    if (!std::isfinite(in.cjor) || in.cjor <= 0.0) return "CJOR must be finite and positive";
    if (!std::isfinite(in.csrh) || in.csrh < 0.0) return "CSRH must be finite and nonnegative";
    if (!std::isfinite(in.ctat) || in.ctat < 0.0) return "CTAT must be finite and nonnegative";
    if (!std::isfinite(in.cbbt) || in.cbbt < 0.0) return "CBBT must be finite and nonnegative";
    if (!std::isfinite(in.vbr) || in.vbr <= 0.0) return "VBR must be finite and positive";
    if (!std::isfinite(in.ftd) || in.ftd < 0.0) return "FTD must be finite and nonnegative";
    return {};
}

inline double juncap2SafeExp(double x) { return std::exp(std::clamp(x, -80.0, 80.0)); }

inline double juncap2IdealMid(double vak, double phi_td, double v_max) {
    if (!(phi_td > 0.0)) return 1.0;
    const double vmax = std::max(0.0, v_max);
    if (vak < vmax) return juncap2SafeExp(vak / phi_td);
    return (1.0 + (vak - vmax) / phi_td) * juncap2SafeExp(vmax / phi_td);
}

inline double juncap2IdealMidDerivative(double vak, double phi_td, double v_max) {
    if (!(phi_td > 0.0)) return 0.0;
    const double vmax = std::max(0.0, v_max);
    return (vak < vmax ? juncap2SafeExp(vak / phi_td) : juncap2SafeExp(vmax / phi_td)) / phi_td;
}

inline double juncap2Hyp2(double value, double lower, double scale) {
    const double s = std::max(std::abs(scale), 1.0e-15);
    const double x = (value - lower) / s;
    return lower + 0.5 * (x + std::sqrt(x * x + 4.0)) * s;
}

inline double juncap2SrhCurrent(const Juncap2Inputs& in, double vak) {
    if (in.csrh == 0.0 || in.ftd == 0.0) return 0.0;
    const double phi = std::max(in.phi_td, 1.0e-15);
    const double mid = std::max(juncap2IdealMid(vak, phi, in.v_max), 1.0e-30);
    const double zinv = std::sqrt(mid), z = 1.0 / zinv, p = std::clamp(in.p, 0.0, 0.999999);
    double psi = 0.0;
    if (vak > 0.0) {
        psi = phi * std::log(std::max(z + 2.0 + p / ((z + 1.0) * (z + 3.0)), 1.0e-30));
    } else {
        psi = -0.5 * vak + phi * std::log(std::max(1.0 + 2.0 * zinv +
            p / ((1.0 + zinv) * (1.0 + 3.0 * zinv)), 1.0e-30));
    }
    const double vj = juncap2Hyp2(vak, in.vbi_min - 2.0 * psi, phi);
    const double denominator = std::max(in.vbi - vj, 1.0e-15);
    const double wstep = 1.0 - std::sqrt(std::max(0.0, 1.0 - 2.0 * psi / denominator));
    const double log_ratio = std::log(std::max(wstep / std::max(1.0 - wstep, 1.0e-15), 1.0e-30));
    const double w = std::max(0.0, wstep + (0.5 * wstep * wstep * log_ratio + wstep) * (1.0 - 2.0 * p));
    const double depletion = std::max(0.0, in.xjun * in.eps_si / std::max(in.cjor, 1.0e-30) *
        std::pow(std::max((in.vbi - vj) / std::max(in.vbir, 1.0e-15), 0.0), p));
    return in.csrh * in.ftd * (zinv - 1.0) * w * depletion;
}

inline double juncap2BbtCurrent(const Juncap2Inputs& in, double vak) {
    if (in.cbbt == 0.0) return 0.0;
    const double p = std::clamp(in.p, 0.0, 0.999999);
    const double phi = std::max(in.phi_tr, 1.0e-15);
    const double vbbt = juncap2Hyp2(vak, in.vbir - in.delta_vbi, phi);
    const double wdep = std::max(in.xjun * in.eps_si / std::max(in.cjor, 1.0e-30) *
        std::pow(std::max((in.vbir - vbbt) / std::max(in.vbir, 1.0e-15), 0.0), p), 1.0e-30);
    const double fmax = std::max((in.vbir - vbbt) / (wdep * std::max(1.0 - p, 1.0e-12)), 1.0e-30);
    const double fbbt = in.fbbtr * (1.0 + in.stfbbt * (in.tkd - in.tkr));
    return in.cbbt * vak * fmax * fmax * juncap2SafeExp(-fbbt / fmax);
}

inline double juncap2TatCurrent(const Juncap2Inputs& in, double vak) {
    if (in.ctat == 0.0 || in.ftd == 0.0) return 0.0;
    const double phi = std::max(in.phi_tr, 1.0e-15);
    const double reverse = juncap2Hyp2(-vak, 0.0, phi);
    const double fieldScale = std::max(in.fbbtr * std::max(in.xjun, 1.0e-15), 1.0e-15);
    const double fieldTerm = reverse / (reverse + fieldScale);
    return -in.ctat * in.ftd * reverse * fieldTerm *
        juncap2SafeExp(-fieldScale / (reverse + phi));
}

inline double juncap2AvalancheCurrent(const Juncap2Inputs& in, double vak) {
    if (in.vbr >= 1.0e9) return 0.0;
    const double phi = std::max(in.phi_tr, 1.0e-15);
    const double over = juncap2Hyp2(-vak - in.vbr, 0.0, phi);
    if (over <= 0.0) return 0.0;
    const double seed = std::max(std::abs(in.idsat) + std::abs(in.csrh) +
                                     std::abs(in.ctat) + std::abs(in.cbbt),
                                 1.0e-30);
    return -seed * (juncap2SafeExp(over / phi) - 1.0);
}

inline Juncap2Result juncap2Evaluate(const Juncap2Inputs& in) {
    Juncap2Result out;
    if (const auto reason = juncap2Validate(in); !reason.empty()) {
        out.supported = false;
        out.reason = reason;
        return out;
    }
    out.current = (juncap2IdealMid(in.vak, in.phi_td, in.v_max) - 1.0) * in.idsat +
        juncap2SrhCurrent(in, in.vak) + juncap2BbtCurrent(in, in.vak) +
        juncap2TatCurrent(in, in.vak) + juncap2AvalancheCurrent(in, in.vak);
    out.dcurrent_dvak = juncap2IdealMidDerivative(in.vak, in.phi_td, in.v_max) * in.idsat;
    if (in.csrh != 0.0 || in.cbbt != 0.0 || in.ctat != 0.0 || in.vbr < 1.0e9) {
        const double h = std::max(1.0e-7, std::abs(in.vak) * 1.0e-6);
        const auto branch = [&in](double vak) {
            return juncap2SrhCurrent(in, vak) + juncap2BbtCurrent(in, vak) +
                juncap2TatCurrent(in, vak) + juncap2AvalancheCurrent(in, vak);
        };
        out.dcurrent_dvak += (branch(in.vak + h) - branch(in.vak - h)) / (2.0 * h);
    }
    return out;
}

} // namespace gspice
