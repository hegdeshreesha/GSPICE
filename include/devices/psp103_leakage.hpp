#ifndef GSPICE_PSP103_LEAKAGE_HPP
#define GSPICE_PSP103_LEAKAGE_HPP

#include "../gsdi.hpp"

#include <cmath>

namespace gspice {

inline double psp103ScalarMinA(double x, double y, double a) {
    return 0.5 * (x + y - std::sqrt((x - y) * (x - y) + a * a));
}

inline double psp103ScalarMaxA(double x, double y, double a) {
    return 0.5 * (x + y + std::sqrt((x - y) * (x - y) + a * a));
}

struct Psp103GidlInputs {
    double vov = 0.0;
    double v = 0.0;
    double cgidl = 0.0;
    double a = 0.0;
    double b = 0.0;
};

inline double psp103GidlCurrent(const Psp103GidlInputs& in) {
    if (in.vov >= 0.0 || in.a == 0.0) return 0.0;
    const double vtov = std::sqrt(in.vov * in.vov + in.cgidl * in.cgidl * in.v * in.v + 1e-6);
    const double t = in.v * vtov * in.vov;
    return -in.a * t * std::exp(-in.b / vtov);
}

struct Psp103GateOverlapInputs {
    double vgx = 0.0;
    double psi_ov = 0.0;
    double vov = 0.0;
    double igov = 0.0;
    double bov = 0.0;
    double gcq = 0.0;
    double chib = 0.0;
    double gco = 0.0;
    double gc2ov = 0.375;
    double gc3ov = 0.063;
};

inline double psp103GateOverlapCurrent(const Psp103GateOverlapInputs& in) {
    const double vstar = std::sqrt(in.vov * in.vov + 1e-6);
    const double zg = in.gc3ov < 0.0
        ? psp103ScalarMinA(psp103ScalarMinA(vstar, in.chib, in.gcq), 0.0, 1e-6)
        : vstar / in.chib;
    const double phi = 1.0;
    const double fs1 = 3.0 * phi + in.psi_ov / phi;
    const double fs2 = -3.0 - in.gco;
    const double fs3 = 30.0 * in.vgx;
    const double fsov = psp103ScalarMaxA(fs2, psp103ScalarMinA(fs1, fs3, 0.9), 0.3);
    return in.igov * fsov * std::exp(in.bov * (-1.5 + zg * (in.gc2ov + in.gc3ov * zg)));
}

struct Psp103GateChannelInputs {
    double xds_dc = 0.0;
    double vdse_dc = 0.0;
    double v_sb_dc = 0.0;
    double phi_t = 0.02585;
    double vgs = 0.0;
    double alpha_b = 0.0;
    double voxm_dc = 0.0;
    double xm_dc = 0.0;
    double xg_dc = 0.0;
    double delta_psi_dc = 0.0;
    double h_dc = 1.0;
    double b = 0.0;
    double gco = 0.0;
    double gcq = 0.0;
    double chib = 0.0;
    double ig_inv = 0.0;
    double gc2 = 0.375;
    double gc3 = 0.063;
};

struct Psp103GateChannelResult {
    double igc = 0.0;
    double igd = 0.0;
    double igs = 0.0;
    double igb = 0.0;
};

inline Psp103GateChannelResult psp103GateChannelCurrent(const Psp103GateChannelInputs& in) {
    if (in.xg_dc <= 0.0 || in.ig_inv == 0.0) return {};
    const double vm = in.v_sb_dc + in.phi_t * (in.xds_dc / 2.0 -
        std::log1p(std::exp((in.xds_dc - in.vdse_dc / in.phi_t) / 2.0)));
    const double psi_t = psp103ScalarMinA(0.0, in.voxm_dc + in.gco * in.phi_t, 0.01);
    const double voxm = std::sqrt(in.voxm_dc * in.voxm_dc + 1e-6);
    const double zg = in.gc3 < 0.0
        ? psp103ScalarMinA(psp103ScalarMinA(voxm, in.chib, in.gcq), 0.0, 1e-6)
        : voxm / in.chib;
    const double delta_si = std::exp(in.xm_dc - (in.alpha_b + vm - psi_t) / in.phi_t);
    const double fs = std::log((1.0 + delta_si) /
        (1.0 + delta_si * std::exp((-in.vgs + in.v_sb_dc - vm) / in.phi_t)));
    const double igco = in.ig_inv * fs *
        std::exp(in.b * (-1.5 + zg * (in.gc2 + in.gc3 * zg)));
    double pgc = 1.0;
    double pgd = 0.5;
    if (in.xg_dc > 0.0) {
        const double slope = in.b * (in.gc2 + 2.0 * in.gc3 * zg);
        if (slope != 0.0) {
            const double u0 = in.chib / slope;
            const double x = in.delta_psi_dc / (2.0 * u0);
            const double b = u0 / in.h_dc;
            const double bg = b * (1.0 - b) / 2.0;
            const double ag = 0.5 - 3.0 * bg;
            const double sinh_over_x = std::abs(x) < 1e-8 ? 1.0 : std::sinh(x) / x;
            const double coth_minus_inv_x = std::abs(x) < 1e-6
                ? x / 3.0 : std::cosh(x) / std::sinh(x) - 1.0 / x;
            const double pgc0 = (1.0 - b) * sinh_over_x + b * std::cosh(x);
            pgd = pgc0 / 2.0 - bg * std::sinh(x) - ag * sinh_over_x * coth_minus_inv_x;
            pgc = pgc0;
        }
    }
    const double sg = 0.5 * (1.0 + in.xg_dc /
        std::sqrt(in.xg_dc * in.xg_dc + 1e-6));
    const double igc = igco * pgc * sg;
    const double igd = igco * pgd * sg;
    return {igc, igd, igc - igd, igco * pgc * (1.0 - sg)};
}

inline void psp103AppendCurrentBranch(
    GsdiEvalResult& result, int positive_equation, int negative_equation,
    int positive_unknown, int negative_unknown, double current,
    double derivative, int conservation_group = -1) {
    result.static_residual.push_back({positive_equation, current, conservation_group});
    result.static_residual.push_back({negative_equation, -current, conservation_group});
    result.static_jacobian.push_back({positive_equation, positive_unknown, derivative, conservation_group});
    result.static_jacobian.push_back({positive_equation, negative_unknown, -derivative, conservation_group});
    result.static_jacobian.push_back({negative_equation, positive_unknown, -derivative, conservation_group});
    result.static_jacobian.push_back({negative_equation, negative_unknown, derivative, conservation_group});
}

} // namespace gspice

#endif
