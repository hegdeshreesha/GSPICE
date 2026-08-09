#ifndef GSPICE_PSP103_CORE_HPP
#define GSPICE_PSP103_CORE_HPP

#include "gmc_dual.hpp"

#include <cmath>

namespace gspice {

inline GmcDual4 psp103MinA(const GmcDual4& x, const GmcDual4& y, double a) {
    return 0.5 * (x + y - gmcSqrt((x - y) * (x - y) + a));
}

inline GmcDual4 psp103MaxA(const GmcDual4& x, const GmcDual4& y, double a) {
    return 0.5 * (x + y + gmcSqrt((x - y) * (x - y) + a));
}

inline GmcDual4 psp103Chi(const GmcDual4& y) {
    const GmcDual4 y2 = y * y;
    return y2 / (2.0 + y2);
}

inline GmcDual4 psp103ChiPrime(const GmcDual4& y) {
    const GmcDual4 y2 = y * y;
    const GmcDual4 denominator = 2.0 + y2;
    return 4.0 * y / (denominator * denominator);
}

inline GmcDual4 psp103ChiSecond(const GmcDual4& y) {
    const GmcDual4 y2 = y * y;
    const GmcDual4 denominator = 2.0 + y2;
    return (8.0 - 12.0 * y2) / (denominator * denominator * denominator);
}

inline GmcDual4 psp103Sigma1(
    const GmcDual4& a,
    const GmcDual4& c,
    const GmcDual4& tau,
    const GmcDual4& eta) {
    const GmcDual4 nu = a + c;
    const GmcDual4 mutau = nu * nu + tau * (0.5 * c * c - a);
    const GmcDual4 correction = (nu / mutau) * tau * tau * c * (c * c / 3.0 - a);
    return eta + a * nu * tau / (mutau + correction);
}

inline GmcDual4 psp103Sigma2(
    const GmcDual4& a,
    const GmcDual4& b,
    const GmcDual4& c,
    const GmcDual4& tau,
    const GmcDual4& eta) {
    const GmcDual4 nu = a + c;
    const GmcDual4 mutau = nu * nu + tau * (0.5 * c * c - a * b);
    const GmcDual4 correction =
        (nu / mutau) * tau * tau * c * (c * c / 3.0 - a * b);
    return eta + a * nu * tau / (mutau + correction);
}

inline GmcDual4 psp103PositiveRoot(const GmcDual4& value) {
    return value.value > 0.0 ? gmcSqrt(value) : GmcDual4::constant(0.0);
}

inline GmcDual4 psp103MinA(
    const GmcDual4& x, const GmcDual4& y, const GmcDual4& z, double a) {
    return psp103MinA(psp103MinA(x, y, a), z, a);
}

struct Psp103TerminalConditioning {
    GmcDual4 vdsx;
    GmcDual4 phi_v;
    GmcDual4 vsb_star;
};

// PSP103 equations 4.95-4.98 for the DC terminal conditioning path.
inline Psp103TerminalConditioning psp103ConditionTerminals(
    const GmcDual4& vds,
    const GmcDual4& vsb,
    double b_phi,
    double phi_x,
    double a_phi,
    double phi_x_star) {
    const GmcDual4 vds_squared = vds * vds;
    Psp103TerminalConditioning result;
    result.vdsx = vds_squared / gmcSqrt(vds_squared + 0.01) + 0.1;
    result.phi_v = psp103MinA(vsb, vsb + vds, b_phi) + phi_x;
    result.vsb_star = vsb - psp103MinA(
        result.phi_v, GmcDual4::constant(0.0), a_phi) + phi_x_star;
    return result;
}

struct Psp103SurfacePotentialState {
    GmcDual4 xi;
    GmcDual4 xns;
    GmcDual4 delta_ns;
    GmcDual4 xg;
};

// PSP103 equations 4.123-4.125 and 4.121. The complete surface-potential
// solve is assembled from these auxiliary functions in the next model block.
inline Psp103SurfacePotentialState psp103SourceSide(
    const GmcDual4& v_gb_star,
    const GmcDual4& v_sb_star,
    double phi_b,
    double phi_t_star,
    double body_factor) {
    Psp103SurfacePotentialState state;
    state.xi = GmcDual4::constant(1.0 + body_factor / std::sqrt(2.0));
    state.xns = (GmcDual4::constant(phi_b) + v_sb_star) / phi_t_star;
    state.delta_ns = gmcExp(-state.xns);
    state.xg = v_gb_star / phi_t_star;
    return state;
}

// Explicit PSP103 source-side surface-potential approximation. The branch
// equations follow the model's Appendix-A auxiliary functions and are kept
// separate so the future drain-side implementation can share them.
inline GmcDual4 psp103SourceSurfacePotential(
    const Psp103SurfacePotentialState& state,
    const GmcDual4& body_factor) {
    const double margin = 1.0e-5 * state.xi.value;
    if (std::abs(state.xg.value) <= margin) {
        const double inv_sqrt_2 = 1.0 / std::sqrt(2.0);
        return state.xg / state.xi *
            (1.0 + state.xg * (1.0 - state.delta_ns) * body_factor /
                (state.xi * state.xi) * (1.0 / 6.0) * inv_sqrt_2);
    }

    const GmcDual4 g2 = body_factor * body_factor;
    if (state.xg.value < -margin) {
        const GmcDual4 yg = -state.xg;
        const GmcDual4 ysub = 1.25 * yg / state.xi;
        const GmcDual4 eta = 0.5 *
            (ysub + 10.0 - gmcSqrt((ysub - 6.0) * (ysub - 6.0) + 64.0));
        const GmcDual4 temp = yg - eta;
        const GmcDual4 a = temp * temp + g2 * (eta + 1.0);
        const GmcDual4 c = 2.0 * temp - g2;
        const GmcDual4 tau = -eta + gmcLog(a / g2);
        const GmcDual4 y0 = psp103Sigma1(a, c, tau, eta);
        const GmcDual4 delta0 = gmcExp(y0);
        const GmcDual4 delta1 = 1.0 / delta0;
        const GmcDual4 chi0 = psp103Chi(y0);
        const GmcDual4 chi1 = psp103ChiPrime(y0);
        const GmcDual4 chi2 = psp103ChiSecond(y0);
        const GmcDual4 channel = yg - y0;
        const GmcDual4 delta_term = state.delta_ns * delta1;
        const GmcDual4 p = 2.0 * channel + g2 *
            (delta0 - 1.0 - delta_term + state.delta_ns * (1.0 - chi1));
        const GmcDual4 q = channel * channel - g2 *
            (delta0 - y0 - 1.0 + delta_term +
             state.delta_ns * (y0 - 1.0 - chi0));
        const GmcDual4 discriminant = p * p - 2.0 * q *
            (2.0 - g2 * (delta0 + delta_term - state.delta_ns * chi2));
        return -y0 - 2.0 * q / (p + psp103PositiveRoot(discriminant));
    }

    const GmcDual4 xg1 = 1.0 / (1.25 + body_factor * 0.7324648775);
    const GmcDual4 xbar = state.xg / state.xi *
        (1.0 + ((state.xi * 1.25 * xg1 - 1.0) * xg1) * state.xg);
    const GmcDual4 w = 1.0 - gmcExp(-xbar);
    const GmcDual4 x1 = state.xg + 0.5 * g2 - body_factor *
        psp103PositiveRoot(state.xg + 0.25 * g2 - w);
    const GmcDual4 bx = state.xns + 3.0;
    const GmcDual4 eta = psp103MinA(x1, bx, 5.0) -
        0.5 * (bx - gmcSqrt(bx * bx + 5.0));
    const GmcDual4 channel = state.xg - eta;
    const GmcDual4 exp_neg_eta = gmcExp(-eta);
    const GmcDual4 chi0 = psp103Chi(eta);
    const GmcDual4 chi1 = psp103ChiPrime(eta);
    const GmcDual4 chi2 = psp103ChiSecond(eta);
    const GmcDual4 a = psp103MaxA(
        channel * channel - g2 * (exp_neg_eta + eta - 1.0 -
            state.delta_ns * (eta + 1.0 + chi0)),
        GmcDual4::constant(1.0e-40), 1.0e-40);
    const GmcDual4 b = 1.0 - 0.5 * g2 *
        (exp_neg_eta - state.delta_ns * chi2);
    const GmcDual4 c = 2.0 * channel + g2 *
        (1.0 - exp_neg_eta - state.delta_ns * (1.0 + chi1));
    const GmcDual4 tau = state.xns - eta + gmcLog(a / g2);
    const GmcDual4 x0 = psp103Sigma2(a, b, c, tau, eta);
    const GmcDual4 delta0 = state.delta_ns * gmcExp(x0);
    const GmcDual4 delta1 = gmcExp(-x0);
    const GmcDual4 chi_x0 = psp103Chi(x0);
    const GmcDual4 chi_x0_prime = psp103ChiPrime(x0);
    const GmcDual4 chi_x0_second = psp103ChiSecond(x0);
    const GmcDual4 channel0 = state.xg - x0;
    const GmcDual4 p = 2.0 * channel0 + g2 *
        (1.0 - delta1 + delta0 - state.delta_ns * (1.0 + chi_x0_prime));
    const GmcDual4 q = channel0 * channel0 - g2 *
        (delta1 + x0 - 1.0 + delta0 -
         state.delta_ns * (x0 + 1.0 + chi_x0));
    const GmcDual4 discriminant = p * p - 2.0 * q *
        (2.0 - g2 * (delta1 + delta0 - state.delta_ns * chi_x0_second));
    return x0 + 2.0 * q / (p + psp103PositiveRoot(discriminant));
}

// PSP103 drain-side approximation. Unlike the source-side solve, the drain
// path uses the common transition/positive branch after x1 is formed.
inline GmcDual4 psp103DrainSurfacePotential(
    const Psp103SurfacePotentialState& state,
    const GmcDual4& body_factor) {
    const double margin = 1.0e-5 * state.xi.value;
    if (std::abs(state.xg.value) <= margin) {
        return state.xg / state.xi *
            (1.0 + state.xg * (1.0 - state.delta_ns) * body_factor /
                (state.xi * state.xi) * (1.0 / 6.0) / std::sqrt(2.0));
    }

    const GmcDual4 g2 = body_factor * body_factor;
    const GmcDual4 xg1 = 1.0 / (1.25 + body_factor * 0.7324648775);
    const GmcDual4 xbar = state.xg / state.xi *
        (1.0 + ((state.xi * 1.25 * xg1 - 1.0) * xg1) * state.xg);
    const GmcDual4 w = 1.0 - gmcExp(-xbar);
    const GmcDual4 x1 = state.xg + 0.5 * g2 - body_factor *
        psp103PositiveRoot(state.xg + 0.25 * g2 - w);
    const GmcDual4 bx = state.xns + 3.0;
    const GmcDual4 eta = psp103MinA(x1, bx, 5.0) -
        0.5 * (bx - gmcSqrt(bx * bx + 5.0));
    const GmcDual4 channel = state.xg - eta;
    const GmcDual4 exp_neg_eta = gmcExp(-eta);
    const GmcDual4 chi0 = psp103Chi(eta);
    const GmcDual4 chi1 = psp103ChiPrime(eta);
    const GmcDual4 chi2 = psp103ChiSecond(eta);
    const GmcDual4 a = psp103MaxA(
        channel * channel - g2 * (exp_neg_eta + eta - 1.0 -
            state.delta_ns * (eta + 1.0 + chi0)),
        GmcDual4::constant(1.0e-40), 1.0e-40);
    const GmcDual4 b = 1.0 - 0.5 * g2 *
        (exp_neg_eta - state.delta_ns * chi2);
    const GmcDual4 c = 2.0 * channel + g2 *
        (1.0 - exp_neg_eta - state.delta_ns * (1.0 + chi1));
    const GmcDual4 tau = state.xns - eta + gmcLog(a / g2);
    const GmcDual4 x0 = psp103Sigma2(a, b, c, tau, eta);
    const GmcDual4 delta0 = state.delta_ns * gmcExp(x0);
    const GmcDual4 delta1 = gmcExp(-x0);
    const GmcDual4 chi_x0 = psp103Chi(x0);
    const GmcDual4 chi_x0_prime = psp103ChiPrime(x0);
    const GmcDual4 chi_x0_second = psp103ChiSecond(x0);
    const GmcDual4 channel0 = state.xg - x0;
    const GmcDual4 p = 2.0 * channel0 + g2 *
        (1.0 - delta1 + delta0 - state.delta_ns * (1.0 + chi_x0_prime));
    const GmcDual4 q = channel0 * channel0 - g2 *
        (delta1 + x0 - 1.0 + delta0 -
         state.delta_ns * (x0 + 1.0 + chi_x0));
    const GmcDual4 discriminant = p * p - 2.0 * q *
        (2.0 - g2 * (delta1 + delta0 - state.delta_ns * chi_x0_second));
    return x0 + 2.0 * q / (p + psp103PositiveRoot(discriminant));
}

struct Psp103IntrinsicCurrentInputs {
    double beta = 0.0;
    double f_delta_l = 1.0;
    double q_im_star = 0.0;
    double g_vsat = 1.0;
    GmcDual4 delta_psi = GmcDual4::constant(0.0);
};

struct Psp103ChannelInputs {
    GmcDual4 xs;
    GmcDual4 xd;
    GmcDual4 ds;
    GmcDual4 dd;
    GmcDual4 vds;
    GmcDual4 vdse;
    GmcDual4 vdsx;
    GmcDual4 delta_psi;
    double phi_t_star = 0.026;
    double body_factor = 1.0;
    double eta_p = 0.0;
    double rsg = 0.0;
    double theta_r = 0.0;
    double rho_b = 0.0;
    double eeff0 = 0.0;
    double eta_mu = 0.0;
    double mu_e = 0.0;
    double theta_mu = 1.0;
    double c_s = 0.0;
    double rho_mu_x = 0.0;
    double vp = 1.0;
    double alp = 0.0;
    double alp1 = 0.0;
    double alp2 = 0.0;
    double xi_tb = 0.0;
    double theta_sat = 1.0;
    double thesatg = 0.0;
    double beta = 0.0;
    bool pmos = false;
};

struct Psp103ChannelState {
    GmcDual4 xm;
    GmcDual4 em;
    GmcDual4 dm;
    GmcDual4 pm;
    GmcDual4 xgm;
    GmcDual4 qim;
    GmcDual4 alpha_m;
    GmcDual4 qim_star;
    GmcDual4 qbm;
    GmcDual4 gmob;
    GmcDual4 g_delta_l;
    GmcDual4 f_delta_l;
    GmcDual4 g_vsat;
    GmcDual4 ids;
};

// PSP103 equations 4.169 and 4.186-4.206. Inputs are the already conditioned
// source/drain surface-potential quantities; no model defaults are invented.
inline Psp103ChannelState psp103Channel(const Psp103ChannelInputs& input) {
    Psp103ChannelState state;
    const GmcDual4 g = GmcDual4::constant(input.body_factor);
    const GmcDual4 g2 = g * g;
    const GmcDual4 phi_t = GmcDual4::constant(input.phi_t_star);
    const GmcDual4 es = gmcExp(-input.xs);
    const GmcDual4 ed = gmcExp(-input.xd);
    state.xm = 0.5 * (input.xs + input.xd);
    state.em = gmcSqrt(es * ed);
    const GmcDual4 d_bar = 0.5 * (input.ds + input.dd);
    const GmcDual4 xds = input.xd - input.xs;
    state.dm = d_bar + 0.125 * xds * xds * (state.em - 2.0 / g2);
    state.pm = state.xm - 1.0 + state.em;
    state.xgm = g * gmcSqrt(state.dm + state.pm);
    state.qim = g2 * phi_t * state.dm /
        (state.xgm + g * gmcSqrt(state.pm));
    state.alpha_m = input.eta_p +
        g * (1.0 - state.em) / (2.0 * gmcSqrt(state.pm));
    state.qim_star = state.qim + phi_t * state.alpha_m;
    state.qbm = phi_t * g * gmcSqrt(state.pm);

    const GmcDual4 rho_g = input.rsg >= 0.0
        ? 1.0 / (1.0 + input.rsg * state.qim)
        : 1.0 / (1.0 - input.rsg * state.qim);
    const GmcDual4 rho_s = input.theta_r * input.rho_b * rho_g * state.qim;
    const GmcDual4 eeff = input.eeff0 * (state.qbm + input.eta_mu * state.qim);
    state.gmob = 1.0 + gmcPow(input.mu_e * eeff, input.theta_mu) +
        input.c_s * gmcPow(state.qbm / (state.qim + state.qbm), input.theta_mu) +
        input.rho_mu_x + rho_s;

    const GmcDual4 r1 = state.qim / state.qim_star;
    const GmcDual4 r2 = phi_t * state.alpha_m / state.qim_star;
    const GmcDual4 t1 = gmcLog(
        (1.0 + (input.vds - input.delta_psi) / input.vp) /
        (1.0 + (input.vdse - input.delta_psi) / input.vp));
    const GmcDual4 t2 = gmcLog(1.0 + input.vdsx / input.vp);
    const GmcDual4 delta_l = input.alp * t1;
    state.g_delta_l = 1.0 / (1.0 + delta_l + delta_l * delta_l);
    const GmcDual4 delta_l1 =
        (input.alp + input.alp1 / state.qim_star * r1) * t1 +
        input.alp2 * state.qbm * r2 * r2 * t2;
    state.f_delta_l = (1.0 + delta_l1 + delta_l1 * delta_l1) * state.g_delta_l;

    const GmcDual4 wsat = 100.0 * state.qim * input.xi_tb /
        (100.0 + state.qim * input.xi_tb);
    const GmcDual4 theta_sat = input.thesatg >= 0.0
        ? input.theta_sat / (state.gmob * state.g_delta_l) * (1.0 + input.thesatg * wsat)
        : input.theta_sat / (state.gmob * state.g_delta_l) /
            (1.0 - input.thesatg * wsat);
    const GmcDual4 zsat = input.pmos
        ? (theta_sat * input.delta_psi) * (theta_sat * input.delta_psi) /
            (1.0 + theta_sat * input.delta_psi)
        : (theta_sat * input.delta_psi) * (theta_sat * input.delta_psi);
    state.g_vsat = state.gmob * state.g_delta_l * 0.5 *
        (1.0 + gmcSqrt(1.0 + 2.0 * zsat));
    state.ids = input.beta * state.f_delta_l * state.qim_star /
        state.g_vsat * input.delta_psi;
    return state;
}

struct Psp103ChargeInputs {
    GmcDual4 voxm;
    GmcDual4 qim;
    GmcDual4 alpha_m;
    GmcDual4 delta_psi;
    GmcDual4 g_delta_l;
    GmcDual4 h;
    double c_ox_qm = 0.0;
    double eta_p = 0.0;
    bool source_partition = false;
};

struct Psp103ChargeState {
    GmcDual4 qg;
    GmcDual4 qd;
    GmcDual4 qi;
    GmcDual4 qb;
    GmcDual4 q_delta_l;
};

// PSP103 equations 4.275-4.284. The returned terminal charges conserve
// charge by construction: QG + QD + QI + QB = 0.
inline Psp103ChargeState psp103Charges(const Psp103ChargeInputs& input) {
    Psp103ChargeState state;
    const GmcDual4 half_delta = 0.5 * input.delta_psi;
    state.q_delta_l = input.source_partition
        ? GmcDual4::constant(0.0)
        : (1.0 - input.g_delta_l) * (input.qim - input.alpha_m * half_delta);
    const GmcDual4 q_delta_l = input.c_ox_qm * state.q_delta_l;
    const GmcDual4 q_delta_l_star = q_delta_l * (1.0 + input.g_delta_l);
    const GmcDual4 fj = input.delta_psi / (2.0 * input.h);
    state.qg = input.c_ox_qm *
        (input.voxm + input.eta_p * half_delta *
            (input.g_delta_l / 3.0 * fj + input.g_delta_l - 1.0));
    if (input.source_partition) {
        state.qd = -input.c_ox_qm * input.g_delta_l * input.g_delta_l * 0.5 *
            (input.qim + input.alpha_m * half_delta * (fj - 2.0));
    } else {
        state.qd = -input.c_ox_qm * input.g_delta_l * input.g_delta_l * 0.5 *
            (input.qim + input.alpha_m * input.delta_psi / 6.0 *
                (fj * fj / 5.0 + fj - 1.0)) - q_delta_l_star;
    }
    state.qi = -input.c_ox_qm * input.g_delta_l *
        (input.qim + input.alpha_m * input.delta_psi / 6.0 * fj) - q_delta_l;
    state.qb = -(state.qg + state.qd + state.qi);
    return state;
}

// PSP103 equation 4.214 for the intrinsic drain-source channel current.
// Only the active-channel branch is evaluated; callers must provide the
// validated auxiliary quantities produced by the preceding PSP blocks.
inline GmcDual4 psp103IntrinsicDrainCurrent(
    const GmcDual4& xg_dc,
    const Psp103IntrinsicCurrentInputs& inputs) {
    if (xg_dc.value <= 0.0) return GmcDual4::constant(0.0);
    const double scale = inputs.beta * inputs.f_delta_l * inputs.q_im_star /
                         inputs.g_vsat;
    return inputs.delta_psi * scale;
}

// ---------------------------------------------------------------------------
// PSP103.4.0 SP calc DC core (PSP103_SPCalculation.include, faithful port).
// The sp_s / sp_s_d surface-potential macros are implemented by
// psp103SourceSurfacePotential / psp103DrainSurfacePotential above; the
// surrounding bias definition, DIBL/SCE, Vdsat/Vdse, drain-side solve, CLM,
// velocity saturation and Voxm/H block is assembled here. Voltages vgs, vds,
// vsb are already polarity-converted and source-drain interchanged by the
// caller, matching the module's SPcalc_dc context.
// ---------------------------------------------------------------------------

// PSP103 constants (Common103_macrodefs.include / PSP103_macrodefs.include)
namespace psp103_const {
constexpr double se = 4.6051701859880916e+02;
constexpr double se05 = 2.3025850929940458e+02;
constexpr double ke = 1.0e-200;
constexpr double ke05 = 1.0e-100;
constexpr double ke05inv = 1.0e100;
constexpr double oneThird = 3.3333333333333333e-01;
constexpr double oneSixth = 1.6666666666666667e-01;
constexpr double invSqrt2 = 7.0710678118654746e-01;
constexpr double vdsat_lim_factor = 3.912023005;
} // namespace psp103_const

// P3: 3rd order polynomial expansion of exp()
inline double psp103P3(double u) {
    return 1.0 + u * (1.0 + 0.5 * (u * (1.0 + u * psp103_const::oneThird)));
}

inline GmcDual4 psp103P3(const GmcDual4& u) {
    return 1.0 + u * (1.0 + 0.5 * (u * (1.0 + u * psp103_const::oneThird)));
}

// expl_low: exp() with 3rd order polynomial extrapolation for very low values
inline GmcDual4 psp103ExplLow(const GmcDual4& x) {
    if (x.value > -psp103_const::se05) return gmcExp(x);
    return psp103_const::ke05 / psp103P3(-psp103_const::se05 - x);
}

// expl_high: exp() with 3rd order polynomial extrapolation for very high values
inline GmcDual4 psp103ExplHigh(const GmcDual4& x) {
    if (x.value < psp103_const::se05) return gmcExp(x);
    return psp103_const::ke05inv * psp103P3(x - psp103_const::se05);
}

struct Psp103DcCoreInputs {
    // Temperature-scaled core quantities (TempScaling)
    double phib = 0.9;       // phib_dc
    double g0 = 1.0;         // G_0_dc
    double phit0 = 0.0259;   // phit0 = phit * (1 + CT * rTn)
    double vfb_t = 0.0;      // VFB_T
    double kp = 0.0;         // polysilicon depletion
    // Local (geometry scaled) parameters
    double cf = 0.0;         // CF_p (DIBL)
    double cfd = 0.0;        // CFD_p
    double cfb = 0.0;        // CFB_p
    double psce = 0.0;       // PSCE_p (SCE)
    double psce_d = 0.0;     // PSCED_p
    double psce_b = 0.0;     // PSCEB_p
    double dnsub = 0.0;      // DNSUB_p
    double vnsub = 0.0;      // VNSUB_p
    double nslp = 0.0;       // NSLP_p
    double rsb = 0.0;        // RSB_p
    double rsg = 0.0;        // RSG_p
    double ther = 0.0;       // THER_i = 2*BET_i*RS_T
    double e_eff0 = 0.0;     // E_eff0
    double eta_mu = 0.5;     // eta_mu
    double eta_mu1 = 0.5;    // eta_mu1
    double mue_t = 0.0;      // MUE_T
    double themu_t = 1.5;    // THEMU_T
    double cs_t = 0.0;       // CS_T
    double xcor_t = 0.0;     // XCOR_T
    double thesat_b = 0.0;   // THESATB_p
    double thesat_g = 0.0;   // THESATG_p
    double theta_sat_t = 1.0; // THESAT_T
    double ax = 18.0;        // AX_p (linear/saturation transition)
    double alp = 0.0;        // ALP_p (CLM)
    double alp1 = 0.0;       // ALP1_p
    double alp2 = 0.0;       // ALP2_p
    double vp = 1.0;         // VP_p
    double bet = 0.0;        // BET_i = FACTUO*BETN_T*CoxPrime
    // Terminal conditioning (SWGEO=1 values)
    double phix = 0.0;       // phix_dc
    double aphi = 0.0;       // aphi_dc == bphi_dc
    double phix1 = 0.0;      // phix1_dc
    // NUD effect (SWNUD)
    double vsbnud = 0.0;
    double dvsbnud = 1.0;
    double gfac_nud = 1.0;
    double us1 = 0.0;
    double us21 = 0.0;
    bool sw_nud = false;
    bool pmos = false;
};

struct Psp103DcCoreResult {
    GmcDual4 xg;        // xg_dc
    GmcDual4 ids;       // BET_i * FdL_dc * qim1_dc * dps_dc * Gvsatinv_dc
    GmcDual4 qim;
    GmcDual4 qim1;
    GmcDual4 qbm;
    GmcDual4 alpha;
    GmcDual4 dps;
    GmcDual4 voxm;
    GmcDual4 g_delta_l; // GdL_dc
    GmcDual4 f_delta_l; // FdL_dc
    GmcDual4 g_vsat;
    GmcDual4 g_vsatinv;
    GmcDual4 h;         // H_dc
    GmcDual4 eta_p;
    GmcDual4 qeff1;
    GmcDual4 vdsat;
    GmcDual4 udse;
    GmcDual4 x_ds;
    GmcDual4 x_m;
    GmcDual4 gf;        // Gf_dc (body factor)
};

inline Psp103DcCoreResult psp103DcCore(
    const Psp103DcCoreInputs& in, GmcDual4 vgs, GmcDual4 vds, GmcDual4 vsb) {
    Psp103DcCoreResult result;

    // Vdsx smoothing (module line 2076)
    const GmcDual4 vds2 = vds * vds;
    const GmcDual4 vdsx = vds2 / (gmcSqrt(vds2 + 0.01) + 0.1);

    // Conditioning of terminal voltages (SPcalc_dc lines 2086-2087)
    GmcDual4 temp = psp103MinA(vds + vsb, vsb, in.aphi) + in.phix;
    GmcDual4 vsbstar = vsb - psp103MinA(temp, GmcDual4::constant(0.0), in.aphi) +
        GmcDual4::constant(in.phix1);
    const GmcDual4 vsbstar_dc_tmp = vsbstar;

    // NUD effect (lines 2091-2098)
    if (in.sw_nud && in.gfac_nud != 1.0) {
        const GmcDual4 vmb = vsbstar + 0.5 * (vds - vdsx);
        const GmcDual4 us = gmcSqrt(vmb + in.phib) -
            GmcDual4::constant(std::sqrt(in.phib));
        temp = 2.0 * (us - in.us1) / in.us21 - 1.0;
        const GmcDual4 usnew = us - 0.25 * (1.0 - in.gfac_nud) * in.us21 *
            (temp + gmcSqrt(temp * temp + 0.4804530139182));
        const GmcDual4 vmbnew = usnew * usnew + 2.0 * std::sqrt(in.phib) * usnew;
        vsbstar = vmbnew - 0.5 * (vds - vdsx);
    }

    // Bias definition (lines 51-65)
    const GmcDual4 vgbstar = vgs + vsbstar;
    GmcDual4 vgb1 = vgbstar - GmcDual4::constant(in.vfb_t);
    const GmcDual4 vsbx = vsbstar + 0.5 * (vds - vdsx);
    GmcDual4 vdsp = vdsx;
    if (in.cfd >= 1.0e-10) {
        vdsp = 2.0 * (gmcSqrt(1.0 + in.cfd * vdsx) - 1.0) / in.cfd;
    }
    const GmcDual4 delvg = in.cf * vdsp * (1.0 + in.cfb * vsbx); // DIBL
    const GmcDual4 dphit1 = in.psce * (1.0 + in.psce_d * vdsx) *
        (1.0 + in.psce_b * vsbx); // SCE on subthreshold slope
    const GmcDual4 phit1 = in.phit0 * (1.0 + dphit1);
    const GmcDual4 inv_phit1 = 1.0 / phit1;
    vgb1 = vgb1 + delvg;
    const GmcDual4 xg = vgb1 * inv_phit1;

    // Bias dependent body factor (lines 67-73)
    GmcDual4 gf = GmcDual4::constant(in.g0);
    if (in.dnsub > 0.0) {
        const GmcDual4 dnsub = in.dnsub * psp103MaxA(
            GmcDual4::constant(0.0), vgs + vsb - in.vnsub, in.nslp);
        gf = in.g0 * gmcSqrt(1.0 + dnsub);
    }
    const GmcDual4 g2 = gf * gf;
    const GmcDual4 inv_g2 = 1.0 / g2;
    const GmcDual4 xi = 1.0 + gf * psp103_const::invSqrt2;
    const GmcDual4 inv_xi = 1.0 / xi;
    const GmcDual4 ux = vsbstar * inv_phit1;
    const GmcDual4 xn_s = GmcDual4::constant(in.phib) * inv_phit1 + ux;
    GmcDual4 delta_ns;
    if (xn_s.value < psp103_const::se) {
        delta_ns = gmcExp(-xn_s);
    } else {
        delta_ns = psp103_const::ke / psp103P3(xn_s - psp103_const::se);
    }
    const double margin = 1.0e-5 * xi.value;

    // Source surface potential (sp_s)
    Psp103SurfacePotentialState state;
    state.xi = xi;
    state.xns = xn_s;
    state.delta_ns = delta_ns;
    state.xg = xg;
    GmcDual4 x_s = psp103SourceSurfacePotential(state, gf);
    GmcDual4 x_d = x_s;
    GmcDual4 x_m = x_s;
    GmcDual4 x_ds = GmcDual4::constant(0.0);
    (void)margin;

    // Core PSP current calculation (lines 93-330)
    const GmcDual4 vdsat_lim = psp103_const::vdsat_lim_factor * phit1;
    GmcDual4 voxm = GmcDual4::constant(0.0);
    GmcDual4 qeff1 = GmcDual4::constant(0.0);
    GmcDual4 vdsat = vdsat_lim;
    GmcDual4 vdse = vds;
    GmcDual4 alpha = GmcDual4::constant(1.0);
    GmcDual4 sqm = GmcDual4::constant(0.0);
    GmcDual4 qim = GmcDual4::constant(0.0);
    GmcDual4 qim1 = GmcDual4::constant(0.0);
    GmcDual4 qbm = GmcDual4::constant(0.0);
    GmcDual4 eta_p = GmcDual4::constant(1.0);
    GmcDual4 g_delta_l = GmcDual4::constant(1.0);
    GmcDual4 g_vsat = GmcDual4::constant(1.0);
    GmcDual4 g_vsatinv = GmcDual4::constant(1.0);
    GmcDual4 h = GmcDual4::constant(1.0);
    GmcDual4 xgm = GmcDual4::constant(0.0);
    GmcDual4 dps = GmcDual4::constant(0.0);
    GmcDual4 s1 = GmcDual4::constant(0.0);
    GmcDual4 d_l = GmcDual4::constant(0.0);
    GmcDual4 udse = GmcDual4::constant(0.0);

    if (xg.value > 0.0) {
        // delta_1s / Es / Ds / Ps / alpha at source (lines 103-134)
        GmcDual4 delta_1s = GmcDual4::constant(0.0);
        temp = 1.0 / (2.0 + x_s * x_s);
        const GmcDual4 xi0s = x_s * x_s * temp;
        const GmcDual4 xi1s = 4.0 * (x_s * temp * temp);
        const GmcDual4 xi2s = (8.0 * temp - 12.0 * xi0s) * temp * temp;
        GmcDual4 es;
        if (x_s.value < psp103_const::se05) {
            delta_1s = gmcExp(x_s);
            es = 1.0 / delta_1s;
            delta_1s = delta_ns * delta_1s;
        } else if (x_s.value > (xn_s.value - psp103_const::se05)) {
            delta_1s = gmcExp(x_s - xn_s);
            es = delta_ns / delta_1s;
        } else {
            delta_1s = psp103_const::ke05 / psp103P3(xn_s - x_s - psp103_const::se05);
            es = psp103_const::ke05 / psp103P3(x_s - psp103_const::se05);
        }
        GmcDual4 ds = delta_1s - delta_ns * (x_s + 1.0 + xi0s);
        GmcDual4 ps;
        if (x_s.value < 1.0e-5) {
            ps = 0.5 * (x_s * x_s * (1.0 - psp103_const::oneThird *
                (x_s * (1.0 - 0.25 * x_s))));
            ds = psp103_const::oneSixth * (delta_ns * x_s * x_s * x_s *
                (1.0 + 1.75 * x_s));
            temp = gmcSqrt(1.0 - psp103_const::oneThird *
                (x_s * (1.0 - 0.25 * x_s)));
            sqm = psp103_const::invSqrt2 * (x_s * temp);
            alpha = 1.0 + gf * psp103_const::invSqrt2 *
                (1.0 - 0.5 * x_s + psp103_const::oneSixth * (x_s * x_s)) / temp;
        } else {
            ps = x_s - 1.0 + es;
            sqm = gmcSqrt(ps);
            alpha = 1.0 + 0.5 * (gf * (1.0 - es) / sqm);
        }
        GmcDual4 em = es;
        const GmcDual4 ed_src = em;
        GmcDual4 dm = ds;
        const GmcDual4 dd_src = dm;

        // Drain saturation voltage (lines 137-185)
        const GmcDual4 rxcor = (1.0 + 0.2 * in.xcor_t * vsbx) /
            (1.0 + in.xcor_t * vsbx);
        GmcDual4 gmob = GmcDual4::constant(1.0);
        GmcDual4 xitsb = GmcDual4::constant(1.0);
        GmcDual4 thesat1 = GmcDual4::constant(1.0);
        GmcDual4 rhob = GmcDual4::constant(0.0);
        if (ds.value > psp103_const::ke05) {
            const GmcDual4 xgs = gf * gmcSqrt(ps + ds);
            const GmcDual4 qis = g2 * ds * phit1 / (xgs + gf * sqm);
            const GmcDual4 qbs = sqm * gf * phit1;
            if (in.rsb < 0.0) {
                rhob = 1.0 / (1.0 - in.rsb * vsbx);
            } else {
                rhob = 1.0 + in.rsb * vsbx;
            }
            GmcDual4 gr_mob;
            if (in.rsg < 0.0) {
                gr_mob = 1.0 - in.rsg * qis;
            } else {
                gr_mob = 1.0 / (1.0 + in.rsg * qis);
            }
            const GmcDual4 gr = in.ther * (rhob * gr_mob * qis);
            const GmcDual4 eeffm = in.e_eff0 * (qbs + in.eta_mu * qis);
            const GmcDual4 mutmp = gmcPow(eeffm * in.mue_t, in.themu_t) +
                in.cs_t * (ps / (ps + ds + 1.0e-14));
            gmob = (1.0 + mutmp + gr) * rxcor;
            if (in.thesat_b < 0.0) {
                xitsb = 1.0 / (1.0 - in.thesat_b * vsbx);
            } else {
                xitsb = 1.0 + in.thesat_b * vsbx;
            }
            const GmcDual4 temp2 = qis * xitsb;
            const GmcDual4 wsat = 100.0 * (temp2 / (100.0 + temp2));
            GmcDual4 vsat_mob;
            if (in.thesat_g < 0.0) {
                vsat_mob = 1.0 / (1.0 - in.thesat_g * wsat);
            } else {
                vsat_mob = 1.0 + in.thesat_g * wsat;
            }
            thesat1 = in.theta_sat_t * (vsat_mob / gmob);
            const GmcDual4 phi_inf = qis / alpha + phit1;
            GmcDual4 ysat = thesat1 * phi_inf * psp103_const::invSqrt2;
            if (in.pmos) {
                ysat = ysat / gmcSqrt(1.0 + ysat);
            }
            const GmcDual4 za = 2.0 / (1.0 + gmcSqrt(1.0 + 4.0 * ysat));
            const GmcDual4 temp1 = za * ysat;
            const GmcDual4 phi_0 = phi_inf * za * (1.0 + 0.86 *
                (temp1 * (1.0 - temp1 * za) / (1.0 + 4.0 * (temp1 * temp1 * za))));
            const GmcDual4 asat = xgs + 0.5 * g2;
            const GmcDual4 phi_2 = 0.98 * (g2 * ds * phit1 /
                (asat + gmcSqrt(asat * asat - g2 * ds * 0.98)));
            const GmcDual4 phi_0_2 = phi_0 + phi_2;
            const GmcDual4 phi0_phi2 = 2.0 * (phi_0 * phi_2);
            const GmcDual4 phi_sat = phi0_phi2 / (phi_0_2 +
                gmcSqrt(phi_0_2 * phi_0_2 - 1.98 * phi0_phi2));
            vdsat = phi_sat - phit1 * gmcLog(1.0 + phi_sat *
                (phi_sat - 2.0 * asat * phit1) * inv_g2 / (phit1 * phit1 * ds));
        }
        temp = gmcPow(vds / vdsat, in.ax);
        vdse = vds * gmcPow(1.0 + temp, -1.0 / in.ax);

        // Surface potential at drain side (lines 190-231)
        udse = vdse * inv_phit1;
        const GmcDual4 xn_d = xn_s + udse;
        GmcDual4 k_ds;
        if (udse.value < psp103_const::se) {
            k_ds = gmcExp(-udse);
        } else {
            k_ds = psp103_const::ke / psp103P3(udse - psp103_const::se);
        }
        const GmcDual4 delta_nd = delta_ns * k_ds;

        Psp103SurfacePotentialState drain_state;
        drain_state.xi = xi;
        drain_state.xns = xn_d;
        drain_state.delta_ns = delta_nd;
        drain_state.xg = xg;
        x_d = psp103DrainSurfacePotential(drain_state, gf);
        x_ds = x_d - x_s;

        // Approximation for extremely small x_ds (lines 203-210)
        if (x_ds.value < 1.0e-10) {
            const GmcDual4 p_c = 2.0 * (xg - x_s) + g2 * (1.0 - es +
                delta_1s * k_ds - delta_nd * (1.0 + xi1s));
            const GmcDual4 q_c = g2 * (1.0 - k_ds) * ds;
            temp = 2.0 - g2 * (es + delta_1s * k_ds - delta_nd * xi2s);
            temp = p_c * p_c - 2.0 * (temp * q_c);
            x_ds = 2.0 * (q_c / (p_c + gmcSqrt(temp)));
            x_d = x_s + x_ds;
        }
        // Dd / Ed at drain (lines 213-231)
        dps = x_ds * phit1;
        const GmcDual4 xi0d = x_d * x_d / (2.0 + x_d * x_d);
        GmcDual4 ed;
        GmcDual4 dd;
        if (x_d.value < psp103_const::se05) {
            ed = gmcExp(-x_d);
            if (x_d.value < 1.0e-5) {
                dd = psp103_const::oneSixth * delta_nd * x_d * x_d * x_d *
                    (1.0 + 1.75 * x_d);
            } else {
                dd = delta_nd * (1.0 / ed - x_d - 1.0 - xi0d);
            }
        } else {
            if (x_d.value > (xn_d.value - psp103_const::se05)) {
                temp = gmcExp(x_d - xn_d);
                ed = delta_nd / temp;
                dd = temp - delta_nd * (x_d + 1.0 + xi0d);
            } else {
                ed = psp103_const::ke05 / psp103P3(x_d - psp103_const::se05);
                temp = psp103_const::ke05 /
                    psp103P3(xn_d - x_d - psp103_const::se05);
                dd = temp - delta_nd * (x_d + 1.0 + xi0d);
            }
        }

        // Mid-point surface potential (lines 234-280)
        x_m = 0.5 * (x_s + x_d);
        em = GmcDual4::constant(0.0);
        temp = ed * es;
        if (temp.value > 0.0) {
            em = gmcSqrt(temp);
        }
        const GmcDual4 d_bar = 0.5 * (ds + dd);
        dm = d_bar + 0.125 * (x_ds * x_ds * (em - 2.0 * inv_g2));
        GmcDual4 pm = GmcDual4::constant(0.0);

        if (x_m.value < 1.0e-5) {
            pm = 0.5 * (x_m * x_m * (1.0 -
                psp103_const::oneThird * (x_m * (1.0 - 0.25 * x_m))));
            xgm = gf * gmcSqrt(dm + pm);
            if (in.kp > 0.0) {
                eta_p = 1.0 / gmcSqrt(1.0 + in.kp * xgm);
            }
            temp = gmcSqrt(1.0 - psp103_const::oneThird *
                (x_m * (1.0 - 0.25 * x_m)));
            sqm = psp103_const::invSqrt2 * (x_m * temp);
            alpha = eta_p + psp103_const::invSqrt2 *
                (gf * (1.0 - 0.5 * x_m +
                 psp103_const::oneSixth * (x_m * x_m)) / temp);
        } else {
            pm = x_m - 1.0 + em;
            xgm = gf * gmcSqrt(dm + pm);
            if (in.kp > 0.0) {
                const GmcDual4 d0 = 1.0 - em + 2.0 * (xgm * inv_g2);
                eta_p = 1.0 / gmcSqrt(1.0 + in.kp * xgm);
                temp = eta_p / (eta_p + 1.0);
                const GmcDual4 x_pm = in.kp * (temp * temp * g2 * dm);
                const GmcDual4 p_pd = 2.0 * (xgm - x_pm) + g2 * (1.0 - em + dm);
                const GmcDual4 q_pd = x_pm * (x_pm - 2.0 * xgm);
                const GmcDual4 xi_pd = 1.0 - 0.5 * (g2 * (em + dm));
                const GmcDual4 u_pd = q_pd * p_pd /
                    (p_pd * p_pd - xi_pd * q_pd);
                x_m = x_m + u_pd;
                const GmcDual4 km = gmcExp(u_pd);
                em = em / km;
                dm = dm * km;
                pm = x_m - 1.0 + em;
                xgm = gf * gmcSqrt(dm + pm);
                const GmcDual4 km0 = 1.0 - em + 2.0 * (xgm * eta_p * inv_g2);
                x_ds = x_ds * km * (d0 + d_bar) / (km0 + km * d_bar);
                dps = x_ds * phit1;
            }
            sqm = gmcSqrt(pm);
            alpha = eta_p + 0.5 * (gf * (1.0 - em) / sqm);
        }

        // Potential midpoint inversion charge (lines 283-285)
        qim = phit1 * (g2 * dm / (xgm + gf * sqm));
        qim1 = qim + phit1 * alpha;
        qbm = sqm * gf * phit1;

        // Series resistance (lines 288-293)
        GmcDual4 gr_mob;
        if (in.rsg < 0.0) {
            gr_mob = 1.0 - in.rsg * qim;
        } else {
            gr_mob = 1.0 / (1.0 + in.rsg * qim);
        }
        const GmcDual4 gr = in.ther * (rhob * gr_mob * qim);

        // Mobility reduction (lines 296-300)
        const GmcDual4 qeff = qbm + in.eta_mu * qim;
        qeff1 = qbm + in.eta_mu1 * qim;
        const GmcDual4 eeffm = in.e_eff0 * qeff;
        const GmcDual4 mutmp = gmcPow(eeffm * in.mue_t, in.themu_t) +
            in.cs_t * (pm / (pm + dm + 1.0e-14));
        gmob = (1.0 + mutmp + gr) * rxcor;

        // Channel length modulation (lines 303-305)
        s1 = gmcLog((1.0 + (vds - dps) / in.vp) /
            (1.0 + (vdse - dps) / in.vp));
        d_l = in.alp * s1;
        g_delta_l = 1.0 / (1.0 + d_l + d_l * d_l);

        // Velocity saturation (lines 308-322)
        const GmcDual4 temp2 = qim * xitsb;
        const GmcDual4 wsat = 100.0 * (temp2 / (100.0 + temp2));
        const GmcDual4 gmob_d_l = gmob * g_delta_l;
        GmcDual4 vsat_mob;
        if (in.thesat_g < 0.0) {
            vsat_mob = 1.0 / (1.0 - in.thesat_g * wsat);
        } else {
            vsat_mob = 1.0 + in.thesat_g * wsat;
        }
        thesat1 = in.theta_sat_t * (vsat_mob / gmob_d_l);
        GmcDual4 zsat = thesat1 * thesat1 * dps * dps;
        if (in.pmos) {
            zsat = zsat / (1.0 + thesat1 * dps);
        }
        g_vsat = 0.5 * (gmob_d_l * (1.0 + gmcSqrt(1.0 + 2.0 * zsat)));
        g_vsatinv = 1.0 / g_vsat;

        // Variables for intrinsic charges and gate current (lines 325-329)
        voxm = xgm * phit1;
        temp = gmob_d_l * g_vsatinv;
        const GmcDual4 alpha1 = alpha * (1.0 + 0.5 * (zsat * temp * temp));
        h = temp * qim1 / alpha1;
    } else {
        // xg <= 0: accumulation/depletion branch (lines 95-101)
        const GmcDual4 xgm_neg = xg - x_s;
        voxm = xgm_neg * phit1;
        qeff1 = voxm;
        vdsat = vdsat_lim;
        vdse = vds;
        g_delta_l = GmcDual4::constant(1.0);
        g_vsat = GmcDual4::constant(1.0);
        g_vsatinv = GmcDual4::constant(1.0);
        h = GmcDual4::constant(1.0);
    }

    result.xg = xg;
    result.gf = gf;
    result.qim = qim;
    result.qim1 = qim1;
    result.qbm = qbm;
    result.alpha = alpha;
    result.dps = dps;
    result.voxm = voxm;
    result.g_delta_l = g_delta_l;
    result.g_vsat = g_vsat;
    result.g_vsatinv = g_vsatinv;
    result.h = h;
    result.eta_p = eta_p;
    result.qeff1 = qeff1;
    result.vdsat = vdsat;
    result.udse = udse;
    result.x_ds = x_ds;
    result.x_m = x_m;

    // FdL (module lines 2108-2115)
    if (xg.value > 0.0) {
        const GmcDual4 qim1_1 = 1.0 / qim1;
        const GmcDual4 r1 = qim * qim1_1;
        const GmcDual4 r2 = phit1 * (alpha * qim1_1);
        const GmcDual4 s2 = gmcLog(1.0 + vdsx * (1.0 / in.vp));
        const GmcDual4 d_l1 = d_l + in.alp1 * (qim1_1 * r1 * s1) +
            in.alp2 * (qbm * r2 * r2 * s2);
        result.f_delta_l = (1.0 + d_l1 + d_l1 * d_l1) * g_delta_l;
        result.ids = in.bet * result.f_delta_l * qim1 * dps * g_vsatinv;
    } else {
        result.f_delta_l = GmcDual4::constant(1.0);
        result.ids = GmcDual4::constant(0.0);
    }
    return result;
}

} // namespace gspice

#endif // GSPICE_PSP103_CORE_HPP
