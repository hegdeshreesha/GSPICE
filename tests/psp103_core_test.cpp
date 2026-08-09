#include "devices/psp103_core.hpp"

#include <cassert>
#include <cmath>

int main() {
    using gspice::GmcDual4;

    const auto vds = GmcDual4::variable(0.5, 2);
    const auto conditioned = gspice::psp103ConditionTerminals(
        vds, GmcDual4::constant(0.1), 0.2, 0.05, 0.1, 0.02);
    assert(conditioned.vdsx.value > 0.1);
    assert(std::isfinite(conditioned.vsb_star.value));
    assert(std::isfinite(conditioned.vdsx.derivative[2]));

    const auto min_a = gspice::psp103MinA(
        GmcDual4::constant(1.0), GmcDual4::constant(2.0), 1e-8);
    const auto max_a = gspice::psp103MaxA(
        GmcDual4::constant(1.0), GmcDual4::constant(2.0), 1e-8);
    assert(std::abs(min_a.value - 1.0) < 1e-4);
    assert(std::abs(max_a.value - 2.0) < 1e-4);

    const auto y = GmcDual4::variable(0.4, 3);
    const auto chi = gspice::psp103Chi(y);
    const auto chi_prime = gspice::psp103ChiPrime(y);
    const auto chi_second = gspice::psp103ChiSecond(y);
    assert(std::isfinite(chi.value));
    assert(std::isfinite(chi_prime.value));
    assert(std::isfinite(chi_second.value));
    assert(chi.derivative[3] > 0.0);

    const auto sigma1 = gspice::psp103Sigma1(
        GmcDual4::constant(2.0), GmcDual4::constant(0.4),
        GmcDual4::constant(0.2), GmcDual4::constant(-0.1));
    const auto sigma2 = gspice::psp103Sigma2(
        GmcDual4::constant(2.0), GmcDual4::constant(0.8),
        GmcDual4::constant(0.4), GmcDual4::constant(0.2),
        GmcDual4::constant(-0.1));
    assert(std::isfinite(sigma1.value));
    assert(std::isfinite(sigma2.value));

    const auto vgb = GmcDual4::variable(0.9, 0);
    const auto vsb = GmcDual4::variable(0.1, 1);
    const auto surface = gspice::psp103SourceSide(vgb, vsb, 0.7, 0.026, 0.8);
    assert(std::abs(surface.xg.value - (0.9 / 0.026)) < 1e-12);
    assert(std::abs(surface.xns.derivative[1] - 1.0 / 0.026) < 1e-12);
    assert(std::abs(surface.xg.derivative[0] - 1.0 / 0.026) < 1e-12);

    const auto source_potential = gspice::psp103SourceSurfacePotential(
        surface, GmcDual4::constant(0.8));
    assert(std::isfinite(source_potential.value));
    for (double derivative : source_potential.derivative) assert(std::isfinite(derivative));

    const auto drain_potential = gspice::psp103DrainSurfacePotential(
        surface, GmcDual4::constant(0.8));
    assert(std::isfinite(drain_potential.value));
    for (double derivative : drain_potential.derivative) assert(std::isfinite(derivative));

    for (double gate_bias : {-0.5, -1.0e-8, 1.0e-8, 0.5}) {
        const auto branch_state = gspice::psp103SourceSide(
            GmcDual4::constant(gate_bias), GmcDual4::constant(0.1),
            0.7, 0.026, 0.8);
        const auto potential = gspice::psp103SourceSurfacePotential(
            branch_state, GmcDual4::constant(0.8));
        assert(std::isfinite(potential.value));
    }

    gspice::Psp103IntrinsicCurrentInputs inputs;
    inputs.beta = 2.0e-3;
    inputs.f_delta_l = 0.9;
    inputs.q_im_star = 0.8;
    inputs.g_vsat = 1.25;
    inputs.delta_psi = GmcDual4::variable(0.3, 2);
    const auto ids = gspice::psp103IntrinsicDrainCurrent(surface.xg, inputs);
    const double scale = inputs.beta * inputs.f_delta_l * inputs.q_im_star / inputs.g_vsat;
    assert(std::abs(ids.value - 0.3 * scale) < 1e-15);
    assert(std::abs(ids.derivative[2] - scale) < 1e-15);

    const auto cutoff = gspice::psp103IntrinsicDrainCurrent(
        GmcDual4::constant(-1.0), inputs);
    assert(cutoff.value == 0.0);

    gspice::Psp103ChannelInputs channel;
    channel.xs = GmcDual4::constant(1.1);
    channel.xd = GmcDual4::constant(1.4);
    channel.ds = GmcDual4::constant(0.2);
    channel.dd = GmcDual4::constant(0.1);
    channel.vds = GmcDual4::variable(0.8, 0);
    channel.vdse = GmcDual4::variable(0.6, 0);
    channel.vdsx = GmcDual4::variable(0.7, 0);
    channel.delta_psi = GmcDual4::variable(0.25, 1);
    channel.body_factor = 0.8;
    channel.phi_t_star = 0.026;
    channel.eta_p = 0.1;
    channel.eeff0 = 1.0;
    channel.mu_e = 0.2;
    channel.theta_mu = 1.0;
    channel.vp = 1.0;
    channel.xi_tb = 0.5;
    channel.theta_sat = 1.0;
    channel.beta = 2.0e-3;
    const auto channel_state = gspice::psp103Channel(channel);
    assert(std::isfinite(channel_state.ids.value));
    assert(std::isfinite(channel_state.ids.derivative[0]));
    assert(std::isfinite(channel_state.ids.derivative[1]));

    gspice::Psp103ChargeInputs charge;
    charge.voxm = GmcDual4::variable(0.2, 2);
    charge.qim = GmcDual4::constant(0.5);
    charge.alpha_m = GmcDual4::constant(0.3);
    charge.delta_psi = GmcDual4::variable(0.25, 3);
    charge.g_delta_l = GmcDual4::constant(0.9);
    charge.h = GmcDual4::constant(1.1);
    charge.c_ox_qm = 1.5;
    charge.eta_p = 0.2;
    const auto charges = gspice::psp103Charges(charge);
    assert(std::isfinite(charges.qg.value));
    assert(std::isfinite(charges.qd.value));
    assert(std::isfinite(charges.qi.value));
    assert(std::abs(charges.qg.value + charges.qd.value +
                    charges.qi.value + charges.qb.value) < 1e-12);

    gspice::Psp103DcCoreInputs dc;
    dc.phib = 0.867;
    dc.g0 = 0.851;
    dc.phit0 = 0.02585;
    dc.vfb_t = -0.35;
    dc.kp = 0.0;
    dc.cf = 0.18;
    dc.cfd = 0.0;
    dc.cfb = 0.0;
    dc.psce = 0.0;
    dc.e_eff0 = 4.5e7;
    dc.eta_mu = 0.5;
    dc.eta_mu1 = 0.5;
    dc.mue_t = 1.0;
    dc.themu_t = 2.0;
    dc.cs_t = 1.0;
    dc.ther = 0.0;
    dc.thesat_b = 0.0;
    dc.thesat_g = 0.0;
    dc.theta_sat_t = 1.0;
    dc.ax = 18.0;
    dc.alp = 1.5;
    dc.alp1 = 0.0;
    dc.alp2 = 0.0;
    dc.vp = 3.0;
    dc.bet = 2.0e-3;
    dc.aphi = 0.0025 * dc.phib * dc.phib;

    const auto dc_off = gspice::psp103DcCore(
        dc, GmcDual4::constant(-1.0), GmcDual4::constant(0.5),
        GmcDual4::constant(0.0));
    assert(dc_off.ids.value == 0.0);
    assert(dc_off.xg.value <= 0.0);

    const auto dc_on = gspice::psp103DcCore(
        dc, GmcDual4::variable(1.2, 0), GmcDual4::variable(0.5, 2),
        GmcDual4::constant(0.0));
    assert(std::isfinite(dc_on.ids.value));
    for (double derivative : dc_on.ids.derivative) assert(std::isfinite(derivative));
    assert(dc_on.ids.value > 0.0);
    assert(dc_on.qim1.value > 0.0);
    assert(dc_on.dps.value > 0.0);
    assert(dc_on.g_vsatinv.value > 0.0);
    assert(std::abs(dc_on.f_delta_l.value - 1.0) < 0.5);

    double prev = -1.0;
    for (double vg : {0.2, 0.4, 0.6, 0.8, 1.0, 1.2}) {
        const auto r = gspice::psp103DcCore(
            dc, GmcDual4::constant(vg), GmcDual4::constant(0.5),
            GmcDual4::constant(0.0));
        assert(std::isfinite(r.ids.value));
        assert(r.ids.value >= prev - 1e-9);
        prev = r.ids.value;
    }
}
