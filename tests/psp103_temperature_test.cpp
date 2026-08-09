#include "devices/psp103_temperature.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

namespace {

// Reference NMOS model card: C:\EDA\_psp_ref_work\psp103_nmos-2.mod (SWGEO=1)
gspice::Psp103ParameterSet refNmosModel() {
    using gspice::Psp103ParameterSet;
    return Psp103ParameterSet::from({
        {"TYPE", 1.0},
        {"TR", 27.0}, {"DTA", 0.0}, {"SWGEO", 1.0}, {"QMC", 1.0},
        {"LVARO", -10e-9}, {"LAP", 10e-9}, {"WVARO", 10e-9}, {"WOT", 0.0},
        {"VFBO", -1.1}, {"STVFBO", 5e-4},
        {"TOXO", 1.5e-9}, {"EPSROXO", 3.9},
        {"NSUBO", 3e23}, {"NSUBW", 0.0}, {"WSEG", 1.5e-10},
        {"NPCK", 1e24}, {"WSEGP", 0.9e-8}, {"LPCK", 5.5e-8},
        {"FOL1", 2e-2}, {"FOL2", 5e-6},
        {"FACNEFFACO", 0.8}, {"GFACNUDO", 0.1},
        {"VSBNUDO", 0.0}, {"DVSBNUDO", 1.0},
        {"VNSUBO", 0.0}, {"NSLPO", 0.05}, {"DNSUBO", 0.0}, {"DPHIBO", 0.0},
        {"NPO", 1.5e26}, {"NPL", 10e-18},
        {"CTO", 5e-15}, {"CTL", 4e-2}, {"CTLEXP", 0.6},
        {"CFL", 3e-4}, {"CFLEXP", 2.0}, {"CFW", 5e-3}, {"CFBO", 0.3},
        {"UO", 3.5e-2},
        {"FBET1", -0.3}, {"FBET1W", 0.15}, {"LP1", 1.5e-7}, {"LP1W", -2.5e-2},
        {"FBET2", 50.0}, {"LP2", 8.5e-10},
        {"BETW1", 5e-2}, {"BETW2", -2e-2}, {"WBET", 5e-10},
        {"STBETO", 1.75}, {"STBETL", -2e-2}, {"STBETW", -2e-3}, {"STBETLW", -3e-3},
        {"MUEO", 0.6}, {"MUEW", -1.2e-2}, {"STMUEO", 0.5},
        {"THEMUO", 2.75}, {"STTHEMUO", -0.1},
        {"CSO", 1e-2}, {"CSL", 0.0}, {"CSLEXP", 1.0}, {"CSW", 0.0}, {"CSLW", 0.0},
        {"STCSO", -5.0},
        {"XCORO", 0.15}, {"XCORL", 2e-3}, {"XCORW", -3e-2}, {"XCORLW", -3.5e-3},
        {"STXCORO", 1.25},
        {"FETAO", 1.0},
        {"RSW1", 50.0}, {"RSW2", 0.0}, {"STRSO", -2.0}, {"RSBO", 0.0}, {"RSGO", 0.0},
        {"THESATO", 1e-6}, {"THESATL", 0.6}, {"THESATLEXP", 0.75},
        {"THESATW", -1e-2}, {"THESATLW", 0.0},
        {"STTHESATO", 1.5}, {"STTHESATL", -2.5e-2}, {"STTHESATW", -2e-2},
        {"STTHESATLW", -5e-3},
        {"THESATBO", 0.15}, {"THESATGO", 0.75},
        {"AXO", 20.0}, {"AXL", 0.2},
        {"ALPL", 7e-3}, {"ALPLEXP", 0.6}, {"ALPW", 0.0},
        {"ALP1L1", 2.5e-2}, {"ALP1LEXP", 0.5}, {"ALP1L2", 0.0}, {"ALP1W", 0.0},
        {"ALP2L1", 0.5}, {"ALP2LEXP", 0.5}, {"ALP2L2", 0.0}, {"ALP2W", 0.0},
        {"VPO", 0.25},
    });
}

void testRefNmosShortChannel() {
    using namespace gspice;
    const auto model = refNmosModel();
    const auto prepared = model.prepare({{"W", 1e-6}, {"L", 0.1e-6}}, 21.0);
    const auto s = psp103PrepareDevice(model, prepared, 21.0);

    // geometry
    assert(std::abs(s.cox_prime - 2.3021e-2) < 0.025e-2);
    assert(s.swgeo == 1 && s.channel_type == 1 && s.qmc == 1.0);
    assert(std::abs(s.mult_i - 1.0) < 1e-12);
    assert(s.neff > 7.3e23 && s.neff < 8.2e23);
    assert(s.neffac > 5.8e23 && s.neffac < 6.7e23);
    assert(std::abs(s.eta_mu - 0.5) < 1e-12);
    assert(std::abs(s.eta_mu1 - 0.5) < 1e-12);
    assert(std::abs(s.vsbnud) < 1e-12 && std::abs(s.dvsbnud - 1.0) < 1e-12);
    assert(s.gfac_nud > 0.05 && s.gfac_nud < 0.15);

    // temperature scaling outputs (hand-checked vs TempScaling)
    assert(s.phib_dc > 0.98 && s.phib_dc < 1.03);
    assert(s.g0_dc > 1.78 && s.g0_dc < 1.90);
    assert(s.kp > 0.010 && s.kp < 0.0115);
    // BET_i = FACTUO*BETN_T*CoxPrime (module line 361): 0.00806 for this model
    assert(s.bet_i > 0.0075 && s.bet_i < 0.0090);
    assert(s.vfb_t > -1.13 && s.vfb_t < -1.10);
    assert(s.ther > 0.7 && s.ther < 0.85);
    assert(s.qq > 0.18 && s.qq < 0.21);
    assert(s.e_eff0 > 2.1 && s.e_eff0 < 2.45);
    assert(s.phit0 > 0.029 && s.phit0 < 0.032);
    assert(s.us1 < 1e-9);
    assert(s.us21 > 0.38 && s.us21 < 0.45);
    assert(s.phib_ac > 0.95 && s.phib_ac < 1.05);
    assert(s.g0_ac > 1.55 && s.g0_ac < 1.75);

    // temperature-scaled mobility / resistance
    assert(s.mue_t > 0.55 && s.mue_t < 0.65);
    assert(s.themu_t > 2.6 && s.themu_t < 2.85);
    assert(s.cs_t > 8.5e-3 && s.cs_t < 9.5e-3);
    assert(s.xcor_t > 0.13 && s.xcor_t < 0.16);
    assert(s.rs_t > 42.0 && s.rs_t < 52.0);

    // geometry-scaled local parameters
    assert(s.vp == 0.25);
    assert(s.ax > 4.9 && s.ax < 5.5);
    assert(s.alp > 0.032 && s.alp < 0.037);
    assert(s.alp1 > 0.088 && s.alp1 < 0.10);
    assert(s.alp2 > 1.76 && s.alp2 < 2.02);
    assert(s.theta_sat_t > 2.5 && s.theta_sat_t < 3.5);
    assert(s.cf > 0.055 && s.cf < 0.067);
    assert(s.cfb > 0.28 && s.cfb < 0.32);
    assert(s.ct > 0.17 && s.ct < 0.22);
}

void testRefNmosTemperatureAndGeometryDependence() {
    using namespace gspice;
    const auto model = refNmosModel();
    const auto prepared21 = model.prepare({{"W", 1e-6}, {"L", 0.1e-6}}, 21.0);
    const auto s21 = psp103PrepareDevice(model, prepared21, 21.0);
    const auto s80 = psp103PrepareDevice(model, prepared21, 80.0);

    // hotter device: lower beta (STBET>0), lower phib, higher vfb_t
    assert(s80.bet_i < s21.bet_i);
    assert(s80.phib_dc < s21.phib_dc);
    assert(s80.vfb_t > s21.vfb_t);
    assert(std::abs(s80.delta_t - 53.0) < 1e-9);
    assert(std::abs(s21.delta_t + 6.0) < 1e-9);

    // longer channel: smaller alp / bet / neff, larger ax
    const auto preparedLong = model.prepare({{"W", 1e-6}, {"L", 2e-6}}, 21.0);
    const auto sLong = psp103PrepareDevice(model, preparedLong, 21.0);
    assert(sLong.alp < s21.alp);
    assert(sLong.bet_i < s21.bet_i);
    assert(sLong.neff < s21.neff);
    assert(sLong.ax > s21.ax);
    assert(sLong.alp2 < s21.alp2);
}

} // namespace

int main() {
    using namespace gspice;
    const auto model = Psp103ParameterSet::from({
        {"TNOM", 27.0}, {"VFB", 0.2}, {"STVFB", 1e-3},
        {"BETN", 2e-3}, {"STBETN", 1.0}, {"RS", 10.0},
        {"STRS", 0.5}, {"FNT", 1.0}});
    const auto prepared = model.prepare({{"W", 1e-6}, {"L", 0.1e-6}}, 27.0);
    const auto scaled = psp103ScaleTemperature(model, prepared, 127.0);

    assert(std::abs(scaled.delta_t - 100.0) < 1e-12);
    assert(std::abs(scaled.vfb - 0.3) < 1e-12);
    assert(std::abs(scaled.beta - 2e-3 * (300.15 / 400.15)) < 1e-12);
    assert(std::abs(scaled.rs - 10.0 * std::pow(300.15 / 400.15, 0.5)) < 1e-12);
    assert(scaled.rs > 0.0);
    assert(scaled.nt > 0.0);

    testRefNmosShortChannel();
    testRefNmosTemperatureAndGeometryDependence();
    std::printf("PREPARE DEVICE TEST PASS\n");
    return 0;
}
