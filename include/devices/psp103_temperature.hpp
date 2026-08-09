#ifndef GSPICE_PSP103_TEMPERATURE_HPP
#define GSPICE_PSP103_TEMPERATURE_HPP

#include "psp103_parameters.hpp"

#include <algorithm>
#include <cmath>

namespace gspice {

struct Psp103TemperatureScaled {
    double tkr = 300.15;
    double tka = 300.15;
    double delta_t = 0.0;
    double ratio_ref_device = 1.0;
    double vfb = 0.0;
    double beta = 0.0;
    double theta_mu = 0.0;
    double mu_e = 0.0;
    double cs = 0.0;
    double theta_cs = 0.0;
    double rs = 0.0;
    double theta_sat = 0.0;
    double theta_sat_ac = 0.0;
    double ig_inv = 0.0;
    double ig_ov = 0.0;
    double ig_ovd = 0.0;
    double agidl = 0.0;
    double bgidl = 0.0;
    double agidld = 0.0;
    double bgidld = 0.0;
    double rth = 0.0;
    double nt = 0.0;
};

struct Psp103SelfHeatingInputs {
    double power_w = 0.0;
    double previous_rise_c = 0.0;
    double dt_s = 0.0;
    double rth_k_per_w = 0.0;
    double cth_j_per_k = 0.0;
};

struct Psp103SelfHeatingResult {
    double residual_w = 0.0;
    double dresidual_drise_w_per_k = 0.0;
};

inline Psp103SelfHeatingResult psp103SelfHeatingResidual(
    const Psp103SelfHeatingInputs& in, double rise_c) {
    const double rth = std::max(in.rth_k_per_w, 0.0);
    const double cth = std::max(in.cth_j_per_k, 0.0);
    const bool transient = in.dt_s > 0.0 && cth > 0.0;
    const double conductance = rth > 0.0 ? 1.0 / rth : 0.0;
    const double capacitance = transient ? cth / in.dt_s : 0.0;
    return {
        capacitance * (rise_c - in.previous_rise_c) + conductance * rise_c - in.power_w,
        capacitance + conductance};
}

// PSP physical constants (Common103_macrodefs.include)
namespace psp103_phys {
constexpr double kelvin_conversion = 273.15;
constexpr double kbol = 1.3806505e-23;
constexpr double qele = 1.6021918e-19;
constexpr double eps0 = 8.8541878176e-12;
constexpr double epsr_si = 11.8;
constexpr double eps_si = eps0 * epsr_si;
constexpr double len = 1.0e-6;
constexpr double wen = 1.0e-6;
constexpr double qmn = 5.951993;
constexpr double qmp = 7.448711;
constexpr double twoThirds = 6.6666666666666667e-01;
} // namespace psp_phys

inline double clipLow(double value, double min) { return value > min ? value : min; }
inline double clipHigh(double value, double max) { return value < max ? value : max; }
inline double clipBoth(double value, double min, double max) {
    return value > min ? (value < max ? value : max) : min;
}

// PSP minimum with smooth rolloff: MINA(x,y,a) (PSP103_macrodefs.include)
inline double pspMinA(double x, double y, double a) {
    return 0.5 * (x + y - std::sqrt((x - y) * (x - y) + a));
}

struct Psp103DeviceSetup {
    // temperatures (K)
    double tkr = 300.15;
    double tkd = 300.15;
    double delta_t = 0.0;
    double r_tn = 1.0;
    double ln_r_tn = 0.0;
    double phit = 0.025852;
    double inv_phit = 1.0 / 0.025852;
    double eg = 1.17;
    double phib_fac = 1.0e-3;
    double phita = 0.025852;
    double inv_phita = 1.0 / 0.025852;

    // TempScaling outputs
    double phit0 = 0.025852;
    double phib_dc = 0.9;
    double g0_dc = 0.9;
    double kp = 0.0;
    double sqrt_phib_dc = std::sqrt(0.9);
    double phix_dc = 0.0;
    double aphi_dc = 0.0;
    double phix1_dc = 0.0;
    double alpha_b = 0.0;
    double us1 = 0.0;
    double us21 = 0.0;
    double vfb_t = 0.0;
    double bet_i = 0.0;
    double themu_t = 1.5;
    double mue_t = 0.0;
    double cs_t = 0.0;
    double xcor_t = 0.0;
    double rs_t = 0.0;
    double ther = 0.0;
    double theta_sat_t = 1.0;
    double a2_t = 0.0;

    // local (geometry scaled) parameters
    double cf = 0.0;
    double cfd = 0.0;
    double cfb = 0.0;
    double psce = 0.0;
    double psce_b = 0.0;
    double psce_d = 0.0;
    double dnsub = 0.0;
    double vnsub = 0.0;
    double nslp = 0.05;
    double rsb = 0.0;
    double rsg = 0.0;
    double thesat_b = 0.0;
    double thesat_g = 0.0;
    double ax = 18.0;
    double alp = 0.0;
    double alp1 = 0.0;
    double alp2 = 0.0;
    double vp = 1.0;
    double vsbnud = 0.0;
    double dvsbnud = 1.0;
    double gfac_nud = 1.0;
    int sw_nud = 0;

    // process / internal
    double neff = 5.0e23;
    double neffac = 5.0e23;
    double cox_prime = 0.0;
    double cox_over_q = 0.0;
    double tox_sq = 0.0;
    double qq = 0.0;
    double qmc = 0.0;
    double e_eff0 = 0.0;
    double eta_mu = 0.5;
    double eta_mu1 = 0.5;
    double npol = 0.0;
    double ct = 0.0;
    double dphib = 0.0;
    int channel_type = 1;
    double mult_i = 1.0;
    int swgeo = 1;

    // charge-model extras (populated later)
    double phib_ac = 0.9;
    double g0_ac = 0.9;
    double phix_ac = 0.0;
    double aphi_ac = 0.0;
    double phix1_ac = 0.0;

    // Charge / overlap geometry (SPcalc_ac extrinsic block)
    double qlim2 = 0.0;        // 100 * phit^2
    double tox_i = 0.0;        // TOX_i (m)
    double epsrox_i = 0.0;     // EPSROX_i
    double cox_i = 0.0;        // total intrinsic COX (F): eps*Wcv*Lcv/tox
    double cgov_i = 0.0;       // source-overlap oxide cap (F)
    double cgovd_i = 0.0;      // drain-overlap oxide cap (F)
    double cgbov_i = 0.0;      // gate-bulk overlap cap (F)
    double cfr_i = 0.0;        // source outer-fringe cap (F)
    double cfrd_i = 0.0;       // drain outer-fringe cap (F)
    // Overlap surface potential constants (sp_ovInit macro)
    double gov_s = 0.0, gov2_s = 0.0, spov_eps2_s = 0.0;
    double spov_a_s = 0.0, spov_delta1_s = 0.0;
    double gov_d = 0.0, gov2_d = 0.0, spov_eps2_d = 0.0;
    double spov_a_d = 0.0, spov_delta1_d = 0.0;
};

// Geometry scaling (SWGEO=1) + PSP103 TempInitialize + TempScaling.
// Implements PSP103_module.include lines ~1100-1468 (geometry) and
// PSP103_macrodefs.include TempInitialize / TempScaling.
inline Psp103DeviceSetup psp103PrepareDevice(
    const Psp103ParameterSet& model,
    const Psp103PreparedModel& prepared,
    double device_c,
    double thermal_rise_c = 0.0) {
    Psp103DeviceSetup s;

    const double eps0 = psp103_phys::eps0;
    const double eps_si = psp103_phys::eps_si;
    const double qe = psp103_phys::qele;
    const double kb = psp103_phys::kbol;
    const double l_ref = psp103_phys::len;
    const double w_ref = psp103_phys::wen;

    s.channel_type = prepared.polarity < 0 ? -1 : 1;
    const bool pmos = s.channel_type == -1;

    const double toxo = clipLow(model.get("TOXO", 2.0e-9), 1.0e-10);
    const double epsroxo = clipLow(model.get("EPSROXO", 3.9), 1.0);
    const double nsub0 = clipLow(model.get("NSUBO", 3.0e23), 1.0e20);
    const double wseg = clipLow(model.get("WSEG", 1.0e-8), 1.0e-10);
    const double npck = clipLow(model.get("NPCK", 1.0e24), 0.0);
    const double wsegp = clipLow(model.get("WSEGP", 1.0e-8), 1.0e-10);
    const double lpck = clipLow(model.get("LPCK", 1.0e-8), 1.0e-10);
    const double nsubw = model.get("NSUBW", 0.0);
    const double npckw = model.get("NPCKW", 0.0);
    const double lpckw = model.get("LPCKW", 0.0);
    const double fol1 = model.get("FOL1", 0.0);
    const double fol2 = model.get("FOL2", 0.0);
    const double facneffaco = model.get("FACNEFFACO", 1.0);
    const double facneffacl = model.get("FACNEFFACL", 0.0);
    const double facneffacw = model.get("FACNEFFACW", 0.0);
    const double facneffac_lw = model.get("FACNEFFACLW", 0.0);
    const double gfacnudo = model.get("GFACNUDO", 1.0);
    const double gfacnudl = model.get("GFACNUDL", 0.0);
    const double gfacnudl_exp = model.get("GFACNUDLEXP", 1.0);
    const double gfacnudw = model.get("GFACNUDW", 0.0);
    const double gfacnud_lw = model.get("GFACNUDLW", 0.0);
    const double vsbnudo = model.get("VSBNUDO", 0.0);
    const double dvsbnudo = model.get("DVSBNUDO", 1.0);
    const double vnsubo = model.get("VNSUBO", 0.0);
    const double nslpo = clipLow(model.get("NSLPO", 0.05), 1.0e-3);
    const double dnsubo = clipBoth(model.get("DNSUBO", 0.0), 0.0, 1.0);
    const double dphibo = model.get("DPHIBO", 0.0);
    const double dphibl = model.get("DPHIBL", 0.0);
    const double dphibl_exp = model.get("DPHIBLEXP", 1.0);
    const double dphibw = model.get("DPHIBW", 0.0);
    const double dphib_lw = model.get("DPHIBLW", 0.0);
    const double npo = model.get("NPO", 1.0e26);
    const double npl = model.get("NPL", 0.0);
    const double cto = model.get("CTO", 0.0);
    const double ctl = model.get("CTL", 0.0);
    const double ctl_exp = model.get("CTLEXP", 1.0);
    const double ctw = model.get("CTW", 0.0);
    const double ctlw = model.get("CTLW", 0.0);
    const double cfl = model.get("CFL", 0.0);
    const double cfl_exp = model.get("CFLEXP", 2.0);
    const double cfw = model.get("CFW", 0.0);
    const double cfbo = model.get("CFBO", 0.0);
    const double u0 = model.get("UO", 5.0e-2);
    const double fbet1 = model.get("FBET1", 0.0);
    const double fbet1w = model.get("FBET1W", 0.0);
    const double lp1 = clipLow(model.get("LP1", 1.0e-8), 1.0e-10);
    const double lp1w = model.get("LP1W", 0.0);
    const double fbet2 = model.get("FBET2", 0.0);
    const double lp2 = clipLow(model.get("LP2", 1.0e-8), 1.0e-10);
    const double betw1 = model.get("BETW1", 0.0);
    const double betw2 = model.get("BETW2", 0.0);
    const double wbet = clipLow(model.get("WBET", 1.0e-9), 1.0e-10);
    const double stbeto = model.get("STBETO", 1.0);
    const double stbetl = model.get("STBETL", 0.0);
    const double stbetw = model.get("STBETW", 0.0);
    const double stbetlw = model.get("STBETLW", 0.0);
    const double mueo = model.get("MUEO", 0.5);
    const double muew = model.get("MUEW", 0.0);
    const double stmueo = model.get("STMUEO", 0.0);
    const double themuo = model.get("THEMUO", 1.5);
    const double stthemuo = model.get("STTHEMUO", 1.5);
    const double cso = model.get("CSO", 0.0);
    const double csl = model.get("CSL", 0.0);
    const double csl_exp = model.get("CSLEXP", 1.0);
    const double csw = model.get("CSW", 0.0);
    const double cslw = model.get("CSLW", 0.0);
    const double stcso = model.get("STCSO", 0.0);
    const double xcoro = model.get("XCORO", 0.0);
    const double xcorl = model.get("XCORL", 0.0);
    const double xcorw = model.get("XCORW", 0.0);
    const double xcor_lw = model.get("XCORLW", 0.0);
    const double stxcoro = model.get("STXCORO", 0.0);
    const double fetao = model.get("FETAO", 1.0);
    const double rsw1 = model.get("RSW1", 50.0);
    const double rsw2 = model.get("RSW2", 0.0);
    const double strso = model.get("STRSO", 1.0);
    const double rsbo = model.get("RSBO", 0.0);
    const double rsgo = model.get("RSGO", 0.0);
    const double thesato = model.get("THESATO", 0.0);
    const double thesatl = model.get("THESATL", 0.0);
    const double thesatl_exp = model.get("THESATLEXP", 1.0);
    const double thesatw = model.get("THESATW", 0.0);
    const double thesatlw = model.get("THESATLW", 0.0);
    const double stthesato = model.get("STTHESATO", 1.0);
    const double stthesatl = model.get("STTHESATL", 0.0);
    const double stthesatw = model.get("STTHESATW", 0.0);
    const double stthesatlw = model.get("STTHESATLW", 0.0);
    const double thesatbo = model.get("THESATBO", 0.0);
    const double thesatgo = model.get("THESATGO", 0.0);
    const double axo = model.get("AXO", 18.0);
    const double axl = clipLow(model.get("AXL", 0.0), 0.0);
    const double alpl = model.get("ALPL", 5.0e-4);
    const double alpl_exp = model.get("ALPLEXP", 1.0);
    const double alpw = model.get("ALPW", 0.0);
    const double alp1l1 = model.get("ALP1L1", 0.0);
    const double alp1l_exp = model.get("ALP1LEXP", 0.5);
    const double alp1l2 = clipLow(model.get("ALP1L2", 0.0), 0.0);
    const double alp1w = model.get("ALP1W", 0.0);
    const double alp2l1 = model.get("ALP2L1", 0.0);
    const double alp2l_exp = model.get("ALP2LEXP", 0.5);
    const double alp2l2 = clipLow(model.get("ALP2L2", 0.0), 0.0);
    const double alp2w = model.get("ALP2W", 0.0);
    const double vpo = model.get("VPO", 0.05);
    const double swgeo = model.get("SWGEO", 1.0);
    const double qmc = clipLow(model.get("QMC", 1.0), 0.0);
    const double factuo = clipLow(model.get("FACTUO", 1.0), 0.0);
    const double delvto = model.get("DELVTO", 0.0);
    const double stvfbo = model.get("STVFBO", 5.0e-4);
    const double stvfbl = model.get("STVFBL", 0.0);
    const double stvfbw = model.get("STVFBW", 0.0);
    const double stvfblw = model.get("STVFBLW", 0.0);
    const double vfbo = model.get("VFBO", -1.0);

    // geometry (initial_instance)
    double nf = clipLow(prepared.geometry.fingers, 1.0);
    nf = std::floor(nf + 0.5);
    const double inv_nf = 1.0 / nf;
    const double L_i = clipLow(prepared.geometry.length_m, 1.0e-9);
    const double W_i = clipLow(prepared.geometry.width_m * inv_nf, 1.0e-9);
    const double iL = l_ref / L_i;
    const double iW = w_ref / W_i;
    const double lvaro = model.get("LVARO", 0.0);
    const double lvarl = model.get("LVARL", 0.0);
    const double lvarw = model.get("LVARW", 0.0);
    const double wvaro = model.get("WVARO", 0.0);
    const double wvarl = model.get("WVARL", 0.0);
    const double wvarw = model.get("WVARW", 0.0);
    const double lap = model.get("LAP", 0.0);
    const double wot = model.get("WOT", 0.0);
    const double toxovo = clipLow(model.get("TOXOVO", 2.0e-9), 1.0e-10);
    const double toxovdo = clipLow(model.get("TOXOVDO", 2.0e-9), 1.0e-10);
    const double novo = clipBoth(model.get("NOVO", 5.0e25), 1.0e20, 1.0e27);
    const double novdo = clipBoth(model.get("NOVDO", 5.0e25), 1.0e20, 1.0e27);
    const double lov_i = clipLow(model.get("LOV", 0.0), 0.0);
    const double lovd_i = clipLow(model.get("LOVD", 0.0), 0.0);
    const double cgbovl = model.get("CGBOVL", 0.0);
    const double cfrw = model.get("CFRW", 0.0);
    const double cfrdw = model.get("CFRDW", 0.0);
    const double del_lps = lvaro * (1.0 + lvarl * iL) * (1.0 + lvarw * iW);
    const double del_wod = wvaro * (1.0 + wvarl * iL) * (1.0 + wvarw * iW);
    const double LE = clipLow(L_i + del_lps - 2.0 * lap, 1.0e-9);
    const double WE = clipLow(W_i + del_wod - 2.0 * wot, 1.0e-9);
    const double dlq = model.get("DLQ", 0.0);
    const double dwq = model.get("DWQ", 0.0);
    const double LEcv = clipLow(L_i + del_lps - 2.0 * lap + dlq, 1.0e-9);
    const double WEcv = clipLow(W_i + del_wod - 2.0 * wot + dwq, 1.0e-9);
    const double Lcv = clipLow(L_i + del_lps + dlq, 1.0e-9);
    const double Wcv = clipLow(W_i + del_wod + dwq, 1.0e-9);
    const double iLE = l_ref / LE;
    const double iWE = w_ref / WE;

    // SWGEO=1 geometry scaling (lines 1311-1468); SWGEO=0 local fallback
    double VFB_p = model.get("VFB", vfbo);
    double STVFB_p = model.get("STVFB", stvfbo);
    double NEFF_p = model.get("NEFF", nsub0);
    double DPHIB_p = model.get("DPHIB", 0.0);
    double DELVTAC_p = model.get("DELVTAC", 0.0);
    double BETN_p = model.get("BETN", 7.0e-2);
    double STBET_p = model.get("STBET", 1.0);
    double MUE_p = model.get("MUE", 0.5);
    double STMUE_p = model.get("STMUE", 0.0);
    double THEMU_p = model.get("THEMU", 1.5);
    double STTHEMU_p = model.get("STTHEMU", 1.5);
    double CS_p = model.get("CS", 0.0);
    double STCS_p = model.get("STCS", 0.0);
    double XCOR_p = model.get("XCOR", 0.0);
    double STXCOR_p = model.get("STXCOR", 0.0);
    double FETA_p = model.get("FETA", fetao);
    double RS_p = model.get("RS", 30.0);
    double STRS_p = model.get("STRS", 1.0);
    double THESAT_p = model.get("THESAT", 1.0);
    double STTHESAT_p = model.get("STTHESAT", 1.0);
    double AX_p = model.get("AX", axo);
    double ALP_p = model.get("ALP", 0.01);
    double ALP1_p = model.get("ALP1", 0.0);
    double ALP2_p = model.get("ALP2", 0.0);
    double VP_p = model.get("VP", vpo);
    double CF_p = model.get("CF", 0.0);
    double CFD_p = model.get("CFD", 0.0);
    double CFB_p = model.get("CFB", 0.0);
    double PSCE_p = model.get("PSCE", 0.0);
    double PSCEB_p = model.get("PSCEB", 0.0);
    double PSCED_p = model.get("PSCED", 0.0);
    double TOXOV_p = model.get("TOXOV", 2.0e-9);
    double TOXOVD_p = model.get("TOXOVD", 2.0e-9);
    double NOV_p = model.get("NOV", 5.0e25);
    double NOVD_p = model.get("NOVD", 5.0e25);
    double CGOV_p = model.get("CGOV", 1.0e-15);
    double CGOVD_p = model.get("CGOVD", 1.0e-15);
    double CGBOV_p = model.get("CGBOV", 0.0);
    double CFR_p = model.get("CFR", 0.0);
    double CFRD_p = model.get("CFRD", 0.0);
    double VSBNUD_p = model.get("VSBNUD", vsbnudo);
    double DVSBNUD_p = model.get("DVSBNUD", dvsbnudo);
    double GFACNUD_p = model.get("GFACNUD", gfacnudo);
    double FACNEFFAC_p = model.get("FACNEFFAC", 1.0);
    double TOX_p = model.get("TOX", toxo);
    double EPSROX_p = model.get("EPSROX", epsroxo);
    double NP_p = model.get("NP", npo);
    double CT_p = model.get("CT", 0.0);

    const bool do_geo = (swgeo == 1.0 || swgeo == 2.0);
    if (do_geo) {
        VFB_p = vfbo + model.get("VFBL", 0.0) * iLE +
            model.get("VFBW", 0.0) * iWE +
            model.get("VFBLW", 0.0) * iLE * iWE;
        STVFB_p = stvfbo + stvfbl * iLE + stvfbw * iWE + stvfblw * iLE * iWE;
        const double nsub0e = nsub0 * std::max(
            1.0 + nsubw * iWE * std::log(1.0 + WE / wseg), 1.0e-3);
        const double npcke = npck * std::max(
            1.0 + npckw * iWE * std::log(1.0 + WE / wsegp), 1.0e-3);
        const double lpcke = lpck * std::max(
            1.0 + lpckw * iWE * std::log(1.0 + WE / wsegp), 1.0e-3);
        double NSUB;
        if (LE > 2.0 * lpcke) {
            const double aa = 7.5e10;
            const double bb = std::sqrt(nsub0e + 0.5 * npcke) - std::sqrt(nsub0e);
            const double sq = std::sqrt(nsub0e) +
                aa * std::log(1.0 + 2.0 * lpcke / LE * (std::exp(bb / aa) - 1.0));
            NSUB = sq * sq;
        } else if (LE >= lpcke) {
            NSUB = nsub0e + npcke * lpcke / LE;
        } else {
            NSUB = nsub0e + npcke * (2.0 - LE / lpcke);
        }
        NEFF_p = NSUB * (1.0 - fol1 * iLE - fol2 * iLE * iLE);
        FACNEFFAC_p = facneffaco + facneffacl * iLE + facneffacw * iWE +
            facneffac_lw * iLE * iWE;
        GFACNUD_p = gfacnudo + gfacnudl * std::pow(iLE, gfacnudl_exp) +
            gfacnudw * iWE + gfacnud_lw * iLE * iWE;
        VSBNUD_p = vsbnudo;
        DVSBNUD_p = dvsbnudo;
        DPHIB_p = dphibo + dphibl * std::pow(iLE, dphibl_exp) +
            dphibw * iWE + dphib_lw * iLE * iWE;
        DELVTAC_p = model.get("DELVTACO", 0.0) +
            model.get("DELVTACL", 0.0) * std::pow(iLE, model.get("DELVTACLEXP", 1.0)) +
            model.get("DELVTACW", 0.0) * iWE +
            model.get("DELVTACLW", 0.0) * iLE * iWE;
        NP_p = clipLow(npo * std::max(1.0e-6, 1.0 + npl * iLE), 0.0);
        CT_p = (cto + ctl * std::pow(iLE, ctl_exp)) * (1.0 + ctw * iWE) *
            (1.0 + ctlw * iLE * iWE);

        CF_p = cfl * std::pow(iLE, cfl_exp) * (1.0 + cfw * iWE);
        CFD_p = model.get("CFDO", 0.0);
        CFB_p = cfbo;
        PSCE_p = model.get("PSCEL", 0.0) *
            std::pow(iLE, model.get("PSCELEXP", 2.0)) *
            (1.0 + model.get("PSCEW", 0.0) * iWE);
        PSCED_p = model.get("PSCEDO", 0.0);
        PSCEB_p = model.get("PSCEBO", 0.0);

        const double fbet1e = fbet1 * (1.0 + fbet1w * iWE);
        const double lp1e = lp1 * std::max(1.0 + lp1w * iWE, 1.0e-3);
        double gpe = 1.0 + fbet1e * lp1e / LE * (1.0 - std::exp(-LE / lp1e)) +
            fbet2 * lp2 / LE * (1.0 - std::exp(-LE / lp2));
        gpe = std::max(gpe, 1.0e-15);
        const double gwe = 1.0 + betw1 * iWE +
            betw2 * iWE * std::log(1.0 + WE / wbet);
        BETN_p = u0 * WE / (gpe * LE) * gwe;
        STBET_p = stbeto + stbetl * iLE + stbetw * iWE + stbetlw * iLE * iWE;
        MUE_p = mueo * (1.0 + muew * iWE);
        STMUE_p = stmueo;
        THEMU_p = themuo;
        STTHEMU_p = stthemuo;
        CS_p = (cso + csl * std::pow(iLE, csl_exp)) * (1.0 + csw * iWE) *
            (1.0 + cslw * iLE * iWE);
        STCS_p = stcso;
        XCOR_p = xcoro * (1.0 + xcorl * iLE) * (1.0 + xcorw * iWE) *
            (1.0 + xcor_lw * iLE * iWE);
        STXCOR_p = stxcoro;
        FETA_p = fetao;
        RS_p = rsw1 * iWE * (1.0 + rsw2 * iWE);
        STRS_p = strso;
        THESAT_p = (thesato + thesatl * (gwe / gpe) *
            std::pow(iLE, thesatl_exp)) * (1.0 + thesatw * iWE) *
            (1.0 + thesatlw * iLE * iWE);
        STTHESAT_p = stthesato + stthesatl * iLE + stthesatw * iWE +
            stthesatlw * iLE * iWE;
        AX_p = axo / (1.0 + axl * iLE);
        ALP_p = alpl * std::pow(iLE, alpl_exp) * (1.0 + alpw * iWE);
        double tmpx = std::pow(iLE, alp1l_exp);
        ALP1_p = alp1l1 * tmpx * (1.0 + alp1w * iWE) /
            (1.0 + alp1l2 * iLE * tmpx);
        tmpx = std::pow(iLE, alp2l_exp);
        ALP2_p = alp2l1 * tmpx * (1.0 + alp2w * iWE) /
            (1.0 + alp2l2 * iLE * tmpx);
        VP_p = vpo;
        TOXOV_p = toxovo;
        TOXOVD_p = toxovdo;
        NOV_p = novo;
        NOVD_p = novdo;
    }

    // charge-model overlap geometry: SWGEO=1 computed caps (line 1428ff),
    // SWGEO=0 falls back to the model-card values already read above
    double COX_p = model.get("COX", 1.0e-15);
    if (do_geo) {
        COX_p = eps0 * epsroxo * WEcv * LEcv / toxo;
        CGOV_p = eps0 * epsroxo * WEcv * lov_i / toxovo;
        CGOVD_p = eps0 * epsroxo * WEcv * lovd_i / toxovdo;
        CGBOV_p = cgbovl * Lcv / l_ref;
        CFR_p = cfrw * Wcv / w_ref;
        CFRD_p = cfrdw * Wcv / w_ref;
    }

    // clipping of local parameters (lines 1562-1677)
    const double vfb_i = VFB_p;
    const double stvfb_i = STVFB_p;
    const double tox_i = clipLow(do_geo ? toxo : TOX_p, 1.0e-10);
    const double epsrox_i = clipLow(do_geo ? epsroxo : EPSROX_p, 1.0);
    const double neff_i = clipBoth(NEFF_p, 1.0e20, 1.0e26);
    const double facneffac_i = clipLow(FACNEFFAC_p, 0.0);
    const double gfacnud_i = clipLow(GFACNUD_p, 0.01);
    const double vsbnud_i = clipLow(VSBNUD_p, 0.0);
    const double dvsbnud_i = clipLow(DVSBNUD_p, 0.1);
    const double vnsub_i = do_geo ? vnsubo : model.get("VNSUB", 0.0);
    const double nslp_i = clipLow(do_geo ? nslpo : model.get("NSLP", 0.05), 1.0e-3);
    const double dnsub_i = clipBoth(do_geo ? dnsubo : model.get("DNSUB", 0.0), 0.0, 1.0);
    const double dphib_i = DPHIB_p;
    const double delvtac_i = DELVTAC_p;
    const double npol_i = clipLow(NP_p, 0.0);
    const double ct_i = clipLow(CT_p, 0.0);
    const double cf_i = clipLow(CF_p, 0.0);
    const double cfd_i = clipLow(CFD_p, 0.0);
    const double cfb_i = clipBoth(CFB_p, 0.0, 1.0);
    const double psce_i = clipLow(PSCE_p, 0.0);
    const double psce_b_i = clipBoth(PSCEB_p, 0.0, 1.0);
    const double psce_d_i = clipLow(PSCED_p, 0.0);
    const double betn_i = clipLow(BETN_p, 0.0);
    const double stbet_i = STBET_p;
    const double mue_i = clipLow(MUE_p, 0.0);
    const double stmue_i = STMUE_p;
    const double themu_i = clipLow(THEMU_p, 0.0);
    const double stthemu_i = STTHEMU_p;
    const double cs_i = clipLow(CS_p, 0.0);
    const double stcs_i = STCS_p;
    const double xcor_i = clipLow(XCOR_p, 0.0);
    const double stxcor_i = STXCOR_p;
    const double feta_i = clipLow(FETA_p, 0.0);
    const double rs_i = clipLow(RS_p, 0.0);
    const double strs_i = STRS_p;
    const double rsb_i = clipBoth(do_geo ? rsbo : model.get("RSB", 0.0), -0.5, 1.0);
    const double rsg_i = clipLow(do_geo ? rsgo : model.get("RSG", 0.0), -0.5);
    const double thesat_i = clipLow(THESAT_p, 0.0);
    const double stthesat_i = STTHESAT_p;
    const double thesatb_i = clipBoth(do_geo ? thesatbo : model.get("THESATB", 0.0), -0.5, 1.0);
    const double thesatg_i = clipLow(do_geo ? thesatgo : model.get("THESATG", 0.0), -0.5);
    const double ax_i = clipLow(AX_p, 2.0);
    const double alp_i = clipLow(ALP_p, 0.0);
    const double alp1_i = clipLow(ALP1_p, 0.0);
    const double alp2_i = clipLow(ALP2_p, 0.0);
    const double vp_i = clipLow(VP_p, 1.0e-10);

    // overlap/local charge caps (clipping + symmetric-connection, lines 1681-1694)
    const double toxov_i = clipLow(TOXOV_p, 1.0e-10);
    const double toxovd_i = clipLow(TOXOVD_p, 1.0e-10);
    const double nov_i = clipBoth(NOV_p, 1.0e20, 1.0e27);
    double novd_i = clipBoth(NOVD_p, 1.0e20, 1.0e27);
    const double cgov_i = clipLow(CGOV_p, 0.0);
    double cgovd_i = clipLow(CGOVD_p, 0.0);
    const double cgbov_i = clipLow(CGBOV_p, 0.0);
    const double cfr_i = clipLow(CFR_p, 0.0);
    double cfrd_i = clipLow(CFRD_p, 0.0);
    // symmetric-junction drain aliasing (SWJUNASYM=0 default)
    const int swjunasym = static_cast<int>(model.get("SWJUNASYM", 0.0));
    const double toxovd = (swjunasym == 0) ? toxov_i : toxovd_i;
    if (swjunasym == 0) {
        novd_i = nov_i;
        cgovd_i = cgov_i;
        cfrd_i = cfr_i;
    }
    const double cox_i = clipLow(COX_p, 0.0);

    // overlap surface-potential coefficients Gov/sp_ovInit computed after the
    // TempScaling block below (needs inv_phit/phit at device temperature)

    // internal process parameters (lines 1697-1727)
    const double eps_ox = eps0 * epsrox_i;
    const double cox_prime = eps_ox / tox_i;
    const double tox_sq = tox_i * tox_i;
    const double cox_over_q = cox_prime / qe;
    const double neffac_i = clipBoth(facneffac_i * neff_i, 1.0e20, 1.0e26);
    double qq = 0.0;
    if (qmc > 0.0) {
        qq = 0.4 * psp103_phys::qmn * qmc * std::pow(cox_prime, psp103_phys::twoThirds);
        if (pmos) qq = psp103_phys::qmp / psp103_phys::qmn * qq;
    }
    const double e_eff0 = 1.0e-8 * cox_prime / eps_si;
    double eta_mu = 0.5 * feta_i;
    double eta_mu1 = 0.5;
    if (pmos) {
        eta_mu = (1.0 / 3.0) * feta_i;
        eta_mu1 = 1.0 / 3.0;
    }

    // TempInitialize + TempScaling
    const double tkr = model.get("TR", prepared.temperature.nominal_c) + 273.15;
    const double tkd = device_c + 273.15 + model.get("DTA", 0.0) + thermal_rise_c;
    const double delta_t = tkd - tkr;
    const double r_tn = tkr / tkd;
    const double ln_r_tn = std::log(r_tn);
    const double phit = tkd * kb / qe;
    const double inv_phit = 1.0 / phit;

    // overlap surface-potential coefficients Gov/sp_ovInit (lines 1464-1488)
    const double coxov_prime = eps0 * epsrox_i / toxov_i;
    const double coxov_prime_d = eps0 * epsrox_i / toxovd;
    const double gov_s = std::sqrt(2.0 * qe * nov_i * eps_si * inv_phit) / coxov_prime;
    const double gov_d = std::sqrt(2.0 * qe * novd_i * eps_si * inv_phit) / coxov_prime_d;
    double spov_eps2_s = 0.0, spov_a_s = 0.0, spov_delta1_s = 0.0;
    double spov_eps2_d = 0.0, spov_a_d = 0.0, spov_delta1_d = 0.0;
    {
        const double inv_gov_s = 1.0 / gov_s;
        const double spov_eps_s = 3.1 * gov_s + 8.5;
        spov_eps2_s = spov_eps_s * spov_eps_s;
        const double spov_delta_s = 0.5 * spov_eps_s;
        double spov_a;
        if (inv_gov_s < 0.06) {
            spov_a = 64.0 * inv_gov_s;
        } else if (inv_gov_s <= 0.45) {
            spov_a = 22.0 * inv_gov_s + 3.0;
        } else if (inv_gov_s <= 1.6) {
            spov_a = -7.2 * inv_gov_s + 15.5;
        } else {
            spov_a = gov_s;
        }
        spov_a_s = spov_a;
        spov_delta1_s = spov_delta_s - gov_s *
            std::sqrt(spov_delta_s + gov_s * gov_s * 0.25 + spov_a);
    }
    {
        const double inv_gov_d = 1.0 / gov_d;
        const double spov_eps_d = 3.1 * gov_d + 8.5;
        spov_eps2_d = spov_eps_d * spov_eps_d;
        const double spov_delta_d = 0.5 * spov_eps_d;
        double spov_a;
        if (inv_gov_d < 0.06) {
            spov_a = 64.0 * inv_gov_d;
        } else if (inv_gov_d <= 0.45) {
            spov_a = 22.0 * inv_gov_d + 3.0;
        } else if (inv_gov_d <= 1.6) {
            spov_a = -7.2 * inv_gov_d + 15.5;
        } else {
            spov_a = gov_d;
        }
        spov_a_d = spov_a;
        spov_delta1_d = spov_delta_d - gov_d *
            std::sqrt(spov_delta_d + gov_d * gov_d * 0.25 + spov_a);
    }
    const double eg = 1.179 - 9.025e-5 * tkd - 3.05e-7 * tkd * tkd;
    const double phib_fac = std::max(
        (1.045 + 4.5e-4 * tkd) * (0.523 + 1.4e-3 * tkd - 1.48e-6 * tkd * tkd) *
            tkd * tkd / 9.0e4,
        1.0e-3);

    const double phit0 = phit * (1.0 + ct_i * r_tn);
    double phib_dc = eg + dphib_i +
        2.0 * phit * std::log(neff_i * std::pow(phib_fac, -0.75) * 4.0e-26);
    phib_dc = std::max(phib_dc, 5.0e-2);
    double g0_dc = std::sqrt(2.0 * qe * neff_i * eps_si * inv_phit) / cox_prime;

    double kp = 0.0;
    if (npol_i > 0.0) {
        const double arg2max = 8.0e7 / tox_sq;
        double npol = std::max(npol_i, arg2max);
        npol = std::max(5.0e24, npol);
        kp = 2.0 * cox_prime * cox_prime * phit / (qe * npol * eps_si);
    }

    if (qmc > 0.0) {
        const double qb0 = std::sqrt(phit * g0_dc * g0_dc * phib_dc);
        const double dphibq = 0.75 * qq * std::pow(qb0, psp103_phys::twoThirds);
        phib_dc = phib_dc + dphibq;
        g0_dc = g0_dc * (1.0 + 2.0 * psp103_phys::twoThirds * dphibq / qb0);
    }
    const double sqrt_phib_dc = std::sqrt(phib_dc);
    const double phix_dc = 0.95 * phib_dc;
    const double aphi_dc = 0.0025 * phib_dc * phib_dc;
    const double phix2 = 0.5 * std::sqrt(aphi_dc);
    const double phix1_dc = pspMinA(phix_dc - phix2, 0.0, aphi_dc);
    const double alpha_b = 0.5 * (phib_dc + eg);
    const double us1 = std::sqrt(vsbnud_i + phib_dc) - sqrt_phib_dc;
    const double us21 = std::sqrt(vsbnud_i + dvsbnud_i + phib_dc) -
        sqrt_phib_dc - us1;

    const double vfb_t = vfb_i + stvfb_i * delta_t + delvto;
    const double tf_bet = std::exp(stbet_i * ln_r_tn);
    const double betn_t = betn_i * tf_bet;
    const double bet_i = factuo * betn_t * cox_prime;
    const double themu_t = themu_i * std::exp(stthemu_i * ln_r_tn);
    const double tf_mue = std::exp(stmue_i * ln_r_tn);
    const double mue_t = mue_i * tf_mue;
    const double tf_cs = std::exp(stcs_i * ln_r_tn);
    const double cs_t = cs_i * tf_cs;
    const double tf_xcor = std::exp(stxcor_i * ln_r_tn);
    const double xcor_t = xcor_i * tf_xcor;
    const double tf_ther = std::exp(strs_i * ln_r_tn);
    const double rs_t = rs_i * tf_ther;
    const double ther = 2.0 * bet_i * rs_t;
    const double tf_thesat = std::exp(stthesat_i * ln_r_tn);
    const double theta_sat_t = thesat_i * tf_thesat;
    const double a2_t = model.get("A2O", 10.0) *
        std::exp(-model.get("STA2O", 0.0) * ln_r_tn);

    // separate capacitance (SWDELVTAC path always populated)
    double phib_ac = eg + delvtac_i + dphib_i +
        2.0 * phit * std::log(neffac_i * std::pow(phib_fac, -0.75) * 4.0e-26);
    phib_ac = std::max(phib_ac, 5.0e-2);
    double g0_ac = std::sqrt(2.0 * qe * neffac_i * eps_si * inv_phit) / cox_prime;
    if (qmc > 0.0) {
        const double qb0 = std::sqrt(phit * g0_ac * g0_ac * phib_ac);
        const double dphibq = 0.75 * qq * std::pow(qb0, psp103_phys::twoThirds);
        phib_ac = phib_ac + dphibq;
        g0_ac = g0_ac * (1.0 + 2.0 * psp103_phys::twoThirds * dphibq / qb0);
    }
    const double phix_ac = 0.95 * phib_ac;
    const double aphi_ac = 0.0025 * phib_ac * phib_ac;
    const double phix2_ac = 0.5 * std::sqrt(aphi_ac);
    const double phix1_ac = pspMinA(phix_ac - phix2_ac, 0.0, aphi_ac);

    s.tkr = tkr; s.tkd = tkd; s.delta_t = delta_t; s.r_tn = r_tn;
    s.ln_r_tn = ln_r_tn; s.phit = phit; s.inv_phit = inv_phit;
    s.eg = eg; s.phib_fac = phib_fac; s.phita = phit; s.inv_phita = inv_phit;
    s.phit0 = phit0; s.phib_dc = phib_dc; s.g0_dc = g0_dc; s.kp = kp;
    s.sqrt_phib_dc = sqrt_phib_dc; s.phix_dc = phix_dc; s.aphi_dc = aphi_dc;
    s.phix1_dc = phix1_dc; s.alpha_b = alpha_b; s.us1 = us1; s.us21 = us21;
    s.vfb_t = vfb_t; s.bet_i = bet_i; s.themu_t = themu_t; s.mue_t = mue_t;
    s.cs_t = cs_t; s.xcor_t = xcor_t; s.rs_t = rs_t; s.ther = ther;
    s.theta_sat_t = theta_sat_t; s.a2_t = a2_t;
    s.cf = cf_i; s.cfd = cfd_i; s.cfb = cfb_i;
    s.psce = psce_i; s.psce_b = psce_b_i; s.psce_d = psce_d_i;
    s.dnsub = dnsub_i; s.vnsub = vnsub_i; s.nslp = nslp_i;
    s.rsb = rsb_i; s.rsg = rsg_i; s.thesat_b = thesatb_i; s.thesat_g = thesatg_i;
    s.ax = ax_i; s.alp = alp_i; s.alp1 = alp1_i; s.alp2 = alp2_i; s.vp = vp_i;
    s.vsbnud = vsbnud_i; s.dvsbnud = dvsbnud_i; s.gfac_nud = gfacnud_i;
    s.sw_nud = static_cast<int>(model.get("SWNUD", 0.0));
    s.neff = neff_i; s.neffac = neffac_i; s.cox_prime = cox_prime;
    s.cox_over_q = cox_over_q; s.tox_sq = tox_sq; s.qq = qq; s.qmc = qmc;
    s.e_eff0 = e_eff0; s.eta_mu = eta_mu; s.eta_mu1 = eta_mu1;
    s.npol = npol_i; s.ct = ct_i; s.dphib = dphib_i;
    s.channel_type = prepared.polarity < 0 ? -1 : 1;
    s.mult_i = prepared.geometry.multiplier * nf;
    s.swgeo = static_cast<int>(swgeo);
    s.phib_ac = phib_ac; s.g0_ac = g0_ac; s.phix_ac = phix_ac;
    s.aphi_ac = aphi_ac; s.phix1_ac = phix1_ac;
    s.qlim2 = 100.0 * phit * phit;
    s.tox_i = tox_i; s.epsrox_i = epsrox_i;
    s.cox_i = cox_i; s.cgov_i = cgov_i; s.cgovd_i = cgovd_i;
    s.cgbov_i = cgbov_i; s.cfr_i = cfr_i; s.cfrd_i = cfrd_i;
    s.gov_s = gov_s; s.gov2_s = gov_s * gov_s;
    s.spov_eps2_s = spov_eps2_s; s.spov_a_s = spov_a_s;
    s.spov_delta1_s = spov_delta1_s;
    s.gov_d = gov_d; s.gov2_d = gov_d * gov_d;
    s.spov_eps2_d = spov_eps2_d; s.spov_a_d = spov_a_d;
    s.spov_delta1_d = spov_delta1_d;
    return s;
}

inline double psp103TemperatureScale(double value, double ratio, double exponent) {
    return value * std::pow(ratio, exponent);
}

// PSP103 internal temperature scaling equations 4.1-4.3, 4.11, 4.52-4.82.
// Missing optional parameters retain PSP's zero/default behavior; no
// electrical behavior is enabled by this preprocessing object alone.
inline Psp103TemperatureScaled psp103ScaleTemperature(
    const Psp103ParameterSet& model,
    const Psp103PreparedModel& prepared,
    double device_c,
    double thermal_rise_c = 0.0) {
    Psp103TemperatureScaled scaled;
    const double t0 = prepared.temperature.nominal_c + 273.15;
    const double tkr = t0 + model.get("TR", 0.0);
    const double tka = t0 + model.get("TA", 0.0) +
        model.get("DTA", 0.0) + thermal_rise_c;
    const double tkd = device_c + 273.15;
    scaled.tkr = tkr;
    scaled.tka = tka;
    scaled.delta_t = tkd - tkr;
    scaled.ratio_ref_device = tkr / std::max(tkd, 1e-12);

    const double ratio = scaled.ratio_ref_device;
    scaled.vfb = model.alias(0.0, {"VFB", "VFB0", "VFBO"}) +
        model.get("STVFB", 0.0) * scaled.delta_t;
    scaled.beta = psp103TemperatureScale(
        model.alias(0.0, {"BETN", "BET", "BETO", "BETA", "KP"}), ratio,
        model.get("STBETN", 0.0));
    scaled.theta_mu = psp103TemperatureScale(
        model.get("THEMU", 0.0), ratio, model.get("STTHEMU", 0.0));
    scaled.mu_e = psp103TemperatureScale(
        model.get("MUE", 0.0), ratio, model.get("STMUE", 0.0));
    scaled.cs = psp103TemperatureScale(
        model.get("CS", 0.0), ratio, model.get("STCS", 0.0));
    scaled.theta_cs = psp103TemperatureScale(
        model.get("THECS", 0.0), ratio, model.get("STTHECS", 0.0));
    scaled.rs = psp103TemperatureScale(
        model.get("RS", 0.0), ratio, model.get("STRS", 0.0));
    scaled.theta_sat = psp103TemperatureScale(
        model.get("THESAT", 0.0), ratio, model.get("STTHESAT", 0.0));
    scaled.theta_sat_ac = psp103TemperatureScale(
        model.get("THESATAC", 0.0), ratio, model.get("STTHESAT", 0.0));
    const double stig = model.alias(0.0, {"STIG", "STIGO"});
    scaled.ig_inv = psp103TemperatureScale(
        model.get("IGINV", 0.0), tka / std::max(tkr, 1e-12), stig);
    scaled.ig_ov = psp103TemperatureScale(
        model.get("IGOV", 0.0), tka / std::max(tkr, 1e-12), stig);
    scaled.ig_ovd = psp103TemperatureScale(
        model.get("IGOVD", 0.0), tka / std::max(tkr, 1e-12), stig);
    scaled.agidl = model.sum({"AGIDL", "AGIDLW"});
    scaled.bgidl = model.sum({"BGIDL", "BGIDLO"}) *
        std::max(1.0 + model.alias(0.0, {"STBGIDL", "STBGIDLO"}) * scaled.delta_t, 0.0);
    scaled.agidld = model.sum({"AGIDLD", "AGIDLDW"});
    scaled.bgidld = model.sum({"BGIDLD", "BGIDLDO"}) *
        std::max(1.0 + model.alias(0.0, {"STBGIDLD", "STBGIDLDO"}) * scaled.delta_t, 0.0);
    scaled.rth = model.get("RTH", 0.0) * std::pow(
        std::max(tka * tkr, 1e-12) / (t0 * t0), model.get("STRTH", 0.0));
    scaled.nt = model.get("FNT", 0.0) * 4.0 * 1.380649e-23 * tkd;
    return scaled;
}

} // namespace gspice

#endif // GSPICE_PSP103_TEMPERATURE_HPP
