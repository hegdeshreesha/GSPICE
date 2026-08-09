// Debug dump: compare GSPICE PSP103 prepare internals @ T=21C against the
// VACASK reference psp103v4 OP capture (Vg=0.75 Vd=0.1; op1.raw, 73 vars).
// Optional bias arguments feed the native DC-core pass and print the SP-stage
// quantities (sp_ x rows) that the Python oracle transcription is diffed
// against by tools/compare_psp103_sp.py.
#include "devices/psp103_model.hpp"
#include "devices/psp103_temperature.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace gspice;

static int failures = 0;

static void cmp(const char* name, double mine, double ref, double tol) {
    double ratio = (ref != 0.0) ? mine / ref : 0.0;
    bool ok = (ref == 0.0) ? (std::abs(mine) <= tol) : (std::abs(ratio - 1.0) <= tol);
    printf("%-14s mine=% .8e  ref=% .8e  ratio=% .6f  %s\n", name, mine, ref, ratio, ok ? "OK" : "DIFF");
    if (!ok) ++failures;
}

namespace {

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
        {"RSW1", 50.0}, {"RSW2", 5.0e-2}, {"STRSO", -2.0}, {"RSBO", 0.0}, {"RSGO", 0.0},
        {"THESATO", 1e-6}, {"THESATL", 0.6}, {"THESATLEXP", 0.75},
        {"THESATW", -1e-2}, {"THESATLW", 0.0},
        {"STTHESATO", 1.5}, {"STTHESATL", -2.5e-2}, {"STTHESATW", -2e-2},
        {"STTHESATLW", -5e-3},
        {"THESATBO", 0.15}, {"THESATGO", 0.75},
        {"AXO", 20.0}, {"AXL", 0.2},
        {"ALPL", 7e-3}, {"ALPLEXP", 0.6}, {"ALPW", 5.0e-2},
        {"ALP1L1", 2.5e-2}, {"ALP1LEXP", 0.4}, {"ALP1L2", 0.1}, {"ALP1W", 8.5e-3},
        {"ALP2L1", 0.5}, {"ALP2LEXP", 0.0}, {"ALP2L2", 0.5}, {"ALP2W", -0.2},
        {"VPO", 0.25},
        {"STA2O", -0.5},
    });
}

} // namespace

int main(int argc, char** argv) {
    double vg = 0.75, vd = 0.1, vsb = 0.0;
    if (argc >= 2) vg = std::atof(argv[1]);
    if (argc >= 3) vd = std::atof(argv[2]);
    if (argc >= 4) vsb = std::atof(argv[3]);

    const auto model = refNmosModel();
    const auto prepared = model.prepare({{"W", 1e-6}, {"L", 0.1e-6}}, 21.0);
    const Psp103DeviceSetup s = psp103PrepareDevice(model, prepared, 21.0);

    printf("== my internals @ T=21C, W=1u L=0.1u (VACASK ref from op1.raw) ==\n");
    cmp("vfb_t", s.vfb_t, -1.103000000000000e+00, 1e-6);
    cmp("neff_i", s.neff, 7.744023323615161e+23, 1e-6);
    cmp("neffac_i", s.neffac, 6.195218658892129e+23, 1e-6);
    cmp("gfac_nud", s.gfac_nud, 1.000000000000000e-01, 1e-6);
    cmp("vsbnud", s.vsbnud, 0.000000000000000e+00, 0.0);
    cmp("dvsbnud", s.dvsbnud, 1.000000000000000e+00, 1e-6);
    cmp("vnsub", s.vnsub, 0.000000000000000e+00, 0.0);
    cmp("nslp", s.nslp, 5.000000000000000e-02, 1e-6);
    cmp("dnsub", s.dnsub, 0.000000000000000e+00, 0.0);
    cmp("dphib", s.dphib, 0.000000000000000e+00, 0.0);
    cmp("npol_i", s.npol, 1.500000000000000e+26, 1e-6);
    cmp("ct_i", s.ct, 1.972428037702950e-01, 1e-6);
    cmp("cf_i", s.cf, 6.152758132956150e-02, 1e-6);
    cmp("cfd_i", s.cfd, 0.000000000000000e+00, 0.0);
    cmp("cfb_i", s.cfb, 3.000000000000000e-01, 1e-6);
    cmp("psce_i", s.psce, 0.000000000000000e+00, 0.0);
    cmp("psce_b", s.psce_b, 0.000000000000000e+00, 0.0);
    cmp("psce_d", s.psce_d, 0.000000000000000e+00, 0.0);
    cmp("bet_i", s.bet_i, 3.503223284354082e-01 * s.cox_prime, 1e-6);
    cmp("mue_t", s.mue_t, 5.988873852719270e-01, 1e-6);
    cmp("themu_t", s.themu_t, 2.744452664310870e+00, 1e-6);
    cmp("cs_t", s.cs_t, 9.039668931005800e-03, 1e-6);
    cmp("xcor_t", s.xcor_t, 1.459291810903350e-01, 1e-6);
    cmp("rs_t", s.rs_t, 4.989926309788670e+01, 1e-6);
    cmp("theta_sat", s.theta_sat_t, 3.005798081035560e+00, 1e-6);
    cmp("thesat_b", s.thesat_b, 1.500000000000000e-01, 1e-6);
    cmp("thesat_g", s.thesat_g, 7.500000000000000e-01, 1e-6);
    cmp("ax_i", s.ax, 5.185185185185185e+00, 1e-6);
    cmp("alp_i", s.alp, 3.622627732612760e-02, 1e-6);
    cmp("alp1_i", s.alp1, 1.421307857822620e-02, 1e-6);
    cmp("alp2_i", s.alp2, 4.924439812402293e-02, 1e-6);
    cmp("vp_i", s.vp, 2.500000000000000e-01, 1e-6);
    cmp("a2_t", s.a2_t, 1.010147393327733e+01, 1e-6);

    printf("\nfailures: %d\n", failures);

    // --- DC-core pass (mirrors psp103_model.hpp input mapping) ----------
    Psp103DcCoreInputs in;
    in.phib = s.phib_dc;
    in.g0 = s.g0_dc;
    in.phit0 = s.phit0;
    in.vfb_t = s.vfb_t;
    in.kp = s.kp;
    in.cf = s.cf;
    in.cfd = s.cfd;
    in.cfb = s.cfb;
    in.psce = s.psce;
    in.psce_d = s.psce_d;
    in.psce_b = s.psce_b;
    in.dnsub = s.dnsub;
    in.vnsub = s.vnsub;
    in.nslp = s.nslp;
    in.rsb = s.rsb;
    in.rsg = s.rsg;
    in.ther = s.ther;
    in.e_eff0 = s.e_eff0;
    in.eta_mu = s.eta_mu;
    in.eta_mu1 = s.eta_mu1;
    in.mue_t = s.mue_t;
    in.themu_t = s.themu_t;
    in.cs_t = s.cs_t;
    in.xcor_t = s.xcor_t;
    in.thesat_b = s.thesat_b;
    in.thesat_g = s.thesat_g;
    in.theta_sat_t = s.theta_sat_t;
    in.ax = s.ax;
    in.alp = s.alp;
    in.alp1 = s.alp1;
    in.alp2 = s.alp2;
    in.vp = s.vp;
    in.bet = s.bet_i;
    in.phix = s.phix_dc;
    in.aphi = s.aphi_dc;
    in.phix1 = s.phix1_dc;
    in.vsbnud = s.vsbnud;
    in.dvsbnud = s.dvsbnud;
    in.gfac_nud = s.gfac_nud;
    in.us1 = s.us1;
    in.us21 = s.us21;
    in.sw_nud = s.sw_nud != 0;
    in.pmos = s.channel_type == -1;

    printf("\n== native SP pass @ Vg=%.4f Vd=%.4f Vsb=%.4f (for oracle diff) ==\n",
           vg, vd, vsb);
    const auto core = psp103DcCore(in, vg, vd, vsb);
    printf("setup phib % .8e\n", s.phib_dc);
    printf("setup g0    % .8e\n", s.g0_dc);
    printf("setup phit0 % .8e\n", s.phit0);
    printf("setup vfb_t % .8e\n", s.vfb_t);
    printf("setup kp    % .8e\n", s.kp);
    printf("setup cf    % .8e\n", s.cf);
    printf("setup cfd   % .8e\n", s.cfd);
    printf("setup cfb   % .8e\n", s.cfb);
    printf("setup psce  % .8e\n", s.psce);
    printf("setup psce_d % .8e\n", s.psce_d);
    printf("setup psce_b % .8e\n", s.psce_b);
    printf("setup dnsub % .8e\n", s.dnsub);
    printf("setup vnsub % .8e\n", s.vnsub);
    printf("setup nslp  % .8e\n", s.nslp);
    printf("setup rsb   % .8e\n", s.rsb);
    printf("setup rsg   % .8e\n", s.rsg);
    printf("setup ther  % .8e\n", s.ther);
    printf("setup e_eff0 % .8e\n", s.e_eff0);
    printf("setup eta_mu % .8e\n", s.eta_mu);
    printf("setup eta_mu1 % .8e\n", s.eta_mu1);
    printf("setup mue_t % .8e\n", s.mue_t);
    printf("setup themu_t % .8e\n", s.themu_t);
    printf("setup cs_t  % .8e\n", s.cs_t);
    printf("setup xcor_t % .8e\n", s.xcor_t);
    printf("setup thesat_b % .8e\n", s.thesat_b);
    printf("setup thesat_g % .8e\n", s.thesat_g);
    printf("setup theta_sat_t % .8e\n", s.theta_sat_t);
    printf("setup ax    % .8e\n", s.ax);
    printf("setup alp   % .8e\n", s.alp);
    printf("setup vp    % .8e\n", s.vp);
    printf("setup bet   % .8e\n", s.bet_i);
    printf("setup phix  % .8e\n", s.phix_dc);
    printf("setup aphi  % .8e\n", s.aphi_dc);
    printf("setup phix1 % .8e\n", s.phix1_dc);

    printf("sp xg      % .8e\n", core.xg.value);
    printf("sp x_ds    % .8e\n", core.x_ds.value);
    printf("sp x_m     % .8e\n", core.x_m.value);
    printf("sp dpsi    % .8e\n", core.dps.value);
    printf("sp voxm    % .8e\n", core.voxm.value);
    printf("sp vdsat   % .8e\n", core.vdsat.value);
    printf("sp vdse    % .8e\n", core.udse.value);
    printf("sp alpha   % .8e\n", core.alpha.value);
    printf("sp eta_p   % .8e\n", core.eta_p.value);
    printf("sp qeff1   % .8e\n", core.qeff1.value);
    printf("sp qim     % .8e\n", core.qim.value);
    printf("sp qbm     % .8e\n", core.qbm.value);
    printf("sp g_vsatinv % .8e\n", core.g_vsatinv.value);
    printf("sp g_delta_l % .8e\n", core.g_delta_l.value);
    printf("sp f_delta_l % .8e\n", core.f_delta_l.value);
    printf("sp h       % .8e\n", core.h.value);
    printf("sp ids     % .8e\n", core.ids.value);
    printf("sp gf      % .8e\n", core.gf.value);

    return failures == 0 ? 0 : 1;
}