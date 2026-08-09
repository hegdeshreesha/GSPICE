#include "devices/psp103_model.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

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
        {"TOXOVO", 1.5e-9}, {"TOXOVDO", 2e-9},
        {"LOV", 10e-9}, {"LOVD", 0.0},
        {"NOVO", 7.5e25}, {"NOVDO", 5e25},
        {"CFL", 3e-4}, {"CFLEXP", 2.0}, {"CFW", 5e-3}, {"CFBO", 0.3},
        {"CFRW", 5e-17}, {"CFRDW", 0.0}, {"CGBOVL", 0.0},
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
    });
}

} // namespace

int main() {
    using namespace gspice;
    const auto model = refNmosModel();
    Psp103Mosfet mos("M1", 0, 1, 2, 3, model, {{"W", 1e-6}, {"L", 0.1e-6}}, 21.0);

    // VACASK .OP reference (op1.raw), bias Vd=0.1, Vg=0.75, Vs=Vb=0, T=21 C.
    // sigVds>0, so the loadDynamic OP block that keeps Qd/Qs as-is applies
    // (PSP103_module.include lines 2965-2982). Derived caps use these combos:
    //   cdg=-dQd/dVg cdb=-dQd/dVb cds=cdd-cdg-cdb
    //   cgd=-dQg/dVd cgb=-dQg/dVb cgs=cgg-cgd-cgb
    //   csd=-dQs/dVd csb=-dQs/dVb css=csg+csd+csb
    //   cbd=-dQb/dVd cbg=-dQb/dVg cbs=cbb-cbd-cbg
    VectorReal x(4);
    x[0] = 0.1;   // Vd
    x[1] = 0.75;  // Vg
    x[2] = 0.0;   // Vs
    x[3] = 0.0;   // Vb
    const auto c = mos.intrinsicAndOverlapCharges(x);

    const double cdd = c.qd.derivative[0];
    const double cdg = -c.qd.derivative[1];
    const double cdb = -c.qd.derivative[3];
    const double cds = cdd - cdg - cdb;

    const double cgd = -c.qg.derivative[0];
    const double cgg = c.qg.derivative[1];
    const double cgb = -c.qg.derivative[3];
    const double cgs = cgg - cgd - cgb;

    const double csd = -c.qs.derivative[0];
    const double csg = -c.qs.derivative[1];
    const double csb = -c.qs.derivative[3];
    const double css = csg + csd + csb;

    const double cbd = -c.qb.derivative[0];
    const double cbg = -c.qb.derivative[1];
    const double cbb = c.qb.derivative[3];
    const double cbs = cbb - cbd - cbg;

    struct Probe {
        const char* name;
        double got;
        double ref;
    };
    const std::vector<Probe> probes = {
        {"cgg", cgg, 1.087866262081805e-15},
        {"cgs", cgs, 6.403404986218363e-16},
        {"cgd", cgd, 4.178879603243965e-16},
        {"cgb", cgb, 2.963780313557200e-17},
        {"cdd", cdd, 3.421843146469608e-16},
        {"cds", cds, -2.618001204265176e-16},
        {"cdg", cdg, 5.208417443440902e-16},
        {"cdb", cdb, 8.314269072938812e-17},
        {"css", css, 4.653195128504072e-16},
        {"csg", csg, 5.258158779534403e-16},
        {"csb", csb, 8.501569260883162e-17},
        {"csd", csd, -1.455120577118647e-16},
{"cbs", cbs, 8.677913465508833e-17},
        {"cbb", cbb, 8.677913465508833e-17 + 6.980841203442896e-17 + 4.120863978427443e-17},
        {"cbg", cbg, 4.120863978427443e-17},
        {"cbd", cbd, 6.980841203442896e-17},
        {"cgsol", c.qfgs.derivative[1], 2.738533302338999e-16},
        {"cgdol", c.qfgd.derivative[1], 2.724679529935561e-16},
    };

    std::printf("%6s  %14s  %14s  %9s\n", "name", "gspice", "ref", "ratio");
    int failures = 0;
    for (const auto& p : probes) {
        const double ratio = std::abs(p.ref) > 0.0 ? p.got / p.ref : 0.0;
        std::printf("%6s  %14.6e  %14.6e  %9.5f\n", p.name, p.got, p.ref, ratio);
        if (std::abs(ratio - 1.0) > 5e-3) ++failures;
    }
    std::printf("FAILURES: %d\n", failures);
    return failures == 0 ? 0 : 1;
}