#include "devices/psp103_model.hpp"

#include <algorithm>
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
    });
}

} // namespace

int main() {
    using namespace gspice;
    const auto model = refNmosModel();
    Psp103Mosfet mos("M1", 0, 1, 2, 3, model, {{"W", 1e-6}, {"L", 0.1e-6}}, 21.0);

    // ngspice reference Id (Vd=0.1, Vb=0, Vs=0, T=21C): psp_ref_idvg.raw
    const std::vector<double> refVg = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                                       0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
                                       0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
                                       1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35,
                                       1.40, 1.45, 1.50};
    const std::vector<double> refId = {2.798872181e-09, 1.093550907e-08, 4.292254664e-08,
                                       1.682544634e-07, 6.454529279e-07, 2.286294201e-06,
                                       6.784016912e-06, 1.586210528e-05, 2.964320107e-05,
                                       4.655416315e-05, 6.475628261e-05, 8.300684720e-05,
                                       1.006196178e-04, 1.172433758e-04, 1.327120092e-04,
                                       1.469624528e-04, 1.599904060e-04, 1.718255770e-04,
                                       1.825172459e-04, 1.921255333e-04, 2.007159691e-04,
                                       2.083560523e-04, 2.151130513e-04, 2.210526015e-04,
                                       2.262378256e-04, 2.307288046e-04, 2.345822852e-04,
                                       2.378515494e-04, 2.405863942e-04, 2.430697660e-04,
                                       2.448505700e-04};

    DaeRequest request;
    request.staticResidual = true;
    request.staticJacobian = true;
    constexpr double maxRelError = 1.0e-3;
    double worstRelError = 0.0;
    double worstVg = 0.0;
    bool failed = false;
    std::printf("     Vg        Id(gspice)     Id(ref)       ratio\n");
    for (size_t i = 0; i < refVg.size(); ++i) {
        VectorReal x(4);
        x[0] = 0.1;
        x[1] = refVg[i];
        x[2] = 0.0;
        x[3] = 0.0;
        DaeEvaluation evaluation;
        if (!mos.evaluateDae(x, request, evaluation)) {
            std::printf("%10.3f  EVAL FAIL\n", refVg[i]);
            failed = true;
            continue;
        }
        double drainResidual = 0.0;
        double sourceResidual = 0.0;
        for (const auto& term : evaluation.staticResidual) {
            if (term.equation == 0) drainResidual += term.value;
            if (term.equation == 2) sourceResidual += term.value;
        }
        const double id = drainResidual;
        const double relError = std::abs(id - refId[i]) / std::max(std::abs(refId[i]), 1.0e-15);
        if (relError > worstRelError) {
            worstRelError = relError;
            worstVg = refVg[i];
        }
        std::printf("%10.3f  %12.5e  %12.5e  %9.4f\n", refVg[i], id, refId[i],
                    id / refId[i]);
        (void)sourceResidual;
    }
    std::printf("PSP103 IdVg probe parity: points=%zu max_rel_err=%.6g at vg=%.3f limit=%.6g\n",
                refVg.size(), worstRelError, worstVg, maxRelError);
    return failed || worstRelError > maxRelError ? 1 : 0;
}
