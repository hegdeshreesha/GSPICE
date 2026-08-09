#include "feature_registry.hpp"
#include <cassert>

int main() {
    const auto& registry = gspice::FeatureRegistry::instance();
    const auto juncap = registry.getFeature("juncap2_ideal_srh_bbt");
    assert(juncap.available && juncap.maturity == gspice::FeatureMaturity::Tested);
    const auto tat = registry.getFeature("juncap2_tat");
    assert(tat.available && tat.maturity == gspice::FeatureMaturity::Tested);
    const auto avalanche = registry.getFeature("juncap2_avalanche");
    assert(avalanche.available && avalanche.maturity == gspice::FeatureMaturity::Tested);
    const auto psp = registry.getFeature("psp103_native");
    assert(psp.available && psp.maturity == gspice::FeatureMaturity::Tested);
    const auto psp_full = registry.getFeature("psp103_full");
    assert(!psp_full.available && psp_full.maturity == gspice::FeatureMaturity::Wired);
    const auto psp_ref = registry.getFeature("psp103_full_reference");
    assert(!psp_ref.available && psp_ref.maturity == gspice::FeatureMaturity::Prototype);
    const auto bsim3 = registry.getFeature("bsim3_gsdi");
    assert(bsim3.available && bsim3.maturity == gspice::FeatureMaturity::Tested);
    const auto bsim4 = registry.getFeature("bsim4_gsdi");
    assert(bsim4.available && bsim4.maturity == gspice::FeatureMaturity::Tested);
    const auto bsim4_ref = registry.getFeature("bsim4_full_reference");
    assert(!bsim4_ref.available && bsim4_ref.maturity == gspice::FeatureMaturity::Prototype);
    const auto generated_noise = registry.getFeature("gmc_generated_noise");
    assert(generated_noise.available && generated_noise.maturity == gspice::FeatureMaturity::Tested);
    const auto rf_smoke = registry.getFeature("ihp_psp_rf_smoke");
    assert(rf_smoke.available && rf_smoke.maturity == gspice::FeatureMaturity::Tested);
    const auto rf_full = registry.getFeature("ihp_psp_rf_full");
    assert(rf_full.available && rf_full.maturity == gspice::FeatureMaturity::Tested);
    const auto passive_smoke = registry.getFeature("ihp_passive_wrapper_smoke");
    assert(passive_smoke.available && passive_smoke.maturity == gspice::FeatureMaturity::Tested);
    const auto mosvar_cv = registry.getFeature("ihp_mosvar_cv_smoke");
    assert(mosvar_cv.available && mosvar_cv.maturity == gspice::FeatureMaturity::Tested);
    const auto passive_full = registry.getFeature("ihp_passive_wrapper_full");
    assert(passive_full.available && passive_full.maturity == gspice::FeatureMaturity::Tested);
}
