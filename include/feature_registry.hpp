#ifndef GSPICE_FEATURE_REGISTRY_HPP
#define GSPICE_FEATURE_REGISTRY_HPP

#include <string>
#include <map>
#include <iostream>

#ifndef GSPICE_VERSION
#define GSPICE_VERSION "development"
#endif

namespace gspice {

enum class FeatureMaturity {
    Prototype,    // Code outline or option exists without full integration
    Wired,        // Parsed and connected, but lacks full numerical test verification
    Tested,       // Numerical verification test exists and passes
    Validated     // Validated against an independent benchmark
};

inline std::string maturityToString(FeatureMaturity maturity) {
    switch (maturity) {
        case FeatureMaturity::Prototype: return "prototype";
        case FeatureMaturity::Wired: return "wired";
        case FeatureMaturity::Tested: return "tested";
        case FeatureMaturity::Validated: return "validated";
    }
    return "unknown";
}

struct FeatureInfo {
    std::string name;
    FeatureMaturity maturity;
    bool available;
    std::string notes;
};

class FeatureRegistry {
public:
    static FeatureRegistry& instance() {
        static FeatureRegistry reg;
        return reg;
    }

    void registerFeature(const std::string& name, FeatureMaturity maturity, bool available, const std::string& notes = "") {
        features_[name] = {name, maturity, available, notes};
    }

    FeatureInfo getFeature(const std::string& name) const {
        auto it = features_.find(name);
        if (it != features_.end()) return it->second;
        return {name, FeatureMaturity::Prototype, false, "unregistered"};
    }

    std::string getCapabilitiesJson() const {
        std::string json = "{\n";
        json += "  \"name\": \"GSPICE\",\n";
        json += "  \"version\": \"" GSPICE_VERSION "\",\n";
        json += "  \"maturity\": \"academic-beta\",\n";
#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
        json += "  \"sparse_backend\": \"SuiteSparse-KLU\",\n";
#else
        json += "  \"sparse_backend\": \"internal-sparse\",\n";
#endif
        json += "  \"features\": {\n";
        bool first = true;
        for (const auto& kv : features_) {
            if (!first) json += ",\n";
            first = false;
            json += "    \"" + kv.first + "\": {\"maturity\": \"" + maturityToString(kv.second.maturity) +
                    "\", \"available\": " + (kv.second.available ? "true" : "false") + "}";
        }
        json += "\n  },\n";
        json += "  \"analyses\": {\n";
        json += "    \"op\": \"tested\", \"dc\": \"tested\", \"tran\": \"tested\", \"ac\": \"tested\",\n";
        json += "    \"noise\": \"tested\", \"stb\": \"prototype\", \"pss\": \"prototype\", \"pac\": \"prototype\", \"psspac\": \"prototype\", \"pssstb\": \"prototype\", \"pnoise\": \"prototype\", \"hb\": \"experimental\"\n";
        json += "  },\n";
        json += "  \"experimental\": [\"stb\", \"pss\", \"pac\", \"psspac\", \"pssstb\", \"pnoise\", \"hb\", \"fastspice\", \"multirate\", \"ticer\", \"binary_raw\", \"c_api\"],\n";
        json += "  \"outputs\": [\"spice-ascii-raw\", \"csv\"]\n";
        json += "}\n";
        return json;
    }

private:
    FeatureRegistry() {
        // Classify features strictly according to Phase 0-8 truth gates
        registerFeature("op", FeatureMaturity::Validated, true, "Operating point analysis");
        registerFeature("dc", FeatureMaturity::Validated, true, "DC sweep analysis");
        registerFeature("tran", FeatureMaturity::Validated, true, "Transient analysis");
        registerFeature("ac", FeatureMaturity::Validated, true, "AC small-signal analysis");
        registerFeature("native_compact_models", FeatureMaturity::Tested, true, "Native GMC/GSDI model registry with fail-closed unsupported models");
        registerFeature("bsim3_gsdi", FeatureMaturity::Tested, true, "Native BSIM3 level-49 OP/DC/AC/tran/noise route through GSDI/GMC with binned geometry, leakage, and verbose compatibility diagnostics");
        registerFeature("bsim4_gsdi", FeatureMaturity::Tested, true, "Native BSIM4 level-54 OP/DC/AC/tran/noise route through GSDI/GMC");
        registerFeature("bsim3_full_reference", FeatureMaturity::Prototype, false, "Full independent benchmark parity is not yet complete");
        registerFeature("bsim4_full_reference", FeatureMaturity::Prototype, false, "Full independent BSIM4.8.3 benchmark parity is not yet complete");
        registerFeature("gmc_generated_noise", FeatureMaturity::Tested, true, "GMC-generated white_noise/flicker_noise sources emit through GSDI noise vectors");
        registerFeature("juncap2_ideal_srh_bbt", FeatureMaturity::Tested, true, "Native JUNCAP2 ideal, SRH, and BBT current branches");
        registerFeature("juncap2_tat", FeatureMaturity::Tested, true, "Native smoothed TAT branch with derivative coverage");
        registerFeature("juncap2_avalanche", FeatureMaturity::Tested, true, "Native smoothed reverse-breakdown branch with derivative coverage");
        registerFeature("psp103_native", FeatureMaturity::Tested, true, "Native PSP103 OP/DC/AC/tran/noise path with parser and IHP alias coverage; not full OpenVAF equation parity");
        registerFeature("psp103_full", FeatureMaturity::Wired, false, "Full PSP103/OpenVAF equation parity is not complete; use psp103_native for the tested native subset");
        registerFeature("psp103_full_reference", FeatureMaturity::Prototype, false, "Full PSP103/OpenVAF equation parity and arbitrary OSDI runtime loading are not yet complete");
        registerFeature("ihp_psp103_matrix", FeatureMaturity::Tested, true, "IHP LV PSP OP/DC/AC/tran/noise, corner, temperature, geometry, and mismatch smoke coverage");
        registerFeature("ihp_psp_rf_smoke", FeatureMaturity::Tested, true, "IHP LV PSP RF model cards route through native PSP and pass small-signal smoke coverage");
        registerFeature("ihp_psp_rf_full", FeatureMaturity::Tested, true, "IHP LV PSP RF wrappers route natively with NG/M/DTA aliasing and RF gate-resistance stamping");
        registerFeature("ihp_passive_wrapper_smoke", FeatureMaturity::Tested, true, "IHP passive wrappers run natively with R3 effective geometry/contact/corner/temperature terms and CMIM area/perimeter/temperature capacitance");
        registerFeature("ihp_mosvar_cv_smoke", FeatureMaturity::Tested, true, "IHP HV svaricap routes through native voltage-dependent MOSVAR charge plus parasitic RC pieces");
        registerFeature("ihp_passive_wrapper_full", FeatureMaturity::Tested, true, "IHP resistor, CMIM, parasitic-cap, tap, and MOSVAR wrappers route natively without user-visible compatibility warnings");
        registerFeature("c_api", FeatureMaturity::Tested, true, "C API interface (SimulatorContext)");
        registerFeature("simulator_core", FeatureMaturity::Tested, true, "Decoupled simulator context");
        registerFeature("adjoint_sensitivity", FeatureMaturity::Tested, true, "Adjoint linear solver sensitivity engine");
        registerFeature("nport_sparam", FeatureMaturity::Tested, true, "Touchstone S-parameter Y-matrix conversion");
        registerFeature("btf_decomposition", FeatureMaturity::Tested, true, "SuiteSparse BTF block triangular decomposition");
        registerFeature("multitone_hb", FeatureMaturity::Tested, true, "Multi-tone Harmonic Balance lattice generator");
        registerFeature("python_zero_copy", FeatureMaturity::Tested, true, "Zero-copy NumPy array buffer exporter");
        registerFeature("provenance_stamping", FeatureMaturity::Tested, true, "Reproducibility provenance JSON metadata");
        registerFeature("hb", FeatureMaturity::Wired, false, "Harmonic balance (experimental)");
        registerFeature("stb", FeatureMaturity::Prototype, true, "Return-ratio stability analysis");
        registerFeature("pss", FeatureMaturity::Prototype, true, "Driven and autonomous variational-Newton shooting periodic steady state");
        registerFeature("pac", FeatureMaturity::Prototype, true, "Small-signal envelope AC");
        registerFeature("psspac", FeatureMaturity::Prototype, true, "PSS-orbit LPTV-coupled sideband AC");
        registerFeature("pssstb", FeatureMaturity::Prototype, true, "PSS-orbit LPTV-coupled stability analysis");
        registerFeature("pnoise", FeatureMaturity::Prototype, true, "PSS-orbit LPTV-coupled sideband periodic noise");
        registerFeature("fastspice", FeatureMaturity::Prototype, false, "FastSPICE partitioning (experimental)");
        registerFeature("multirate", FeatureMaturity::Prototype, false, "Multirate transient (experimental)");
        registerFeature("ticer", FeatureMaturity::Prototype, false, "TICER reduction (experimental)");
        registerFeature("binary_raw", FeatureMaturity::Wired, false, "Binary RAW exporter (experimental)");
    }

    std::map<std::string, FeatureInfo> features_;
};

} // namespace gspice

#endif // GSPICE_FEATURE_REGISTRY_HPP
