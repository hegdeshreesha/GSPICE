#ifndef GSPICE_GMC_JUNCAP_EXPRESS_HPP
#define GSPICE_GMC_JUNCAP_EXPRESS_HPP

#include "../gdi.hpp"
#include "juncap_express.hpp"
#include "juncap2_core.hpp"

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

namespace gspice {

inline double juncapParam(
    const std::unordered_map<std::string, double>& params,
    const char* name, double fallback) {
    const auto it = params.find(name);
    return it == params.end() ? fallback : it->second;
}

class GmcJuncapExpressInstance final : public GdiInstanceBase {
public:
    JuncapExpressInputs model;

    bool evaluate(const double* voltages, const GdiEvaluationRequest& request,
                  GdiEvaluationResult& result) override {
        result.clear();
        if (voltages == nullptr) return false;
        model.vak = voltages[0] - voltages[1];
        const auto evaluated = juncapExpressEvaluate(model);

        if (request.static_residual) {
            result.static_residual.push_back({0, evaluated.current});
            result.static_residual.push_back({1, -evaluated.current});
        }
        if (request.static_jacobian) {
            result.static_jacobian.push_back({0, 0, evaluated.dcurrent_dvak});
            result.static_jacobian.push_back({0, 1, -evaluated.dcurrent_dvak});
            result.static_jacobian.push_back({1, 0, -evaluated.dcurrent_dvak});
            result.static_jacobian.push_back({1, 1, evaluated.dcurrent_dvak});
        }
        if (request.dynamic_residual) {
            result.dynamic_residual.push_back({0, evaluated.charge});
            result.dynamic_residual.push_back({1, -evaluated.charge});
        }
        if (request.dynamic_jacobian) {
            result.dynamic_jacobian.push_back({0, 0, evaluated.capacitance});
            result.dynamic_jacobian.push_back({0, 1, -evaluated.capacitance});
            result.dynamic_jacobian.push_back({1, 0, -evaluated.capacitance});
            result.dynamic_jacobian.push_back({1, 1, evaluated.capacitance});
        }
        if (request.noise) {
            result.noise.push_back({0, 1, evaluated.noise_psd, "juncap-shot"});
        }
        return true;
    }

    std::size_t numTerminals() const override { return 2; }
    std::size_t numInternalNodes() const override { return 0; }
    std::size_t numStates() const override { return 0; }
};

class GmcJuncapExpressModule final : public GdiModuleBase {
public:
    std::string modelName() const override { return "juncap-express"; }

    std::vector<GdiParameter> listParameters() const override {
        return {
            {"SWJUNEXP", 1.0, "", "must be 1 for this native module", true},
            {"PHITD", 0.02585, "V", "thermal voltage", true},
            {"ISATFOR1", 0.0, "A", "first forward saturation current", true},
            {"MFOR1", 1.0, "", "first forward ideality factor", true},
            {"ISATFOR2", 0.0, "A", "second forward saturation current", true},
            {"MFOR2", 1.0, "", "second forward ideality factor", true},
            {"ISATREV", 0.0, "A", "reverse saturation current", true},
            {"MREV", 1.0, "", "reverse ideality factor", true},
            {"TYPE", 1.0, "", "junction polarity", true},
            {"MULT", 1.0, "", "multiplicity", true},
            {"FJUNQ", 0.0, "", "charge component threshold", true},
            {"AB", 1.0, "m2", "bottom junction area", true},
            {"LS", 0.0, "m", "STI-edge junction length", true},
            {"LG", 0.0, "m", "gate-edge junction length", true},
            {"CJORBOT", 0.0, "F/m2", "bottom zero-bias capacitance", true},
            {"VBIBOT", 1.0, "V", "bottom built-in voltage", true},
            {"PBOT", 0.5, "", "bottom grading coefficient", true},
            {"VJ", 0.0, "V", "conditioned junction voltage", false},
            {"DVJ", 1.0, "", "junction voltage derivative", false},
        };
    }

    std::unique_ptr<GdiInstanceBase> createInstance(
        const std::unordered_map<std::string, double>& modelParams,
        const std::unordered_map<std::string, double>& instanceParams) override {
        auto instance = std::make_unique<GmcJuncapExpressInstance>();
        auto& p = instance->model;
        if (juncapParam(modelParams, "SWJUNEXP", 1.0) != 1.0) return nullptr;
        p.phi_td = juncapParam(modelParams, "PHITD", p.phi_td);
        p.isat_for1 = juncapParam(modelParams, "ISATFOR1", p.isat_for1);
        p.m_for1 = juncapParam(modelParams, "MFOR1", p.m_for1);
        p.isat_for2 = juncapParam(modelParams, "ISATFOR2", p.isat_for2);
        p.m_for2 = juncapParam(modelParams, "MFOR2", p.m_for2);
        p.isat_rev = juncapParam(modelParams, "ISATREV", p.isat_rev);
        p.m_rev = juncapParam(modelParams, "MREV", p.m_rev);
        p.type = juncapParam(modelParams, "TYPE", p.type);
        p.mult = juncapParam(modelParams, "MULT", p.mult);
        p.fjunq = juncapParam(modelParams, "FJUNQ", p.fjunq);
        p.bottom.area = juncapParam(modelParams, "AB", p.bottom.area);
        p.bottom.cjo = juncapParam(modelParams, "CJORBOT", p.bottom.cjo);
        p.bottom.vbi = juncapParam(modelParams, "VBIBOT", p.bottom.vbi);
        p.bottom.p = juncapParam(modelParams, "PBOT", p.bottom.p);
        p.vj = juncapParam(instanceParams, "VJ", p.vj);
        p.dvj_dvak = juncapParam(instanceParams, "DVJ", p.dvj_dvak);
        return instance;
    }
};

class GmcJuncap2Instance final : public GdiInstanceBase {
public:
    Juncap2Inputs model;
    double type = 1.0;
    double mult = 1.0;

    bool evaluate(const double* voltages, const GdiEvaluationRequest& request,
                  GdiEvaluationResult& result) override {
        result.clear();
        if (voltages == nullptr) return false;
        model.vak = type * (voltages[0] - voltages[1]);
        const auto evaluated = juncap2Evaluate(model);
        if (!evaluated.supported) return false;
        const double scale = type * mult;
        const double current = scale * evaluated.current;
        const double conductance = mult * evaluated.dcurrent_dvak;
        if (request.static_residual) {
            result.static_residual.push_back({0, current});
            result.static_residual.push_back({1, -current});
        }
        if (request.static_jacobian) {
            result.static_jacobian.push_back({0, 0, conductance});
            result.static_jacobian.push_back({0, 1, -conductance});
            result.static_jacobian.push_back({1, 0, -conductance});
            result.static_jacobian.push_back({1, 1, conductance});
        }
        if (request.noise) result.noise.push_back({0, 1, 2.0 * 1.602176634e-19 * std::abs(current), "juncap2-shot"});
        return true;
    }

    std::size_t numTerminals() const override { return 2; }
    std::size_t numInternalNodes() const override { return 0; }
    std::size_t numStates() const override { return 0; }
};

class GmcJuncap2Module final : public GdiModuleBase {
public:
    std::string modelName() const override { return "juncap2"; }
    std::vector<GdiParameter> listParameters() const override {
        return {
            {"PHITD", 0.02585, "V", "thermal voltage", true},
            {"IDSAT", 0.0, "A", "ideal diode saturation current", true},
            {"CSRH", 0.0, "A/V", "SRH coefficient", true},
            {"CTAT", 0.0, "A/V", "trap-assisted tunneling coefficient", true},
            {"CBBT", 0.0, "A/V", "band-to-band tunneling coefficient", true},
            {"PHITR", 0.02585, "V", "BBT smoothing voltage", true},
            {"DELTAVBI", 0.0, "V", "BBT built-in-voltage shift", true},
            {"FBBTR", 1.0, "V/m", "BBT reference field", true},
            {"STFBBT", 0.0, "1/K", "BBT temperature coefficient", true},
            {"VBR", 1.0e9, "V", "reverse breakdown voltage", true},
            {"VMAX", 1.034, "V", "exponential continuation voltage", true},
            {"VBIBOT", 0.7, "V", "minimum built-in voltage", true},
            {"VBIRBOT", 0.7, "V", "junction built-in reference", true},
            {"PBOT", 0.5, "", "grading coefficient", true},
            {"XJUN", 1.0, "m", "junction depth", true},
            {"CJORBOT", 1.0e-3, "F/m2", "zero-bias junction capacitance", true},
            {"FTD", 1.0, "", "temperature scaling", true},
            {"TYPE", 1.0, "", "junction polarity", true},
            {"MULT", 1.0, "", "multiplicity", true}
        };
    }
    std::unique_ptr<GdiInstanceBase> createInstance(
        const std::unordered_map<std::string, double>& modelParams,
        const std::unordered_map<std::string, double>&) override {
        auto instance = std::make_unique<GmcJuncap2Instance>();
        auto& p = instance->model;
        p.phi_td = juncapParam(modelParams, "PHITD", p.phi_td);
        p.idsat = juncapParam(modelParams, "IDSAT", p.idsat);
        p.csrh = juncapParam(modelParams, "CSRH", p.csrh);
        p.ctat = juncapParam(modelParams, "CTAT", p.ctat);
        p.cbbt = juncapParam(modelParams, "CBBT", p.cbbt);
        p.phi_tr = juncapParam(modelParams, "PHITR", p.phi_tr);
        p.delta_vbi = juncapParam(modelParams, "DELTAVBI", p.delta_vbi);
        p.fbbtr = juncapParam(modelParams, "FBBTR", p.fbbtr);
        p.stfbbt = juncapParam(modelParams, "STFBBT", p.stfbbt);
        p.vbr = juncapParam(modelParams, "VBR", p.vbr);
        p.v_max = juncapParam(modelParams, "VMAX", p.v_max);
        p.vbi = juncapParam(modelParams, "VBIBOT", p.vbi);
        p.vbi_min = juncapParam(modelParams, "VBIBOT", p.vbi_min);
        p.vbir = juncapParam(modelParams, "VBIRBOT", p.vbir);
        p.p = juncapParam(modelParams, "PBOT", p.p);
        p.xjun = juncapParam(modelParams, "XJUN", p.xjun);
        p.cjor = juncapParam(modelParams, "CJORBOT", p.cjor);
        p.ftd = juncapParam(modelParams, "FTD", p.ftd);
        instance->type = juncapParam(modelParams, "TYPE", instance->type);
        instance->mult = juncapParam(modelParams, "MULT", instance->mult);
        if (!juncap2Evaluate(p).supported) return nullptr;
        return instance;
    }
};

} // namespace gspice

#endif
