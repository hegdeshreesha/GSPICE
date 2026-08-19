#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace gspice {

struct Bsim3Validation {
    bool valid = false;
    std::string reason;
    explicit operator bool() const { return valid; }
};

struct Bsim3PreparedModel {
    Bsim3Validation validation;
    double level = 49.0;
    double temperature_k = 300.15;
    double nominal_temperature_k = 300.15;
    double width = 1.0e-6;
    double length = 1.0e-6;
    double leff = 1.0e-6;
    double weff = 1.0e-6;
    double vth0 = 0.4;
    double kp = 120.0e-6;
    double u0 = 0.05;
    double vsat = 1.0e5;
    double gamma1 = 0.0;
    double phi = 0.7;
    double nfactor = 1.0;
    double k1 = 0.5;
    double k2 = 0.0;
    double lambda = 0.0;
    double cox = 3.453133e-3;
    double toxe = 1.0e-8;
    double beta = 120.0e-6;
    double cdep0 = 1.0e-3;
    double ua = 2.25e-9;
    double ub = 5.9e-19;
    double uc = -4.7e-11;
    double kf = 0.0;
    double af = 1.0;
    double ef = 1.0;
    double is = 0.0;
    double js = 0.0;
    double jsw = 0.0;
    double xti = 3.0;
    double eg = 1.11;
    double nj = 1.0;
    double aigc = 0.0;
    double bigc = 0.0;
    double cigc = 0.0;
    double agidl = 0.0;
    double bgidl = 2.3e9;
    double cgidl = 0.5;
    double egidl = 0.8;
    double agisl = 0.0;
    double bgisl = 2.3e9;
    double cgisl = 0.5;
    double egisl = 0.8;
};

class Bsim3ParameterSet {
public:
    static Bsim3ParameterSet from(const std::unordered_map<std::string, double>& values) {
        Bsim3ParameterSet result;
        for (const auto& [name, value] : values) result.values_[normalize(name)] = value;
        return result;
    }

    Bsim3PreparedModel prepare(double width, double length, double temperature_c = 27.0) const {
        Bsim3PreparedModel out;
        out.width = width;
        out.length = length;
        out.temperature_k = temperature_c + 273.15;
        out.nominal_temperature_k = get({"TNOM"}, 27.0) + 273.15;
        out.level = get({"LEVEL"}, 49.0);
        out.vth0 = get({"VTH0", "VT0", "VTO"}, 0.4);
        const bool hasExplicitBeta = values_.count("KP") > 0 || values_.count("BETA") > 0;
        out.kp = get({"KP", "BETA"}, std::numeric_limits<double>::quiet_NaN());
        out.u0 = get({"U0", "MOBMOD"}, 0.05);
        // BSIM cards conventionally use cm^2/(V*s); the native evaluator uses SI.
        if (out.u0 > 1.0) out.u0 *= 1.0e-4;
        out.vsat = get({"VSAT", "VSATT"}, 1.0e5);
        out.gamma1 = get({"GAMMA1", "GAMMA"}, 0.0);
        out.phi = get({"PHIN", "PHI"}, 0.7);
        out.nfactor = get({"NFACTOR", "N"}, 1.0);
        // The DC core uses a single K1-style body coefficient.  Cards that set
        // only the legacy GAMMA (e.g. the ngspice BSIM3 reference decks) are
        // mapped onto that coefficient; an explicit K1 always wins.
        out.k1 = values_.count("K1") > 0 ? get({"K1"}, 0.5)
            : (out.gamma1 != 0.0 ? out.gamma1 : 0.5);
        out.k2 = get({"K2"}, 0.0);
        out.lambda = get({"PCLM", "LAMBDA"}, 0.0);
        out.ua = get({"UA"}, 2.25e-9);
        out.ub = get({"UB"}, 5.9e-19);
        out.uc = get({"UC"}, -4.7e-11);
        const double dl = get({"XL", "DL"}, 0.0);
        const double dw = get({"XW", "DW"}, 0.0);
        const double toxe = get({"TOXE", "TOX"}, 1.0e-8);
        out.toxe = toxe;
        out.kf = get({"KF"}, 0.0);
        out.af = get({"AF"}, 1.0);
        out.ef = get({"EF"}, 1.0);
        out.leff = length - 2.0 * dl;
        out.weff = width - 2.0 * dw;
        out.vth0 = getBinned("VTH0", out.vth0, out.leff, out.weff);
        out.k1 = getBinned("K1", out.k1, out.leff, out.weff);
        out.k2 = getBinned("K2", out.k2, out.leff, out.weff);
        out.lambda = getBinned("PCLM", out.lambda, out.leff, out.weff);
        out.nfactor = getBinned("NFACTOR", out.nfactor, out.leff, out.weff);
        out.ua = getBinned("UA", out.ua, out.leff, out.weff);
        out.ub = getBinned("UB", out.ub, out.leff, out.weff);
        out.uc = getBinned("UC", out.uc, out.leff, out.weff);
        out.vsat = getBinned("VSAT", out.vsat, out.leff, out.weff);
        out.cox = 3.453133e-11 / std::max(toxe, 1.0e-12);
        out.cdep0 = std::sqrt(1.602176634e-19 * 1.03594e-10 *
                              1.7e17 * 1.0e6 /
                              (2.0 * std::max(out.phi, 1.0e-6)));
        if (!hasExplicitBeta || !std::isfinite(out.kp)) {
            out.kp = out.u0 * out.cox;
        }
        out.beta = out.kp * out.weff / std::max(out.leff, 1.0e-15);
        out.is = get({"IS"}, 0.0);
        out.js = get({"JS"}, 0.0);
        out.jsw = get({"JSW"}, 0.0);
        out.xti = get({"XTI"}, 3.0);
        out.eg = get({"EG"}, 1.11);
        out.nj = get({"NJ"}, 1.0);
        out.aigc = get({"AIGC"}, 0.0);
        out.bigc = get({"BIGC"}, 0.0);
        out.cigc = get({"CIGC"}, 0.0);
        out.agidl = get({"AGIDL"}, 0.0);
        out.bgidl = get({"BGIDL"}, 2.3e9);
        out.cgidl = get({"CGIDL"}, 0.5);
        out.egidl = get({"EGIDL"}, 0.8);
        out.agisl = get({"AGISL"}, out.agidl);
        out.bgisl = get({"BGISL"}, out.bgidl);
        out.cgisl = get({"CGISL"}, out.cgidl);
        out.egisl = get({"EGISL"}, out.egidl);
        out.is = getBinned("IS", out.is, out.leff, out.weff);
        out.js = getBinned("JS", out.js, out.leff, out.weff);
        out.jsw = getBinned("JSW", out.jsw, out.leff, out.weff);
        out.nj = getBinned("NJ", out.nj, out.leff, out.weff);
        out.aigc = getBinned("AIGC", out.aigc, out.leff, out.weff);
        out.bigc = getBinned("BIGC", out.bigc, out.leff, out.weff);
        out.cigc = getBinned("CIGC", out.cigc, out.leff, out.weff);
        out.agidl = getBinned("AGIDL", out.agidl, out.leff, out.weff);
        out.bgidl = getBinned("BGIDL", out.bgidl, out.leff, out.weff);
        out.cgidl = getBinned("CGIDL", out.cgidl, out.leff, out.weff);
        out.egidl = getBinned("EGIDL", out.egidl, out.leff, out.weff);
        out.agisl = getBinned("AGISL", out.agisl, out.leff, out.weff);
        out.bgisl = getBinned("BGISL", out.bgisl, out.leff, out.weff);
        out.cgisl = getBinned("CGISL", out.cgisl, out.leff, out.weff);
        out.egisl = getBinned("EGISL", out.egisl, out.leff, out.weff);
        {
            const double vt_nominal = 8.617333262e-5 * out.nominal_temperature_k;
            const double vt_temperature = 8.617333262e-5 * out.temperature_k;
            const double junctionTemperatureScale =
                std::pow(std::max(out.temperature_k / out.nominal_temperature_k, 1.0e-12), out.xti) *
                std::exp(std::clamp(out.eg / vt_nominal - out.eg / vt_temperature, -700.0, 80.0));
            out.is *= junctionTemperatureScale;
            out.js *= junctionTemperatureScale;
            out.jsw *= junctionTemperatureScale;
        }

        if (out.level < 49.0 || out.level > 49.99) out.validation.reason = "only BSIM3 level 49 is supported";
        else if (!std::isfinite(out.temperature_k) || out.temperature_k <= 0.0 ||
                 !std::isfinite(out.nominal_temperature_k) || out.nominal_temperature_k <= 0.0)
            out.validation.reason = "temperature and TNOM must be finite and above absolute zero";
        else if (!std::isfinite(width) || width <= 0.0 || !std::isfinite(length) || length <= 0.0) out.validation.reason = "W and L must be finite and positive";
        else if (!(out.leff > 0.0) || !(out.weff > 0.0)) out.validation.reason = "effective W/L became non-positive";
        else if (!std::isfinite(out.kp) || out.kp < 0.0 || !std::isfinite(out.vsat) || out.vsat <= 0.0) out.validation.reason = "KP must be nonnegative and VSAT positive";
        else if (!std::isfinite(out.toxe) || out.toxe <= 0.0) out.validation.reason = "TOXE must be finite and positive";
        else if (!std::isfinite(out.kf) || out.kf < 0.0 || !std::isfinite(out.af) || out.af <= 0.0 || !std::isfinite(out.ef) || out.ef <= 0.0) out.validation.reason = "KF must be nonnegative and AF/EF positive";
        else if (!std::isfinite(out.is) || out.is < 0.0 || !std::isfinite(out.js) ||
                 out.js < 0.0 || !std::isfinite(out.jsw) || out.jsw < 0.0 ||
                 !std::isfinite(out.xti) || !std::isfinite(out.eg) || out.eg < 0.0 ||
                 !std::isfinite(out.nj) || out.nj <= 0.0)
            out.validation.reason = "BSIM3 junction leakage parameters are invalid";
        else if (!std::isfinite(out.aigc) || out.aigc < 0.0 || !std::isfinite(out.bigc) ||
                 !std::isfinite(out.cigc) || !std::isfinite(out.agidl) || out.agidl < 0.0 ||
                 !std::isfinite(out.bgidl) || out.bgidl < 0.0 || !std::isfinite(out.cgidl) ||
                 out.cgidl < 0.0 || !std::isfinite(out.egidl) || out.egidl < 0.0 ||
                 !std::isfinite(out.agisl) || out.agisl < 0.0 || !std::isfinite(out.bgisl) ||
                 out.bgisl < 0.0 || !std::isfinite(out.cgisl) || out.cgisl < 0.0 ||
                 !std::isfinite(out.egisl) || out.egisl < 0.0)
            out.validation.reason = "BSIM3 gate/GIDL leakage parameters are invalid";
        else out.validation.valid = true;
        return out;
    }

private:
    double getBinned(const char* name, double fallback, double leff, double weff) const {
        const std::string base(name);
        const auto value = [&](const std::string& key, double defaultValue) {
            const auto it = values_.find(key);
            return it == values_.end() ? defaultValue : it->second;
        };
        const double p0 = value(base, fallback);
        const double pl = value("L" + base, 0.0);
        const double pw = value("W" + base, 0.0);
        const double pp = value("P" + base, 0.0);
        return p0 + pl / leff + pw / weff + pp / (leff * weff);
    }

    static std::string normalize(std::string name) {
        for (char& c : name) if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
        return name;
    }
    double get(std::initializer_list<const char*> names, double fallback) const {
        for (const char* name : names) {
            const auto it = values_.find(name);
            if (it != values_.end()) return it->second;
        }
        return fallback;
    }
    double toxe() const { return get({"TOXE", "TOX"}, 1.0e-8); }
    std::unordered_map<std::string, double> values_;
};

struct Bsim3ImplementedParameters {
    static bool isSupported(std::string_view name) {
        std::string upper(name);
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char c) {
                           return static_cast<char>(std::toupper(c));
                       });
        if (supported().find(upper) != supported().end()) return true;
        if (upper.size() > 1 && (upper[0] == 'L' || upper[0] == 'W' || upper[0] == 'P')) {
            return supported().find(upper.substr(1)) != supported().end();
        }
        return false;
    }

    static const std::unordered_set<std::string>& supported() {
        static const std::unordered_set<std::string> supported = {
            "AF", "AGIDL", "AGISL", "AIGC", "BETA", "BGIDL", "BGISL",
            "BIGC", "CGIDL", "CGISL", "CIGC", "CJSW", "CJ", "DL", "DW",
            "EF", "EG", "EGIDL", "EGISL", "GAMMA", "GAMMA1", "IS",
            "JS", "JSW", "K1", "K2", "KF", "KP", "LAMBDA", "LEVEL",
            "N", "NFACTOR", "NJ", "PCLM", "PHI", "PHIN", "TNOM",
            "TOX", "TOXE", "U0", "UA", "UB", "UC", "VSAT", "VSATT",
            "VT0", "VTH0", "VTO", "XL", "XTI", "XW",
        };
        return supported;
    }
};

} // namespace gspice
