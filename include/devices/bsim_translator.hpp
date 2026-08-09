#pragma once

#include <cmath>
#include <initializer_list>
#include <string>
#include <unordered_map>
#include <vector>

namespace gspice {

enum class BsimFamily { Unsupported, Bsim3, Bsim4 };

struct BsimTranslation {
    BsimFamily family = BsimFamily::Unsupported;
    int level = 0;
    std::unordered_map<std::string, double> parameters;
    std::vector<std::string> unknown_parameters;
    std::string reason;
    explicit operator bool() const { return family != BsimFamily::Unsupported && reason.empty(); }
};

class BsimModelTranslator {
public:
    static BsimTranslation translate(const std::unordered_map<std::string, double>& raw) {
        BsimTranslation out;
        for (const auto& [name, value] : raw) {
            const std::string key = normalize(name);
            if (!std::isfinite(value)) {
                out.reason = "parameter '" + key + "' is not finite";
                return out;
            }
            out.parameters[key] = value;
        }
        out.level = static_cast<int>(get(out.parameters, {"LEVEL"}, 0.0));
        if (out.level >= 49 && out.level <= 53) out.family = BsimFamily::Bsim3;
        else if (out.level == 54) out.family = BsimFamily::Bsim4;
        else {
            out.reason = "unsupported BSIM level; expected 49-53 for BSIM3 or 54 for BSIM4";
            return out;
        }
        canonicalize(out.parameters, "VTH0", {"VTH0", "VT0", "VTO"});
        canonicalize(out.parameters, "KP", {"KP", "BETA"});
        canonicalize(out.parameters, "U0", {"U0", "MOBMOD"});
        canonicalize(out.parameters, "VSAT", {"VSAT", "VSATT"});
        canonicalize(out.parameters, "TOXE", {"TOXE", "TOX"});
        canonicalize(out.parameters, "PCLM", {"PCLM", "LAMBDA"});
        canonicalize(out.parameters, "NFACTOR", {"NFACTOR", "N"});
        canonicalize(out.parameters, "GAMMA1", {"GAMMA1", "GAMMA"});
        canonicalize(out.parameters, "PHI", {"PHI", "PHIN"});
        canonicalize(out.parameters, "KF", {"KF"});
        canonicalize(out.parameters, "AF", {"AF"});
        canonicalize(out.parameters, "EF", {"EF"});
        for (const char* name : {"VTH0", "KP", "U0", "VSAT", "TOXE", "KF", "AF", "EF"}) {
            if (!out.parameters.count(name)) out.parameters[name] = defaults(name);
        }
        if (out.parameters["TOXE"] <= 0.0 || out.parameters["VSAT"] <= 0.0 ||
            out.parameters["KP"] < 0.0 || out.parameters["U0"] < 0.0 ||
            out.parameters["KF"] < 0.0 || out.parameters["AF"] <= 0.0 || out.parameters["EF"] <= 0.0) {
            out.reason = "BSIM core parameters have invalid signs";
        }
        return out;
    }

private:
    static std::string normalize(std::string name) {
        for (char& c : name) if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
        return name;
    }
    static double get(const std::unordered_map<std::string, double>& values,
                      std::initializer_list<const char*> names, double fallback) {
        for (const char* name : names) {
            const auto it = values.find(name);
            if (it != values.end()) return it->second;
        }
        return fallback;
    }
    static void canonicalize(std::unordered_map<std::string, double>& values,
                             const char* canonical, std::initializer_list<const char*> names) {
        if (values.count(canonical)) return;
        for (const char* name : names) {
            const auto it = values.find(name);
            if (it != values.end()) {
                values[canonical] = it->second;
                return;
            }
        }
    }
    static double defaults(const char* name) {
        if (std::string(name) == "VTH0") return 0.4;
        if (std::string(name) == "KP") return 120.0e-6;
        if (std::string(name) == "U0") return 0.05;
        if (std::string(name) == "VSAT") return 1.0e5;
        if (std::string(name) == "KF") return 0.0;
        if (std::string(name) == "AF" || std::string(name) == "EF") return 1.0;
        return 1.0e-8;
    }
};

} // namespace gspice
