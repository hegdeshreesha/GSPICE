#ifndef GSPICE_PSP103_PARAMETERS_HPP
#define GSPICE_PSP103_PARAMETERS_HPP

#include "gsdi.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <string>
#include <unordered_set>
#include <vector>

namespace gspice {

struct Psp103Geometry {
    double width_m = 1e-6;
    double length_m = 1e-6;
    double fingers = 1.0;
    double multiplier = 1.0;
};

struct Psp103Temperature {
    double ambient_c = 27.0;
    double nominal_c = 27.0;
    double kelvin = 300.15;
    double thermal_voltage = 0.025852;
};

struct Psp103PreparedModel {
    int polarity = 1;
    Psp103Geometry geometry;
    Psp103Temperature temperature;
    double effective_width_m = 1e-6;
    double effective_length_m = 1e-6;
    double width_length_ratio = 1.0;

    bool valid() const {
        return (polarity == 1 || polarity == -1) &&
               std::isfinite(effective_width_m) && effective_width_m > 0.0 &&
               std::isfinite(effective_length_m) && effective_length_m > 0.0 &&
               std::isfinite(geometry.fingers) && geometry.fingers > 0.0 &&
               std::isfinite(geometry.multiplier) && geometry.multiplier > 0.0 &&
               std::isfinite(temperature.kelvin) && temperature.kelvin > 0.0 &&
               std::isfinite(temperature.thermal_voltage) &&
               temperature.thermal_voltage > 0.0;
    }
};

struct Psp103Validation {
    bool valid = false;
    std::vector<std::string> missing;

    explicit operator bool() const { return valid; }
};

class Psp103ParameterSet {
public:
    static Psp103ParameterSet from(const GsdiParamMap& raw) {
        Psp103ParameterSet result;
        for (const auto& [name, value] : raw) {
            result.values_[normalize(name)] = value;
        }
        return result;
    }

    bool has(const std::string& name) const {
        return values_.find(normalize(name)) != values_.end();
    }

    double get(const std::string& name, double fallback) const {
        const auto it = values_.find(normalize(name));
        return it == values_.end() ? fallback : it->second;
    }

    double alias(double fallback, std::initializer_list<const char*> names) const {
        for (const char* name : names) {
            const auto it = values_.find(normalize(name));
            if (it != values_.end()) return it->second;
        }
        return fallback;
    }

    double sum(std::initializer_list<const char*> names) const {
        double total = 0.0;
        for (const char* name : names) {
            const auto it = values_.find(normalize(name));
            if (it != values_.end() && std::isfinite(it->second)) total += it->second;
        }
        return total;
    }

    Psp103PreparedModel prepare(
        const GsdiParamMap& instance,
        double ambient_c) const {
        Psp103PreparedModel prepared;
        prepared.polarity = alias(1.0, {"TYPE"}) < 0.0 ? -1 : 1;
        const double device_rise_c = instanceValue(instance, {"DTA", "TRISE"}, 0.0);
        prepared.temperature.ambient_c = ambient_c + device_rise_c;
        prepared.temperature.nominal_c = alias(27.0, {"TNOM", "TNOM_C"});
        prepared.temperature.kelvin = prepared.temperature.ambient_c + 273.15;
        prepared.temperature.thermal_voltage =
            8.617333262145e-5 * prepared.temperature.kelvin;

        prepared.geometry.width_m = instanceValue(instance, {"W", "WIDTH"}, 1e-6);
        prepared.geometry.length_m = instanceValue(instance, {"L", "LENGTH"}, 1e-6);
        prepared.geometry.fingers = instanceValue(instance, {"NF", "NFIN", "NG"}, 1.0);
        prepared.geometry.multiplier = instanceValue(instance, {"M", "MULT"}, 1.0);
        prepared.effective_width_m = prepared.geometry.width_m *
            prepared.geometry.fingers * prepared.geometry.multiplier;
        prepared.effective_length_m = prepared.geometry.length_m;
        prepared.width_length_ratio = prepared.effective_width_m /
            std::max(prepared.effective_length_m, 1e-30);
        return prepared;
    }

    Psp103Validation validateIntrinsicCore() const {
        Psp103Validation result;
        const std::initializer_list<std::initializer_list<const char*>> groups = {
            {"VFB", "VFB0", "VFBO"}, {"TOX", "TOXO"},
            {"EPSROX", "EPSROXO"}, {"NDEP", "NDEPO"},
            {"BET", "BETO"}, {"PHIB", "PHIBO"}};
        for (const auto names : groups) {
            bool found = false;
            for (const char* name : names) {
                if (has(name)) {
                    found = true;
                    break;
                }
            }
            if (!found) result.missing.emplace_back(*names.begin());
        }
        result.valid = result.missing.empty();
        return result;
    }

    bool featureRequested(std::initializer_list<const char*> switches) const {
        for (const char* name : switches) {
            if (get(name, 0.0) != 0.0) return true;
        }
        return false;
    }

    static bool handlesModelParameter(std::string name) {
        static const std::unordered_set<std::string> handled = {
            "A1L", "A1O", "A1W", "A2O", "A3L", "A3O", "A3W", "A4L", "A4O", "A4W",
            "AF", "AF0", "AGIDL", "AGIDLD", "ALP", "ALP1",
            "AGIDLDW", "AGIDLW", "ALPNOI", "ALP1L1", "ALP1L2", "ALP1LEXP", "ALP1W", "ALP2", "ALP2L1",
            "ALP2L2", "ALP2LEXP", "ALP2W", "ALPL", "ALPLEXP", "ALPW",
            "AX", "AXL", "AXO", "BET", "BETA", "BETEDGEW", "BETN", "BETO", "BETW1",
            "BETW2", "BGIDL", "BGIDLD", "BGIDLDO", "BGIDLO", "CBBTBOT",
            "CBBTBOTD", "CBBTGAT", "CBBTGATD", "CBBTSTI", "CBBTSTID",
            "CF", "CFB", "CFBEDGEO", "CFBO", "CFD",
            "CFDO", "CFL", "CFLEXP", "CFR", "CFRD", "CFRDW", "CFRW",
            "CFW", "CFDEDGEO", "CFEDGEL", "CFEDGELEXP", "CFEDGEW", "CGBOV", "CGBOVL",
            "CGIDLDO", "CGIDLO", "CGOV", "CGOVD", "CHIBO", "CLM", "COX",
            "CJORBOT", "CJORBOTD", "CJORGAT", "CJORGATD", "CJORSTI",
            "CJORSTID", "CS", "CSL", "CSLEXP", "CSLW", "CSO", "CSRHBOT",
            "CSRHBOTD", "CSRHGAT", "CSRHGATD", "CSRHSTI", "CSRHSTID",
            "CSW", "CT", "CTATBOT", "CTATBOTD", "CTATGAT", "CTATGATD",
            "CTATSTI", "CTATSTID", "CTBO", "CTEDGEL", "CTEDGELEXP", "CTEDGEO", "CTGO", "CTL",
            "CTLEXP", "CTLW", "CTO", "CTW", "DELVTAC", "DELVTACL",
            "DELVTACLEXP", "DELVTACLW", "DELVTACO", "DELVTACW", "DELVTO",
            "DLSIL", "DLQ", "DNSUB", "DNSUBO", "DPHIB", "DPHIBEDGEL",
            "DPHIBEDGELEXP", "DPHIBEDGELW", "DPHIBEDGEO", "DPHIBEDGEW", "DPHIBL", "DPHIBLEXP",
            "DPHIBLW", "DPHIBO", "DPHIBW", "DTA", "DVSBNUD", "DVSBNUDO",
            "DWQ", "EF", "EF0", "EFEDGEO", "EFO", "EPSROX", "EPSROXO", "FACNEFFAC",
            "FACNEFFACL", "FACNEFFACLW", "FACNEFFACO", "FACNEFFACW",
            "FACTUO", "FBBTRBOT", "FBBTRBOTD", "FBBTRGAT", "FBBTRGATD",
            "FBBTRSTI", "FBBTRSTID", "FBET1", "FBET1W", "FBET2", "FBETEDGE", "FETA", "FETAO",
            "FJUNQ", "FJUNQD", "FNT", "FNTEDGEO", "FNTEXCL", "FNTO", "FREV",
            "FOL1", "FOL2", "GFACNUD", "GFACNUDL", "GFACNUDLEXP",
            "GFACNUDLW", "GFACNUDO", "GFACNUDW", "GC2O", "GC3O", "GCOO",
            "IGINV", "IGINVLW", "IGOV", "IGOVD", "IGOVDW", "IGOVW",
            "IDSATRBOT", "IDSATRBOTD", "IDSATRGAT", "IDSATRGATD", "IDSATRSTI", "IDSATRSTID",
            "IMAX", "KF", "KF0", "KP", "KUO", "KUOWEL", "KUOWELW", "KUOWEO",
            "KUOWEW", "KVSAT", "KVTHO", "KVTHOWEL", "KVTHOWELW", "KVTHOWEO",
            "KVTHOWEW", "L", "LAM", "LAMBDA", "LAMDA", "LAP", "LEVEL",
            "LINTNOI", "LKUO", "LKVTHO", "LLODKUO", "LLODVTH", "LODETAO",
            "LOV", "LOVD", "LP1", "LP1W", "LP2", "LPCK", "LPCKW", "LPEDGE", "LVARL",
            "LVARO", "LVARW", "MOBILITY", "MUE", "MUEO", "MUEW", "N", "N0",
            "MEFFTATBOT", "MEFFTATBOTD", "MEFFTATGAT", "MEFFTATGATD",
            "MEFFTATSTI", "MEFFTATSTID", "MUNQSO",
            "NFACTOR", "NEFF", "NOV", "NOVD", "NOVDO", "NOVO", "NP",
            "NPCK", "NPCKW", "NPL", "NPO", "NFAEDGELW", "NFALW", "NFBEDGELW",
            "NFBLW", "NFCEDGELW", "NFCLW", "NSLP", "NSLPO", "NSUBEDGEL",
            "NSUBEDGELEXP", "NSUBEDGELW", "NSUBEDGEO", "NSUBEDGEW", "NSUBO",
            "NSUBW", "PBOT", "PBOTD", "PBRBOT", "PBRBOTD", "PBRGAT", "PBRGATD",
            "PBRSTI", "PBRSTID", "PGAT", "PGATD", "PHIB", "PHIBO", "PHIGBOT",
            "PHIGBOTD", "PHIGGAT", "PHIGGATD", "PHIGSTI", "PHIGSTID", "PKUO",
            "PKVTHO", "PSCE", "PSCEB", "PSCEBEDGEO", "PSCEBO", "PSCED",
            "PSCEDEDGEO", "PSCEEDGEL", "PSCEEDGELEXP", "PSCEEDGEW",
            "PSCEDO", "PSCEL", "PSCELEXP", "PSCEW", "PSTI", "PSTID", "QMC", "RS", "RSB",
            "RSBO", "RSG", "RSGO", "RSH", "RSHD", "RSHG", "RBULKO", "RGO",
            "RJUNDO", "RJUNSO", "RINT", "RSW1", "RSW2", "RTH", "RVPOLY",
            "RWELLO", "SAREF", "SBREF", "SCREF", "SFL", "SFLN",
            "ST2VFBO", "STA2O", "STBET", "STBETEDGEL", "STBETEDGELW", "STBETEDGEO",
            "STBETEDGEW", "STBETL", "STBETLW", "STBETN", "STBETO",
            "STBETW", "STBGIDL", "STBGIDLD", "STBGIDLDO", "STBGIDLO",
            "STCS", "STCSO", "STCTO", "STETAO", "STFBBTBOT", "STFBBTBOTD",
            "STFBBTGAT", "STFBBTGATD", "STFBBTSTI", "STFBBTSTID", "STIG", "STIGO",
            "STMUE", "STMUEO", "STRS", "STRSO", "STRTH", "STTHECS",
            "STTHEMU", "STTHEMUO", "STTHESAT", "STTHESATL", "STTHESATLW",
            "STTHESATO", "STTHESATW", "STTHECSO", "STVFB", "STVFB0",
            "STVFBEDGEL", "STVFBEDGELW", "STVFBEDGEO", "STVFBEDGEW", "STVFBL",
            "STVFBLW", "STVFBO", "STVFBW", "STXCOR", "STXCORO", "SWDELVTAC",
            "SWEDGE", "SWGEO", "SWGIDL", "SWIGATE", "SWIGN", "SWIMPACT",
            "SWJUNCAP", "SWJUNEXP", "SWJUNASYM", "SWNQS", "SWNUD", "SWSOA",
            "TA", "THECS", "THECSO", "THEMU", "THEMUO",
            "THESAT", "THESATAC", "THESATB", "THESATBO", "THESATG",
            "THESATGO", "THESATL", "THESATLEXP", "THESATLW", "THESATO",
            "THESATW", "TKUO", "TNOM", "TNOM_C", "TOX", "TOXO", "TOXOV", "TOXOVD",
            "TOXOVDO", "TOXOVO", "TR", "TRJ", "TYPE", "U0", "UO", "VFB", "VFB0",
            "VFBEDGEO", "VFBL", "VFBLW", "VFBO", "VFBW", "VJUNREF", "VJUNREFD",
            "VDB_MAX", "VDS_MAX", "VGB_MAX", "VGD_MAX", "VGS_MAX",
            "VNSUB", "VNSUBO", "VP", "VBIRBOT", "VBIRBOTD", "VBIRGAT",
            "VBIRGATD", "VBIRSTI", "VBIRSTID", "VBRBOT", "VBRBOTD", "VBRGAT",
            "VBRGATD", "VBRSTI", "VBRSTID",
            "VPO", "VSB_MAX", "VSBNUD", "VSBNUDO", "VT0", "VTH", "VTH0", "VTO",
            "W", "WBET", "WEB", "WEC", "WEDGE", "WEDGEW", "WKUO", "WKVTHO",
            "WLOD", "WLODKUO", "WLODVTH", "WOT", "WSEG", "WSEGP", "WVARL",
            "WVARO", "WVARW", "X", "XCOR", "XCORL", "XCORLW", "XCORO", "XCORW",
            "XJUNGAT", "XJUNGATD", "XJUNSTI", "XJUNSTID", "Y"
        };
        return handled.find(normalize(name)) != handled.end();
    }

private:
    static std::string normalize(std::string value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        return value;
    }

    static double instanceValue(
        const GsdiParamMap& instance,
        std::initializer_list<const char*> names,
        double fallback) {
        for (const auto& [key, value] : instance) {
            for (const char* name : names) {
                if (normalize(key) == normalize(name)) return value;
            }
        }
        return fallback;
    }

    GsdiParamMap values_;
};

} // namespace gspice

#endif // GSPICE_PSP103_PARAMETERS_HPP
