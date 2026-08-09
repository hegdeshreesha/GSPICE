#pragma once

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace gspice {

struct Bsim4Validation {
    bool valid = false;
    std::string reason;
    explicit operator bool() const { return valid; }
};

struct Bsim4PreparedModel {
    Bsim4Validation validation;
    double level = 54.0;
    double temperature_k = 300.15;
    double nominal_temperature_k = 300.15;
    double width = 1.0e-6;
    double length = 1.0e-6;
    double leff = 1.0e-6;
    double weff = 1.0e-6;
    double toxe = 1.0e-8;
    double cox = 3.453133e-3;
    double toxm = 1.0e-8;
    double toxp = 1.0e-8;
    double coxp = 3.453133e-3;
    double vth0 = 0.4;
    double phi = 0.7;
    double sqrtPhi = std::sqrt(0.7);
    double vbi = 0.0;
    double vtfbphi2 = 0.0;
    double vfb = -0.3;
    double ndep = 1.0e17;
    double minv = 0.0;
    double voff = -0.08;
    double cdep0 = 1.0e-3;
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
    double kf = 0.0;
    double af = 1.0;
    double ef = 1.0;
    double kp = 120.0e-6;
    double u0 = 0.05;
    double vsat = 1.0e5;
    double ute = -1.5;
    double at = 3.3e4;
    double k1 = 0.5;
    double k2 = 0.0;
    double k1ox = 0.0;
    double k2ox = 0.0;
    double k3 = 80.0;
    double k3b = 0.0;
    double w0 = 2.5e-6;
    double lpe0 = 1.74e-7;
    double lpeb = 0.0;
    double phin = 0.0;
    double kt1 = -0.11;
    double kt1l = 0.0;
    double kt2 = 0.022;
    double cdsc = 2.4e-4;
    double cdscb = 0.0;
    double cdscd = 0.0;
    double cit = 0.0;
    double nsub = 6.0e16;
    double nsd = 1.0e20;
    double ados = 1.0;
    double bdos = 1.0;
    int mobmod = 0;
    double a0 = 1.0;
    double ags = 0.0;
    double b0 = 0.0;
    double b1 = 0.0;
    double a1 = 0.0;
    double a2 = 1.0;
    double keta = 0.0;
    double xj = 1.5e-7;
    double xdep0 = 1.0e-7;
    double nfactor = 1.0;
    double ua = 0.0;
    double ub = 0.0;
    double uc = 0.0;
    double ua1 = 1.0e-9;
    double ub1 = -1.0e-18;
    double uc1 = -0.056e-9;
    double dvt0 = 0.0;
    double dvt1 = 0.0;
    double dvt2 = 0.0;
    double dvt0w = 0.0;
    double dvt1w = 0.0;
    double dvt2w = 0.0;
    double pclm = 1.3;
    double pdibl1 = 0.39;
    double pdibl2 = 0.0086;
    double dsub = 0.56;
    double drout = 0.56;
    double thetaRout = 0.0086;
    double pvag = 0.0;
    double pdiblb = 0.0;
    double fprout = 0.0;
    double lambda = 0.0;
    int rdsmod = 0;
    double litl = 0.0;
    double eta0 = 0.08;
    double etab = 0.0;
    double rdsw = 200.0;
    double rdswmin = 0.0;
    double rds0 = 200.0;
    double prwg = 1.0;
    double prwb = 0.0;
    double wr = 1.0;
    double prt = 0.0;
    double delta = 0.01;
    double cgso = 0.0;
    double cgdo = 0.0;
    double cgbo = 0.0;
    double cj = 0.0;
    double cjsw = 0.0;
    double pb = 0.8;
    double pbsw = 0.8;
    double mj = 0.5;
    double mjsw = 0.33;
    double tcj = 0.0;
    double tpb = 0.0;
    double tcjsw = 0.0;
    double tpbsw = 0.0;
    double xpart = 0.5;
    double beta = 120.0e-6;
};

class Bsim4ParameterSet {
public:
    static Bsim4ParameterSet from(const std::unordered_map<std::string, double>& values) {
        Bsim4ParameterSet result;
        for (const auto& [name, value] : values) result.values_[normalize(name)] = value;
        return result;
    }

    Bsim4PreparedModel prepare(double width, double length, double temperature_c = 27.0) const {
        Bsim4PreparedModel out;
        out.width = width;
        out.length = length;
        out.temperature_k = temperature_c + 273.15;
        out.nominal_temperature_k = get({"TNOM"}, 27.0) + 273.15;
        out.level = get({"LEVEL"}, 54.0);
        out.vth0 = get({"VTH0", "VT0", "VTH"}, 0.4);
        out.ndep = get({"NDEP"}, 1.7e17);
        out.phin = get({"PHIN"}, 0.0);
        {
            // BSIM4.8.3 phi from NDEP and intrinsic density (b4temp.c:1323).
            const double phiGiven = get({"PHI"}, std::numeric_limits<double>::quiet_NaN());
            const double tnom = out.nominal_temperature_k;
            const double vtm0 = (1.380649e-23 / 1.602176634e-19) * tnom;
            const double eg0 = 1.16 - 7.02e-4 * tnom * tnom / (tnom + 1108.0);
            const double ni = 1.45e10 * (tnom / 300.15) * std::sqrt(tnom / 300.15) *
                std::exp(std::clamp(21.5565981 - eg0 / (2.0 * vtm0), -700.0, 700.0));
            out.phi = std::isfinite(phiGiven)
                ? phiGiven
                : vtm0 * std::log(std::max(out.ndep, 1.0e10) / ni) + out.phin + 0.4;
            out.sqrtPhi = std::sqrt(std::max(out.phi, 1.0e-6));
            out.nsub = get({"NSUB"}, 6.0e16);
            out.nsd = get({"NSD"}, 1.0e20);
            out.vbi = vtm0 * std::log(std::max(out.nsd, 1.0e10) *
                std::max(out.ndep, 1.0e10) / (ni * ni));
        }
        out.minv = get({"MINV"}, 0.0);
        out.voff = get({"VOFF"}, -0.08);
        out.kp = get({"KP", "BETA"}, 0.0);
        out.u0 = get({"U0", "MOBMOD"}, 0.05);
        if (out.u0 > 1.0) out.u0 *= 1.0e-4;
        out.vsat = get({"VSAT", "VSATT"}, 1.0e5);
        out.ute = get({"UTE"}, -1.5);
        out.at = get({"AT"}, 3.3e4);
        out.k1 = get({"K1"}, 0.5);
        out.k2 = get({"K2"}, 0.0);
        out.k3 = get({"K3"}, 80.0);
        out.k3b = get({"K3B"}, 0.0);
        out.w0 = get({"W0"}, 2.5e-6);
        out.lpe0 = get({"LPE0"}, 1.74e-7);
        out.lpeb = get({"LPEB"}, 0.0);
        out.kt1 = get({"KT1"}, -0.11);
        out.kt1l = get({"KT1L"}, 0.0);
        out.kt2 = get({"KT2"}, 0.022);
        out.cdsc = get({"CDSC"}, 2.4e-4);
        out.cdscb = get({"CDSCB"}, 0.0);
        out.cdscd = get({"CDSCD"}, 0.0);
        out.cit = get({"CIT"}, 0.0);
        out.ados = get({"ADOS"}, 1.0);
        out.bdos = get({"BDOS"}, 1.0);
        out.mobmod = static_cast<int>(get({"MOBMOD"}, 0.0));
        out.a0 = get({"A0"}, 1.0);
        out.ags = get({"AGS"}, 0.0);
        out.b0 = get({"B0"}, 0.0);
        out.b1 = get({"B1"}, 0.0);
        out.a1 = get({"A1"}, 0.0);
        out.a2 = get({"A2"}, 1.0);
        out.keta = get({"KETA"}, 0.0);
        out.xj = get({"XJ"}, 1.5e-7);
        out.xdep0 = get({"XDEP0"}, 1.0e-7);
        out.nfactor = get({"NFACTOR", "N"}, 1.0);
        out.ua = get({"UA"}, 1.0e-9);
        out.ub = get({"UB"}, 1.0e-19);
        out.uc = get({"UC"}, -0.0465e-9);
        out.ua1 = get({"UA1"}, 1.0e-9);
        out.ub1 = get({"UB1"}, -1.0e-18);
        out.uc1 = get({"UC1"}, -0.056e-9);
        out.dvt0 = get({"DVT0"}, 2.2);
        out.dvt1 = get({"DVT1"}, 0.53);
        out.dvt2 = get({"DVT2"}, -0.032);
        out.dvt0w = get({"DVT0W"}, 0.0);
        out.dvt1w = get({"DVT1W"}, 5.3e6);
        out.dvt2w = get({"DVT2W"}, 0.0);
        out.pclm = get({"PCLM"}, 1.3);
        out.pdibl1 = get({"PDIBL1"}, 0.39);
        out.pdibl2 = get({"PDIBL2"}, 0.0086);
        out.dsub = get({"DSUB"}, 0.56);
        out.drout = get({"DROUT"}, 0.56);
        out.pvag = get({"PVAG"}, 0.0);
        out.pdiblb = get({"PDIBLB"}, 0.0);
        out.fprout = get({"FPROUT"}, 0.0);
        out.lambda = get({"LAMBDA"}, 0.0);
        out.rdsmod = static_cast<int>(get({"RDSMOD"}, 0.0));
        out.eta0 = get({"ETA0"}, 0.08);
        out.etab = get({"ETAB"}, -0.07);
        out.rdsw = get({"RDSW"}, 200.0);
        out.rdswmin = get({"RDSWMIN"}, 0.0);
        out.prwg = get({"PRWG"}, 1.0);
        out.prwb = get({"PRWB"}, 0.0);
        out.wr = get({"WR"}, 1.0);
        out.prt = get({"PRT"}, 0.0);
        out.delta = get({"DELTA"}, 0.01);
        out.cgso = get({"CGSO"}, 0.0);
        out.cgdo = get({"CGDO"}, 0.0);
        out.cgbo = get({"CGBO"}, 0.0);
        out.cj = get({"CJ"}, 0.0);
        out.cjsw = get({"CJSW"}, 0.0);
        out.pb = get({"PB"}, 0.8);
        out.pbsw = get({"PBSW"}, 0.8);
        out.mj = get({"MJ"}, 0.5);
        out.mjsw = get({"MJSW"}, 0.33);
        out.tcj = get({"TCJ"}, 0.0);
        out.tpb = get({"TPB"}, 0.0);
        out.tcjsw = get({"TCJSW"}, 0.0);
        out.tpbsw = get({"TPBSW"}, 0.0);
        out.xpart = get({"XPART"}, 0.5);
        const double dl = get({"DL"}, 0.0);
        const double dw = get({"DW"}, 0.0);
        out.toxe = get({"TOXE", "TOX"}, 1.0e-8);
        out.toxm = get({"TOXM"}, out.toxe);
        out.toxp = get({"TOXP"}, out.toxe);
        out.coxp = 3.453133e-11 / std::max(out.toxp, 1.0e-12);
        out.leff = length - 2.0 * dl;
        out.weff = width - 2.0 * dw;
        out.xj = getBinned("XJ", out.xj, out.leff, out.weff);
        out.xdep0 = getBinned("XDEP0", out.xdep0, out.leff, out.weff);
        out.phi = getBinned("PHI", out.phi, out.leff, out.weff);
        out.ndep = getBinned("NDEP", out.ndep, out.leff, out.weff);
        out.voff = getBinned("VOFF", out.voff, out.leff, out.weff);
        out.k1 = getBinned("K1", out.k1, out.leff, out.weff);
        out.k2 = getBinned("K2", out.k2, out.leff, out.weff);
        out.a0 = getBinned("A0", out.a0, out.leff, out.weff);
        out.ags = getBinned("AGS", out.ags, out.leff, out.weff);
        out.b0 = getBinned("B0", out.b0, out.leff, out.weff);
        out.b1 = getBinned("B1", out.b1, out.leff, out.weff);
        out.a1 = getBinned("A1", out.a1, out.leff, out.weff);
        out.a2 = getBinned("A2", out.a2, out.leff, out.weff);
        out.keta = getBinned("KETA", out.keta, out.leff, out.weff);
        out.nfactor = getBinned("NFACTOR", out.nfactor, out.leff, out.weff);
        out.vth0 = getBinned("VTH0", out.vth0, out.leff, out.weff);
        out.vsat = getBinned("VSAT", out.vsat, out.leff, out.weff);
        out.ua = getBinned("UA", out.ua, out.leff, out.weff);
        out.ub = getBinned("UB", out.ub, out.leff, out.weff);
        out.uc = getBinned("UC", out.uc, out.leff, out.weff);
        out.dvt0 = getBinned("DVT0", out.dvt0, out.leff, out.weff);
        out.dvt1 = getBinned("DVT1", out.dvt1, out.leff, out.weff);
        out.dvt2 = getBinned("DVT2", out.dvt2, out.leff, out.weff);
        out.dvt0w = getBinned("DVT0W", out.dvt0w, out.leff, out.weff);
        out.dvt1w = getBinned("DVT1W", out.dvt1w, out.leff, out.weff);
        out.dvt2w = getBinned("DVT2W", out.dvt2w, out.leff, out.weff);
        out.pclm = getBinned("PCLM", out.pclm, out.leff, out.weff);
        out.pdibl1 = getBinned("PDIBL1", out.pdibl1, out.leff, out.weff);
        out.pdibl2 = getBinned("PDIBL2", out.pdibl2, out.leff, out.weff);
        out.dsub = getBinned("DSUB", out.dsub, out.leff, out.weff);
        out.drout = getBinned("DROUT", out.drout, out.leff, out.weff);
        out.pvag = getBinned("PVAG", out.pvag, out.leff, out.weff);
        out.pdiblb = getBinned("PDIBLB", out.pdiblb, out.leff, out.weff);
        out.fprout = getBinned("FPROUT", out.fprout, out.leff, out.weff);
        out.lambda = getBinned("LAMBDA", out.lambda, out.leff, out.weff);
        out.eta0 = getBinned("ETA0", out.eta0, out.leff, out.weff);
        out.etab = getBinned("ETAB", out.etab, out.leff, out.weff);
        out.rdsw = getBinned("RDSW", out.rdsw, out.leff, out.weff);
        out.rdswmin = getBinned("RDSWMIN", out.rdswmin, out.leff, out.weff);
        out.prwg = getBinned("PRWG", out.prwg, out.leff, out.weff);
        out.prwb = getBinned("PRWB", out.prwb, out.leff, out.weff);
        out.wr = getBinned("WR", out.wr, out.leff, out.weff);
        out.prt = getBinned("PRT", out.prt, out.leff, out.weff);
        out.delta = getBinned("DELTA", out.delta, out.leff, out.weff);
        out.cgso = getBinned("CGSO", out.cgso, out.leff, out.weff);
        out.cgdo = getBinned("CGDO", out.cgdo, out.leff, out.weff);
        out.cgbo = getBinned("CGBO", out.cgbo, out.leff, out.weff);
        out.cj = getBinned("CJ", out.cj, out.leff, out.weff);
        out.cjsw = getBinned("CJSW", out.cjsw, out.leff, out.weff);
        out.pb = getBinned("PB", out.pb, out.leff, out.weff);
        out.pbsw = getBinned("PBSW", out.pbsw, out.leff, out.weff);
        out.mj = getBinned("MJ", out.mj, out.leff, out.weff);
        out.mjsw = getBinned("MJSW", out.mjsw, out.leff, out.weff);
        out.cox = 3.453133e-11 / std::max(out.toxe, 1.0e-12);
        const double delta_temperature = out.temperature_k - out.nominal_temperature_k;
        out.ua *= 1.0 + out.ua1 * delta_temperature;
        out.ub *= 1.0 + out.ub1 * delta_temperature;
        out.uc *= 1.0 + out.uc1 * delta_temperature;
        out.u0 *= std::pow(std::max(out.temperature_k / out.nominal_temperature_k, 1.0e-12), out.ute);
        out.vsat *= std::max(0.0, 1.0 - out.at * delta_temperature);
        out.cj *= std::max(0.0, 1.0 + out.tcj * delta_temperature);
        out.cjsw *= std::max(0.0, 1.0 + out.tcjsw * delta_temperature);
        out.pb = std::max(1.0e-6, out.pb + out.tpb * delta_temperature);
        out.pbsw = std::max(1.0e-6, out.pbsw + out.tpbsw * delta_temperature);
        // BSIM4 does not require the legacy Level-1 KP card.  When it is
        // absent, the low-field coefficient is mobility times oxide
        // capacitance; retain an explicit KP/BETA override when supplied.
        if (out.kp == 0.0) out.kp = out.u0 * out.cox;
        out.cdep0 = std::sqrt(1.602176634e-19 * 1.03594e-10 *
                              std::max(out.ndep, 1.0e10) * 1.0e6 /
                              (2.0 * std::max(out.phi, 1.0e-6)));
        // BSIM4.8.3 depletion depth at zero body bias (b4temp.c:1341).
        out.xdep0 = std::sqrt(2.0 * 1.03594e-10 /
                              (1.602176634e-19 * std::max(out.ndep, 1.0e10) * 1.0e6)) *
                    out.sqrtPhi;
        // BSIM4.8.3 body factors scaled by toxm (b4temp.c:1516,1802).
        out.k1ox = out.k1 * out.toxe / std::max(out.toxm, 1.0e-12);
        out.k2ox = out.k2 * out.toxe / std::max(out.toxm, 1.0e-12);
        // vtfbphi2 for the quantum-affected oxide capacitance (b4temp.c:1786).
        out.vtfbphi2 = std::max(0.0, 4.0 * out.k1 * out.sqrtPhi);
        out.vfb = out.vth0 - out.phi - out.k1 * out.sqrtPhi;
        out.litl = std::sqrt(std::max(3.0 * out.xj * out.toxe, 1.0e-30));
        const double epsRatioToxe = (1.03594e-10 / (3.9 * 8.8541878128e-12)) *
                                    std::max(out.toxe, 1.0e-12);
        const double routExponent = out.drout * out.leff /
            std::sqrt(std::max(epsRatioToxe * out.xdep0, 1.0e-30));
        const double routExp = std::exp(std::clamp(routExponent, -700.0, 80.0));
        const double routDelta = routExp - 1.0;
        out.thetaRout = out.pdibl1 * routExp /
            (routDelta * routDelta + 2.0 * routExp * 1.0e-300) + out.pdibl2;
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
        out.kf = get({"KF"}, 0.0);
        out.af = get({"AF"}, 1.0);
        out.ef = get({"EF"}, 1.0);
        out.is = getBinned("IS", out.is, out.leff, out.weff);
        out.js = getBinned("JS", out.js, out.leff, out.weff);
        out.jsw = getBinned("JSW", out.jsw, out.leff, out.weff);
        // Junction leakage scales with temperature after the geometry binned
        // values are resolved (XTI/EG are read above).
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
        out.kf = getBinned("KF", out.kf, out.leff, out.weff);
        out.af = getBinned("AF", out.af, out.leff, out.weff);
        out.ef = getBinned("EF", out.ef, out.leff, out.weff);
        out.beta = out.kp * out.weff / std::max(out.leff, 1.0e-15);
        const double width_um = std::max(out.weff * 1.0e6, 1.0e-12);
        const double rdsTemperature = std::max(
            0.0, 1.0 + out.prt * (out.temperature_k - out.nominal_temperature_k));
        const double widthScale = std::pow(width_um, std::max(out.wr, 1.0e-12));
        out.rds0 = out.rdsw * rdsTemperature / widthScale;
        out.rdswmin = out.rdswmin * rdsTemperature / widthScale;

        if (!std::isfinite(out.level) || out.level != 54.0) out.validation.reason = "only BSIM4 level 54 is supported";
        else if (!std::isfinite(out.temperature_k) || out.temperature_k <= 0.0 ||
                 !std::isfinite(out.nominal_temperature_k) || out.nominal_temperature_k <= 0.0)
            out.validation.reason = "temperature and TNOM must be above absolute zero";
        else if (!std::isfinite(width) || width <= 0.0 || !std::isfinite(length) || length <= 0.0)
            out.validation.reason = "W and L must be finite and positive";
        else if (!(out.leff > 0.0) || !(out.weff > 0.0)) out.validation.reason = "effective W/L became non-positive";
        else if (!std::isfinite(out.kp) || out.kp < 0.0 || !std::isfinite(out.u0) || out.u0 < 0.0 ||
                 !std::isfinite(out.vsat) || out.vsat <= 0.0 || !std::isfinite(out.toxe) || out.toxe <= 0.0)
            out.validation.reason = "KP/U0 must be nonnegative and VSAT/TOXE positive";
        else if (!std::isfinite(out.nfactor) || out.nfactor <= 0.0 || !std::isfinite(out.xpart) ||
                 out.xpart < 0.0 || out.xpart > 1.0)
            out.validation.reason = "NFACTOR must be positive and XPART must be in [0,1]";
        else if (!std::isfinite(out.a2) || out.a2 <= 0.0 || !std::isfinite(out.fprout) ||
                 out.fprout < 0.0 || !std::isfinite(out.lambda) || out.lambda < 0.0)
            out.validation.reason = "A2 must be positive and FPROUT/LAMBDA nonnegative";
        else if (!std::isfinite(out.is) || out.is < 0.0 || !std::isfinite(out.js) ||
                 out.js < 0.0 || !std::isfinite(out.jsw) || out.jsw < 0.0)
            out.validation.reason = "IS/JS/JSW must be finite and nonnegative";
        else if (!std::isfinite(out.xti) || !std::isfinite(out.eg) || out.eg < 0.0)
            out.validation.reason = "XTI/EG junction temperature parameters are invalid";
        else if (!std::isfinite(out.agidl) || out.agidl < 0.0 || !std::isfinite(out.bgidl) ||
                 out.bgidl < 0.0 || !std::isfinite(out.cgidl) || out.cgidl < 0.0 ||
                 !std::isfinite(out.egidl) || out.egidl < 0.0 || !std::isfinite(out.agisl) ||
                 out.agisl < 0.0 || !std::isfinite(out.bgisl) || out.bgisl < 0.0 ||
                 !std::isfinite(out.cgisl) || out.cgisl < 0.0 || !std::isfinite(out.egisl) ||
                 out.egisl < 0.0)
            out.validation.reason = "GIDL/GISL parameters must be finite and nonnegative";
        else if (!std::isfinite(out.cj) || out.cj < 0.0 || !std::isfinite(out.cjsw) ||
                 out.cjsw < 0.0 || !std::isfinite(out.pb) || out.pb <= 0.0 ||
                 !std::isfinite(out.pbsw) || out.pbsw <= 0.0 || !std::isfinite(out.mj) ||
                 out.mj < 0.0 || out.mj >= 1.0 || !std::isfinite(out.mjsw) ||
                 out.mjsw < 0.0 || out.mjsw >= 1.0)
            out.validation.reason = "junction capacitance parameters are invalid";
        else if (!std::isfinite(out.tcj) || !std::isfinite(out.tpb) ||
                 !std::isfinite(out.tcjsw) || !std::isfinite(out.tpbsw))
            out.validation.reason = "junction capacitance temperature coefficients are invalid";
        else if (out.rdsmod != 0 && out.rdsmod != 1)
            out.validation.reason = "RDSMOD must be 0 or 1";
        else if (!std::isfinite(out.wr) || out.wr <= 0.0 || !std::isfinite(out.prt) ||
                 !std::isfinite(out.ute) || !std::isfinite(out.at) || out.vsat <= 0.0)
            out.validation.reason = "WR/UTE/AT/PRT temperature parameters are invalid";
        else out.validation.valid = true;
        return out;
    }

private:
    double getBinned(const char* name, double fallback, double leff, double weff) const {
        const std::string base(name);
        const auto value = [&](const std::string& key, double default_value) {
            const auto it = values_.find(key);
            return it == values_.end() ? default_value : it->second;
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

    std::unordered_map<std::string, double> values_;
};

// Allow-list of model-card parameter names consumed by the native BSIM4
// evaluator. Kept in sync by tools/audit_bsim4_parameters.py; the audit tool
// fails when a native read-set name is missing from this list.
struct Bsim4ImplementedParameters {
    static bool isSupported(std::string_view name) {
        std::string upper(name);
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char c) {
                           return static_cast<char>(std::toupper(c));
                       });
        return supported().find(upper) != supported().end();
    }

    static const std::unordered_set<std::string>& supported() {
        static const std::unordered_set<std::string> supported = {
            "A0", "A1", "A2", "ADOS", "AF", "AGIDL", "AGISL", "AGS",
            "AIGC", "AT", "B0", "B1", "BDOS", "BETA", "BGIDL", "BGISL",
            "BIGC", "CDSC", "CDSCB", "CDSCD", "CGBO", "CGDO", "CGIDL",
            "CGISL", "CGSO", "CIGC", "CIT", "CJ", "CJSW", "DELTA", "DL",
            "DROUT", "DSUB", "DVT0", "DVT0W", "DVT1", "DVT1W", "DVT2",
            "DVT2W", "DW", "EF", "EG", "EGIDL", "EGISL", "ETA0", "ETAB",
            "FPROUT", "IS", "JS", "JSW", "K1", "K2", "K3", "K3B", "KETA",
            "KF", "KP", "KT1", "KT1L", "KT2", "LAMBDA", "LEVEL",
            "LPE0", "LPEB",
            "MINV", "MJ", "MJSW", "MOBMOD", "N", "NDEP", "NFACTOR", "NJ",
            "NSD", "NSUB", "PB", "PBSW", "PCLM", "PDIBL1", "PDIBL2",
            "PDIBLB", "PHI", "PHIN", "PRT", "PRWB", "PRWG", "PVAG", "RDSMOD",
            "RDSW", "RDSWMIN", "TCJ", "TCJSW", "TNOM", "TOX", "TOXE", "TOXM",
            "TOXP", "TPB", "TPBSW", "U0", "UA", "UA1", "UB", "UB1", "UC",
            "UC1", "UTE", "VOFF", "VSAT", "VSATT", "VT0", "VTH", "VTH0",
            "W0", "WR", "XDEP0", "XJ", "XPART", "XTI",
        };
        return supported;
    }
};

} // namespace gspice
