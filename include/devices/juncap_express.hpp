#ifndef GSPICE_JUNCAP_EXPRESS_HPP
#define GSPICE_JUNCAP_EXPRESS_HPP

#include "../gsdi.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace gspice {

struct JuncapExpressValidation {
    bool supported = false;
    std::string reason;
};

inline JuncapExpressValidation validateJuncapExpress(int swjunexp) {
    if (swjunexp == 1) return {true, "JUNCAP Express enabled"};
    return {false, "full JUNCAP2 is not implemented; refusing Express fallback"};
}

struct JuncapExpressComponent {
    double area = 0.0;
    double side = 0.0;
    double edge = 0.0;
    double isat = 0.0;
    double m = 1.0;
    double cjo = 0.0;
    double vbi = 1.0;
    double p = 0.5;
    double charge_linear = 0.0;
};

struct JuncapExpressInputs {
    double vak = 0.0;
    double vj = 0.0;
    double dvj_dvak = 1.0;
    double phi_td = 0.02585;
    double isat_for1 = 0.0;
    double m_for1 = 1.0;
    double isat_for2 = 0.0;
    double m_for2 = 1.0;
    double isat_rev = 0.0;
    double m_rev = 1.0;
    double type = 1.0;
    double mult = 1.0;
    double fjunq = 0.0;
    JuncapExpressComponent bottom;
    JuncapExpressComponent sti;
    JuncapExpressComponent gate;
};

struct JuncapExpressResult {
    double current = 0.0;
    double dcurrent_dvak = 0.0;
    double charge = 0.0;
    double capacitance = 0.0;
    double noise_psd = 0.0;
};

inline double juncapExpressG(double vak, double isat, double m, double phi_td) {
    if (isat == 0.0) return 0.0;
    const double exponent = std::clamp(vak * m / phi_td, -80.0, 80.0);
    return isat * std::expm1(exponent);
}

inline double juncapExpressGDerivative(double vak, double isat, double m, double phi_td) {
    if (isat == 0.0) return 0.0;
    const double exponent = std::clamp(vak * m / phi_td, -80.0, 80.0);
    return isat * m / phi_td * std::exp(exponent);
}

inline double juncapExpressComponentCharge(
    const JuncapExpressComponent& c, double vak, double vj, double dvj_dvak,
    double fjunq, double ztot, double& dqdv) {
    const double z = c.area * c.cjo + c.side * c.cjo + c.edge * c.cjo;
    if (z <= fjunq * ztot || c.cjo == 0.0) {
        dqdv = 0.0;
        return 0.0;
    }
    const double x = std::max(1.0 - vj / c.vbi, 1e-12);
    const double p = c.p;
    const double q = std::abs(1.0 - p) < 1e-12
        ? -c.cjo * c.vbi * std::log(x)
        : c.cjo * c.vbi / (1.0 - p) * (1.0 - std::pow(x, 1.0 - p));
    dqdv = c.cjo * std::pow(x, -p) * dvj_dvak + c.charge_linear * (1.0 - dvj_dvak);
    return q + c.charge_linear * (vak - vj);
}

inline JuncapExpressResult juncapExpressEvaluate(const JuncapExpressInputs& in) {
    const double forward1 = juncapExpressG(in.vak, in.isat_for1, in.m_for1, in.phi_td);
    const double forward2 = juncapExpressG(in.vak, in.isat_for2, in.m_for2, in.phi_td);
    const double reverse = -juncapExpressG(-in.vak, in.isat_rev, in.m_rev, in.phi_td);
    const double dcurrent = juncapExpressGDerivative(in.vak, in.isat_for1, in.m_for1, in.phi_td) +
        juncapExpressGDerivative(in.vak, in.isat_for2, in.m_for2, in.phi_td) +
        juncapExpressGDerivative(-in.vak, in.isat_rev, in.m_rev, in.phi_td);
    const double ztot = in.bottom.area * in.bottom.cjo + in.sti.side * in.sti.cjo +
        in.gate.edge * in.gate.cjo;
    double cbot = 0.0, csti = 0.0, cgate = 0.0;
    const double qbot = juncapExpressComponentCharge(in.bottom, in.vak, in.vj, in.dvj_dvak, in.fjunq, ztot, cbot);
    const double qsti = juncapExpressComponentCharge(in.sti, in.vak, in.vj, in.dvj_dvak, in.fjunq, ztot, csti);
    const double qgate = juncapExpressComponentCharge(in.gate, in.vak, in.vj, in.dvj_dvak, in.fjunq, ztot, cgate);
    return {
        in.type * in.mult * (forward1 + forward2 + reverse),
        in.type * in.mult * dcurrent,
        in.type * in.mult * (qbot + qsti + qgate),
        in.type * in.mult * (cbot + csti + cgate),
        2.0 * 1.602176634e-19 * std::abs(in.type * in.mult * (forward1 + forward2 + reverse))};
}

inline void juncapExpressAppendChargeBranch(
    GsdiEvalResult& result, int positive_equation, int negative_equation,
    int positive_unknown, int negative_unknown, double charge,
    double capacitance, int conservation_group = -1) {
    result.dynamic_residual.push_back({positive_equation, charge, conservation_group});
    result.dynamic_residual.push_back({negative_equation, -charge, conservation_group});
    result.dynamic_jacobian.push_back({positive_equation, positive_unknown, capacitance, conservation_group});
    result.dynamic_jacobian.push_back({positive_equation, negative_unknown, -capacitance, conservation_group});
    result.dynamic_jacobian.push_back({negative_equation, positive_unknown, -capacitance, conservation_group});
    result.dynamic_jacobian.push_back({negative_equation, negative_unknown, capacitance, conservation_group});
}

} // namespace gspice

#endif
