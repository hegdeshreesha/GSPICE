#include "devices/juncap_express.hpp"

#include <cassert>
#include <cmath>

int main() {
    assert(gspice::validateJuncapExpress(1).supported);
    assert(!gspice::validateJuncapExpress(0).supported);

    gspice::JuncapExpressInputs in;
    in.vak = 0.6;
    in.vj = 0.5;
    in.phi_td = 0.026;
    in.isat_for1 = 1e-14;
    in.m_for1 = 1.0;
    in.bottom.area = 1.0;
    in.bottom.cjo = 1e-12;
    in.bottom.vbi = 0.8;
    in.bottom.p = 0.5;
    const auto result = gspice::juncapExpressEvaluate(in);
    assert(result.current > 0.0);
    assert(result.dcurrent_dvak > 0.0);
    assert(result.charge > 0.0);
    assert(result.capacitance > 0.0);
    assert(std::abs(result.noise_psd - 2.0 * 1.602176634e-19 * result.current) < 1e-30);

    gspice::GsdiEvalResult gsdi;
    gspice::juncapExpressAppendChargeBranch(gsdi, 1, 2, 1, 2, result.charge, result.capacitance, 3);
    assert(gsdi.dynamic_residual.size() == 2);
    assert(gsdi.dynamic_jacobian.size() == 4);
    assert(gsdi.finite());
    return 0;
}
