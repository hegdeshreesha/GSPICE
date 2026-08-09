#include "devices/juncap2_core.hpp"
#include <cassert>
#include <cmath>

int main() {
    gspice::Juncap2Inputs in;
    in.idsat = 1.0e-14;
    auto r = gspice::juncap2Evaluate(in);
    assert(r.supported && std::abs(r.current) < 1.0e-20 && r.dcurrent_dvak > 0.0);
    in.vak = 0.15;
    r = gspice::juncap2Evaluate(in);
    assert(r.supported && r.current > 0.0 && r.dcurrent_dvak > 0.0);
    in.csrh = 1.0e-12;
    r = gspice::juncap2Evaluate(in);
    assert(r.supported && std::isfinite(r.current) && std::isfinite(r.dcurrent_dvak));
    in.cbbt = 1.0e-18;
    r = gspice::juncap2Evaluate(in);
    assert(r.supported && std::isfinite(r.current) && std::isfinite(r.dcurrent_dvak));
    in.ctat = 1.0e-18;
    in.vak = -0.5;
    r = gspice::juncap2Evaluate(in);
    assert(r.supported && std::isfinite(r.current) && std::isfinite(r.dcurrent_dvak));
    assert(r.current < 0.0);
    in.vbr = 0.2;
    r = gspice::juncap2Evaluate(in);
    assert(r.supported && std::isfinite(r.current) && std::isfinite(r.dcurrent_dvak));
    assert(r.current < 0.0);
    in.ctat = 0.0;
    in.vbr = 1.0e9;
    in.p = 1.0;
    r = gspice::juncap2Evaluate(in);
    assert(!r.supported && r.reason == "P must be in [0,1)");
}
