#include "devices/bsim4_parameters.hpp"

#include <cassert>
#include <cmath>

int main() {
    const auto prepared = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"U0", 500.0}, {"TOXE", 1.0e-8},
        {"DL", 1.0e-8}, {"DW", 2.0e-8}, {"XPART", 0.5}
    }).prepare(1.0e-6, 1.0e-6, 27.0);
    assert(prepared.validation);
    assert(prepared.u0 == 0.05);
    assert(prepared.kp > 1.7e-4 && prepared.kp < 1.8e-4);
    assert(prepared.delta == 0.01);
    assert(prepared.dvt0 == 2.2);
    assert(prepared.dvt1 == 0.53);
    assert(prepared.dvt2 == -0.032);
    assert(prepared.dvt1w == 5.3e6);
    assert(prepared.dsub == 0.56);
    assert(prepared.drout == 0.56);
    assert(prepared.pvag == 0.0);
    assert(prepared.pdiblb == 0.0);
    assert(prepared.a1 == 0.0);
    assert(prepared.a2 == 1.0);
    assert(prepared.fprout == 0.0);
    assert(prepared.lambda == 0.0);
    assert(prepared.rdsmod == 0);
    // rds0 = RDS / (Weff/Wnom)^WR: with DW = 2e-8 the effective width is
    // 0.96 um and WR = 1.0, so RDS scales by 1/0.96.
    assert(std::abs(prepared.rds0 - 200.0 / 0.96) < 1.0e-9);
    assert(prepared.prwg == 1.0);
    assert(prepared.wr == 1.0);
    assert(prepared.prt == 0.0);
    assert(prepared.ute == -1.5);
    assert(prepared.at == 3.3e4);
    assert(prepared.ua1 == 1.0e-9);
    assert(prepared.ub1 == -1.0e-18);
    assert(prepared.uc1 == -0.056e-9);
    assert(prepared.litl > 0.0);
    assert(prepared.thetaRout > 0.0);
    assert(prepared.leff > 0.0 && prepared.weff > 0.0);
    assert(prepared.beta > 0.0 && prepared.cox > 0.0);

    const auto binned = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"VTH0", 0.4}, {"LVTH0", 1.0e-7},
        {"WVTH0", 2.0e-7}, {"PVTH0", 3.0e-13}
    }).prepare(1.0e-6, 1.0e-6);
    assert(binned.validation);
    assert(std::abs(binned.vth0 - 1.0) < 1.0e-12);

    const auto passiveBinned = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"CGSO", 1.0e-10}, {"LCGSO", 1.0e-16},
        {"IS", 1.0e-14}, {"LIS", 1.0e-20}, {"XJ", 1.5e-7},
        {"LXJ", 1.0e-14}
    }).prepare(1.0e-6, 1.0e-6);
    assert(passiveBinned.validation);
    // Geometry bins divide the L/W/P coefficients by Leff/Weff (both 1e-6 m
    // here, no DL/DW): cgso = 1e-10 + 1e-16/1e-6 = 2e-10.
    assert(std::abs(passiveBinned.cgso - 2e-10) < 1.0e-16);
    assert(std::abs(passiveBinned.is - 2e-14) < 1.0e-20);
    assert(std::abs(passiveBinned.xj - 1.6e-7) < 1.0e-14);

    const auto bodyBinned = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"K1", 0.5}, {"LK1", 1.0e-7},
        {"ETA0", 0.08}, {"LETA0", 1.0e-8}, {"RDSWMIN", 10.0},
        {"LRDSWMIN", 1.0e-14}
    }).prepare(1.0e-6, 1.0e-6);
    assert(bodyBinned.validation);
    assert(std::abs(bodyBinned.k1 - 0.6) < 1.0e-12);
    assert(std::abs(bodyBinned.eta0 - 0.09) < 1.0e-12);
    assert(std::abs(bodyBinned.rdswmin - 10.00000001) < 1.0e-10);

    const auto hot = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"U0", 500.0}, {"VSAT", 1.0e5},
        {"TOXE", 1.0e-8}, {"UA", 2.0e-9}, {"UB", 5.0e-19},
        {"UC", 1.0e-10}, {"AT", 0.005}
    }).prepare(1.0e-6, 1.0e-6, 100.0);
    assert(hot.validation);
    assert(hot.u0 < 0.05);
    assert(hot.vsat < 1.0e5);
    assert(hot.ua > 2.0e-9);
    assert(hot.ub < 5.0e-19);
    assert(hot.uc < 1.0e-10);

    const auto hotJunction = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"JS", 1.0e-14}, {"CJ", 1.0e-3}, {"CJSW", 2.0e-10},
        {"PB", 0.8}, {"PBSW", 0.8}, {"TCJ", 1.0e-3},
        {"TCJSW", 2.0e-3}, {"TPB", -1.0e-3}, {"TPBSW", -2.0e-3},
        {"AT", 0.005}
    }).prepare(1.0e-6, 1.0e-6, 100.0);
    assert(hotJunction.validation);
    assert(hotJunction.cj > 1.0e-3);
    assert(hotJunction.cjsw > 2.0e-10);
    assert(hotJunction.js > 1.0e-14);
    assert(hotJunction.pb < 0.8);
    assert(hotJunction.pbsw < 0.8);

    const auto invalid = gspice::Bsim4ParameterSet::from({{"LEVEL", 49.0}})
                             .prepare(1.0e-6, 1.0e-6);
    assert(!invalid.validation);
}
