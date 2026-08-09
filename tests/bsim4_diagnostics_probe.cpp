#include "devices/bsim4_dc_core.hpp"

#include <iomanip>
#include <iostream>

int main() {
    const auto model = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"VTH0", 0.40}, {"U0", 500.0},
        {"TOXE", 1.0e-8}, {"XJ", 0.15e-6}, {"K1", 0.0}, {"K2", 0.0},
        {"DVT0", 0.0}, {"DVT1", 0.53}, {"DVT2", -0.032}, {"UA", 2.0e-9},
        {"UB", 5.0e-19}, {"VSAT", 1.0e5}
    }).prepare(1.0e-6, 1.0e-6);

    std::cout << std::scientific << std::setprecision(12);
    for (int i = 0; i <= 10; ++i) {
        const double gate = 0.1 * i;
        const auto result = gspice::bsim4EvaluateDc(model, {0.8, gate, 0.0, 0.0});
        if (!result.valid) return 2;
        std::cout << gate << ' ' << result.current[0] << ' '
                  << result.state.vth << ' ' << result.state.vgtEff << ' '
                  << result.state.vdsat << ' ' << result.state.abulk << ' '
                  << result.state.vasat << ' ' << result.jacobian[0][1]
                  << ' ' << result.jacobian[0][0] << '\n';
    }
}
