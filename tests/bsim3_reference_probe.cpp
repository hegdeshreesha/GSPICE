#include "devices/bsim3_dc_core.hpp"

#include <iomanip>
#include <iostream>

int main() {
    const auto model = gspice::Bsim3ParameterSet::from({
        {"LEVEL", 49.0}, {"VTH0", 0.40}, {"KP", 120.0e-6},
        {"GAMMA", 0.50}, {"PHI", 0.60}, {"U0", 500.0},
        {"TOXE", 1.0e-8}, {"XJ", 0.15e-6}
    }).prepare(1.0e-6, 1.0e-6);
    std::cout << std::scientific << std::setprecision(12);
    for (int i = 0; i <= 10; ++i) {
        const auto result = gspice::bsim3EvaluateDc(model, {0.8, 0.1 * i, 0.0, 0.0});
        if (!result.valid) return 2;
        std::cout << 0.1 * i << " " << result.current[0] << "\n";
    }
}
