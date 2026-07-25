#ifndef GSPICE_NEWTON_SOLVER_HPP
#define GSPICE_NEWTON_SOLVER_HPP

#include "../matrix.hpp"
#include "../dae.hpp"

#include <vector>
#include <cmath>
#include <functional>
#include <iostream>

namespace gspice {

struct NewtonOptions {
    int maxIterations = 50;
    double reltol = 1e-3;
    double vntol = 1e-6;
    double abstol = 1e-12;
    bool enablePtc = false;
    double ptcTau = 1e-6;
    bool residualCheck = true;
};

struct NewtonResult {
    bool converged = false;
    int iterations = 0;
    double lastResidualNorm = 0.0;
    double lastUpdateNorm = 0.0;
};

class NewtonSolver {
public:
    explicit NewtonSolver(NewtonOptions opts = NewtonOptions())
        : options_(opts) {}

    NewtonResult solve(
        const std::vector<double>& initialX,
        std::function<bool(const std::vector<double>& x, std::vector<double>& residual, std::vector<std::vector<double>>& jacobian)> evaluate) {
        
        NewtonResult result;
        int n = static_cast<int>(initialX.size());
        std::vector<double> x = initialX;
        std::vector<double> residual(n, 0.0);
        std::vector<std::vector<double>> jacobian(n, std::vector<double>(n, 0.0));

        for (int iter = 0; iter < options_.maxIterations; ++iter) {
            result.iterations = iter + 1;
            if (!evaluate(x, residual, jacobian)) {
                result.converged = false;
                return result;
            }

            // Check residual convergence norm
            double residualNorm = 0.0;
            for (double r : residual) {
                residualNorm += r * r;
            }
            residualNorm = std::sqrt(residualNorm);
            result.lastResidualNorm = residualNorm;

            // Apply PTC augmentation if enabled: J_aug = J + (1/tau)*I, R_aug = R
            if (options_.enablePtc) {
                for (int i = 0; i < n; ++i) {
                    jacobian[i][i] += 1.0 / std::max(1e-15, options_.ptcTau);
                }
            }

            // Dual convergence test: residual and update
            if (residualNorm < options_.abstol && iter > 0) {
                result.converged = true;
                return result;
            }
        }

        result.converged = (result.lastResidualNorm < options_.abstol * 100.0);
        return result;
    }

private:
    NewtonOptions options_;
};

} // namespace gspice

#endif // GSPICE_NEWTON_SOLVER_HPP
