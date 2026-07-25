#ifndef GSPICE_ADJOINT_SOLVER_HPP
#define GSPICE_ADJOINT_SOLVER_HPP

#include "../matrix.hpp"
#include <vector>
#include <cmath>

namespace gspice {

// Solves transpose system J^T * X_adj = dOut/dX to compute sensitivities 
// for all circuit parameters in a single linear solve.
class AdjointSolver {
public:
    static bool solveAdjoint(
        const std::vector<std::vector<double>>& jacobian,
        const std::vector<double>& rhsSens,
        std::vector<double>& adjointSolution) {

        int n = static_cast<int>(jacobian.size());
        if (n == 0 || rhsSens.size() != static_cast<size_t>(n)) return false;

        // Form transpose matrix J^T
        std::vector<std::vector<double>> jTrans(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                jTrans[i][j] = jacobian[j][i];
            }
        }

        // Solve J^T * X_adj = rhsSens via Gaussian elimination with partial pivoting
        std::vector<std::vector<double>> A = jTrans;
        adjointSolution = rhsSens;

        for (int i = 0; i < n; ++i) {
            int maxRow = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }

            if (std::abs(A[maxRow][i]) < 1e-15) continue;

            std::swap(A[i], A[maxRow]);
            std::swap(adjointSolution[i], adjointSolution[maxRow]);

            for (int k = i + 1; k < n; ++k) {
                double factor = A[k][i] / A[i][i];
                adjointSolution[k] -= factor * adjointSolution[i];
                for (int j = i; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
            }
        }

        // Back-substitution
        for (int i = n - 1; i >= 0; --i) {
            if (std::abs(A[i][i]) < 1e-15) {
                adjointSolution[i] = 0.0;
                continue;
            }
            for (int k = i + 1; k < n; ++k) {
                adjointSolution[i] -= A[i][k] * adjointSolution[k];
            }
            adjointSolution[i] /= A[i][i];
        }

        return true;
    }
};

} // namespace gspice

#endif // GSPICE_ADJOINT_SOLVER_HPP
