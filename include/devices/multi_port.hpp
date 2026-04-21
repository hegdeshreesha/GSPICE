#ifndef GSPICE_MULTIPORT_HPP
#define GSPICE_MULTIPORT_HPP

#include "device.hpp"
#include "touchstone.hpp"
#include <string>
#include <vector>
#include <complex>

namespace gspice {

class MultiPort : public Device {
public:
    MultiPort(const std::string& name, const std::vector<int>& nodes, int num_ports, const std::string& filename)
        : Device(name), nodes_(nodes), num_ports_(num_ports), ts_(filename, num_ports) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        // DC: For simplicity, treat as a direct connection (short) for DC with a small conductance.
        // Alternatively, use the S-parameters at f=0 to compute a DC conductance matrix.
        // For robustness, just add a small G to ground to prevent floating nodes.
        double G = 1e-6; 
        for (size_t i = 0; i < nodes_.size(); ++i) {
            J.add(nodes_[i], nodes_[i], G);
        }
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        double freq = omega / (2.0 * M_PI);
        MatrixComplex S = ts_.getS(freq);
        double z0 = ts_.getZ0();

        // Convert S to Y: (I + S) * (Y * Z0) = (I - S)
        MatrixComplex A(num_ports_);
        MatrixComplex B(num_ports_);
        for (int i = 0; i < num_ports_; ++i) {
            for (int j = 0; j < num_ports_; ++j) {
                std::complex<double> delta = (i == j) ? 1.0 : 0.0;
                A(i, j) = delta + S(i, j);
                B(i, j) = delta - S(i, j);
            }
        }

        MatrixComplex Y(num_ports_);
        for (int col = 0; col < num_ports_; ++col) {
            VectorComplex B_col(num_ports_);
            for (int row = 0; row < num_ports_; ++row) B_col[row] = B(row, col);
            
            // Solve A * X_col = B_col
            VectorComplex X_col = A.solve(B_col);
            for (int row = 0; row < num_ports_; ++row) {
                Y(row, col) = X_col[row] / z0;
            }
        }

        // Stamp Y into MNA
        for (int r = 0; r < num_ports_; ++r) {
            for (int c = 0; c < num_ports_; ++c) {
                J.add(nodes_[r], nodes_[c], Y(r, c));
            }
        }
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        int K = 2 * n_harms + 1;
        for (int h = 1; h <= n_harms; ++h) {
            double freq_h = f_fund * h;
            MatrixComplex S = ts_.getS(freq_h);
            double z0 = ts_.getZ0();

            MatrixComplex A(num_ports_);
            MatrixComplex B(num_ports_);
            for (int i = 0; i < num_ports_; ++i) {
                for (int j = 0; j < num_ports_; ++j) {
                    std::complex<double> delta = (i == j) ? 1.0 : 0.0;
                    A(i, j) = delta + S(i, j);
                    B(i, j) = delta - S(i, j);
                }
            }

            MatrixComplex Y(num_ports_);
            for (int col = 0; col < num_ports_; ++col) {
                VectorComplex B_col(num_ports_);
                for (int row = 0; row < num_ports_; ++row) B_col[row] = B(row, col);
                VectorComplex X_col = A.solve(B_col);
                for (int row = 0; row < num_ports_; ++row) {
                    Y(row, col) = X_col[row] / z0;
                }
            }

            int cos_idx = 2 * h - 1; 
            int sin_idx = 2 * h;

            for (int r = 0; r < num_ports_; ++r) {
                for (int c = 0; c < num_ports_; ++c) {
                    double g = Y(r, c).real();
                    double b_sus = Y(r, c).imag();

                    // Standard Harmonic Balance linear admittance stamp
                    J.add(nodes_[r] * K + cos_idx, nodes_[c] * K + cos_idx, g);
                    J.add(nodes_[r] * K + cos_idx, nodes_[c] * K + sin_idx, -b_sus);
                    J.add(nodes_[r] * K + sin_idx, nodes_[c] * K + cos_idx, b_sus);
                    J.add(nodes_[r] * K + sin_idx, nodes_[c] * K + sin_idx, g);
                }
            }
        }
        
        // DC Stamp for HB (h=0)
        double G = 1e-6;
        for (int r = 0; r < num_ports_; ++r) {
            J.add(nodes_[r] * K, nodes_[r] * K, G);
        }
    }

private:
    std::vector<int> nodes_;
    int num_ports_;
    Touchstone ts_;
};

} // namespace gspice

#endif // GSPICE_MULTIPORT_HPP
