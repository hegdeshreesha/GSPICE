#ifndef GSPICE_MATRIX_HPP
#define GSPICE_MATRIX_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

namespace gspice {

/**
 * Templated Vector class.
 * Support for both double (DC/Tran) and std::complex<double> (AC/HB).
 */
template <typename T>
class Vector {
public:
    Vector(int size = 0) : data_(size, T(0)) {}

    T& operator[](int index) {
        return data_[index];
    }

    const T& operator[](int index) const {
        return data_[index];
    }

    void add(int index, T value) {
        if (index >= 0 && index < (int)data_.size()) {
            data_[index] += value;
        }
    }

    int getSize() const { return (int)data_.size(); }

    void print() const {
        for (const T& val : data_) {
            if constexpr (std::is_same_v<T, std::complex<double>>) {
                std::cout << std::setw(15) << val.real() << " + " << val.imag() << "j" << std::endl;
            } else {
                std::cout << std::setw(10) << val << std::endl;
            }
        }
    }

private:
    std::vector<T> data_;
};

/**
 * Templated Matrix class.
 */
template <typename T>
class Matrix {
public:
    Matrix(int size) : size_(size), data_(size * size, T(0)) {}

    T& operator()(int row, int col) {
        return data_[row * size_ + col];
    }

    void add(int row, int col, T value) {
        if (row >= 0 && col >= 0 && row < size_ && col < size_) {
            data_[row * size_ + col] += value;
        }
    }

    void print() const {
        for (int i = 0; i < size_; ++i) {
            for (int j = 0; j < size_; ++j) {
                if constexpr (std::is_same_v<T, std::complex<double>>) {
                     std::cout << "(" << data_[i * size_ + j].real() << "," << data_[i * size_ + j].imag() << ") ";
                } else {
                    std::cout << std::setw(10) << data_[i * size_ + j] << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    int getSize() const { return size_; }

    /**
     * Solves Ax = b using Gaussian Elimination with partial pivoting.
     */
    Vector<T> solve(const Vector<T>& b_in) {
        int n = size_;
        std::vector<T> A = data_;
        Vector<T> x(n);
        Vector<T> b = b_in;

        for (int i = 0; i < n; i++) {
            // 1. Pivot search
            double maxVal = std::abs(A[i * n + i]);
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (std::abs(A[k * n + i]) > maxVal) {
                    maxVal = std::abs(A[k * n + i]);
                    maxRow = k;
                }
            }

            // 2. Swap
            for (int k = i; k < n; k++) std::swap(A[maxRow * n + k], A[i * n + k]);
            std::swap(b[maxRow], b[i]);

            // 3. Elimination
            if (std::abs(A[i * n + i]) < 1e-25) continue; 
            for (int k = i + 1; k < n; k++) {
                T c = A[k * n + i] / A[i * n + i];
                for (int j = i; j < n; j++) {
                    A[k * n + j] -= c * A[i * n + j];
                }
                b[k] -= c * b[i];
            }
        }

        // 4. Back substitution
        for (int i = n - 1; i >= 0; i--) {
            if (std::abs(A[i * n + i]) < 1e-25) {
                x[i] = T(0);
                continue;
            }
            T sum = b[i];
            for (int k = i + 1; k < n; k++) {
                sum -= A[i * n + k] * x[k];
            }
            x[i] = sum / A[i * n + i];
        }
        return x;
    }

    /**
     * Solves A^T * x = b (The Adjoint System).
     * Used for Noise and Sensitivity analysis.
     */
    Vector<T> solveTranspose(const Vector<T>& b_in) {
        int n = size_;
        Matrix<T> At(n);
        // Transpose the matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                At(j, i) = (*this)(i, j);
            }
        }
        return At.solve(b_in);
    }

private:
    int size_;
    std::vector<T> data_;
};

// Typedefs for convenience
using MatrixReal = Matrix<double>;
using VectorReal = Vector<double>;
using MatrixComplex = Matrix<std::complex<double>>;
using VectorComplex = Vector<std::complex<double>>;

} // namespace gspice

#endif // GSPICE_MATRIX_HPP
