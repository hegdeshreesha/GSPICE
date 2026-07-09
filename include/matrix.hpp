#ifndef GSPICE_MATRIX_HPP
#define GSPICE_MATRIX_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <atomic>
#include <omp.h>

namespace gspice {

/**
 * Templated Vector class.
 * Support for both double (DC/Tran) and std::complex<double> (AC/HB).
 * Thread-safe for parallel stamping.
 */
template <typename T>
class Vector {
public:
    Vector(int size = 0) : data_(size, T(0)) {
        initializeThreadBuffers();
    }

    Vector(const Vector& other) : data_(other.snapshotData()) {
        initializeThreadBuffers();
    }

    Vector& operator=(const Vector& other) {
        if (this != &other) {
            data_ = other.snapshotData();
            initializeThreadBuffers();
        }
        return *this;
    }

    T& operator[](int index) {
        collapseThreadData();
        return data_[index];
    }

    const T& operator[](int index) const {
        collapseThreadData();
        return data_[index];
    }

    void add(int index, T value) {
        if (index >= 0 && index < (int)data_.size()) {
            if (omp_in_parallel()) {
                thread_data_[omp_get_thread_num()][index] += value;
                has_thread_data_.store(true, std::memory_order_relaxed);
            } else {
                data_[index] += value;
            }
        }
    }

    int getSize() const { return (int)data_.size(); }

    void clear() {
        collapseThreadData();
        std::fill(data_.begin(), data_.end(), T(0));
        clearThreadBuffers();
    }

    std::vector<T> snapshotData() const {
        collapseThreadData();
        return data_;
    }

    void print() const {
        collapseThreadData();
        for (const T& val : data_) {
            if constexpr (std::is_same_v<T, std::complex<double>>) {
                std::cout << std::setw(15) << val.real() << " + " << val.imag() << "j" << std::endl;
            } else {
                std::cout << std::setw(10) << val << std::endl;
            }
        }
    }

private:
    void initializeThreadBuffers() {
        const int threads = std::max(1, omp_get_max_threads());
        thread_data_.assign(static_cast<size_t>(threads), std::vector<T>(data_.size(), T(0)));
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    void clearThreadBuffers() {
        for (auto& buffer : thread_data_) {
            std::fill(buffer.begin(), buffer.end(), T(0));
        }
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    void collapseThreadData() const {
        if (!has_thread_data_.load(std::memory_order_relaxed)) return;
        for (auto& buffer : thread_data_) {
            for (size_t i = 0; i < data_.size(); ++i) {
                data_[i] += buffer[i];
                buffer[i] = T(0);
            }
        }
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    mutable std::vector<T> data_;
    mutable std::vector<std::vector<T>> thread_data_;
    mutable std::atomic<bool> has_thread_data_{false};
};

/**
 * Templated Matrix class.
 */
template <typename T>
class Matrix {
public:
    Matrix() : size_(0), data_() {}
    Matrix(int size) : size_(size), data_(size * size, T(0)) {}

    T& operator()(int row, int col) {
        return data_[row * size_ + col];
    }

    const T& operator()(int row, int col) const {
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

    Vector<T> solve(const Vector<T>& b_in) {
        int n = size_;
        std::vector<T> A = data_;
        Vector<T> x(n);
        Vector<T> b = b_in;

        for (int i = 0; i < n; i++) {
            double maxVal = std::abs(A[i * n + i]);
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (std::abs(A[k * n + i]) > maxVal) {
                    maxVal = std::abs(A[k * n + i]);
                    maxRow = k;
                }
            }
            for (int k = i; k < n; k++) std::swap(A[maxRow * n + k], A[i * n + k]);
            std::swap(b[maxRow], b[i]);

            if (std::abs(A[i * n + i]) < 1e-25) {
                throw std::runtime_error("Singular matrix: zero pivot at row " + std::to_string(i));
            }
            for (int k = i + 1; k < n; k++) {
                T c = A[k * n + i] / A[i * n + i];
                for (int j = i; j < n; j++) {
                    A[k * n + j] -= c * A[i * n + j];
                }
                b[k] -= c * b[i];
            }
        }

        for (int i = n - 1; i >= 0; i--) {
            if (std::abs(A[i * n + i]) < 1e-25) {
                throw std::runtime_error("Singular matrix during back substitution at row " + std::to_string(i));
            }
            T sum = b[i];
            for (int k = i + 1; k < n; k++) {
                sum -= A[i * n + k] * x[k];
            }
            x[i] = sum / A[i * n + i];
        }
        return x;
    }

    Vector<T> solveTranspose(const Vector<T>& b_in) {
        int n = size_;
        Matrix<T> At(n);
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

using MatrixReal = Matrix<double>;
using VectorReal = Vector<double>;
using MatrixComplex = Matrix<std::complex<double>>;
using VectorComplex = Vector<std::complex<double>>;

} // namespace gspice

#endif // GSPICE_MATRIX_HPP
