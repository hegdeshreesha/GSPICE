#ifndef GSPICE_SPARSE_MATRIX_HPP
#define GSPICE_SPARSE_MATRIX_HPP

#include <vector>
#include <complex>
#include <algorithm>
#include <atomic>
#include <omp.h>
#include "matrix.hpp"

namespace gspice {

/**
 * A Sparse Matrix representation using a Coordinate Map (Triplet).
 * Thread-safe for parallel stamping.
 */
template <typename T>
class SparseMatrix {
public:
    struct Entry {
        int row;
        int col;
        T value;
    };

    SparseMatrix(int size) : size_(size) {
        initializeThreadBuffers();
    }

    void add(int row, int col, T value) {
        if (row >= 0 && col >= 0 && row < size_ && col < size_) {
            if (omp_in_parallel()) {
                thread_triplets_[static_cast<size_t>(omp_get_thread_num())].push_back({row, col, value});
                has_thread_data_.store(true, std::memory_order_relaxed);
            } else {
                triplets_.push_back({row, col, value});
            }
        }
    }

    int getSize() const { return size_; }

    /**
     * Returns a dense version for debugging or small solves.
     */
    Matrix<T> toDense() const {
        collapseThreadData();
        Matrix<T> dense(size_);
        for (const auto& entry : combinedEntries()) {
            dense(entry.row, entry.col) = entry.value;
        }
        return dense;
    }

    /**
     * Clear the matrix for a new iteration.
     */
    void clear() {
        collapseThreadData();
        triplets_.clear();
        clearThreadBuffers();
    }

    std::vector<Entry> getEntries() const {
        collapseThreadData();
        return combinedEntries();
    }

    /**
     * Converts the matrix to Compressed Column Storage (CCS) format.
     * Required for KLU and other high-performance sparse solvers.
     */
    void toCCS(std::vector<int>& Ap, std::vector<int>& Ai, std::vector<T>& Ax) const {
        const auto entries = getEntries();
        Ap.assign(size_ + 1, 0);
        Ai.clear();
        Ax.clear();
        std::vector<std::vector<std::pair<int, T>>> columns(static_cast<size_t>(size_));
        for (const auto& entry : entries) {
            columns[static_cast<size_t>(entry.col)].push_back({entry.row, entry.value});
        }

        int nnz = 0;
        for (int j = 0; j < size_; ++j) {
            Ap[j] = nnz;
            for (const auto& [row, value] : columns[static_cast<size_t>(j)]) {
                Ai.push_back(row);
                Ax.push_back(value);
                ++nnz;
            }
        }
        Ap[size_] = nnz;
    }

private:
    void initializeThreadBuffers() {
        const int threads = std::max(1, omp_get_max_threads());
        thread_triplets_.assign(static_cast<size_t>(threads), {});
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    void clearThreadBuffers() {
        for (auto& bucket : thread_triplets_) {
            bucket.clear();
        }
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    void collapseThreadData() const {
        if (!has_thread_data_.load(std::memory_order_relaxed)) return;
        for (auto& bucket : thread_triplets_) {
            triplets_.insert(triplets_.end(), bucket.begin(), bucket.end());
            bucket.clear();
        }
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    std::vector<Entry> combinedEntries() const {
        std::vector<Entry> entries = triplets_;
        std::sort(entries.begin(), entries.end(), [](const Entry& lhs, const Entry& rhs) {
            if (lhs.row != rhs.row) return lhs.row < rhs.row;
            return lhs.col < rhs.col;
        });

        std::vector<Entry> combined;
        combined.reserve(entries.size());
        for (const auto& entry : entries) {
            if (combined.empty() || combined.back().row != entry.row || combined.back().col != entry.col) {
                combined.push_back(entry);
            } else {
                combined.back().value += entry.value;
            }
        }
        return combined;
    }

    int size_;
    mutable std::vector<Entry> triplets_;
    mutable std::vector<std::vector<Entry>> thread_triplets_;
    mutable std::atomic<bool> has_thread_data_{false};
};

using SparseMatrixReal = SparseMatrix<double>;
using SparseMatrixComplex = SparseMatrix<std::complex<double>>;

} // namespace gspice

#endif // GSPICE_SPARSE_MATRIX_HPP
