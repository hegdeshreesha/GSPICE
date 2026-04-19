#ifndef GSPICE_SPARSE_MATRIX_HPP
#define GSPICE_SPARSE_MATRIX_HPP

#include <vector>
#include <map>
#include <complex>
#include "matrix.hpp"

namespace gspice {

/**
 * A Sparse Matrix representation using a Coordinate Map (Triplet) 
 * for easy assembly, which can then be converted to CCS for KLU/Eigen.
 */
template <typename T>
class SparseMatrix {
public:
    SparseMatrix(int size) : size_(size) {}

    /**
     * Add a value to the sparse matrix.
     * Coordinate format is easiest for 'stamping' logic.
     */
    void add(int row, int col, T value) {
        if (row >= 0 && col >= 0 && row < size_ && col < size_) {
            triplets_[{row, col}] += value;
        }
    }

    int getSize() const { return size_; }

    /**
     * Returns a dense version for debugging or small solves.
     */
    Matrix<T> toDense() const {
        Matrix<T> dense(size_);
        for (auto const& [pos, val] : triplets_) {
            dense(pos.first, pos.second) = val;
        }
        return dense;
    }

    /**
     * Clear the matrix for a new iteration.
     */
    void clear() {
        triplets_.clear();
    }

    // Access to internal map for solver bridges
    const std::map<std::pair<int, int>, T>& getTriplets() const {
        return triplets_;
    }

    /**
     * Converts the matrix to Compressed Column Storage (CCS) format.
     * Required for KLU and other high-performance sparse solvers.
     */
    void toCCS(std::vector<int>& Ap, std::vector<int>& Ai, std::vector<T>& Ax) const {
        Ap.assign(size_ + 1, 0);
        Ai.clear();
        Ax.clear();

        // Sort by column then row (map already sorts by key which is {row, col})
        // We need {col, row} order for CCS
        std::map<int, std::map<int, T>> col_major;
        for (auto const& [pos, val] : triplets_) {
            col_major[pos.second][pos.first] = val;
        }

        int nnz = 0;
        for (int j = 0; j < size_; j++) {
            Ap[j] = nnz;
            if (col_major.count(j)) {
                for (auto const& [i, val] : col_major[j]) {
                    Ai.push_back(i);
                    Ax.push_back(val);
                    nnz++;
                }
            }
        }
        Ap[size_] = nnz;
    }

private:
    int size_;
    // map of {row, col} -> value
    std::map<std::pair<int, int>, T> triplets_;
};

using SparseMatrixReal = SparseMatrix<double>;
using SparseMatrixComplex = SparseMatrix<std::complex<double>>;

} // namespace gspice

#endif // GSPICE_SPARSE_MATRIX_HPP
