#ifndef GSPICE_SPARSE_MATRIX_HPP
#define GSPICE_SPARSE_MATRIX_HPP

#include <vector>
#include <complex>
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <functional>
#include <omp.h>
#include <unordered_map>
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

    struct EntryKey {
        int row = 0;
        int col = 0;

        bool operator==(const EntryKey& other) const {
            return row == other.row && col == other.col;
        }
    };

    struct EntryKeyHash {
        std::size_t operator()(const EntryKey& key) const noexcept {
            std::size_t seed = static_cast<std::size_t>(key.row + 1) * 1469598103934665603ull;
            seed ^= static_cast<std::size_t>(key.col + 1) * 1099511628211ull;
            seed += 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    SparseMatrix(int size) : size_(size) {
        initializeThreadBuffers();
    }

    void add(int row, int col, T value) {
        if (row >= 0 && col >= 0 && row < size_ && col < size_) {
            if (addToCachedSlot(row, col, value)) {
                return;
            }
            if (omp_in_parallel()) {
                const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());
                if (tid < thread_triplets_.size()) {
                    thread_triplets_[tid].push_back({row, col, value});
                    has_thread_data_.store(true, std::memory_order_relaxed);
                } else {
                    triplets_.push_back({row, col, value});
                }
            } else {
                triplets_.push_back({row, col, value});
            }
        }
    }

    int getSize() const { return size_; }

    void setStructureCacheEnabled(bool enabled) {
        structure_cache_enabled_ = enabled;
        if (!enabled) {
            clearStructureCache();
        }
    }

    bool structureCacheEnabled() const { return structure_cache_enabled_; }
    bool structureCacheReady() const { return structure_cache_valid_; }

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
        triplets_.clear();
        clearThreadBuffers();
        if (structure_cache_enabled_ && structure_cache_valid_) {
            for (auto& entry : cached_entries_) {
                entry.value = T(0);
            }
            clearCachedThreadBuffers();
            cache_pattern_miss_ = false;
        } else {
            clearStructureCache();
        }
    }

    std::vector<Entry> getEntries() const {
        collapseThreadData();
        collapseCachedThreadData();
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
        thread_cached_values_.assign(static_cast<size_t>(threads), {});
        thread_cached_dirty_.assign(static_cast<size_t>(threads), false);
        has_thread_data_.store(false, std::memory_order_relaxed);
        has_cached_thread_data_.store(false, std::memory_order_relaxed);
    }

    void clearThreadBuffers() {
        for (auto& bucket : thread_triplets_) {
            bucket.clear();
        }
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    void clearCachedThreadBuffers() const {
        for (auto& buffer : thread_cached_values_) {
            std::fill(buffer.begin(), buffer.end(), T(0));
        }
        std::fill(thread_cached_dirty_.begin(), thread_cached_dirty_.end(), false);
        has_cached_thread_data_.store(false, std::memory_order_relaxed);
    }

    void collapseThreadData() const {
        if (!has_thread_data_.load(std::memory_order_relaxed)) return;
        for (auto& bucket : thread_triplets_) {
            triplets_.insert(triplets_.end(), bucket.begin(), bucket.end());
            bucket.clear();
        }
        has_thread_data_.store(false, std::memory_order_relaxed);
    }

    void collapseCachedThreadData() const {
        if (!has_cached_thread_data_.load(std::memory_order_relaxed)) return;
        for (std::size_t t = 0; t < thread_cached_values_.size(); ++t) {
            if (!thread_cached_dirty_[t]) continue;
            auto& buffer = thread_cached_values_[t];
            const std::size_t limit = std::min(buffer.size(), cached_entries_.size());
            for (std::size_t i = 0; i < limit; ++i) {
                cached_entries_[i].value += buffer[i];
                buffer[i] = T(0);
            }
            thread_cached_dirty_[t] = false;
        }
        has_cached_thread_data_.store(false, std::memory_order_relaxed);
    }

    bool addToCachedSlot(int row, int col, T value) {
        if (!structure_cache_enabled_ || !structure_cache_valid_ || cache_pattern_miss_) {
            return false;
        }
        auto it = cached_slot_map_.find({row, col});
        if (it == cached_slot_map_.end()) {
            cache_pattern_miss_ = true;
            return false;
        }
        const std::size_t slot = it->second;
        if (omp_in_parallel()) {
            const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());
            if (tid < thread_cached_values_.size() && slot < thread_cached_values_[tid].size()) {
                thread_cached_values_[tid][slot] += value;
                thread_cached_dirty_[tid] = true;
                has_cached_thread_data_.store(true, std::memory_order_relaxed);
                return true;
            }
        }
        cached_entries_[slot].value += value;
        return true;
    }

    void rebuildStructureCache(const std::vector<Entry>& entries) const {
        cached_entries_ = entries;
        cached_slot_map_.clear();
        cached_slot_map_.reserve(cached_entries_.size() * 2);
        for (std::size_t i = 0; i < cached_entries_.size(); ++i) {
            cached_slot_map_[{cached_entries_[i].row, cached_entries_[i].col}] = i;
        }
        for (auto& buffer : thread_cached_values_) {
            buffer.assign(cached_entries_.size(), T(0));
        }
        std::fill(thread_cached_dirty_.begin(), thread_cached_dirty_.end(), false);
        structure_cache_valid_ = true;
        cache_pattern_miss_ = false;
        has_cached_thread_data_.store(false, std::memory_order_relaxed);
        triplets_.clear();
    }

    void clearStructureCache() const {
        cached_entries_.clear();
        cached_slot_map_.clear();
        for (auto& buffer : thread_cached_values_) {
            buffer.clear();
        }
        std::fill(thread_cached_dirty_.begin(), thread_cached_dirty_.end(), false);
        structure_cache_valid_ = false;
        cache_pattern_miss_ = false;
        has_cached_thread_data_.store(false, std::memory_order_relaxed);
    }

    std::vector<Entry> combinedEntries() const {
        if (structure_cache_enabled_ && structure_cache_valid_ && !cache_pattern_miss_ && triplets_.empty()) {
            return cached_entries_;
        }
        std::vector<Entry> entries = triplets_;
        if (structure_cache_enabled_ && structure_cache_valid_) {
            entries.insert(entries.end(), cached_entries_.begin(), cached_entries_.end());
        }
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
        if (structure_cache_enabled_) {
            rebuildStructureCache(combined);
        }
        return combined;
    }

    int size_;
    mutable std::vector<Entry> triplets_;
    mutable std::vector<std::vector<Entry>> thread_triplets_;
    mutable std::atomic<bool> has_thread_data_{false};
    bool structure_cache_enabled_ = false;
    mutable bool structure_cache_valid_ = false;
    mutable bool cache_pattern_miss_ = false;
    mutable std::vector<Entry> cached_entries_;
    mutable std::unordered_map<EntryKey, std::size_t, EntryKeyHash> cached_slot_map_;
    mutable std::vector<std::vector<T>> thread_cached_values_;
    mutable std::vector<bool> thread_cached_dirty_;
    mutable std::atomic<bool> has_cached_thread_data_{false};
};

using SparseMatrixReal = SparseMatrix<double>;
using SparseMatrixComplex = SparseMatrix<std::complex<double>>;

} // namespace gspice

#endif // GSPICE_SPARSE_MATRIX_HPP
