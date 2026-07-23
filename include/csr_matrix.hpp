#ifndef GSPICE_CSR_MATRIX_HPP
#define GSPICE_CSR_MATRIX_HPP

// ---------------------------------------------------------------------------
// FixedCsrMatrix — Zero-copy, sort-free sparse matrix for Newton iteration.
//
// Design:
//   Phase A (symbolic, once per circuit topology):
//     1. Call beginSymbolic() to enter registration mode.
//     2. Every device calls registerEntry(row, col) for each nonzero it stamps.
//     3. Call finalizeSymbolic(). This deduplicates and sorts the pattern once,
//        builds the CSC Ap[]/Ai[] arrays for KLU, and creates a hash map from
//        (row,col) to the index in Ax[] for O(1) device writes.
//
//   Phase B (numeric, every Newton iteration):
//     1. Call zeroValues() — single memset on Ax[].
//     2. Each device calls add(row, col, value) — resolves to Ax[offset] += value
//        via a direct index lookup; no sort, no dynamic allocation.
//     3. Pass getAp(), getAi(), getAx() directly to KLU klu_refactor().
//
// Thread safety:
//   zeroValues() and add() are thread-safe under the same condition as the
//   original SparseMatrix: device stamps issued from OpenMP parallel regions
//   use per-thread scratch buffers that are flushed with flushThreadBuffers().
//   Serial stamps (outside parallel regions) write directly to Ax[].
// ---------------------------------------------------------------------------

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstring>
#include <omp.h>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace gspice {

// ---- Slot handle returned by registerEntry() --------------------------------
// Devices store these during elaboration and pass them to addViaSlot() during
// Newton iterations. Using the handle avoids hash-map lookups in the hot path.
using CsrSlot = int; // index into Ax_[]; -1 means out-of-range (ground)

// ---- Key for the (row, col) → index map ------------------------------------
struct CsrKey {
    int row, col;
    bool operator==(const CsrKey& o) const { return row == o.row && col == o.col; }
};

struct CsrKeyHash {
    std::size_t operator()(const CsrKey& k) const noexcept {
        // FNV-like mix of two ints
        std::size_t h = static_cast<std::size_t>(k.row) * 2654435761ull;
        h ^= static_cast<std::size_t>(k.col) * 2246822519ull;
        return h;
    }
};

// ---------------------------------------------------------------------------
template <typename T>
class FixedCsrMatrix {
public:
    explicit FixedCsrMatrix(int size)
        : n_(size), symbolic_done_(false) {
        const int threads = std::max(1, omp_get_max_threads());
        thread_scratch_.resize(static_cast<std::size_t>(threads));
    }

    int getSize() const { return n_; }

    // -----------------------------------------------------------------------
    // Phase A — Symbolic registration
    // -----------------------------------------------------------------------

    /// Enter symbolic-build mode. Must be called before registerEntry().
    void beginSymbolic() {
        symbolic_done_ = false;
        pending_entries_.clear();
        slot_map_.clear();
        Ap_.clear();
        Ai_.clear();
        Ax_.clear();
    }

    /// Register a nonzero at (row, col). Out-of-range entries are silently
    /// ignored (like ground connections). Returns the slot handle that can be
    /// stored by the device for fast numeric writes.
    CsrSlot registerEntry(int row, int col) {
        if (row < 0 || col < 0 || row >= n_ || col >= n_) {
            return -1; // ground / out-of-range
        }
        pending_entries_.push_back({row, col});
        return -2; // placeholder; real slot assigned in finalizeSymbolic()
    }

    /// Finalize the symbolic phase. Must be called after all devices have
    /// called registerEntry(). Builds the deduplicated CSC pattern and the
    /// (row,col)→index map. O(nnz log nnz).
    void finalizeSymbolic() {
        // Sort and deduplicate pending (row, col) pairs (column-major for KLU).
        auto& pe = pending_entries_;
        std::sort(pe.begin(), pe.end(), [](const CsrKey& a, const CsrKey& b) {
            return a.col != b.col ? a.col < b.col : a.row < b.row;
        });
        pe.erase(std::unique(pe.begin(), pe.end()), pe.end());

        const int nnz = static_cast<int>(pe.size());
        Ap_.assign(static_cast<std::size_t>(n_) + 1, 0);
        Ai_.assign(static_cast<std::size_t>(nnz), 0);
        Ax_.assign(static_cast<std::size_t>(nnz), T(0));

        // Build Ap (column pointer array).
        for (const auto& e : pe) {
            ++Ap_[static_cast<std::size_t>(e.col) + 1];
        }
        for (int j = 0; j < n_; ++j) {
            Ap_[static_cast<std::size_t>(j) + 1] += Ap_[static_cast<std::size_t>(j)];
        }
        // Build Ai and fill slot_map_.
        slot_map_.reserve(static_cast<std::size_t>(nnz) * 2);
        std::vector<int> col_cursor = Ap_;
        for (const auto& e : pe) {
            const int idx = col_cursor[static_cast<std::size_t>(e.col)]++;
            Ai_[static_cast<std::size_t>(idx)] = e.row;
            slot_map_[{e.row, e.col}] = idx;
        }

        pending_entries_.clear();
        pending_entries_.shrink_to_fit();

        // Pre-allocate per-thread scratch (same size as Ax_).
        for (auto& scratch : thread_scratch_) {
            scratch.assign(static_cast<std::size_t>(nnz), T(0));
        }
        thread_dirty_.assign(thread_scratch_.size(), false);
        symbolic_done_ = true;
    }

    /// Resolve a (row, col) pair to its CsrSlot after finalizeSymbolic().
    /// Devices should call this once during elaboration and cache the result.
    CsrSlot resolveSlot(int row, int col) const {
        if (row < 0 || col < 0 || row >= n_ || col >= n_) return -1;
        auto it = slot_map_.find({row, col});
        if (it == slot_map_.end()) return -1;
        return it->second;
    }

    // -----------------------------------------------------------------------
    // Phase B — Numeric fill (called every Newton iteration)
    // -----------------------------------------------------------------------

    /// Zero all values — single cache-line sweep over Ax_.
    void zeroValues() {
        std::memset(Ax_.data(), 0, Ax_.size() * sizeof(T));
        for (std::size_t t = 0; t < thread_scratch_.size(); ++t) {
            if (thread_dirty_[t]) {
                std::memset(thread_scratch_[t].data(), 0,
                            thread_scratch_[t].size() * sizeof(T));
                thread_dirty_[t] = false;
            }
        }
    }

    /// Add value to (row, col) via a pre-resolved slot handle.
    /// Use this in the hot Newton path — zero overhead, direct index write.
    void addViaSlot(CsrSlot slot, T value) {
        if (slot < 0) return;
        const std::size_t idx = static_cast<std::size_t>(slot);
        if (omp_in_parallel()) {
            const int tid = omp_get_thread_num();
            thread_scratch_[static_cast<std::size_t>(tid)][idx] += value;
            thread_dirty_[static_cast<std::size_t>(tid)] = true;
        } else {
            Ax_[idx] += value;
        }
    }

    /// Add value to (row, col) by hash-map lookup. Slightly slower than
    /// addViaSlot but safe when the slot is not pre-cached (legacy devices).
    void add(int row, int col, T value) {
        addViaSlot(resolveSlot(row, col), value);
    }

    /// Flush per-thread scratch buffers into Ax_. Call after the parallel
    /// stamp region ends and before passing Ax_ to KLU.
    void flushThreadBuffers() {
        for (std::size_t t = 0; t < thread_scratch_.size(); ++t) {
            if (!thread_dirty_[t]) continue;
            const std::size_t sz = thread_scratch_[t].size();
            for (std::size_t i = 0; i < sz; ++i) {
                Ax_[i] += thread_scratch_[t][i];
                thread_scratch_[t][i] = T(0);
            }
            thread_dirty_[t] = false;
        }
    }

    // -----------------------------------------------------------------------
    // CCS accessors for KLU
    // -----------------------------------------------------------------------
    const std::vector<int>& getAp() const { return Ap_; }
    const std::vector<int>& getAi() const { return Ai_; }
    std::vector<T>& getAx() { return Ax_; }
    const std::vector<T>& getAx() const { return Ax_; }
    int getNnz() const { return static_cast<int>(Ax_.size()); }
    bool isSymbolicDone() const { return symbolic_done_; }

    // -----------------------------------------------------------------------
    // Compatibility helpers for old-style `getEntries()` consumers
    // -----------------------------------------------------------------------
    struct Entry { int row; int col; T value; };
    std::vector<Entry> getEntries() const {
        std::vector<Entry> out;
        out.reserve(Ax_.size());
        for (int j = 0; j < n_; ++j) {
            for (int k = Ap_[static_cast<std::size_t>(j)];
                 k < Ap_[static_cast<std::size_t>(j) + 1]; ++k) {
                out.push_back({Ai_[static_cast<std::size_t>(k)], j,
                               Ax_[static_cast<std::size_t>(k)]});
            }
        }
        return out;
    }

private:
    int n_;
    bool symbolic_done_;

    // Symbolic phase staging area.
    std::vector<CsrKey>                             pending_entries_;
    std::unordered_map<CsrKey, int, CsrKeyHash>     slot_map_;

    // CCS arrays — handed directly to KLU.
    std::vector<int>  Ap_; // column pointers [n+1]
    std::vector<int>  Ai_; // row indices     [nnz]
    std::vector<T>    Ax_; // values          [nnz]

    // Per-thread scratch for parallel stamping (avoids atomic conflicts on Ax_).
    std::vector<std::vector<T>>  thread_scratch_;
    std::vector<bool>            thread_dirty_;
};

using FixedCsrMatrixReal    = FixedCsrMatrix<double>;
using FixedCsrMatrixComplex = FixedCsrMatrix<std::complex<double>>;

} // namespace gspice

#endif // GSPICE_CSR_MATRIX_HPP
