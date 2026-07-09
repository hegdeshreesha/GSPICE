#ifndef GSPICE_KLU_SOLVER_HPP
#define GSPICE_KLU_SOLVER_HPP

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <omp.h>
#include "sparse_matrix.hpp"

#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
#include <klu.h>
#endif

namespace gspice {

template <typename T>
struct LinearSolveContext {
    int size = 0;
    std::size_t pattern_hash = 0;
    std::vector<std::vector<int>> row_columns;
    long long pattern_reuse_hits = 0;
    long long pattern_rebuilds = 0;
    std::string backend = "AUTO";
    std::string ordering = "AUTO";
    bool use_singletons = true;
#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
    klu_symbolic* external_symbolic = nullptr;
    std::size_t external_symbolic_hash = 0;
    int external_symbolic_size = 0;
#endif

    ~LinearSolveContext() {
        freeExternalSymbolic();
    }

    void clear() {
        freeExternalSymbolic();
        size = 0;
        pattern_hash = 0;
        row_columns.clear();
        pattern_reuse_hits = 0;
        pattern_rebuilds = 0;
        backend = "AUTO";
        ordering = "AUTO";
        use_singletons = true;
    }

    void freeExternalSymbolic() {
#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
        if (external_symbolic) {
            klu_common common;
            klu_defaults(&common);
            klu_free_symbolic(&external_symbolic, &common);
            external_symbolic = nullptr;
            external_symbolic_hash = 0;
            external_symbolic_size = 0;
        }
#endif
    }
};

template <typename T>
struct LinearSolverStats {
    long long solve_calls = 0;
    long long pattern_reuse_hits = 0;
    long long pattern_rebuilds = 0;
    long long parallel_pivot_scans = 0;
    long long parallel_elimination_passes = 0;
    long long singleton_eliminations = 0;
    long long singleton_core_solves = 0;
    long long singleton_zero_core_solves = 0;
    long long max_singleton_core_size = 0;
    long long external_klu_calls = 0;
    long long external_klu_symbolic_reuse = 0;
    long long external_klu_symbolic_rebuilds = 0;
    long long external_klu_failures = 0;
    double structure_seconds = 0.0;
    double numeric_fill_seconds = 0.0;
    double pivot_seconds = 0.0;
    double elimination_seconds = 0.0;
    double factor_seconds = 0.0;
    double backsolve_seconds = 0.0;
    double external_klu_seconds = 0.0;
    double total_seconds = 0.0;
};

/**
 * Internal sparse direct solver used as the current KLU-compatible path.
 * The interface is intentionally shaped so a real SuiteSparse/KLU backend can
 * replace the numeric core later without touching the simulation engine.
 */
template <typename T>
class KluSolver {
public:
    using Context = LinearSolveContext<T>;
    using Stats = LinearSolverStats<T>;

    static Vector<T> solve(const SparseMatrix<T>& matrix, const Vector<T>& b, Context* context = nullptr) {
        auto total_start = std::chrono::steady_clock::now();
        ++stats().solve_calls;

        const int n = matrix.getSize();
        const auto entries = matrix.getEntries();
        const auto rhs_values = b.snapshotData();

        if (shouldUseExternalKlu(context)) {
            Vector<T> external_solution;
            if (trySuiteSparseKluSolve(entries, n, rhs_values, context, external_solution)) {
                stats().total_seconds += secondsBetween(total_start, std::chrono::steady_clock::now());
                return external_solution;
            }
            if (externalKluExplicitlyRequested(context)) {
                throw std::runtime_error("SuiteSparse/KLU solver requested but unavailable or failed for this matrix");
            }
        }

        auto structure_start = std::chrono::steady_clock::now();
        const std::size_t signature = patternHash(entries, n);
        std::vector<Row> rows = instantiateRows(entries, n, signature, context);
        auto structure_end = std::chrono::steady_clock::now();

        auto fill_start = std::chrono::steady_clock::now();
        populateValues(entries, rows);
        auto fill_end = std::chrono::steady_clock::now();

        auto factor_start = std::chrono::steady_clock::now();
        std::vector<T> rhs = rhs_values;
        const bool use_singletons = !context || context->use_singletons;
        Vector<T> solution = solveRows(rows, rhs, use_singletons);
        auto factor_end = std::chrono::steady_clock::now();

        auto backsolve_start = std::chrono::steady_clock::now();
        auto backsolve_end = std::chrono::steady_clock::now();

        stats().structure_seconds += secondsBetween(structure_start, structure_end);
        stats().numeric_fill_seconds += secondsBetween(fill_start, fill_end);
        stats().factor_seconds += secondsBetween(factor_start, factor_end);
        stats().backsolve_seconds += secondsBetween(backsolve_start, backsolve_end);
        stats().total_seconds += secondsBetween(total_start, backsolve_end);
        return solution;
    }

    static void resetStats() {
        stats() = Stats{};
    }

    static Stats getStats() {
        return stats();
    }

private:
    struct RowEntry {
        int col = 0;
        T value = T(0);
    };

    using Row = std::vector<RowEntry>;

    static constexpr double kPivotTolerance = 1e-25;
    static constexpr double kPruneTolerance = 1e-30;
    static constexpr int kPivotParallelThreshold = 8192;
    static constexpr int kEliminationParallelThreshold = 256;
    static constexpr std::size_t kEliminationParallelWorkThreshold = 65536;

    static Stats& stats() {
        static Stats instance;
        return instance;
    }

    static std::string normalizeOption(std::string value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        value.erase(std::remove(value.begin(), value.end(), '_'), value.end());
        value.erase(std::remove(value.begin(), value.end(), '-'), value.end());
        return value;
    }

    static std::string backendName(const Context* context) {
        return normalizeOption(context ? context->backend : "AUTO");
    }

    static bool externalKluExplicitlyRequested(const Context* context) {
        const std::string backend = backendName(context);
        return backend == "KLU" || backend == "SUITESPARSE" || backend == "SUITESPARSEKLU";
    }

    static bool shouldUseExternalKlu(const Context* context) {
        const std::string backend = backendName(context);
        if (backend == "INTERNAL" || backend == "GSPICE" || backend == "NATIVE") return false;
        if (externalKluExplicitlyRequested(context)) return true;
#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
        return backend == "AUTO";
#else
        return false;
#endif
    }

    static void buildCcs(
        const std::vector<typename SparseMatrix<T>::Entry>& entries,
        int n,
        std::vector<int>& Ap,
        std::vector<int>& Ai,
        std::vector<T>& Ax) {
        Ap.assign(static_cast<size_t>(n) + 1, 0);
        Ai.assign(entries.size(), 0);
        Ax.assign(entries.size(), T(0));
        for (const auto& entry : entries) {
            if (entry.col >= 0 && entry.col < n) {
                ++Ap[static_cast<size_t>(entry.col) + 1];
            }
        }
        for (int col = 0; col < n; ++col) {
            Ap[static_cast<size_t>(col) + 1] += Ap[static_cast<size_t>(col)];
        }
        std::vector<int> next = Ap;
        for (const auto& entry : entries) {
            if (entry.col < 0 || entry.col >= n) continue;
            const int offset = next[static_cast<size_t>(entry.col)]++;
            Ai[static_cast<size_t>(offset)] = entry.row;
            Ax[static_cast<size_t>(offset)] = entry.value;
        }
    }

#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
    static void configureKluCommon(klu_common& common, const Context* context) {
        klu_defaults(&common);
        const std::string ordering = normalizeOption(context ? context->ordering : "AUTO");
        if (ordering == "AMD") {
            common.ordering = 0;
        } else if (ordering == "COLAMD") {
            common.ordering = 1;
        } else if (ordering == "NOBTF") {
            common.btf = 0;
        }
    }

    static klu_symbolic* getOrCreateKluSymbolic(
        const std::vector<int>& Ap,
        const std::vector<int>& Ai,
        int n,
        std::size_t signature,
        Context* context,
        klu_common& common) {
        const std::string ordering = normalizeOption(context ? context->ordering : "AUTO");
        const bool can_cache = context != nullptr;
        if (can_cache &&
            context->external_symbolic &&
            context->external_symbolic_size == n &&
            context->external_symbolic_hash == signature) {
            ++stats().external_klu_symbolic_reuse;
            return context->external_symbolic;
        }

        if (can_cache) {
            context->freeExternalSymbolic();
        }

        klu_symbolic* symbolic = nullptr;
        if (ordering == "NATURAL") {
            symbolic = klu_analyze_given(n, const_cast<int*>(Ap.data()), const_cast<int*>(Ai.data()), nullptr, nullptr, &common);
        } else {
            symbolic = klu_analyze(n, const_cast<int*>(Ap.data()), const_cast<int*>(Ai.data()), &common);
        }
        if (!symbolic) {
            ++stats().external_klu_failures;
            return nullptr;
        }

        ++stats().external_klu_symbolic_rebuilds;
        if (can_cache) {
            context->external_symbolic = symbolic;
            context->external_symbolic_hash = signature;
            context->external_symbolic_size = n;
        }
        return symbolic;
    }
#endif

    static bool trySuiteSparseKluSolve(
        const std::vector<typename SparseMatrix<T>::Entry>& entries,
        int n,
        const std::vector<T>& rhs_values,
        Context* context,
        Vector<T>& solution) {
#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
        const auto external_start = std::chrono::steady_clock::now();
        ++stats().external_klu_calls;

        std::vector<int> Ap;
        std::vector<int> Ai;
        std::vector<T> Ax;
        buildCcs(entries, n, Ap, Ai, Ax);

        klu_common common;
        configureKluCommon(common, context);
        const std::size_t signature = patternHash(entries, n);
        klu_symbolic* symbolic = getOrCreateKluSymbolic(Ap, Ai, n, signature, context, common);
        if (!symbolic) return false;

        bool ok = false;
        if constexpr (std::is_same_v<T, double>) {
            std::vector<double> rhs = rhs_values;
            klu_numeric* numeric = klu_factor(Ap.data(), Ai.data(), Ax.data(), symbolic, &common);
            if (numeric && klu_solve(symbolic, numeric, n, 1, rhs.data(), &common)) {
                solution = Vector<T>(n);
                for (int i = 0; i < n; ++i) solution[i] = rhs[static_cast<size_t>(i)];
                ok = true;
            }
            if (numeric) {
                klu_free_numeric(&numeric, &common);
            }
        } else if constexpr (std::is_same_v<T, std::complex<double>>) {
            std::vector<double> complexAx(entries.size() * 2, 0.0);
            for (size_t i = 0; i < Ax.size(); ++i) {
                complexAx[2 * i] = Ax[i].real();
                complexAx[2 * i + 1] = Ax[i].imag();
            }
            std::vector<double> rhs(static_cast<size_t>(n) * 2, 0.0);
            for (int i = 0; i < n; ++i) {
                rhs[2 * static_cast<size_t>(i)] = rhs_values[static_cast<size_t>(i)].real();
                rhs[2 * static_cast<size_t>(i) + 1] = rhs_values[static_cast<size_t>(i)].imag();
            }
            klu_numeric* numeric = klu_z_factor(Ap.data(), Ai.data(), complexAx.data(), symbolic, &common);
            if (numeric && klu_z_solve(symbolic, numeric, n, 1, rhs.data(), &common)) {
                solution = Vector<T>(n);
                for (int i = 0; i < n; ++i) {
                    solution[i] = T(rhs[2 * static_cast<size_t>(i)], rhs[2 * static_cast<size_t>(i) + 1]);
                }
                ok = true;
            }
            if (numeric) {
                klu_z_free_numeric(&numeric, &common);
            }
        }

        if (!ok) {
            ++stats().external_klu_failures;
        }
        stats().external_klu_seconds += secondsBetween(external_start, std::chrono::steady_clock::now());
        if (!context && symbolic) {
            klu_free_symbolic(&symbolic, &common);
        }
        return ok;
#else
        (void)entries;
        (void)n;
        (void)rhs_values;
        (void)context;
        (void)solution;
        if (externalKluExplicitlyRequested(context)) {
            ++stats().external_klu_failures;
        }
        return false;
#endif
    }

    static double magnitude(const T& value) {
        return std::abs(value);
    }

    template <typename TimePoint>
    static double secondsBetween(TimePoint start, TimePoint end) {
        return std::chrono::duration<double>(end - start).count();
    }

    static std::size_t patternHash(const std::vector<typename SparseMatrix<T>::Entry>& entries, int n) {
        std::size_t seed = static_cast<std::size_t>(n) * 1469598103934665603ull;
        for (const auto& entry : entries) {
            seed ^= static_cast<std::size_t>(entry.row + 1) * 1099511628211ull;
            seed ^= static_cast<std::size_t>(entry.col + 1) * 1469598103934665603ull;
            seed += 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    static std::vector<std::vector<int>> buildTemplateColumns(
        const std::vector<typename SparseMatrix<T>::Entry>& entries,
        int n) {
        std::vector<std::vector<int>> row_columns(static_cast<size_t>(n));
        for (const auto& entry : entries) {
            row_columns[static_cast<size_t>(entry.row)].push_back(entry.col);
        }
        return row_columns;
    }

    static std::vector<Row> instantiateTemplateRows(const std::vector<std::vector<int>>& row_columns) {
        std::vector<Row> rows(row_columns.size());
        for (size_t row = 0; row < row_columns.size(); ++row) {
            rows[row].reserve(row_columns[row].size());
            for (int col : row_columns[row]) {
                rows[row].push_back({col, T(0)});
            }
        }
        return rows;
    }

    static std::vector<Row> instantiateRows(
        const std::vector<typename SparseMatrix<T>::Entry>& entries,
        int n,
        std::size_t signature,
        Context* context) {
        if (context && context->size == n && context->pattern_hash == signature &&
            context->row_columns.size() == static_cast<size_t>(n)) {
            ++context->pattern_reuse_hits;
            ++stats().pattern_reuse_hits;
            return instantiateTemplateRows(context->row_columns);
        }

        auto row_columns = buildTemplateColumns(entries, n);
        if (context) {
            context->size = n;
            context->pattern_hash = signature;
            context->row_columns = row_columns;
            ++context->pattern_rebuilds;
        }
        ++stats().pattern_rebuilds;
        return instantiateTemplateRows(row_columns);
    }

    static void populateValues(
        const std::vector<typename SparseMatrix<T>::Entry>& entries,
        std::vector<Row>& rows) {
        std::vector<size_t> row_offsets(rows.size(), 0);
        for (const auto& entry : entries) {
            Row& row = rows[static_cast<size_t>(entry.row)];
            size_t& offset = row_offsets[static_cast<size_t>(entry.row)];
            if (offset < row.size() && row[offset].col == entry.col) {
                row[offset].value = entry.value;
                ++offset;
                continue;
            }
            auto it = std::lower_bound(row.begin(), row.end(), entry.col, [](const RowEntry& lhs, int rhs_col) {
                return lhs.col < rhs_col;
            });
            if (it == row.end() || it->col != entry.col) {
                row.insert(it, {entry.col, entry.value});
            } else {
                it->value = entry.value;
            }
        }
    }

    static typename Row::iterator findEntry(Row& row, int col) {
        return std::lower_bound(row.begin(), row.end(), col, [](const RowEntry& lhs, int rhs_col) {
            return lhs.col < rhs_col;
        });
    }

    static typename Row::const_iterator findEntry(const Row& row, int col) {
        return std::lower_bound(row.begin(), row.end(), col, [](const RowEntry& lhs, int rhs_col) {
            return lhs.col < rhs_col;
        });
    }

    static void addToRow(Row& row, int col, const T& delta) {
        if (magnitude(delta) <= kPruneTolerance) return;
        auto it = findEntry(row, col);
        if (it == row.end() || it->col != col) {
            row.insert(it, {col, delta});
            return;
        }
        it->value += delta;
        if (magnitude(it->value) <= kPruneTolerance) {
            row.erase(it);
        }
    }

    static void eraseColumn(Row& row, int col) {
        auto it = findEntry(row, col);
        if (it != row.end() && it->col == col) {
            row.erase(it);
        }
    }

    static T pivotValue(const Row& row, int col) {
        auto it = findEntry(row, col);
        if (it == row.end() || it->col != col) {
            return T(0);
        }
        return it->value;
    }

    static int findPivotRow(const std::vector<Row>& rows, int k) {
        const int n = static_cast<int>(rows.size());
        int pivot_row = -1;
        double pivot_mag = 0.0;
        const bool use_parallel =
            omp_get_max_threads() > 1 && (n - k) >= kPivotParallelThreshold;

        if (use_parallel) {
            const int slots = std::max(1, omp_get_max_threads());
            std::vector<int> local_rows(static_cast<size_t>(slots), -1);
            std::vector<double> local_mags(static_cast<size_t>(slots), 0.0);
            #pragma omp parallel
            {
                const int tid = omp_get_thread_num();
                int best_row = -1;
                double best_mag = 0.0;
                #pragma omp for nowait schedule(static)
                for (int r = k; r < n; ++r) {
                    const double mag = magnitude(pivotValue(rows[static_cast<size_t>(r)], k));
                    if (mag > best_mag || (mag == best_mag && best_row >= 0 && r < best_row)) {
                        best_mag = mag;
                        best_row = r;
                    } else if (best_row < 0 && mag > 0.0) {
                        best_mag = mag;
                        best_row = r;
                    }
                }
                local_rows[static_cast<size_t>(tid)] = best_row;
                local_mags[static_cast<size_t>(tid)] = best_mag;
            }

            for (size_t i = 0; i < local_rows.size(); ++i) {
                const int row = local_rows[i];
                const double mag = local_mags[i];
                if (row < 0) continue;
                if (mag > pivot_mag || (mag == pivot_mag && pivot_row >= 0 && row < pivot_row)) {
                    pivot_mag = mag;
                    pivot_row = row;
                } else if (pivot_row < 0) {
                    pivot_mag = mag;
                    pivot_row = row;
                }
            }
            if (pivot_row >= 0) {
                ++stats().parallel_pivot_scans;
            }
            return pivot_row;
        }

        for (int r = k; r < n; ++r) {
            const double mag = magnitude(pivotValue(rows[static_cast<size_t>(r)], k));
            if (mag > pivot_mag) {
                pivot_mag = mag;
                pivot_row = r;
            }
        }
        return pivot_row;
    }

    struct SingletonInfo {
        int row = -1;
        int col = -1;
        T value = T(0);
    };

    static bool findRowSingleton(
        const std::vector<Row>& rows,
        const std::vector<char>& active_rows,
        const std::vector<char>& active_cols,
        SingletonInfo& out) {
        for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
            if (!active_rows[static_cast<size_t>(r)]) continue;
            int count = 0;
            int col = -1;
            T value = T(0);
            for (const auto& entry : rows[static_cast<size_t>(r)]) {
                if (!active_cols[static_cast<size_t>(entry.col)] || magnitude(entry.value) <= kPruneTolerance) {
                    continue;
                }
                ++count;
                col = entry.col;
                value = entry.value;
                if (count > 1) break;
            }
            if (count == 1) {
                out = {r, col, value};
                return true;
            }
        }
        return false;
    }

    static Vector<T> solveRowsWithSingletons(std::vector<Row>& rows, std::vector<T>& rhs) {
        const int n = static_cast<int>(rows.size());
        std::vector<char> active_rows(static_cast<size_t>(n), 1);
        std::vector<char> active_cols(static_cast<size_t>(n), 1);
        std::vector<T> full_solution(static_cast<size_t>(n), T(0));
        int eliminated = 0;

        SingletonInfo singleton;
        while (findRowSingleton(rows, active_rows, active_cols, singleton)) {
            if (magnitude(singleton.value) <= kPivotTolerance) {
                break;
            }
            const T value = rhs[static_cast<size_t>(singleton.row)] / singleton.value;
            full_solution[static_cast<size_t>(singleton.col)] = value;
            active_rows[static_cast<size_t>(singleton.row)] = 0;
            active_cols[static_cast<size_t>(singleton.col)] = 0;
            ++eliminated;

            for (int r = 0; r < n; ++r) {
                if (!active_rows[static_cast<size_t>(r)]) continue;
                auto it = findEntry(rows[static_cast<size_t>(r)], singleton.col);
                if (it == rows[static_cast<size_t>(r)].end() || it->col != singleton.col) continue;
                rhs[static_cast<size_t>(r)] -= it->value * value;
                rows[static_cast<size_t>(r)].erase(it);
            }
        }

        stats().singleton_eliminations += eliminated;
        if (eliminated == 0) {
            factorize(rows, rhs);
            const auto backsolve_start = std::chrono::steady_clock::now();
            Vector<T> direct = backsolve(rows, rhs);
            const auto backsolve_end = std::chrono::steady_clock::now();
            stats().backsolve_seconds += secondsBetween(backsolve_start, backsolve_end);
            return direct;
        }

        std::vector<int> row_map;
        std::vector<int> col_map;
        row_map.reserve(static_cast<size_t>(n - eliminated));
        col_map.reserve(static_cast<size_t>(n - eliminated));
        for (int i = 0; i < n; ++i) {
            if (active_rows[static_cast<size_t>(i)]) row_map.push_back(i);
            if (active_cols[static_cast<size_t>(i)]) col_map.push_back(i);
        }

        if (row_map.size() != col_map.size()) {
            throw std::runtime_error("Singular matrix after singleton reduction");
        }

        const int core_size = static_cast<int>(row_map.size());
        stats().max_singleton_core_size = std::max(stats().max_singleton_core_size, static_cast<long long>(core_size));
        if (core_size == 0) {
            ++stats().singleton_zero_core_solves;
            Vector<T> result(n);
            for (int i = 0; i < n; ++i) result[i] = full_solution[static_cast<size_t>(i)];
            return result;
        }

        ++stats().singleton_core_solves;
        std::vector<int> col_to_core(static_cast<size_t>(n), -1);
        for (int i = 0; i < core_size; ++i) {
            col_to_core[static_cast<size_t>(col_map[static_cast<size_t>(i)])] = i;
        }

        std::vector<Row> core_rows(static_cast<size_t>(core_size));
        std::vector<T> core_rhs(static_cast<size_t>(core_size), T(0));
        for (int r = 0; r < core_size; ++r) {
            const int original_row = row_map[static_cast<size_t>(r)];
            core_rhs[static_cast<size_t>(r)] = rhs[static_cast<size_t>(original_row)];
            for (const auto& entry : rows[static_cast<size_t>(original_row)]) {
                const int mapped_col = entry.col >= 0 && entry.col < n
                    ? col_to_core[static_cast<size_t>(entry.col)]
                    : -1;
                if (mapped_col >= 0 && magnitude(entry.value) > kPruneTolerance) {
                    core_rows[static_cast<size_t>(r)].push_back({mapped_col, entry.value});
                }
            }
        }

        factorize(core_rows, core_rhs);
        const auto backsolve_start = std::chrono::steady_clock::now();
        Vector<T> core_solution = backsolve(core_rows, core_rhs);
        const auto backsolve_end = std::chrono::steady_clock::now();
        stats().backsolve_seconds += secondsBetween(backsolve_start, backsolve_end);

        for (int i = 0; i < core_size; ++i) {
            full_solution[static_cast<size_t>(col_map[static_cast<size_t>(i)])] = core_solution[i];
        }

        Vector<T> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = full_solution[static_cast<size_t>(i)];
        }
        return result;
    }

    static Vector<T> solveRows(std::vector<Row>& rows, std::vector<T>& rhs, bool use_singletons) {
        if (use_singletons) {
            return solveRowsWithSingletons(rows, rhs);
        }
        factorize(rows, rhs);
        const auto backsolve_start = std::chrono::steady_clock::now();
        Vector<T> result = backsolve(rows, rhs);
        const auto backsolve_end = std::chrono::steady_clock::now();
        stats().backsolve_seconds += secondsBetween(backsolve_start, backsolve_end);
        return result;
    }

    static void factorize(std::vector<Row>& rows, std::vector<T>& rhs) {
        const int n = static_cast<int>(rows.size());
        for (int k = 0; k < n; ++k) {
            const auto pivot_start = std::chrono::steady_clock::now();
            const int pivot_row = findPivotRow(rows, k);
            const double pivot_mag = pivot_row >= 0
                ? magnitude(pivotValue(rows[static_cast<size_t>(pivot_row)], k))
                : 0.0;
            const auto pivot_end = std::chrono::steady_clock::now();
            stats().pivot_seconds += secondsBetween(pivot_start, pivot_end);

            if (pivot_row < 0 || pivot_mag <= kPivotTolerance) {
                throw std::runtime_error("Singular matrix: zero pivot at row " + std::to_string(k));
            }

            if (pivot_row != k) {
                std::swap(rows[static_cast<size_t>(pivot_row)], rows[static_cast<size_t>(k)]);
                std::swap(rhs[static_cast<size_t>(pivot_row)], rhs[static_cast<size_t>(k)]);
            }

            const T pivot = pivotValue(rows[static_cast<size_t>(k)], k);
            if (magnitude(pivot) <= kPivotTolerance) {
                throw std::runtime_error("Singular matrix: zero pivot at row " + std::to_string(k));
            }

            const auto pivot_it = findEntry(rows[static_cast<size_t>(k)], k);
            const std::vector<RowEntry> pivot_tail(std::next(pivot_it), rows[static_cast<size_t>(k)].end());
            const T pivot_rhs = rhs[static_cast<size_t>(k)];
            const std::size_t elimination_work =
                static_cast<std::size_t>(std::max(0, n - (k + 1))) * pivot_tail.size();
            const bool use_parallel =
                omp_get_max_threads() > 1 &&
                (n - (k + 1)) >= kEliminationParallelThreshold &&
                elimination_work >= kEliminationParallelWorkThreshold;
            if (use_parallel) {
                ++stats().parallel_elimination_passes;
            }

            const auto elimination_start = std::chrono::steady_clock::now();
            #pragma omp parallel for schedule(static) if(use_parallel)
            for (int r = k + 1; r < n; ++r) {
                Row& target_row = rows[static_cast<size_t>(r)];
                auto target_it = findEntry(target_row, k);
                if (target_it == target_row.end() || target_it->col != k) {
                    continue;
                }

                const T factor = target_it->value / pivot;
                eraseColumn(target_row, k);
                rhs[static_cast<size_t>(r)] -= factor * pivot_rhs;

                for (const auto& pivot_entry : pivot_tail) {
                    addToRow(target_row, pivot_entry.col, -factor * pivot_entry.value);
                }
            }
            const auto elimination_end = std::chrono::steady_clock::now();
            stats().elimination_seconds += secondsBetween(elimination_start, elimination_end);
        }
    }

    static Vector<T> backsolve(const std::vector<Row>& rows, const std::vector<T>& rhs) {
        const int n = static_cast<int>(rows.size());
        std::vector<T> x(static_cast<size_t>(n), T(0));
        for (int i = n - 1; i >= 0; --i) {
            const T pivot = pivotValue(rows[static_cast<size_t>(i)], i);
            if (magnitude(pivot) <= kPivotTolerance) {
                throw std::runtime_error("Singular matrix during back substitution at row " + std::to_string(i));
            }
            T sum = rhs[static_cast<size_t>(i)];
            for (const auto& entry : rows[static_cast<size_t>(i)]) {
                if (entry.col > i) {
                    sum -= entry.value * x[static_cast<size_t>(entry.col)];
                }
            }
            x[static_cast<size_t>(i)] = sum / pivot;
        }

        Vector<T> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = x[static_cast<size_t>(i)];
        }
        return result;
    }
};

using KluSolverReal = KluSolver<double>;
using KluSolverComplex = KluSolver<std::complex<double>>;
using LinearSolveContextReal = LinearSolveContext<double>;
using LinearSolveContextComplex = LinearSolveContext<std::complex<double>>;
using LinearSolverStatsReal = LinearSolverStats<double>;
using LinearSolverStatsComplex = LinearSolverStats<std::complex<double>>;

} // namespace gspice

#endif // GSPICE_KLU_SOLVER_HPP
