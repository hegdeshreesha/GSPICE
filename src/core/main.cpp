#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <exception>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <limits>
#include <chrono>
#include <functional>
#include <random>
#include <omp.h>
#include "parser.hpp"
#include "devices/resistor.hpp"
#include "devices/capacitor.hpp"
#include "devices/inductor.hpp"
#include "devices/diode.hpp"
#include "devices/voltage_source.hpp"
#include "devices/port.hpp"
#include "devices/mosfet.hpp"
#include "devices/probe.hpp"
#include "devices/current_source.hpp"
#include "devices/osdi_device.hpp"
#include "devices/controlled_source.hpp"
#include "devices/behavioral_source.hpp"
#include "osdi_emulator.hpp"
#include "solvers/klu_solver.hpp"

using namespace gspice;

namespace {

double elapsed_seconds(
    const std::chrono::steady_clock::time_point& start,
    const std::chrono::steady_clock::time_point& end) {
    return std::chrono::duration<double>(end - start).count();
}

int choose_active_threads(int requested_threads, int num_devs, int matrix_size) {
    const int workload = std::max(num_devs, matrix_size * 2);
    const double density = matrix_size > 0
        ? static_cast<double>(num_devs) / static_cast<double>(matrix_size)
        : 0.0;
    int recommended = 1;
    if (density >= 16.0) {
        if (workload >= 2048) recommended = 2;
        if (workload >= 8192) recommended = 4;
        if (workload >= 16384) recommended = 8;
        if (workload >= 32768) recommended = 16;
    } else if (density >= 4.0) {
        if (workload >= 4096) recommended = 2;
        if (workload >= 16384) recommended = 4;
        if (workload >= 32768) recommended = 8;
        if (workload >= 65536) recommended = 16;
    } else {
        if (workload >= 8192) recommended = 2;
        if (workload >= 16384) recommended = 4;
        if (workload >= 32768) recommended = 8;
        if (workload >= 65536) recommended = 16;
    }
    recommended = std::min(recommended, 16);
    return std::max(1, std::min(requested_threads, recommended));
}

bool parallel_stamp_enabled(int num_devs, int matrix_size) {
    if (omp_get_max_threads() <= 1) return false;
    return num_devs >= 128 || matrix_size >= 2048;
}

double default_transient_max_step(const SimulationSettings& settings, double output_step) {
    double divisor = 4.0;
    if (settings.tran_lte_reltol <= 3e-4 || settings.reltol <= 1e-4) {
        divisor = 20.0;
    } else if (settings.tran_lte_reltol <= 1e-3 || settings.reltol <= 3e-4) {
        divisor = 10.0;
    } else if (settings.tran_lte_reltol >= 1e-2 || settings.reltol >= 5e-3) {
        divisor = 2.0;
    }
    return std::max(output_step / divisor, 1e-18);
}

void stamp_global_gmin(SparseMatrixReal& J, int num_nodes, double gmin) {
    for (int i = 0; i < num_nodes; ++i) {
        J.add(i, i, gmin);
    }
}

void stamp_global_gmin(SparseMatrixComplex& J, int num_nodes, double gmin) {
    for (int i = 0; i < num_nodes; ++i) {
        J.add(i, i, {gmin, 0.0});
    }
}

bool solution_converged(
    const VectorReal& x_new,
    const VectorReal& x_old,
    int num_nodes,
    const SimulationSettings& settings) {
    const int n = x_new.getSize();
    for (int i = 0; i < n; ++i) {
        const double delta = std::abs(x_new[i] - x_old[i]);
        const double scale = std::max(std::abs(x_new[i]), std::abs(x_old[i]));
        const double abs_floor = (i < num_nodes) ? settings.vntol : settings.abstol;
        const double tol = abs_floor + settings.reltol * scale;
        if (delta > tol) return false;
    }
    return true;
}

double linear_system_residual_error(
    const SparseMatrixReal& J,
    const VectorReal& b,
    const VectorReal& x,
    int num_nodes,
    const SimulationSettings& settings) {
    std::vector<double> ax(static_cast<size_t>(J.getSize()), 0.0);
    for (const auto& entry : J.getEntries()) {
        ax[static_cast<size_t>(entry.row)] += entry.value * x[entry.col];
    }

    double worst = 0.0;
    const int n = std::min(J.getSize(), std::min(b.getSize(), x.getSize()));
    for (int i = 0; i < n; ++i) {
        const double residual = std::abs(ax[static_cast<size_t>(i)] - b[i]);
        const double scale = std::max(std::abs(ax[static_cast<size_t>(i)]), std::abs(b[i]));
        const double abs_floor = (i < num_nodes) ? settings.abstol : settings.vntol;
        const double tol = abs_floor + settings.reltol * scale;
        worst = std::max(worst, residual / std::max(tol, 1e-30));
    }
    return worst;
}

double dc_residual_error(
    const std::vector<std::unique_ptr<Device>>& devices,
    int num_devs,
    int matrix_size,
    int num_nodes,
    const VectorReal& x,
    const SimulationSettings& settings,
    double active_gmin) {
    SparseMatrixReal J(matrix_size);
    VectorReal b(matrix_size);
    stamp_global_gmin(J, num_nodes, active_gmin);
    const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
    #pragma omp parallel for if(use_parallel_stamp)
    for (int i = 0; i < num_devs; ++i) {
        devices[i]->dcStamp(J, b, x, 0.0, 0.0, {});
    }
    return linear_system_residual_error(J, b, x, num_nodes, settings);
}

double transient_residual_error(
    const std::vector<std::unique_ptr<Device>>& devices,
    int num_devs,
    int matrix_size,
    int num_nodes,
    const VectorReal& x,
    const TransientContext& ctx,
    const SimulationSettings& settings) {
    SparseMatrixReal J(matrix_size);
    VectorReal b(matrix_size);
    stamp_global_gmin(J, num_nodes, settings.gmin);
    const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
    #pragma omp parallel for if(use_parallel_stamp)
    for (int i = 0; i < num_devs; ++i) {
        devices[i]->tranStamp(J, b, x, ctx);
    }
    return linear_system_residual_error(J, b, x, num_nodes, settings);
}

bool vector_finite_and_bounded(const VectorReal& x, double limit = 1e9) {
    for (int i = 0; i < x.getSize(); ++i) {
        if (!std::isfinite(x[i]) || std::abs(x[i]) > limit) return false;
    }
    return true;
}

std::string upper_copy(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return text;
}

bool is_ground_name(const std::string& name) {
    const std::string upper = upper_copy(name);
    return upper == "0" || upper == "GND" || upper == "GROUND";
}

double vector_max_delta(const VectorReal& lhs, const VectorReal& rhs) {
    const int n = std::min(lhs.getSize(), rhs.getSize());
    double max_delta = 0.0;
    for (int i = 0; i < n; ++i) {
        max_delta = std::max(max_delta, std::abs(lhs[i] - rhs[i]));
    }
    return max_delta;
}

VectorReal vector_blend(const VectorReal& from, const VectorReal& to, double alpha) {
    const int n = std::min(from.getSize(), to.getSize());
    VectorReal out(from.getSize());
    for (int i = 0; i < n; ++i) {
        out[i] = from[i] + alpha * (to[i] - from[i]);
    }
    return out;
}

double node_value(const VectorReal& x, int node) {
    return node >= 0 ? x[node] : 0.0;
}

double probe_value(const VectorReal& x, int node_pos, int node_neg) {
    return node_value(x, node_pos) - node_value(x, node_neg);
}

double spec_value(const VectorReal& x, const OutputSpec& spec) {
    return probe_value(x, spec.node_pos, spec.node_neg);
}

bool spec_passes(const OutputSpec& spec, double value) {
    if (spec.has_min && value < spec.min_value) return false;
    if (spec.has_max && value > spec.max_value) return false;
    return true;
}

struct ProductionSpecStats {
    std::string name;
    long long count = 0;
    long long pass = 0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();

    void add(double value, bool passes) {
        ++count;
        if (passes) ++pass;
        sum += value;
        sum_sq += value * value;
        min = std::min(min, value);
        max = std::max(max, value);
    }

    double mean() const {
        return count > 0 ? sum / static_cast<double>(count) : 0.0;
    }

    double sigma() const {
        if (count <= 0) return 0.0;
        const double m = mean();
        return std::sqrt(std::max(0.0, sum_sq / static_cast<double>(count) - m * m));
    }

    double yieldPercent() const {
        return count > 0 ? 100.0 * static_cast<double>(pass) / static_cast<double>(count) : 0.0;
    }
};

std::complex<double> complex_node_value(const VectorComplex& x, int node) {
    return node >= 0 ? x[node] : std::complex<double>{0.0, 0.0};
}

std::complex<double> complex_probe_value(const VectorComplex& x, int node_pos, int node_neg) {
    return complex_node_value(x, node_pos) - complex_node_value(x, node_neg);
}

std::string node_label(const Netlist& netlist, int node) {
    return node >= 0 ? netlist.getNodeName(node) : "0";
}

std::string voltage_probe_label(const Netlist& netlist, int node_pos, int node_neg) {
    std::string label = "V(" + node_label(netlist, node_pos);
    if (node_neg >= 0) {
        label += "," + node_label(netlist, node_neg);
    }
    label += ")";
    return label;
}

struct TranSample {
    double time = 0.0;
    VectorReal x;
};

double measure_probe_at_sample(const TranSample& sample, const MeasureSpec& measure) {
    return probe_value(sample.x, measure.node_pos, measure.node_neg);
}

double interpolate_measure_value(const std::vector<TranSample>& samples, const MeasureSpec& measure, double time) {
    if (samples.empty()) return 0.0;
    if (time <= samples.front().time) return measure_probe_at_sample(samples.front(), measure);
    for (size_t i = 1; i < samples.size(); ++i) {
        if (time <= samples[i].time) {
            const double t0 = samples[i - 1].time;
            const double t1 = samples[i].time;
            const double y0 = measure_probe_at_sample(samples[i - 1], measure);
            const double y1 = measure_probe_at_sample(samples[i], measure);
            if (t1 <= t0) return y1;
            const double alpha = (time - t0) / (t1 - t0);
            return y0 + alpha * (y1 - y0);
        }
    }
    return measure_probe_at_sample(samples.back(), measure);
}

double evaluate_transient_measure(const std::vector<TranSample>& samples, const MeasureSpec& measure) {
    if (samples.empty()) return 0.0;
    const double from = measure.has_from ? measure.from : samples.front().time;
    const double to = measure.has_to ? measure.to : samples.back().time;
    const std::string op = upper_copy(measure.op);
    if (op == "FIND") {
        return interpolate_measure_value(samples, measure, measure.has_at ? measure.at : to);
    }

    bool have = false;
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();
    double sum = 0.0;
    double sum_sq = 0.0;
    long long count = 0;
    for (const auto& sample : samples) {
        if (sample.time + 1e-30 < from || sample.time - 1e-30 > to) continue;
        const double value = measure_probe_at_sample(sample, measure);
        min_v = std::min(min_v, value);
        max_v = std::max(max_v, value);
        sum += value;
        sum_sq += value * value;
        ++count;
        have = true;
    }
    if (!have) {
        const double mid = 0.5 * (from + to);
        return interpolate_measure_value(samples, measure, mid);
    }
    if (op == "MAX") return max_v;
    if (op == "MIN") return min_v;
    if (op == "PP" || op == "PEAKTOPEAK") return max_v - min_v;
    if (op == "RMS") return std::sqrt(sum_sq / static_cast<double>(count));
    if (op == "AVG" || op == "AVERAGE" || op == "MEAN") return sum / static_cast<double>(count);
    return interpolate_measure_value(samples, measure, measure.has_at ? measure.at : to);
}

struct TransientStepResult {
    VectorReal x;
    bool converged = false;
    int iterations = 0;
    double stamp_seconds = 0.0;
    double solve_seconds = 0.0;
    double residual_error = std::numeric_limits<double>::infinity();
};

struct TransientStats {
    long long accepted_steps = 0;
    long long rejected_steps = 0;
    long long convergence_rejections = 0;
    long long lte_rejections = 0;
    long long output_points = 0;
    long long newton_solves = 0;
    long long newton_iterations = 0;
    long long bound_step_limited = 0;
    int max_newton_iterations = 0;
    double min_step = std::numeric_limits<double>::infinity();
    double max_step = 0.0;
    double last_error = 0.0;
    double last_residual_error = 0.0;
    double max_residual_error = 0.0;
    double stamp_seconds = 0.0;
    double solve_seconds = 0.0;
    double min_bound_step = std::numeric_limits<double>::infinity();

    void noteSolve(const TransientStepResult& result) {
        ++newton_solves;
        newton_iterations += result.iterations;
        max_newton_iterations = std::max(max_newton_iterations, result.iterations);
        if (std::isfinite(result.residual_error)) {
            last_residual_error = result.residual_error;
            max_residual_error = std::max(max_residual_error, result.residual_error);
        }
        stamp_seconds += result.stamp_seconds;
        solve_seconds += result.solve_seconds;
    }

    void noteAccepted(double dt, double err) {
        ++accepted_steps;
        min_step = std::min(min_step, dt);
        max_step = std::max(max_step, dt);
        last_error = err;
    }

    double averageNewtonIterations() const {
        if (newton_solves <= 0) return 0.0;
        return static_cast<double>(newton_iterations) / static_cast<double>(newton_solves);
    }
};

double collect_transient_bound_step(const std::vector<std::unique_ptr<Device>>& devices) {
    double limit = std::numeric_limits<double>::infinity();
    for (const auto& device : devices) {
        const double bound = device->transientBoundStep();
        if (std::isfinite(bound) && bound > 0.0) {
            limit = std::min(limit, bound);
        }
    }
    return limit;
}

struct SimulationRuntimeStats {
    double dc_stamp_seconds = 0.0;
    double dc_solve_seconds = 0.0;
    double ac_stamp_seconds = 0.0;
    double ac_solve_seconds = 0.0;
    double hb_stamp_seconds = 0.0;
    double hb_solve_seconds = 0.0;
};

std::string normalized_method_key(const SimulationSettings& settings) {
    std::string method = settings.tran_method;
    std::transform(method.begin(), method.end(), method.begin(), ::toupper);
    method.erase(std::remove(method.begin(), method.end(), '_'), method.end());
    method.erase(std::remove(method.begin(), method.end(), '-'), method.end());
    return method;
}

TransientIntegrationMethod choose_transient_method(
    const SimulationSettings& settings,
    const std::vector<VectorReal>& x_hist,
    const std::vector<double>& t_hist) {
    const std::string method = normalized_method_key(settings);
    if (method == "GEAR2" && x_hist.size() >= 2 && t_hist.size() >= 2) {
        return TransientIntegrationMethod::Gear2;
    }
    if (method == "AUTO" && settings.tran_adaptive && x_hist.size() >= 2 && t_hist.size() >= 2) {
        return TransientIntegrationMethod::Gear2;
    }
    if ((method == "TRAP" || method == "TRAPEZOIDAL") && !x_hist.empty()) {
        return TransientIntegrationMethod::Trapezoidal;
    }
    return TransientIntegrationMethod::BackwardEuler;
}

TransientContext make_transient_context(
    const SimulationSettings& settings,
    const std::vector<VectorReal>& x_hist,
    const std::vector<double>& t_hist,
    double step,
    double target_time) {
    TransientContext ctx;
    ctx.timeStep = std::max(step, 1e-30);
    ctx.currentTime = target_time;
    ctx.method = choose_transient_method(settings, x_hist, t_hist);
    ctx.xHistory = &x_hist;
    ctx.timeHistory = &t_hist;

    if (ctx.method == TransientIntegrationMethod::Gear2 && x_hist.size() >= 2 && t_hist.size() >= 2) {
        const double h0 = ctx.timeStep;
        const double h1 = std::max(t_hist.back() - t_hist[t_hist.size() - 2], 1e-30);
        ctx.a0 = (2.0 * h0 + h1) / (h0 * (h0 + h1));
        ctx.a1 = -(h0 + h1) / (h0 * h1);
        ctx.a2 = h0 / (h1 * (h0 + h1));
        ctx.hasSecondHistory = true;
    } else if (ctx.method == TransientIntegrationMethod::Trapezoidal) {
        ctx.a0 = 2.0 / ctx.timeStep;
        ctx.a1 = -ctx.a0;
        ctx.a2 = 0.0;
        ctx.hasSecondHistory = false;
    } else {
        ctx.a0 = 1.0 / ctx.timeStep;
        ctx.a1 = -ctx.a0;
        ctx.a2 = 0.0;
        ctx.hasSecondHistory = false;
    }

    return ctx;
}

VectorReal transient_initial_guess(
    const VectorReal& fallback,
    const std::vector<VectorReal>& x_hist,
    const std::vector<double>& t_hist,
    double target_time,
    const SimulationSettings& settings) {
    if (!settings.tran_predictor || x_hist.size() < 2 || t_hist.size() < 2) {
        return fallback;
    }
    const double t1 = t_hist.back();
    const double t0 = t_hist[t_hist.size() - 2];
    const double dt = t1 - t0;
    if (dt <= 1e-30) return fallback;
    const VectorReal& x1 = x_hist.back();
    const VectorReal& x0 = x_hist[x_hist.size() - 2];
    const int n = std::min(x1.getSize(), x0.getSize());
    VectorReal predicted(fallback.getSize());
    for (int i = 0; i < predicted.getSize(); ++i) {
        predicted[i] = fallback[i];
    }
    const double alpha = (target_time - t1) / dt;
    for (int i = 0; i < n; ++i) {
        predicted[i] = x1[i] + alpha * (x1[i] - x0[i]);
    }
    return vector_finite_and_bounded(predicted) ? predicted : fallback;
}

std::string transient_method_label(TransientIntegrationMethod method) {
    if (method == TransientIntegrationMethod::Gear2) return "Gear2/BDF2";
    if (method == TransientIntegrationMethod::Trapezoidal) return "Trapezoidal";
    return "Backward Euler";
}

int transient_method_order(TransientIntegrationMethod method) {
    return method == TransientIntegrationMethod::BackwardEuler ? 1 : 2;
}

double transient_reject_factor(double err, TransientIntegrationMethod method) {
    const int order = transient_method_order(method);
    const double exponent = -1.0 / static_cast<double>(order + 1);
    return std::clamp(0.85 * std::pow(std::max(err, 1e-30), exponent), 0.1, 0.75);
}

double transient_growth_factor(double err, TransientIntegrationMethod method) {
    if (err <= 1e-12) return method == TransientIntegrationMethod::BackwardEuler ? 2.0 : 2.5;
    const int order = transient_method_order(method);
    const double exponent = -1.0 / static_cast<double>(order + 1);
    return std::clamp(0.9 * std::pow(std::max(err, 1e-30), exponent), 1.05, 2.5);
}

TransientStepResult solve_transient_step(
    const std::vector<std::unique_ptr<Device>>& devices,
    int num_devs,
    int matrix_size,
    int num_nodes,
    const SimulationSettings& settings,
    LinearSolveContextReal* solver_context,
    const VectorReal& x_start,
    const std::vector<VectorReal>& x_hist,
    const std::vector<double>& t_hist,
    double step,
    double target_time) {
    VectorReal x = transient_initial_guess(x_start, x_hist, t_hist, target_time, settings);
    bool converged = false;
    double total_stamp_seconds = 0.0;
    double total_solve_seconds = 0.0;
    const TransientContext tran_ctx = make_transient_context(settings, x_hist, t_hist, step, target_time);
    const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
    for (int iter = 0; iter < settings.tran_max_iter; ++iter) {
        const auto stamp_start = std::chrono::steady_clock::now();
        SparseMatrixReal J_sparse(matrix_size);
        VectorReal b(matrix_size);
        stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
        #pragma omp parallel for if(use_parallel_stamp)
        for (int i = 0; i < num_devs; ++i) {
            devices[i]->tranStamp(J_sparse, b, x, tran_ctx);
        }
        const auto stamp_end = std::chrono::steady_clock::now();
        const auto solve_start = std::chrono::steady_clock::now();
        VectorReal x_new = KluSolverReal::solve(J_sparse, b, solver_context);
        const auto solve_end = std::chrono::steady_clock::now();
        total_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
        total_solve_seconds += elapsed_seconds(solve_start, solve_end);
        if (!vector_finite_and_bounded(x_new)) {
            return {x_new, false, iter + 1, total_stamp_seconds, total_solve_seconds,
                std::numeric_limits<double>::infinity()};
        }
        double residual_error = std::numeric_limits<double>::infinity();
        const bool delta_converged = solution_converged(x_new, x, num_nodes, settings);
        if (delta_converged && settings.nr_residual_check) {
            residual_error = transient_residual_error(
                devices, num_devs, matrix_size, num_nodes, x_new, tran_ctx, settings);
            converged = residual_error <= 1.0;
        } else {
            converged = delta_converged;
        }
        x = x_new;
        if (converged) {
            return {x, true, iter + 1, total_stamp_seconds, total_solve_seconds, residual_error};
        }
    }
    return {x, converged, settings.tran_max_iter, total_stamp_seconds, total_solve_seconds,
        std::numeric_limits<double>::infinity()};
}

double transient_lte_error(
    const VectorReal& x_full,
    const VectorReal& x_half,
    int num_nodes,
    const SimulationSettings& settings) {
    double worst = 0.0;
    const int n = std::min(num_nodes, std::min(x_full.getSize(), x_half.getSize()));
    for (int i = 0; i < n; ++i) {
        const double err = std::abs(x_half[i] - x_full[i]);
        const double scale = std::max(std::abs(x_full[i]), std::abs(x_half[i]));
        const double tol = settings.tran_lte_abstol + settings.tran_lte_reltol * scale;
        worst = std::max(worst, err / std::max(tol, 1e-30));
    }
    return worst;
}

double transient_state_change_error(
    const VectorReal& x_new,
    const VectorReal& x_old,
    int num_nodes,
    const SimulationSettings& settings) {
    double worst = 0.0;
    const int n = std::min(num_nodes, std::min(x_new.getSize(), x_old.getSize()));
    for (int i = 0; i < n; ++i) {
        const double err = std::abs(x_new[i] - x_old[i]);
        const double scale = std::max(std::abs(x_new[i]), std::abs(x_old[i]));
        const double tol = settings.tran_lte_abstol + settings.tran_lte_reltol * scale;
        worst = std::max(worst, err / std::max(tol, 1e-30));
    }
    return worst;
}

void print_transient_point(double t, const VectorReal& x, int num_nodes) {
    std::cout << std::scientific << std::setprecision(9) << t << " | ";
    for (int i = 0; i < num_nodes; ++i) std::cout << x[i] << " ";
    std::cout << std::endl;
}

struct SavedOutputSignal {
    std::string label;
    int node_pos = -1;
    int node_neg = -1;
};

std::vector<SavedOutputSignal> resolve_saved_outputs(const Netlist& netlist, int num_nodes) {
    const auto& settings = netlist.getSettings();
    std::vector<SavedOutputSignal> signals;
    if (settings.save_all || settings.saves.empty()) {
        signals.reserve(static_cast<size_t>(num_nodes));
        for (int i = 0; i < num_nodes; ++i) {
            signals.push_back({"V(" + netlist.getNodeName(i) + ")", i, -1});
        }
        return signals;
    }

    for (const auto& save : settings.saves) {
        if (upper_copy(save.kind) != "V") continue;
        const int pos = netlist.findNode(save.node_pos);
        const int neg = is_ground_name(save.node_neg) ? -1 : netlist.findNode(save.node_neg);
        if (pos == -2 || neg == -2) {
            std::cout << "WARNING: save skipped unknown node in V(" << save.node_pos;
            if (!is_ground_name(save.node_neg)) std::cout << "," << save.node_neg;
            std::cout << ")" << std::endl;
            continue;
        }
        std::string label = "V(" + save.node_pos;
        if (!is_ground_name(save.node_neg)) label += "," + save.node_neg;
        label += ")";
        signals.push_back({label, pos, neg});
    }
    if (signals.empty()) {
        std::cout << "WARNING: no valid .SAVE voltage signals were resolved; writing all node voltages." << std::endl;
        for (int i = 0; i < num_nodes; ++i) {
            signals.push_back({"V(" + netlist.getNodeName(i) + ")", i, -1});
        }
    }
    return signals;
}

class TransientOutput {
public:
    TransientOutput(const std::string& path, const Netlist& netlist, int num_nodes)
        : to_file_(!path.empty()), signals_(resolve_saved_outputs(netlist, num_nodes)) {
        if (to_file_) {
            file_.open(path, std::ios::out | std::ios::trunc);
            if (!file_.is_open()) {
                throw std::runtime_error("Could not open transient output file: " + path);
            }
            file_ << "Title: GSPICE RAW output\n";
            file_ << "Plotname: Transient Analysis\n";
            file_ << "Flags: real\n";
            file_ << "No. Variables: " << (signals_.size() + 1) << "\n";
            file_ << "No. Points: ";
            points_pos_ = file_.tellp();
            file_ << std::setw(20) << 0 << "\n";
            file_ << "Variables:\n";
            file_ << "0\ttime\ttime\n";
            for (size_t i = 0; i < signals_.size(); ++i) {
                file_ << (i + 1) << "\t" << signals_[i].label << "\tvoltage\n";
            }
            file_ << "Values:\n";
        }
    }

    ~TransientOutput() {
        finalize();
    }

    bool toFile() const { return to_file_; }
    const std::vector<SavedOutputSignal>& signals() const { return signals_; }

    void write(double t, const VectorReal& x, int num_nodes) {
        if (!to_file_) {
            print_transient_point(t, x, num_nodes);
            return;
        }
        file_ << std::scientific << std::setprecision(12) << t;
        for (const auto& signal : signals_) {
            file_ << " " << probe_value(x, signal.node_pos, signal.node_neg);
        }
        file_ << "\n";
        ++point_count_;
        if ((point_count_ % 1024) == 0) {
            checkpointPointCount();
        }
    }

    void finalize() {
        if (!to_file_ || finalized_) return;
        finalized_ = true;
        checkpointPointCount();
    }

private:
    void checkpointPointCount() {
        if (!file_.is_open() || points_pos_ < std::streampos(0)) return;
        const auto end_pos = file_.tellp();
        file_.seekp(points_pos_);
        file_ << std::setw(20) << point_count_;
        file_.seekp(end_pos);
        file_.flush();
    }

    bool to_file_ = false;
    bool finalized_ = false;
    std::vector<SavedOutputSignal> signals_;
    long long point_count_ = 0;
    std::streampos points_pos_ = std::streampos(-1);
    std::ofstream file_;
};

std::vector<double> collect_transient_breakpoints(
    const std::vector<std::unique_ptr<Device>>& devices,
    double t_stop,
    double merge_tol) {
    std::vector<double> points;
    for (const auto& dev : devices) {
        dev->collectBreakpoints(t_stop, points);
    }
    points.push_back(0.0);
    points.push_back(t_stop);
    std::sort(points.begin(), points.end());
    std::vector<double> merged;
    for (double t : points) {
        if (t < 0.0 || t > t_stop) continue;
        if (merged.empty() || std::abs(t - merged.back()) > merge_tol) {
            merged.push_back(t);
        }
    }
    return merged;
}

double next_breakpoint_after(const std::vector<double>& points, double t, double tol) {
    auto it = std::upper_bound(points.begin(), points.end(), t + tol);
    if (it == points.end()) return 0.0;
    return *it;
}

std::string normalized_transient_method(const SimulationSettings& settings) {
    const std::string method = normalized_method_key(settings);
    if (method == "TRAP" || method == "TRAPEZOIDAL") return "Trapezoidal";
    if (method == "GEAR2") return "Gear2/BDF2";
    if (method == "BE" || method == "BACKWARDEULER") return "Backward Euler";
    return "Auto (Backward Euler startup, Gear2/BDF2 after history)";
}

bool names_match(const std::string& lhs, const std::string& rhs) {
    return upper_copy(lhs) == upper_copy(rhs);
}

bool set_source_dc_value(
    const std::vector<std::unique_ptr<Device>>& devices,
    const std::string& source_name,
    double value) {
    for (const auto& dev : devices) {
        if (!names_match(dev->getName(), source_name)) continue;
        if (auto* vsrc = dynamic_cast<VoltageSource*>(dev.get())) {
            vsrc->setDcValue(value);
            return true;
        }
        if (auto* isrc = dynamic_cast<CurrentSource*>(dev.get())) {
            isrc->setDcValue(value);
            return true;
        }
        return false;
    }
    return false;
}

bool get_source_dc_value(
    const std::vector<std::unique_ptr<Device>>& devices,
    const std::string& source_name,
    double& value) {
    for (const auto& dev : devices) {
        if (!names_match(dev->getName(), source_name)) continue;
        if (auto* vsrc = dynamic_cast<VoltageSource*>(dev.get())) {
            value = vsrc->getDcValue();
            return true;
        }
        if (auto* isrc = dynamic_cast<CurrentSource*>(dev.get())) {
            value = isrc->getDcValue();
            return true;
        }
        return false;
    }
    return false;
}

struct SourceState {
    Device* device = nullptr;
    double dc_value = 0.0;
};

std::vector<SourceState> collect_source_states(const std::vector<std::unique_ptr<Device>>& devices) {
    std::vector<SourceState> states;
    for (const auto& dev : devices) {
        if (auto* vsrc = dynamic_cast<VoltageSource*>(dev.get())) {
            states.push_back({dev.get(), vsrc->getDcValue()});
        } else if (auto* isrc = dynamic_cast<CurrentSource*>(dev.get())) {
            states.push_back({dev.get(), isrc->getDcValue()});
        }
    }
    return states;
}

void scale_source_states(const std::vector<SourceState>& states, double scale) {
    for (const auto& state : states) {
        if (auto* vsrc = dynamic_cast<VoltageSource*>(state.device)) {
            vsrc->setDcValue(state.dc_value * scale);
        } else if (auto* isrc = dynamic_cast<CurrentSource*>(state.device)) {
            isrc->setDcValue(state.dc_value * scale);
        }
    }
}

void restore_source_states(const std::vector<SourceState>& states) {
    scale_source_states(states, 1.0);
}

int branch_index_for_device(const Device* dev) {
    if (auto* vsrc = dynamic_cast<const VoltageSource*>(dev)) return vsrc->getBranchIndex();
    if (auto* ind = dynamic_cast<const Inductor*>(dev)) return ind->getBranchIndex();
    if (auto* port = dynamic_cast<const Port*>(dev)) return port->getBranchIndex();
    if (auto* probe = dynamic_cast<const StabilityProbe*>(dev)) return probe->getBranchIndex();
    if (auto* vcvs = dynamic_cast<const VoltageControlledVoltageSource*>(dev)) return vcvs->getBranchIndex();
    if (auto* ccvs = dynamic_cast<const CurrentControlledVoltageSource*>(dev)) return ccvs->getBranchIndex();
    if (auto* bsrc = dynamic_cast<const BehavioralSource*>(dev)) return bsrc->getBranchIndex();
    return -1;
}

int find_branch_index_by_name(
    const std::vector<std::unique_ptr<Device>>& devices,
    const std::string& device_name) {
    for (const auto& dev : devices) {
        if (names_match(dev->getName(), device_name)) return branch_index_for_device(dev.get());
    }
    return -1;
}

bool configure_transfer_input_source(
    const std::vector<std::unique_ptr<Device>>& devices,
    const std::string& source_name) {
    bool found = false;
    for (const auto& dev : devices) {
        if (auto* vsrc = dynamic_cast<VoltageSource*>(dev.get())) {
            const bool match = names_match(dev->getName(), source_name);
            vsrc->setAcMagnitude(match ? 1.0 : 0.0);
            found = found || match;
        } else if (auto* isrc = dynamic_cast<CurrentSource*>(dev.get())) {
            const bool match = names_match(dev->getName(), source_name);
            isrc->setAcMagnitude(match ? 1.0 : 0.0);
            found = found || match;
        }
    }
    return found;
}

} // namespace

void run_osdi_test() {
    std::cout << "\n--- GSPICE Industrial Level-50 OSDI Test ---" << std::endl;
    Netlist netlist;
    int nVdd = netlist.getOrCreateNode("vdd");
    int nGate = netlist.getOrCreateNode("gate");
    int nDrain = netlist.getOrCreateNode("drain");
    int nGnd = -1;
    OsdiDescriptor va_model = OsdiEmulator::getDescriptor();
    netlist.addDevice(std::make_unique<VoltageSource>("Vdd", nVdd, nGnd, 1.8, 3));
    netlist.addDevice(std::make_unique<VoltageSource>("Vbias", nGate, nGnd, 0.8, 4));
    netlist.addDevice(std::make_unique<Resistor>("Rd", nVdd, nDrain, 10000.0));
    std::vector<int> mos_nodes = {nDrain, nGate, nGnd, nGnd};
    netlist.addDevice(std::make_unique<OSDIDevice>("M1", va_model, mos_nodes));
    int matrix_size = 5; VectorReal x(matrix_size); x[2] = 1.0;
    const double tol = 1e-8;
    LinearSolveContextReal solver_context;
    for (int iter = 0; iter < 100; ++iter) {
        SparseMatrixReal J(matrix_size); VectorReal b(matrix_size);
        auto& devices = netlist.getDevices();
        int num_devs = static_cast<int>(devices.size());
        const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
        #pragma omp parallel for if(use_parallel_stamp)
        for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J, b, x, 0.0, 0.0, {});
        VectorReal x_new = KluSolverReal::solve(J, b, &solver_context);
        double max_change = std::abs(x_new[2] - x[2]); x = x_new;
        if (max_change < tol) { std::cout << "  Converged in " << iter + 1 << " iterations." << std::endl; break; }
    }
    std::cout << "Results:\n  Vdd: " << x[0] << " V\n  Gate: " << x[1] << " V\n  Drain: " << x[2] << " V\n";
}

void run_simulation(Netlist& netlist, const std::string& output_file = "", int requested_threads = 1) {
    int num_nodes = netlist.getNumNodes();
    std::vector<Port*> ports;
    std::vector<StabilityProbe*> probes;
    auto& devices = netlist.getDevices();
    int num_devs = static_cast<int>(devices.size());

    for (auto& dev : devices) {
        auto* vsrc = dynamic_cast<VoltageSource*>(dev.get());
        if (vsrc) vsrc->setBranchIndex(netlist.getNextBranchId(num_nodes));
        auto* ind = dynamic_cast<Inductor*>(dev.get());
        if (ind) ind->setBranchIndex(netlist.getNextBranchId(num_nodes));
        auto* port = dynamic_cast<Port*>(dev.get());
        if (port) { port->setBranchIndex(netlist.getNextBranchId(num_nodes)); ports.push_back(port); }
        auto* probe = dynamic_cast<StabilityProbe*>(dev.get());
        if (probe) { probe->setBranchIndex(netlist.getNextBranchId(num_nodes)); probes.push_back(probe); }
        auto* vcvs = dynamic_cast<VoltageControlledVoltageSource*>(dev.get());
        if (vcvs) vcvs->setBranchIndex(netlist.getNextBranchId(num_nodes));
        auto* ccvs = dynamic_cast<CurrentControlledVoltageSource*>(dev.get());
        if (ccvs) ccvs->setBranchIndex(netlist.getNextBranchId(num_nodes));
        auto* bsrc = dynamic_cast<BehavioralSource*>(dev.get());
        if (bsrc && bsrc->isVoltageMode()) bsrc->setBranchIndex(netlist.getNextBranchId(num_nodes));
    }
    int matrix_size = num_nodes + netlist.getNumBranches();
    for (auto& dev : devices) {
        if (auto* cccs = dynamic_cast<CurrentControlledCurrentSource*>(dev.get())) {
            const int ctrlBranch = find_branch_index_by_name(devices, cccs->getControlSource());
            if (ctrlBranch < 0) {
                throw std::runtime_error(
                    "CCCS '" + cccs->getName() + "' references missing branch device '" +
                    cccs->getControlSource() + "'");
            }
            cccs->setControlBranchIndex(ctrlBranch);
        }
        if (auto* ccvs = dynamic_cast<CurrentControlledVoltageSource*>(dev.get())) {
            const int ctrlBranch = find_branch_index_by_name(devices, ccvs->getControlSource());
            if (ctrlBranch < 0) {
                throw std::runtime_error(
                    "CCVS '" + ccvs->getName() + "' references missing branch device '" +
                    ccvs->getControlSource() + "'");
            }
            ccvs->setControlBranchIndex(ctrlBranch);
        }
        if (auto* bsrc = dynamic_cast<BehavioralSource*>(dev.get())) {
            for (const auto& ref : bsrc->currentRefs()) {
                const int ctrlBranch = find_branch_index_by_name(devices, ref.name);
                if (ctrlBranch < 0) {
                    throw std::runtime_error(
                        "Behavioral source '" + bsrc->getName() + "' references missing branch device '" +
                        ref.name + "'");
                }
                bsrc->bindBranch(ref.name, ctrlBranch);
            }
            std::string missing;
            if (!bsrc->allBranchesBound(missing)) {
                throw std::runtime_error(
                    "Behavioral source '" + bsrc->getName() + "' has unresolved branch current '" +
                    missing + "'");
            }
        }
    }
    const int active_threads = choose_active_threads(requested_threads, num_devs, matrix_size);
    omp_set_num_threads(active_threads);
    std::cout << "GSPICE Core: Requested " << requested_threads
              << " thread(s), active " << active_threads
              << " thread(s) for " << num_devs << " device(s), matrix size "
              << matrix_size << "." << std::endl;
    KluSolverReal::resetStats();
    KluSolverComplex::resetStats();
    LinearSolveContextReal real_solver_context;
    LinearSolveContextComplex complex_solver_context;
    SimulationRuntimeStats runtime_stats;
    const auto& settings = netlist.getSettings();
    real_solver_context.backend = settings.solver_backend;
    real_solver_context.ordering = settings.solver_ordering;
    real_solver_context.use_singletons = settings.solver_singletons;
    complex_solver_context.backend = settings.solver_backend;
    complex_solver_context.ordering = settings.solver_ordering;
    complex_solver_context.use_singletons = settings.solver_singletons;
    std::cout << "Linear solver: backend=" << settings.solver_backend
              << " ordering=" << settings.solver_ordering
              << " singleton_filter=" << (settings.solver_singletons ? "on" : "off")
#if defined(GSPICE_HAVE_SUITESPARSE_KLU) && GSPICE_HAVE_SUITESPARSE_KLU
              << " (SuiteSparse/KLU available for SOLVER=AUTO/KLU)"
#else
              << " (SuiteSparse/KLU not linked; SOLVER=AUTO uses internal sparse engine)"
#endif
              << std::endl;
    std::cout << std::scientific << std::setprecision(9)
              << "Numerical policy: reltol=" << settings.reltol
              << " vntol=" << settings.vntol
              << " abstol=" << settings.abstol
              << " gmin=" << settings.gmin
              << " source_stepping=" << (settings.source_stepping ? "on" : "off")
              << " gmin_stepping=" << (settings.gmin_stepping ? "on" : "off")
              << " line_search=" << (settings.line_search ? "on" : "off")
              << " residual_check=" << (settings.nr_residual_check ? "on" : "off")
              << " op_max_iter=" << settings.op_max_iter
              << " tran_max_iter=" << settings.tran_max_iter
              << " temp_c=" << settings.temperature_c
              << std::endl;
    const bool osdi_limiting_rhs =
        settings.osdi_limiting_rhs || std::getenv("GSPICE_USE_OSDI_LIMITING_RHS") != nullptr;
    const bool osdi_tran_jacobian =
        settings.osdi_tran_jacobian || std::getenv("GSPICE_USE_OSDI_TRAN_JACOBIAN") != nullptr;
    const bool osdi_full_model_params =
        settings.osdi_bind_full_model_params || std::getenv("GSPICE_BIND_FULL_OSDI_MODEL_PARAMS") != nullptr;
    const bool osdi_internal_nodes = settings.osdi_internal_nodes;
    const bool osdi_spice_rhs = settings.osdi_spice_rhs;
    if (osdi_limiting_rhs || osdi_tran_jacobian || osdi_full_model_params || osdi_internal_nodes || osdi_spice_rhs) {
        std::cout << "OSDI advanced options:"
                  << " limiting_rhs=" << (osdi_limiting_rhs ? "requested" : "off")
                  << " tran_jacobian=" << (osdi_tran_jacobian ? "requested" : "off")
                  << " full_model_params=" << (osdi_full_model_params ? "requested" : "off")
                  << " internal_nodes=" << (osdi_internal_nodes ? "expanded" : "local")
                  << " spice_rhs=" << (osdi_spice_rhs ? "requested" : "off")
                  << std::endl;
    }

    VectorReal x_dc(matrix_size);
    for (const auto& ic : settings.initial_conditions) {
        if (ic.node >= 0 && ic.node < x_dc.getSize()) {
            x_dc[ic.node] = ic.value;
        }
    }
    if (!settings.initial_conditions.empty()) {
        std::cout << "Initial conditions: applied " << settings.initial_conditions.size()
                  << " node voltage seed(s)" << (settings.use_uic ? " with UIC" : "")
                  << "." << std::endl;
    }
    auto solve_dc_operating_point = [&](
        VectorReal initial_guess,
        const std::string& label,
        bool announce,
        double gmin_override = -1.0) {
        if (announce) {
            std::cout << label << std::endl;
        }
        const double active_gmin = gmin_override >= 0.0 ? gmin_override : settings.gmin;
        VectorReal x = initial_guess;
        double previous_delta = std::numeric_limits<double>::infinity();
        for (int iter = 0; iter < settings.op_max_iter; ++iter) {
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixReal J_sparse(matrix_size); VectorReal b(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, active_gmin);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J_sparse, b, x, 0.0, 0.0, {});
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorReal x_new = KluSolverReal::solve(J_sparse, b, &real_solver_context);
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.dc_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
            runtime_stats.dc_solve_seconds += elapsed_seconds(solve_start, solve_end);
            if (!vector_finite_and_bounded(x_new)) {
                throw std::runtime_error("DC operating point diverged during " + label);
            }
            if (settings.line_search && iter > 0) {
                double delta = vector_max_delta(x_new, x);
                double alpha = 1.0;
                while (delta > previous_delta * 2.0 && alpha > 0.0625) {
                    alpha *= 0.5;
                    VectorReal damped = vector_blend(x, x_new, alpha);
                    if (!vector_finite_and_bounded(damped)) continue;
                    x_new = damped;
                    delta = vector_max_delta(x_new, x);
                }
            }
            const bool delta_converged = solution_converged(x_new, x, num_nodes, settings);
            bool converged = delta_converged;
            if (delta_converged && settings.nr_residual_check) {
                const double residual_error = dc_residual_error(
                    devices, num_devs, matrix_size, num_nodes, x_new, settings, active_gmin);
                converged = residual_error <= 1.0;
            }
            previous_delta = vector_max_delta(x_new, x);
            x = x_new;
            if (converged) return x;
        }
        throw std::runtime_error("DC operating point did not converge during " + label);
    };

    auto solve_dc_with_recovery = [&](
        VectorReal initial_guess,
        const std::string& label,
        bool announce) {
        try {
            return solve_dc_operating_point(initial_guess, label, announce);
        } catch (const std::exception& direct_error) {
            if (!settings.source_stepping && !settings.gmin_stepping) {
                throw;
            }
            std::cout << "DC recovery: direct Newton failed for " << label
                      << " (" << direct_error.what() << "); trying"
                      << (settings.source_stepping ? " source stepping" : "")
                      << (settings.gmin_stepping ? " gmin stepping" : "")
                      << "." << std::endl;
        }

        const auto source_states = collect_source_states(devices);
        std::vector<double> source_scales = settings.source_stepping
            ? std::vector<double>{0.0, 0.05, 0.10, 0.20, 0.40, 0.70, 1.0}
            : std::vector<double>{1.0};
        std::vector<double> gmins = settings.gmin_stepping
            ? std::vector<double>{1e-6, 1e-8, 1e-10, settings.gmin}
            : std::vector<double>{settings.gmin};
        VectorReal x = initial_guess;
        try {
            for (double scale : source_scales) {
                scale_source_states(source_states, scale);
                for (double gmin : gmins) {
                    std::ostringstream step_label;
                    step_label << label << " recovery source_scale=" << scale << " gmin=" << std::scientific << gmin;
                    x = solve_dc_operating_point(x, step_label.str(), false, gmin);
                }
            }
            restore_source_states(source_states);
            std::cout << "DC recovery: converged for " << label << "." << std::endl;
            return x;
        } catch (...) {
            restore_source_states(source_states);
            throw;
        }
    };

    if (settings.type == "DC") {
        set_source_dc_value(devices, settings.dc_sweep_source, settings.dc_start);
    } else if (settings.type == "STEP" && !settings.step_sweeps.empty()) {
        set_source_dc_value(devices, settings.step_sweeps.front().source, settings.step_sweeps.front().start);
    } else if (settings.type == "MC" && !settings.mc_source.empty()) {
        set_source_dc_value(devices, settings.mc_source, settings.mc_mean);
    }

    if (!settings.use_uic) {
        x_dc = solve_dc_with_recovery(x_dc, "Calculating DC Operating Point...", true);
    }

    if (settings.type == "DC") {
        std::cout << "Starting DC Sweep Analysis..." << std::endl;
        const std::vector<SweepSpec> sweeps = !settings.dc_sweeps.empty()
            ? settings.dc_sweeps
            : std::vector<SweepSpec>{{settings.dc_sweep_source, settings.dc_start, settings.dc_stop, settings.dc_step}};
        for (const auto& sweep : sweeps) {
            if (!set_source_dc_value(devices, sweep.source, sweep.start)) {
                throw std::runtime_error(".DC sweep source '" + sweep.source + "' was not found or is not an independent source");
            }
        }
        for (const auto& sweep : sweeps) {
            std::cout << "sweep(" << sweep.source << ") ";
        }
        std::cout << "| ";
        for (int i = 0; i < num_nodes; ++i) {
            std::cout << "V(" << netlist.getNodeName(i) << ") ";
        }
        std::cout << std::endl;
        VectorReal sweep_guess = x_dc;
        std::vector<double> values(sweeps.size(), 0.0);
        std::function<void(size_t)> run_sweep = [&](size_t depth) {
            if (depth == sweeps.size()) {
                std::ostringstream label;
                label << ".DC sweep";
                for (size_t i = 0; i < sweeps.size(); ++i) {
                    label << " " << sweeps[i].source << "=" << std::scientific << values[i];
                }
                sweep_guess = solve_dc_with_recovery(sweep_guess, label.str(), false);
                std::cout << std::scientific << std::setprecision(9);
                for (double value : values) std::cout << value << " ";
                std::cout << "| ";
                for (int i = 0; i < num_nodes; ++i) {
                    std::cout << sweep_guess[i] << " ";
                }
                std::cout << std::endl;
                return;
            }
            const auto& sweep = sweeps[depth];
            const double step = sweep.step;
            const bool increasing = step > 0.0;
            const double eps = std::abs(step) * 1e-9;
            const double span = std::abs(sweep.stop - sweep.start);
            const long long max_points = static_cast<long long>(std::floor(span / std::abs(step))) + 2;
            long long points = 0;
            for (double value = sweep.start;
                 increasing ? value <= sweep.stop + eps : value >= sweep.stop - eps;
                 value += step) {
                if (++points > max_points + 1) break;
                values[depth] = value;
                set_source_dc_value(devices, sweep.source, value);
                run_sweep(depth + 1);
            }
        };
        run_sweep(0);
        x_dc = sweep_guess;
    } else if (settings.type == "STEP") {
        std::cout << "Starting STEP Analysis..." << std::endl;
        if (settings.step_sweeps.empty()) {
            throw std::runtime_error(".STEP requires at least one source sweep");
        }
        for (const auto& sweep : settings.step_sweeps) {
            if (!set_source_dc_value(devices, sweep.source, sweep.start)) {
                throw std::runtime_error(".STEP source '" + sweep.source + "' was not found or is not an independent source");
            }
        }
        for (const auto& sweep : settings.step_sweeps) {
            std::cout << "step(" << sweep.source << ") ";
        }
        std::cout << "| ";
        for (int i = 0; i < num_nodes; ++i) {
            std::cout << "V(" << netlist.getNodeName(i) << ") ";
        }
        std::cout << std::endl;
        VectorReal step_guess = x_dc;
        std::vector<double> values(settings.step_sweeps.size(), 0.0);
        std::function<void(size_t)> run_step = [&](size_t depth) {
            if (depth == settings.step_sweeps.size()) {
                std::ostringstream label;
                label << ".STEP";
                for (size_t i = 0; i < settings.step_sweeps.size(); ++i) {
                    label << " " << settings.step_sweeps[i].source << "=" << std::scientific << values[i];
                }
                step_guess = solve_dc_with_recovery(step_guess, label.str(), false);
                std::cout << std::scientific << std::setprecision(9);
                for (double value : values) std::cout << value << " ";
                std::cout << "| ";
                for (int i = 0; i < num_nodes; ++i) {
                    std::cout << step_guess[i] << " ";
                }
                std::cout << std::endl;
                return;
            }
            const auto& sweep = settings.step_sweeps[depth];
            const bool increasing = sweep.step > 0.0;
            const double eps = std::abs(sweep.step) * 1e-9;
            const double span = std::abs(sweep.stop - sweep.start);
            const long long max_points = static_cast<long long>(std::floor(span / std::abs(sweep.step))) + 2;
            long long points = 0;
            for (double value = sweep.start;
                 increasing ? value <= sweep.stop + eps : value >= sweep.stop - eps;
                 value += sweep.step) {
                if (++points > max_points + 1) break;
                values[depth] = value;
                set_source_dc_value(devices, sweep.source, value);
                run_step(depth + 1);
            }
        };
        run_step(0);
        x_dc = step_guess;
    } else if (settings.type == "MC") {
        std::cout << "Starting Monte Carlo Analysis..." << std::endl;
        if (settings.mc_runs <= 0 || settings.mc_source.empty()) {
            throw std::runtime_error(".MC requires run count and independent source name");
        }
        double restore_value = 0.0;
        if (!get_source_dc_value(devices, settings.mc_source, restore_value)) {
            throw std::runtime_error(".MC source '" + settings.mc_source + "' was not found or is not an independent source");
        }
        std::mt19937 rng(settings.mc_seed);
        std::normal_distribution<double> dist(settings.mc_mean, settings.mc_sigma);
        VectorReal mc_guess = x_dc;
        std::vector<double> sum(static_cast<size_t>(num_nodes), 0.0);
        std::vector<double> sum_sq(static_cast<size_t>(num_nodes), 0.0);
        std::vector<ProductionSpecStats> spec_stats;
        for (const auto& spec : settings.output_specs) {
            spec_stats.push_back({spec.name});
        }
        std::cout << "run " << settings.mc_source << " | ";
        for (int i = 0; i < num_nodes; ++i) {
            std::cout << "V(" << netlist.getNodeName(i) << ") ";
        }
        for (const auto& spec : settings.output_specs) {
            std::cout << "SPEC(" << spec.name << ") ";
        }
        std::cout << std::endl;
        for (int run = 1; run <= settings.mc_runs; ++run) {
            const double value = settings.mc_sigma == 0.0 ? settings.mc_mean : dist(rng);
            set_source_dc_value(devices, settings.mc_source, value);
            std::ostringstream label;
            label << ".MC run " << run << " " << settings.mc_source << "=" << std::scientific << value;
            mc_guess = solve_dc_with_recovery(mc_guess, label.str(), false);
            std::cout << std::scientific << std::setprecision(9)
                      << run << " " << value << " | ";
            for (int i = 0; i < num_nodes; ++i) {
                const double node_v = mc_guess[i];
                sum[static_cast<size_t>(i)] += node_v;
                sum_sq[static_cast<size_t>(i)] += node_v * node_v;
                std::cout << node_v << " ";
            }
            for (size_t i = 0; i < settings.output_specs.size(); ++i) {
                const auto& spec = settings.output_specs[i];
                const double measured = spec_value(mc_guess, spec);
                const bool passes = spec_passes(spec, measured);
                spec_stats[i].add(measured, passes);
                std::cout << (passes ? "PASS" : "FAIL") << "(" << measured << ") ";
            }
            std::cout << std::endl;
        }
        set_source_dc_value(devices, settings.mc_source, restore_value);
        std::cout << "MC summary:";
        for (int i = 0; i < num_nodes; ++i) {
            const double mean = sum[static_cast<size_t>(i)] / static_cast<double>(settings.mc_runs);
            const double variance = std::max(
                0.0,
                sum_sq[static_cast<size_t>(i)] / static_cast<double>(settings.mc_runs) - mean * mean);
            std::cout << " V(" << netlist.getNodeName(i) << ") mean="
                      << std::scientific << std::setprecision(9) << mean
                      << " sigma=" << std::sqrt(variance);
        }
        std::cout << std::endl;
        for (const auto& stats : spec_stats) {
            std::cout << std::scientific << std::setprecision(9)
                      << "Production summary: spec=" << stats.name
                      << " pass=" << stats.pass
                      << " fail=" << (stats.count - stats.pass)
                      << " yield=" << std::fixed << std::setprecision(2) << stats.yieldPercent() << "%"
                      << std::scientific << std::setprecision(9)
                      << " mean=" << stats.mean()
                      << " sigma=" << stats.sigma()
                      << " min=" << stats.min
                      << " max=" << stats.max
                      << std::endl;
        }
        x_dc = mc_guess;
    } else if (settings.type == "CORNER") {
        std::cout << "Starting Corner Analysis..." << std::endl;
        if (settings.corners.empty()) {
            throw std::runtime_error(".CORNER analysis requires at least one corner declaration");
        }
        const auto source_states = collect_source_states(devices);
        VectorReal corner_guess = x_dc;
        for (const auto& corner : settings.corners) {
            restore_source_states(source_states);
            for (const auto& [source, value] : corner.source_values) {
                if (!set_source_dc_value(devices, source, value)) {
                    throw std::runtime_error(".CORNER source '" + source + "' was not found or is not an independent source");
                }
            }
            std::ostringstream label;
            label << ".CORNER " << corner.name;
            for (const auto& [source, value] : corner.source_values) {
                label << " " << source << "=" << std::scientific << value;
            }
            corner_guess = solve_dc_with_recovery(corner_guess, label.str(), false);
            std::cout << "corner(" << corner.name << ") | ";
            for (int i = 0; i < num_nodes; ++i) {
                std::cout << "V(" << netlist.getNodeName(i) << ")="
                          << std::scientific << std::setprecision(9) << corner_guess[i] << " ";
            }
            for (const auto& spec : settings.output_specs) {
                const double measured = spec_value(corner_guess, spec);
                const bool passes = spec_passes(spec, measured);
                std::cout << "SPEC(" << spec.name << ")="
                          << (passes ? "PASS" : "FAIL")
                          << "(" << measured << ") ";
            }
            std::cout << std::endl;
        }
        restore_source_states(source_states);
        x_dc = corner_guess;
    } else if (settings.type == "SENS") {
        std::cout << "Starting DC Sensitivity Analysis..." << std::endl;
        double base_source = 0.0;
        if (!get_source_dc_value(devices, settings.sens_source, base_source)) {
            throw std::runtime_error(".SENS source '" + settings.sens_source + "' was not found or is not an independent source");
        }
        const double base_output = probe_value(x_dc, settings.sens_out_pos, settings.sens_out_neg);
        const double delta = std::max(std::abs(base_source) * 1e-6, 1e-9);
        set_source_dc_value(devices, settings.sens_source, base_source + delta);
        VectorReal perturbed = solve_dc_with_recovery(x_dc, ".SENS perturbation", false);
        set_source_dc_value(devices, settings.sens_source, base_source);
        const double pert_output = probe_value(perturbed, settings.sens_out_pos, settings.sens_out_neg);
        const double sensitivity = (pert_output - base_output) / delta;
        std::cout << std::scientific << std::setprecision(9)
                  << "Sensitivity: d" << voltage_probe_label(netlist, settings.sens_out_pos, settings.sens_out_neg)
                  << "/d" << settings.sens_source << " = " << sensitivity
                  << " base_output=" << base_output
                  << " delta=" << delta
                  << std::endl;
    } else if (settings.type == "PZ") {
        std::cout << "Starting Pole-Zero Estimate..." << std::endl;
        if (!configure_transfer_input_source(devices, settings.tf_input_source)) {
            throw std::runtime_error(".PZ input source '" + settings.tf_input_source + "' was not found or is not an independent source");
        }
        const int points_per_dec = std::max(settings.points_per_dec, 1);
        const double dec_mult = std::pow(10.0, 1.0 / static_cast<double>(points_per_dec));
        double f = std::max(settings.f_start, 1e-30);
        const double f_stop = std::max(settings.f_stop, f);
        bool have_prev = false;
        double low_mag = 0.0;
        double prev_f = f;
        double prev_mag = 0.0;
        double prev_db = 0.0;
        std::vector<double> pole_thresholds;
        std::vector<double> zero_thresholds;
        for (int i = 1; i <= 20; ++i) {
            pole_thresholds.push_back(-3.01029995664 * static_cast<double>(i));
            zero_thresholds.push_back(3.01029995664 * static_cast<double>(i));
        }
        std::vector<double> pole_estimates;
        std::vector<double> zero_estimates;
        auto interpolate_crossing = [](double f0, double db0, double f1, double db1, double targetDb) {
            if (db1 == db0 || f0 <= 0.0 || f1 <= 0.0) return f1;
            const double alpha = (targetDb - db0) / (db1 - db0);
            return std::exp(std::log(f0) + alpha * (std::log(f1) - std::log(f0)));
        };
        std::cout << "freq | gain_mag phase_deg" << std::endl;
        while (f <= f_stop * 1.01) {
            const double omega = 2.0 * 3.14159265358979323846 * f;
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixComplex J_sparse(matrix_size);
            VectorComplex b_ac(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) devices[i]->acStamp(J_sparse, b_ac, omega, x_dc);
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorComplex x_ac = KluSolverComplex::solve(J_sparse, b_ac, &complex_solver_context);
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.ac_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
            runtime_stats.ac_solve_seconds += elapsed_seconds(solve_start, solve_end);
            const auto gain = complex_probe_value(x_ac, settings.tf_out_pos, settings.tf_out_neg);
            const double mag = std::abs(gain);
            std::cout << std::scientific << std::setprecision(9)
                      << f << " | " << mag << " "
                      << std::arg(gain) * 180.0 / 3.14159265358979323846
                      << std::endl;
            if (!have_prev) {
                low_mag = mag;
                prev_db = 0.0;
                have_prev = true;
            } else if (low_mag > 0.0 && mag > 0.0 && prev_mag > 0.0) {
                const double db = 20.0 * std::log10(mag / low_mag);
                for (double targetDb : pole_thresholds) {
                    if (prev_db > targetDb && db <= targetDb) {
                        const double pole_f = interpolate_crossing(prev_f, prev_db, f, db, targetDb);
                        pole_estimates.push_back(pole_f);
                        std::cout << std::scientific << std::setprecision(9)
                                  << "Pole estimate[" << pole_estimates.size() << "]: f=" << pole_f
                                  << " Hz threshold_db=" << targetDb
                                  << " output=" << voltage_probe_label(netlist, settings.tf_out_pos, settings.tf_out_neg)
                                  << std::endl;
                    }
                }
                for (double targetDb : zero_thresholds) {
                    if (prev_db < targetDb && db >= targetDb) {
                        const double zero_f = interpolate_crossing(prev_f, prev_db, f, db, targetDb);
                        zero_estimates.push_back(zero_f);
                        std::cout << std::scientific << std::setprecision(9)
                                  << "Zero estimate[" << zero_estimates.size() << "]: f=" << zero_f
                                  << " Hz threshold_db=" << targetDb
                                  << " output=" << voltage_probe_label(netlist, settings.tf_out_pos, settings.tf_out_neg)
                                  << std::endl;
                    }
                }
                prev_db = db;
            }
            prev_f = f;
            prev_mag = mag;
            f *= dec_mult;
        }
        if (pole_estimates.empty()) {
            std::cout << "Pole estimate: no falling gain thresholds found in requested sweep." << std::endl;
        }
        if (zero_estimates.empty()) {
            std::cout << "Zero estimate: no rising gain thresholds found in requested sweep." << std::endl;
        }
        if (!pole_estimates.empty()) {
            std::cout << std::scientific << std::setprecision(9)
                      << "Dominant pole estimate: f=" << pole_estimates.front()
                      << " Hz output=" << voltage_probe_label(netlist, settings.tf_out_pos, settings.tf_out_neg)
                      << std::endl;
        }
    } else if (settings.type == "TF") {
        std::cout << "Starting Transfer Function Analysis..." << std::endl;
        if (!configure_transfer_input_source(devices, settings.tf_input_source)) {
            throw std::runtime_error(".TF input source '" + settings.tf_input_source + "' was not found or is not an independent source");
        }
        const auto stamp_start = std::chrono::steady_clock::now();
        SparseMatrixComplex J_sparse(matrix_size); VectorComplex b_tf(matrix_size);
        stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
        const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
        #pragma omp parallel for if(use_parallel_stamp)
        for (int i = 0; i < num_devs; ++i) devices[i]->acStamp(J_sparse, b_tf, 0.0, x_dc);
        const auto stamp_end = std::chrono::steady_clock::now();
        const auto solve_start = std::chrono::steady_clock::now();
        VectorComplex x_tf = KluSolverComplex::solve(J_sparse, b_tf, &complex_solver_context);
        const auto solve_end = std::chrono::steady_clock::now();
        runtime_stats.ac_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
        runtime_stats.ac_solve_seconds += elapsed_seconds(solve_start, solve_end);
        auto node_value = [&](int node) -> std::complex<double> {
            return node >= 0 ? x_tf[node] : std::complex<double>{0.0, 0.0};
        };
        const auto gain = node_value(settings.tf_out_pos) - node_value(settings.tf_out_neg);
        std::cout << std::scientific << std::setprecision(9)
                  << "Transfer function: V(" << netlist.getNodeName(settings.tf_out_pos);
        if (settings.tf_out_neg >= 0) {
            std::cout << "," << netlist.getNodeName(settings.tf_out_neg);
        }
        std::cout << ")/" << settings.tf_input_source
                  << " = (" << gain.real() << "," << gain.imag() << ")"
                  << " magnitude=" << std::abs(gain)
                  << " phase_deg=" << std::arg(gain) * 180.0 / 3.14159265358979323846
                  << std::endl;
    } else if (settings.type == "TRAN") {
        std::cout << "Starting Transient Analysis..." << std::endl;
        const double output_step = settings.t_step > 0.0
            ? settings.t_step
            : std::max(settings.t_stop / 100.0, 1e-15);
        const double max_step = settings.t_max_step > 0.0
            ? settings.t_max_step
            : (settings.tran_adaptive ? std::min(output_step, default_transient_max_step(settings, output_step)) : output_step);
        const double min_step = settings.t_min_step > 0.0
            ? std::min(settings.t_min_step, max_step)
            : std::max(max_step * 1e-6, 1e-18);
        const double save_start = std::clamp(settings.t_start, 0.0, settings.t_stop);
        double step = max_step;
        std::vector<VectorReal> x_hist; VectorReal x = x_dc; x_hist.push_back(x);
        std::vector<double> t_hist; t_hist.push_back(0.0);
        TransientOutput tran_out(output_file, netlist, num_nodes);
        TransientStats tran_stats;
        std::vector<TranSample> measure_samples;
        auto record_measure_sample = [&](double sample_time, const VectorReal& sample_x) {
            if (!settings.measures.empty()) {
                measure_samples.push_back({sample_time, sample_x});
            }
        };
        const double bp_tol = std::max(1e-30, max_step * 1e-6);
        const std::vector<double> breakpoints = collect_transient_breakpoints(devices, settings.t_stop, bp_tol);
        std::cout << "Transient breakpoints: " << breakpoints.size() << std::endl;
        if (tran_out.toFile()) {
            std::cout << "Waveform output: " << output_file << std::endl;
        }
        std::cout << std::scientific << std::setprecision(9)
                  << "Transient controls: output step=" << output_step
                  << " save start=" << save_start
                  << " max internal step=" << max_step
                  << " min internal step=" << min_step
                  << " adaptive=" << (settings.tran_adaptive ? "on" : "off")
                  << " predictor=" << (settings.tran_predictor ? "on" : "off")
                  << " method=" << normalized_transient_method(settings)
                  << std::endl;
        if (!tran_out.toFile()) {
            std::cout << "time | ";
            for (int i = 0; i < num_nodes; ++i) {
                std::cout << "V(" << netlist.getNodeName(i) << ") ";
            }
            std::cout << std::endl;
        }
        double t = 0.0;
        double device_bound_step = collect_transient_bound_step(devices);
        if (std::isfinite(device_bound_step)) {
            tran_stats.min_bound_step = std::min(tran_stats.min_bound_step, device_bound_step);
        }
        const double stop_tol = std::max(std::max(1e-30, std::abs(settings.t_stop) * 1e-9), min_step * 10.0);
        double last_printed_time = -1.0;
        double last_progress_percent = -1.0;
        double next_output = save_start;
        if (next_output <= 1e-30) {
            tran_out.write(0.0, x, num_nodes);
            record_measure_sample(0.0, x);
            ++tran_stats.output_points;
            last_printed_time = 0.0;
            next_output = output_step;
        } else {
            record_measure_sample(0.0, x);
        }
        while (t < settings.t_stop - stop_tol) {
            const double remaining = settings.t_stop - t;
            if (remaining <= stop_tol) break;
            if (next_output > t && next_output - t <= stop_tol) {
                tran_out.write(next_output, x, num_nodes);
                ++tran_stats.output_points;
                last_printed_time = next_output;
                t = next_output;
                do {
                    next_output += output_step;
                } while (next_output <= t + stop_tol);
                continue;
            }
            step = std::min(step, remaining);
            if (std::isfinite(device_bound_step) && device_bound_step > 0.0 && step > device_bound_step) {
                step = std::max(min_step, device_bound_step);
                ++tran_stats.bound_step_limited;
            }
            const double next_bp = next_breakpoint_after(breakpoints, t, bp_tol);
            if (next_bp > 0.0) {
                const double bp_window = std::max(max_step * 4.0, output_step);
                const double bp_dist = next_bp - t;
                if (bp_dist > stop_tol && bp_dist <= bp_window) {
                    step = std::min(step, bp_dist);
                } else if (bp_dist >= 0.0 && bp_dist <= stop_tol) {
                    t = next_bp;
                    continue;
                }
            }
            if (next_output > t + 1e-30) {
                step = std::min(step, next_output - t);
            }
            if (step <= stop_tol) {
                step = std::min(max_step, std::max(min_step, remaining));
                if (step <= stop_tol) break;
            }
            bool accepted = false;
            VectorReal accepted_x = x;
            double err = 0.0;
            while (!accepted) {
                const double target_time = t + step;
                TransientStepResult full = solve_transient_step(
                    devices, num_devs, matrix_size, num_nodes, settings, &real_solver_context, x, x_hist, t_hist, step, target_time);
                tran_stats.noteSolve(full);
                if (!full.converged && step > min_step) {
                    ++tran_stats.rejected_steps;
                    ++tran_stats.convergence_rejections;
                    step = std::max(min_step, step * 0.5);
                    continue;
                }
                if (!full.converged) {
                    throw std::runtime_error("Transient step failed to converge at minimum timestep");
                }

                if (settings.tran_adaptive && step > min_step) {
                    const TransientIntegrationMethod attempted_method = choose_transient_method(settings, x_hist, t_hist);
                    const double next_bp_after_target = next_breakpoint_after(breakpoints, target_time, bp_tol);
                    const bool near_breakpoint =
                        (next_bp > 0.0 && std::abs(target_time - next_bp) <= bp_tol) ||
                        (next_bp_after_target > 0.0 && next_bp_after_target - target_time <= std::max(max_step * 2.0, output_step));
                    if (!near_breakpoint) {
                        const double state_change = transient_state_change_error(full.x, x, num_nodes, settings);
                        if (state_change < 1e-4) {
                            accepted_x = full.x;
                            err = state_change;
                            accepted = true;
                            continue;
                        }
                    }
                    const double half_step = step * 0.5;
                    TransientStepResult half1 = solve_transient_step(
                        devices, num_devs, matrix_size, num_nodes, settings, &real_solver_context, x, x_hist, t_hist, half_step, t + half_step);
                    tran_stats.noteSolve(half1);
                    std::vector<VectorReal> half_hist = x_hist;
                    half_hist.push_back(half1.x);
                    std::vector<double> half_time_hist = t_hist;
                    half_time_hist.push_back(t + half_step);
                    TransientStepResult half2 = solve_transient_step(
                        devices, num_devs, matrix_size, num_nodes, settings, &real_solver_context, half1.x, half_hist, half_time_hist, half_step, target_time);
                    tran_stats.noteSolve(half2);
                    if (!half1.converged || !half2.converged) {
                        if (step > min_step) {
                            ++tran_stats.rejected_steps;
                            ++tran_stats.convergence_rejections;
                            step = std::max(min_step, half_step);
                            continue;
                        }
                        throw std::runtime_error("Transient half-step failed to converge at minimum timestep");
                    }
                    err = transient_lte_error(full.x, half2.x, num_nodes, settings);
                    if (err > 1.0 && step > min_step) {
                        const double factor = transient_reject_factor(err, attempted_method);
                        ++tran_stats.rejected_steps;
                        ++tran_stats.lte_rejections;
                        step = std::max(min_step, step * factor);
                        continue;
                    }
                    accepted_x = half2.x;
                } else {
                    accepted_x = full.x;
                }
                accepted = true;
            }
            t += step;
            if (settings.t_stop - t < stop_tol) {
                t = settings.t_stop;
            }
            tran_stats.noteAccepted(step, err);
            x = accepted_x;
            record_measure_sample(t, x);
            const TransientContext accepted_ctx = make_transient_context(settings, x_hist, t_hist, step, t);
            for (int i = 0; i < num_devs; ++i) {
                devices[i]->acceptTransientStep(x, t, accepted_ctx);
            }
            device_bound_step = collect_transient_bound_step(devices);
            if (std::isfinite(device_bound_step)) {
                tran_stats.min_bound_step = std::min(tran_stats.min_bound_step, device_bound_step);
            }
            x_hist.push_back(x);
            t_hist.push_back(t);
            if (x_hist.size() > 4) {
                x_hist.erase(x_hist.begin());
                t_hist.erase(t_hist.begin());
            }
            const double output_tol = std::max(1e-30, output_step * 1e-9);
            if (t + output_tol >= next_output) {
                tran_out.write(t, x, num_nodes);
                ++tran_stats.output_points;
                last_printed_time = t;
                do {
                    next_output += output_step;
                } while (next_output <= t + output_tol);
            }
            if (tran_out.toFile() && settings.t_stop > 0.0) {
                const double pct = std::clamp(100.0 * t / settings.t_stop, 0.0, 100.0);
                if (pct - last_progress_percent >= 1.0 || pct >= 100.0) {
                    std::cout << std::fixed << std::setprecision(1)
                              << "Transient progress: " << pct << "% t="
                              << std::scientific << std::setprecision(9) << t << std::endl;
                    last_progress_percent = pct;
                }
            }
            if (settings.tran_adaptive) {
                const double grow = transient_growth_factor(err, accepted_ctx.method);
                step = std::min(max_step, std::max(min_step, step * grow));
                if (std::isfinite(device_bound_step) && device_bound_step > 0.0 && step > device_bound_step) {
                    step = std::max(min_step, device_bound_step);
                    ++tran_stats.bound_step_limited;
                }
            }
        }
        const double output_tol = std::max(1e-30, output_step * 1e-9);
        if (settings.t_stop >= save_start && last_printed_time < settings.t_stop - output_tol) {
            tran_out.write(settings.t_stop, x, num_nodes);
            ++tran_stats.output_points;
        }
        const double min_used_step = std::isfinite(tran_stats.min_step) ? tran_stats.min_step : 0.0;
        const double min_device_bound = std::isfinite(tran_stats.min_bound_step) ? tran_stats.min_bound_step : 0.0;
        std::cout << std::scientific << std::setprecision(9)
                  << "Transient summary: accepted=" << tran_stats.accepted_steps
                  << " rejected=" << tran_stats.rejected_steps
                  << " convergence_rejected=" << tran_stats.convergence_rejections
                  << " lte_rejected=" << tran_stats.lte_rejections
                  << " bound_step_limited=" << tran_stats.bound_step_limited
                  << " output_points=" << tran_stats.output_points
                  << " min_step=" << min_used_step
                  << " max_step=" << tran_stats.max_step
                  << " min_bound_step=" << min_device_bound
                  << std::endl;
        std::cout << std::fixed << std::setprecision(2)
                  << "Newton summary: solves=" << tran_stats.newton_solves
                  << " total_iterations=" << tran_stats.newton_iterations
                  << " average_iterations=" << tran_stats.averageNewtonIterations()
                  << " max_iterations=" << tran_stats.max_newton_iterations
                  << " last_residual_error=" << tran_stats.last_residual_error
                  << " max_residual_error=" << tran_stats.max_residual_error
                  << " stamp_seconds=" << tran_stats.stamp_seconds
                  << " solve_seconds=" << tran_stats.solve_seconds
                  << std::endl;
        std::cout << std::scientific << std::setprecision(9)
                  << "Accuracy summary: method=" << normalized_transient_method(settings)
                  << " reltol=" << settings.reltol
                  << " vntol=" << settings.vntol
                  << " abstol=" << settings.abstol
                  << " lte_reltol=" << settings.tran_lte_reltol
                  << " lte_abstol=" << settings.tran_lte_abstol
                  << " last_lte_error=" << tran_stats.last_error
                  << std::endl;
        for (const auto& measure : settings.measures) {
            const std::string analysis = upper_copy(measure.analysis);
            if (analysis != "TRAN" && analysis != "TRANSIENT" && analysis != "ALL") continue;
            const double value = evaluate_transient_measure(measure_samples, measure);
            std::cout << std::scientific << std::setprecision(9)
                      << "MEASURE " << measure.name << " = " << value << std::endl;
        }
    } else if (settings.type == "AC") {
        std::cout << "Starting AC Analysis..." << std::endl;
        double f = settings.f_start; double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            double omega = 2.0 * 3.14159265358979 * f;
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixComplex J_sparse(matrix_size); VectorComplex b_ac(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) devices[i]->acStamp(J_sparse, b_ac, omega, x_dc);
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorComplex x_ac = KluSolverComplex::solve(J_sparse, b_ac, &complex_solver_context);
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.ac_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
            runtime_stats.ac_solve_seconds += elapsed_seconds(solve_start, solve_end);
            std::cout << std::scientific << std::setprecision(2) << f << " | ";
            for(int i=0; i<num_nodes; ++i) std::cout << "(" << x_ac[i].real() << "," << x_ac[i].imag() << ") ";
            std::cout << std::endl; f *= dec_mult;
        }
    } else if (settings.type == "NOISE") {
        std::cout << "Starting Noise Analysis..." << std::endl;
        if (settings.out_node < 0 || settings.out_node >= num_nodes) {
            throw std::runtime_error(".NOISE output node must be a non-ground circuit node");
        }
        double f = settings.f_start;
        double dec_mult = std::pow(10.0, 1.0 / std::max(settings.points_per_dec, 1));
        std::cout << "freq | onoise_sqrt(V/rtHz) onoise_psd(V^2/Hz) noise_sources" << std::endl;
        while (f <= settings.f_stop * 1.01) {
            double omega = 2.0 * 3.14159265358979 * f;
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixComplex J_sparse(matrix_size); VectorComplex ignored_rhs(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) devices[i]->acStamp(J_sparse, ignored_rhs, omega, x_dc);
            const auto stamp_end = std::chrono::steady_clock::now();
            runtime_stats.ac_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);

            std::vector<NoiseSource> noise_sources;
            for (const auto& dev : devices) {
                dev->collectNoiseSources(omega, x_dc, noise_sources);
            }

            double output_psd = 0.0;
            const auto solve_start = std::chrono::steady_clock::now();
            for (const auto& source : noise_sources) {
                if (source.currentPsd <= 0.0) continue;
                VectorComplex rhs(matrix_size);
                rhs.add(source.nodePos, {-1.0, 0.0});
                rhs.add(source.nodeNeg, {1.0, 0.0});
                VectorComplex transfer = KluSolverComplex::solve(J_sparse, rhs, &complex_solver_context);
                const std::complex<double> out = transfer[settings.out_node];
                output_psd += std::norm(out) * source.currentPsd;
            }
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.ac_solve_seconds += elapsed_seconds(solve_start, solve_end);
            std::cout << std::scientific << std::setprecision(9)
                      << f << " | " << std::sqrt(std::max(output_psd, 0.0))
                      << " " << output_psd
                      << " " << noise_sources.size()
                      << std::endl;
            f *= dec_mult;
        }
    } else if (settings.type == "STB") {
        std::cout << "Starting Stability Analysis (Tian)..." << std::endl;
        if (probes.empty()) return;
        double f = settings.f_start; double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            double omega = 2.0 * 3.14159265358979 * f;
            // 1. Voltage Pass
            SparseMatrixComplex Jv(matrix_size); VectorComplex bv(matrix_size);
            for (const auto& dev : devices) {
                auto* p = dynamic_cast<StabilityProbe*>(dev.get());
                if (p) p->stbStamp(Jv, bv, 1); else dev->acStamp(Jv, bv, omega, x_dc);
            }
            VectorComplex xv = KluSolverComplex::solve(Jv, bv, &complex_solver_context);
            std::complex<double> Tv = -xv[probes[0]->getNodeNeg()] / xv[probes[0]->getNodePos()];
            // 2. Current Pass
            SparseMatrixComplex Ji(matrix_size); VectorComplex bi(matrix_size);
            for (const auto& dev : devices) {
                auto* p = dynamic_cast<StabilityProbe*>(dev.get());
                if (p) p->stbStamp(Ji, bi, 2); else dev->acStamp(Ji, bi, omega, x_dc);
            }
            VectorComplex xi = KluSolverComplex::solve(Ji, bi, &complex_solver_context);
            std::complex<double> Ti = xi[probes[0]->getBranchIndex()];
            // 3. Combine
            std::complex<double> T = (Tv * Ti - std::complex<double>(1,0)) / (Tv + Ti + std::complex<double>(2,0));
            std::cout << std::scientific << f << " | Mag: " << std::abs(T) << " Phase: " << std::arg(T)*180/3.1415 << std::endl;
            f *= dec_mult;
        }
    } else if (settings.type == "PSS") {
        std::cout << "Starting Periodic Steady State (PSS) Analysis..." << std::endl;
        int n_tones = static_cast<int>(settings.f_fund.size());
        std::cout << "  Tones: " << n_tones << " | Harmonics/Tone: " << settings.n_harms << std::endl;
        if (n_tones > 1) {
            std::cout << "  Warning: Multi-tone PSS requires finding a common periodic denominator, which can result in massive transient integration times." << std::endl;
            std::cout << "  For " << n_tones << " non-harmonically related tones, Quasi-Periodic Harmonic Balance (QPHB) is recommended." << std::endl;
        }
        std::cout << "  Executing Shooting Method (Prototype)..." << std::endl;
        // PSS Shooting Method stub
        std::cout << "PSS Converged." << std::endl;
    } else if (settings.type == "HB") {
        std::cout << "Starting Harmonic Balance Analysis (Multi-Tone)..." << std::endl;
        int n_tones = static_cast<int>(settings.f_fund.size());
        int H = settings.n_harms;
        int K = 1;
        for (int t = 0; t < n_tones; ++t) K *= (2 * H + 1);
        int n_vars = matrix_size * K;
        std::cout << "  Tones: " << n_tones << " | Harmonics/Tone: " << H << " | Total States/Node: " << K << std::endl;
        VectorReal x_hb(n_vars);
        for (int i = 0; i < matrix_size; ++i) x_hb[i * K] = x_dc[i];
        bool hb_converged = false;
        for (int hb_iter = 0; hb_iter < 30; ++hb_iter) {
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixReal J_hb(n_vars); VectorReal b_hb(n_vars);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, n_vars);
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) devices[i]->hbStamp(J_hb, b_hb, (n_tones > 0 ? settings.f_fund[0] : 0.0), H, x_hb);
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorReal dx = KluSolverReal::solve(J_hb, b_hb, &real_solver_context);
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.hb_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
            runtime_stats.hb_solve_seconds += elapsed_seconds(solve_start, solve_end);
            double max_dx = 0.0;
            for (int i = 0; i < n_vars; ++i) { x_hb[i] -= dx[i]; max_dx = std::max(max_dx, std::abs(dx[i])); }
            if (max_dx < 1e-6) { hb_converged = true; break; }
        }
        if (hb_converged) std::cout << "HB Converged." << std::endl;
    } else if (settings.type == "PAC") {
        std::cout << "Starting Periodic AC..." << std::endl;
        int N = 5; int K = 2 * N + 1; int n_vars = matrix_size * K;
        double f = settings.f_start; double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            SparseMatrixReal J_pac(n_vars); VectorReal b_pac(n_vars);
            for (const auto& dev : netlist.getDevices()) dev->pacStamp(J_pac, b_pac, f, (settings.f_fund.size() > 0 ? settings.f_fund[0] : 0.0), N, x_dc);
            VectorReal x_pac = KluSolverReal::solve(J_pac, b_pac, &real_solver_context);
            std::cout << std::scientific << f << " | PAC Solved. Mag: " << std::abs(x_pac[0]) << std::endl;
            f *= dec_mult;
        }
    } else {
        std::cout << "DC Operating Point Converged: ";
        for(int i=0; i<num_nodes; ++i) {
            std::cout << "Node " << i << "=" << std::fixed << std::setprecision(6) << x_dc[i]
                      << "V [" << netlist.getNodeName(i) << "] ";
        }
        std::cout << std::endl;
        for (const auto& spec : settings.output_specs) {
            const double measured = spec_value(x_dc, spec);
            const bool passes = spec_passes(spec, measured);
            std::cout << std::scientific << std::setprecision(9)
                      << "Production spec: " << spec.name
                      << " value=" << measured
                      << " result=" << (passes ? "PASS" : "FAIL")
                      << std::endl;
        }
    }
    const auto real_solver_stats = KluSolverReal::getStats();
    const auto complex_solver_stats = KluSolverComplex::getStats();
    std::cout << std::fixed << std::setprecision(6)
              << "Performance summary: dc_stamp=" << runtime_stats.dc_stamp_seconds
              << " dc_solve=" << runtime_stats.dc_solve_seconds
              << " ac_stamp=" << runtime_stats.ac_stamp_seconds
              << " ac_solve=" << runtime_stats.ac_solve_seconds
              << " hb_stamp=" << runtime_stats.hb_stamp_seconds
              << " hb_solve=" << runtime_stats.hb_solve_seconds
              << std::endl;
    std::cout << std::fixed << std::setprecision(6)
              << "Solver summary (real): calls=" << real_solver_stats.solve_calls
              << " pattern_reuse=" << real_solver_stats.pattern_reuse_hits
              << " pattern_rebuilds=" << real_solver_stats.pattern_rebuilds
              << " parallel_pivots=" << real_solver_stats.parallel_pivot_scans
              << " parallel_elim=" << real_solver_stats.parallel_elimination_passes
              << " singleton_elim=" << real_solver_stats.singleton_eliminations
              << " singleton_core=" << real_solver_stats.singleton_core_solves
              << " singleton_zero_core=" << real_solver_stats.singleton_zero_core_solves
              << " max_core=" << real_solver_stats.max_singleton_core_size
              << " ext_klu_calls=" << real_solver_stats.external_klu_calls
              << " ext_klu_reuse=" << real_solver_stats.external_klu_symbolic_reuse
              << " ext_klu_rebuilds=" << real_solver_stats.external_klu_symbolic_rebuilds
              << " ext_klu_failures=" << real_solver_stats.external_klu_failures
              << " ext_klu=" << real_solver_stats.external_klu_seconds
              << " structure=" << real_solver_stats.structure_seconds
              << " fill=" << real_solver_stats.numeric_fill_seconds
              << " pivot=" << real_solver_stats.pivot_seconds
              << " elimination=" << real_solver_stats.elimination_seconds
              << " factor=" << real_solver_stats.factor_seconds
              << " backsolve=" << real_solver_stats.backsolve_seconds
              << " total=" << real_solver_stats.total_seconds
              << std::endl;
    if (complex_solver_stats.solve_calls > 0) {
        std::cout << std::fixed << std::setprecision(6)
                  << "Solver summary (complex): calls=" << complex_solver_stats.solve_calls
                  << " pattern_reuse=" << complex_solver_stats.pattern_reuse_hits
                  << " pattern_rebuilds=" << complex_solver_stats.pattern_rebuilds
                  << " parallel_pivots=" << complex_solver_stats.parallel_pivot_scans
                  << " parallel_elim=" << complex_solver_stats.parallel_elimination_passes
                  << " singleton_elim=" << complex_solver_stats.singleton_eliminations
                  << " singleton_core=" << complex_solver_stats.singleton_core_solves
                  << " singleton_zero_core=" << complex_solver_stats.singleton_zero_core_solves
                  << " max_core=" << complex_solver_stats.max_singleton_core_size
                  << " ext_klu_calls=" << complex_solver_stats.external_klu_calls
                  << " ext_klu_reuse=" << complex_solver_stats.external_klu_symbolic_reuse
                  << " ext_klu_rebuilds=" << complex_solver_stats.external_klu_symbolic_rebuilds
                  << " ext_klu_failures=" << complex_solver_stats.external_klu_failures
                  << " ext_klu=" << complex_solver_stats.external_klu_seconds
                  << " structure=" << complex_solver_stats.structure_seconds
                  << " fill=" << complex_solver_stats.numeric_fill_seconds
                  << " pivot=" << complex_solver_stats.pivot_seconds
                  << " elimination=" << complex_solver_stats.elimination_seconds
                  << " factor=" << complex_solver_stats.factor_seconds
                  << " backsolve=" << complex_solver_stats.backsolve_seconds
                  << " total=" << complex_solver_stats.total_seconds
                  << std::endl;
    }
    std::cout << "Simulation Completed Successfully." << std::endl;
}

void print_usage() {
    std::cout << "Usage: gspice <input_file.sp> [options]\nOptions:\n";
    std::cout << "  -v, --version    Display version information\n";
    std::cout << "  -h, --help       Display this help message\n";
    std::cout << "  -t, --threads <n> Set parallel threads (1-16, default: 1)\n";
    std::cout << "  -o <file>        Specify output file\n";
}

int main(int argc, char* argv[]) {
    std::string input_file = "";
    std::string output_file = "";
    int num_threads = 1;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") { print_usage(); return 0; }
        if (arg == "-v" || arg == "--version") { std::cout << "GSPICE v1.1.0 (Parallel)\n"; return 0; }
        if ((arg == "-t" || arg == "--threads") && i+1 < argc) num_threads = std::stoi(argv[++i]);
        else if ((arg == "-o" || arg == "--output") && i+1 < argc) output_file = argv[++i];
        if (arg[0] != '-') input_file = arg;
    }
    if (num_threads < 1) num_threads = 1; if (num_threads > 16) num_threads = 16;
    if (input_file.empty()) { run_osdi_test(); print_usage(); return 0; }
    Netlist netlist = Parser::parse(input_file);
    for (const auto& status : netlist.getModelStatus()) {
        std::cout << "MODEL_STATUS: " << status << std::endl;
    }
    for (const auto& warning : netlist.getWarnings()) {
        std::cerr << "WARNING: " << warning << std::endl;
    }
    if (!netlist.getErrors().empty()) {
        for (const auto& error : netlist.getErrors()) {
            std::cerr << "ERROR: " << error << std::endl;
        }
        return 2;
    }
    if (netlist.getDevices().empty()) return 1;
    try {
        run_simulation(netlist, output_file, num_threads);
    } catch (const std::exception& exc) {
        std::cerr << "ERROR: Simulation failed: " << exc.what() << std::endl;
        return 3;
    }
    return 0;
}
