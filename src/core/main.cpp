// Ensure M_PI is available on MSVC before any standard math headers.
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
#include <map>
#include <filesystem>
#include <atomic>
#include <omp.h>
#include "parser.hpp"
#include "feature_registry.hpp"
#include "event_queue.hpp"
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
#include "transient_state_store.hpp"
#include "integration_formula.hpp"
#include "transient_control.hpp"
#include "dae_audit.hpp"
#include "csr_matrix.hpp"
#include "homotopy.hpp"
#include "device_arena.hpp"
#include "hb_engine.hpp"
#include "multirate.hpp"
#include "fastspice_engine.hpp"
#include "solvers/parallel_sparse_solver.hpp"
#include "ticer.hpp"

using namespace gspice;

#ifndef GSPICE_VERSION
#define GSPICE_VERSION "development"
#endif

namespace {

double elapsed_seconds(
    const std::chrono::steady_clock::time_point& start,
    const std::chrono::steady_clock::time_point& end) {
    return std::chrono::duration<double>(end - start).count();
}

int choose_active_threads(int requested_threads, int num_devs, int matrix_size) {
    (void)matrix_size;
    if (num_devs < 2) return 1;
    const int max_cpus = omp_get_max_threads();
    const int target = (requested_threads > 0) ? requested_threads : max_cpus;
    return std::clamp(target, 1, std::min(max_cpus, num_devs));
}

bool parallel_stamp_enabled(int num_devs, int matrix_size) {
    if (omp_get_max_threads() <= 1) return false;
    return num_devs >= 256 && matrix_size >= 256;
}

double default_transient_max_step(const SimulationSettings& settings, double output_step) {
    if (settings.t_max_step > 0.0) {
        return std::max(settings.t_max_step, 1e-18);
    }
// SPICE convention: max_step is an UPPER BOUND limit for adaptive timestepping.
    // It allows step size to grow during quiescent flat regions while LTE shrinks step
    // size during active transitions/oscillations.
    const double base = (settings.t_stop > 0.0) ? (settings.t_stop / 50.0) : (output_step * 50.0);
    return std::max({base, output_step * 5.0, 1e-18});
}

VectorReal interpolate_transient_state(
    const VectorReal& previous,
    const VectorReal& current,
    double previous_time,
    double current_time,
    double sample_time) {
    if (current_time <= previous_time) return current;
    const double alpha = std::clamp(
        (sample_time - previous_time) / (current_time - previous_time), 0.0, 1.0);
    VectorReal sample = previous;
    for (int i = 0; i < sample.getSize(); ++i) {
        sample[i] = previous[i] + alpha * (current[i] - previous[i]);
    }
    return sample;
}


void stamp_global_gmin(SparseMatrixReal& J, int num_nodes, double gmin) {
    (void)num_nodes;
    const int n = J.getSize();
    for (int i = 0; i < n; ++i) {
        J.add(i, i, gmin);
    }
}

void stamp_global_gmin(SparseMatrixComplex& J, int num_nodes, double gmin) {
    (void)num_nodes;
    const int n = J.getSize();
    for (int i = 0; i < n; ++i) {
        J.add(i, i, {gmin, 0.0});
    }
}

struct DaeStampStatus {
    bool limitingApplied = false;
    bool bypassed = false;
};

std::uint64_t next_evaluation_epoch() {
    static std::atomic<std::uint64_t> epoch{1};
    return epoch.fetch_add(1, std::memory_order_relaxed);
}

DaeStampStatus stamp_device_dc(
    Device& device,
    SparseMatrixReal& jacobian,
    VectorReal& rhs,
    const VectorReal& x,
    const SimulationSettings* settings = nullptr,
    std::uint64_t evaluationEpoch = 0,
    bool allowBypass = false,
    bool nodeset = false,
    bool highPrecision = false) {
    DaeRequest request;
    request.analysis = DaeAnalysis::OperatingPoint;
    request.staticResidual = true;
    request.dynamicResidual = false;
    request.staticJacobian = true;
    request.dynamicJacobian = false;
    request.highPrecision = highPrecision;
    request.nodeset = nodeset;
    request.allowBypass = allowBypass && settings && settings->nr_bypass;
    request.bypassRelativeTolerance = settings
        ? settings->reltol * settings->nr_bypass_tolerance
        : 0.0;
    request.bypassAbsoluteTolerance = settings
        ? std::max(settings->vntol, settings->abstol) * settings->nr_bypass_tolerance
        : 0.0;
    request.evaluationEpoch = evaluationEpoch;
    DaeEvaluation evaluation;
    if (device.evaluateDae(x, request, evaluation)) {
        stampDaeStatic(evaluation, x, jacobian, rhs);
        return {evaluation.limitingApplied, evaluation.bypassed};
    }
    device.dcStamp(jacobian, rhs, x, 0.0, 0.0, {});
    return {};
}

void stamp_device_ac(
    Device& device,
    SparseMatrixComplex& jacobian,
    VectorComplex& rhs,
    double omega,
    const VectorReal& operatingPoint) {
    DaeRequest request;
    request.analysis = DaeAnalysis::SmallSignal;
    request.staticResidual = false;
    request.dynamicResidual = false;
    request.staticJacobian = true;
    request.dynamicJacobian = true;
    DaeEvaluation evaluation;
    if (device.evaluateDae(operatingPoint, request, evaluation)) {
        stampDaeSmallSignal(evaluation, omega, jacobian);
        return;
    }
    device.acStamp(jacobian, rhs, omega, operatingPoint);
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

std::pair<double, int> worst_solution_update(
    const VectorReal& x_new,
    const VectorReal& x_old,
    int num_nodes,
    const SimulationSettings& settings) {
    double worst = 0.0;
    int worst_index = -1;
    const int n = std::min(x_new.getSize(), x_old.getSize());
    for (int i = 0; i < n; ++i) {
        const double scale = std::max(std::abs(x_new[i]), std::abs(x_old[i]));
        const double absolute = i < num_nodes ? settings.vntol : settings.abstol;
        const double normalized = std::abs(x_new[i] - x_old[i]) /
            std::max(absolute + settings.reltol * scale, 1e-30);
        if (normalized > worst) {
            worst = normalized;
            worst_index = i;
        }
    }
    return {worst, worst_index};
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

double linear_relative_residual(
    const SparseMatrixReal& J,
    const VectorReal& b,
    const VectorReal& x,
    VectorReal* residual = nullptr) {
    const int n = std::min(J.getSize(), std::min(b.getSize(), x.getSize()));
    std::vector<double> ax(static_cast<size_t>(n), 0.0);
    std::vector<double> magnitude(static_cast<size_t>(n), 0.0);
    for (const auto& entry : J.getEntries()) {
        if (entry.row < 0 || entry.row >= n || entry.col < 0 || entry.col >= n) continue;
        const double product = entry.value * x[entry.col];
        ax[static_cast<size_t>(entry.row)] += product;
        magnitude[static_cast<size_t>(entry.row)] += std::abs(product);
    }
    double worst = 0.0;
    for (int i = 0; i < n; ++i) {
        const double r = b[i] - ax[static_cast<size_t>(i)];
        if (residual) (*residual)[i] = r;
        const double denominator = std::abs(b[i]) + magnitude[static_cast<size_t>(i)] + 1e-300;
        worst = std::max(worst, std::abs(r) / denominator);
    }
    return worst;
}

VectorReal solve_with_refinement(
    const SparseMatrixReal& J,
    const VectorReal& b,
    LinearSolveContextReal* context,
    int refinement_steps) {
    VectorReal solution = KluSolverReal::solve(J, b, context);
    for (int pass = 0; pass < refinement_steps; ++pass) {
        VectorReal residual(b.getSize());
        const double relative = linear_relative_residual(J, b, solution, &residual);
        if (!std::isfinite(relative) || relative <= 1e-12) break;
        VectorReal correction = KluSolverReal::solve(J, residual, context);
        bool correction_finite = true;
        for (int i = 0; i < correction.getSize(); ++i) {
            if (!std::isfinite(correction[i]) || std::abs(correction[i]) > 1e9) {
                correction_finite = false;
                break;
            }
        }
        if (!correction_finite) break;
        for (int i = 0; i < solution.getSize(); ++i) solution[i] += correction[i];
    }
    return solution;
}

struct NonlinearResidualCheck {
    double error = std::numeric_limits<double>::infinity();
    bool limitingApplied = false;
    long long bypassedDevices = 0;
};

NonlinearResidualCheck dc_residual_error(
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
    const std::uint64_t epoch = next_evaluation_epoch();
    std::vector<DaeStampStatus> status(static_cast<std::size_t>(num_devs));
    #pragma omp parallel for if(use_parallel_stamp)
    for (int i = 0; i < num_devs; ++i) {
        status[static_cast<std::size_t>(i)] = stamp_device_dc(
            *devices[i], J, b, x, &settings, epoch, false, false, true);
    }
    NonlinearResidualCheck result;
    result.error = linear_system_residual_error(J, b, x, num_nodes, settings);
    for (const auto& item : status) {
        result.limitingApplied = result.limitingApplied || item.limitingApplied;
        result.bypassedDevices += item.bypassed ? 1 : 0;
    }
    return result;
}

class DaeTransientHistoryBank;
DaeStampStatus stamp_device_transient(
    Device& device,
    std::size_t deviceIndex,
    DaeTransientHistoryBank& history,
    SparseMatrixReal& jacobian,
    VectorReal& rhs,
    const VectorReal& x,
    const TransientContext& context,
    const SimulationSettings* settings = nullptr,
    std::uint64_t evaluationEpoch = 0,
    bool allowBypass = false,
    bool highPrecision = false);

NonlinearResidualCheck transient_residual_error(
    const std::vector<std::unique_ptr<Device>>& devices,
    DaeTransientHistoryBank& daeHistory,
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
    const std::uint64_t epoch = next_evaluation_epoch();
    std::vector<DaeStampStatus> status(static_cast<std::size_t>(num_devs));
    #pragma omp parallel for if(use_parallel_stamp)
    for (int i = 0; i < num_devs; ++i) {
        status[static_cast<std::size_t>(i)] = stamp_device_transient(
            *devices[i], static_cast<std::size_t>(i), daeHistory, J, b, x, ctx,
            &settings, epoch, false, true);
    }
    NonlinearResidualCheck result;
    result.error = linear_system_residual_error(J, b, x, num_nodes, settings);
    for (const auto& item : status) {
        result.limitingApplied = result.limitingApplied || item.limitingApplied;
        result.bypassedDevices += item.bypassed ? 1 : 0;
    }
    return result;
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

// Inverse-normal approximation by Peter J. Acklam. This is sufficiently
// accurate for deterministic sampling and keeps the CLI free of a statistics
// library dependency.
double standard_normal_quantile(double probability) {
    const double p = std::clamp(probability, 1e-12, 1.0 - 1e-12);
    const double a1 = -3.969683028665376e+01, a2 = 2.209460984245205e+02;
    const double a3 = -2.759285104469687e+02, a4 = 1.383577518672690e+02;
    const double a5 = -3.066479806614716e+01, a6 = 2.506628277459239e+00;
    const double b1 = -5.447609879822406e+01, b2 = 1.615858368580409e+02;
    const double b3 = -1.556989798598866e+02, b4 = 6.680131188771972e+01;
    const double b5 = -1.328068155288572e+01;
    const double c1 = -7.784894002430293e-03, c2 = -3.223964580411365e-01;
    const double c3 = -2.400758277161838e+00, c4 = -2.549732539343734e+00;
    const double c5 = 4.374664141464968e+00, c6 = 2.938163982698783e+00;
    const double d1 = 7.784695709041462e-03, d2 = 3.224671290700398e-01;
    const double d3 = 2.445134137142996e+00, d4 = 3.754408661907416e+00;
    const double pLow = 0.02425;
    if (p < pLow) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
               ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
    }
    if (p <= 1.0 - pLow) {
        const double q = p - 0.5;
        const double r = q * q;
        return (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
               (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
    }
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
            ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
}

std::vector<double> generate_mc_samples(const SimulationSettings& settings) {
    std::mt19937 rng(settings.mc_seed);
    std::vector<double> probabilities;
    probabilities.reserve(static_cast<size_t>(settings.mc_runs));
    if (settings.mc_latin_hypercube) {
        std::uniform_real_distribution<double> jitter(0.0, 1.0);
        for (int i = 0; i < settings.mc_runs; ++i) {
            probabilities.push_back((static_cast<double>(i) + jitter(rng)) /
                                    static_cast<double>(settings.mc_runs));
        }
        std::shuffle(probabilities.begin(), probabilities.end(), rng);
    }

    std::vector<double> samples;
    samples.reserve(static_cast<size_t>(settings.mc_runs));
    if (upper_copy(settings.mc_distribution) == "UNIFORM") {
        std::uniform_real_distribution<double> distribution(settings.mc_lower, settings.mc_upper);
        for (int i = 0; i < settings.mc_runs; ++i) {
            if (settings.mc_latin_hypercube) {
                const double p = probabilities[static_cast<size_t>(i)];
                samples.push_back(settings.mc_lower + p * (settings.mc_upper - settings.mc_lower));
            } else {
                samples.push_back(distribution(rng));
            }
        }
    } else {
        std::normal_distribution<double> distribution(settings.mc_mean, settings.mc_sigma);
        for (int i = 0; i < settings.mc_runs; ++i) {
            if (settings.mc_latin_hypercube) {
                samples.push_back(settings.mc_mean + settings.mc_sigma *
                    standard_normal_quantile(probabilities[static_cast<size_t>(i)]));
            } else {
                samples.push_back(settings.mc_sigma == 0.0 ? settings.mc_mean : distribution(rng));
            }
        }
    }
    return samples;
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

struct NoiseTransferResult {
    double output_psd = 0.0;
    std::size_t source_count = 0;
};

NoiseTransferResult solve_output_noise_psd(
    const std::vector<std::unique_ptr<Device>>& devices,
    int num_devs,
    int matrix_size,
    int num_nodes,
    int out_node,
    double omega,
    const VectorReal& operating_point,
    const SimulationSettings& settings,
    LinearSolveContextComplex* complex_solver_context) {
    SparseMatrixComplex J_sparse(matrix_size);
    VectorComplex ignored_rhs(matrix_size);
    stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
    const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
    #pragma omp parallel for if(use_parallel_stamp)
    for (int i = 0; i < num_devs; ++i) {
        stamp_device_ac(*devices[i], J_sparse, ignored_rhs, omega, operating_point);
    }

    std::vector<NoiseSource> noise_sources;
    for (const auto& dev : devices) {
        dev->collectNoiseSources(omega, operating_point, noise_sources);
    }

    double output_psd = 0.0;
    for (const auto& source : noise_sources) {
        if (source.currentPsd <= 0.0) continue;
        VectorComplex rhs(matrix_size);
        rhs.add(source.nodePos, {-1.0, 0.0});
        rhs.add(source.nodeNeg, {1.0, 0.0});
        VectorComplex transfer = KluSolverComplex::solve(J_sparse, rhs, complex_solver_context);
        output_psd += std::norm(transfer[out_node]) * source.currentPsd;
    }
    return {output_psd, noise_sources.size()};
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
    double update_error = std::numeric_limits<double>::infinity();
    int update_index = -1;
    int integration_order = 1;
    long long bypassed_devices = 0;
    bool limiting_prevented_convergence = false;
};

struct TransientMatrixWorkspace {
    explicit TransientMatrixWorkspace(int matrix_size)
        : jacobian(matrix_size), rhs(matrix_size) {}

    void setStructureCacheEnabled(bool enabled) {
        jacobian.setStructureCacheEnabled(enabled);
    }

    SparseMatrixReal jacobian;
    VectorReal rhs;
};

std::string transient_failure_message(
    const std::string& label,
    double time,
    double step,
    const TransientStepResult& result,
    int num_nodes) {
    std::ostringstream stream;
    stream << label << " at time=" << std::scientific << std::setprecision(9) << time
           << " step=" << step
           << " iterations=" << result.iterations
           << " update_error=" << result.update_error
           << " residual_error=" << result.residual_error;
    if (result.update_index >= 0) {
        stream << " worst_unknown="
               << (result.update_index < num_nodes ? "node[" : "branch[")
               << (result.update_index < num_nodes ? result.update_index : result.update_index - num_nodes)
               << "]";
    }
    return stream.str();
}

struct TransientStats {
    long long accepted_steps = 0;
    long long rejected_steps = 0;
    long long convergence_rejections = 0;
    long long lte_rejections = 0;
    long long output_points = 0;
    long long newton_solves = 0;
    long long newton_iterations = 0;
    long long bound_step_limited = 0;
    long long bypassed_devices = 0;
    long long limiting_rechecks = 0;
    long long predictor_lte_steps = 0;
    long long step_doubling_audits = 0;
    long long method_switches = 0;
    int max_newton_iterations = 0;
    int max_integration_order = 1;
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
        max_integration_order = std::max(max_integration_order, result.integration_order);
        bypassed_devices += result.bypassed_devices;
        limiting_rechecks += result.limiting_prevented_convergence ? 1 : 0;
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

class DeviceTransientStateArena {
public:
    using Checkpoint = OpaqueTransientStateStore::Checkpoint;

    explicit DeviceTransientStateArena(const std::vector<std::unique_ptr<Device>>& devices)
        : devices_(devices) {
        OpaqueTransientStateLayout layout;
        blocks_.reserve(devices_.size());
        for (const auto& device : devices_) {
            blocks_.push_back(layout.allocate(device->transientStateBytes()));
        }
        layout.seal();
        store_ = std::make_unique<OpaqueTransientStateStore>(layout, 3, 2);
        scratch_.resize(layout.totalSize());
        for (std::size_t i = 0; i < devices_.size(); ++i) {
            devices_[i]->saveTransientStateBytes(store_->initial(blocks_[i]), blocks_[i].size);
        }
    }

    Checkpoint checkpoint() const { return store_->checkpoint(); }

    void restore(const Checkpoint& checkpoint) {
        store_->rollback(checkpoint);
        restoreCurrent();
    }

    void restoreCurrent() {
        for (std::size_t i = 0; i < devices_.size(); ++i) {
            devices_[i]->restoreTransientStateBytes(store_->accepted(blocks_[i]), blocks_[i].size);
        }
    }

    void commitDeviceState() {
        store_->prepareCandidate();
        for (std::size_t i = 0; i < devices_.size(); ++i) {
            devices_[i]->saveTransientStateBytes(store_->candidate(blocks_[i]), blocks_[i].size);
        }
        store_->acceptCandidate();
    }

    void preserveWorkingState() {
        std::byte* scratchBase = scratch_.empty() ? nullptr : scratch_.data();
        for (std::size_t i = 0; i < devices_.size(); ++i) {
            devices_[i]->saveTransientStateBytes(
                scratchBase ? scratchBase + blocks_[i].offset : nullptr, blocks_[i].size);
        }
        scratchValid_ = true;
    }

    void restoreWorkingState() {
        if (!scratchValid_) throw std::logic_error("no preserved transient working state");
        const std::byte* scratchBase = scratch_.empty() ? nullptr : scratch_.data();
        for (std::size_t i = 0; i < devices_.size(); ++i) {
            devices_[i]->restoreTransientStateBytes(
                scratchBase ? scratchBase + blocks_[i].offset : nullptr, blocks_[i].size);
        }
    }

private:
    const std::vector<std::unique_ptr<Device>>& devices_;
    std::vector<OpaqueTransientStateBlock> blocks_;
    std::unique_ptr<OpaqueTransientStateStore> store_;
    std::vector<std::byte> scratch_;
    bool scratchValid_ = false;
};

struct DaeDerivativeAssembly {
    double leading = 0.0;
    DaeHistory known;
};

class DaeTransientHistoryBank {
public:
    struct Checkpoint {
        std::size_t head = 0;
        std::size_t acceptedCount = 0;
        std::uint64_t serial = 0;
    };

    DaeTransientHistoryBank(
        const std::vector<std::unique_ptr<Device>>& devices,
        const VectorReal& initialSolution,
        double initialTime,
        std::size_t historyDepth = 6,
        std::size_t speculativeDepth = 2)
        : historyDepth_(historyDepth), speculativeDepth_(speculativeDepth),
          frames_(historyDepth + speculativeDepth) {
        if (historyDepth == 0 || speculativeDepth == 0) {
            throw std::invalid_argument("DAE history depths must be positive");
        }
        for (auto& frame : frames_) {
            frame.q.resize(devices.size());
            frame.qdot.resize(devices.size());
            frame.active.assign(devices.size(), false);
            frame.qdotValid.assign(devices.size(), false);
        }
        Frame& initial = frames_[head_];
        for (std::size_t i = 0; i < devices.size(); ++i) {
            DaeRequest request;
            request.analysis = DaeAnalysis::Transient;
            request.time = initialTime;
            request.staticResidual = false;
            request.staticJacobian = false;
            request.dynamicResidual = true;
            request.dynamicJacobian = false;
            request.readOnlyState = true;
            request.enableLimiting = false;
            DaeEvaluation evaluation;
            initial.active[i] = devices[i]->evaluateDae(initialSolution, request, evaluation);
            if (initial.active[i]) {
                initial.q[i] = std::move(evaluation.dynamicResidual);
                initial.qdot[i] = zeroLike(initial.q[i]);
            }
        }
    }

    Checkpoint checkpoint() const { return {head_, acceptedCount_, serial_}; }

    void rollback(const Checkpoint& checkpoint) {
        if (checkpoint.serial > serial_ || serial_ - checkpoint.serial > speculativeDepth_) {
            throw std::logic_error("DAE history checkpoint is outside the rollback window");
        }
        head_ = checkpoint.head;
        acceptedCount_ = checkpoint.acceptedCount;
        serial_ = checkpoint.serial;
    }

    bool active(std::size_t deviceIndex) const {
        return deviceIndex < frames_[head_].active.size() && frames_[head_].active[deviceIndex];
    }

    DaeDerivativeAssembly derivativeAssembly(
        std::size_t deviceIndex,
        const TransientContext& context) const {
        DaeDerivativeAssembly assembly;
        const Frame& current = frames_[head_];
        bool derivativeHistoryAvailable = true;
        for (std::size_t age = 0; age < context.derivativeWeights.size(); ++age) {
            if (age > acceptedCount_ || !frameAtAge(age).qdotValid[deviceIndex]) {
                derivativeHistoryAvailable = false;
                break;
            }
        }
        if (!derivativeHistoryAvailable) {
            assembly.leading = 1.0 / std::max(context.timeStep, 1e-30);
            appendScaledDaeResidual(assembly.known, current.q[deviceIndex], -assembly.leading);
            return assembly;
        }
        if (!context.qWeights.empty()) {
            assembly.leading = context.qWeights.front();
            for (std::size_t weight = 1; weight < context.qWeights.size(); ++weight) {
                const std::size_t age = weight - 1;
                if (age > acceptedCount_) break;
                appendScaledDaeResidual(
                    assembly.known, frameAtAge(age).q[deviceIndex], context.qWeights[weight]);
            }
            for (std::size_t age = 0; age < context.derivativeWeights.size(); ++age) {
                appendScaledDaeResidual(
                    assembly.known, frameAtAge(age).qdot[deviceIndex],
                    context.derivativeWeights[age]);
            }
            return assembly;
        }
        assembly.leading = context.a0;
        appendScaledDaeResidual(assembly.known, current.q[deviceIndex], context.a1);
        return assembly;
    }

    void commitAccepted(
        const std::vector<std::unique_ptr<Device>>& devices,
        const VectorReal& solution,
        double time,
        const TransientContext& context) {
        const std::size_t candidate = candidateIndex();
        frames_[candidate] = frames_[head_];
        Frame& next = frames_[candidate];
        for (std::size_t i = 0; i < devices.size(); ++i) {
            if (!frames_[head_].active[i]) continue;
            DaeRequest request;
            request.analysis = DaeAnalysis::Transient;
            request.time = time;
            request.staticResidual = false;
            request.staticJacobian = false;
            request.dynamicResidual = true;
            request.dynamicJacobian = false;
            request.readOnlyState = true;
            request.enableLimiting = false;
            DaeEvaluation evaluation;
            if (!devices[i]->evaluateDae(solution, request, evaluation)) {
                next.active[i] = false;
                continue;
            }
            const auto derivative = derivativeAssembly(i, context);
            next.q[i] = std::move(evaluation.dynamicResidual);
            next.qdot[i] = combineDerivative(next.q[i], derivative.leading, derivative.known);
            next.qdotValid[i] = true;
        }
        head_ = candidate;
        acceptedCount_ = std::min(acceptedCount_ + 1, historyDepth_);
        ++serial_;
    }

    void restart(
        const std::vector<std::unique_ptr<Device>>& devices,
        const VectorReal& solution,
        double time) {
        Frame& current = frames_[head_];
        for (std::size_t i = 0; i < devices.size(); ++i) {
            DaeRequest request;
            request.analysis = DaeAnalysis::Transient;
            request.time = time;
            request.staticResidual = false;
            request.staticJacobian = false;
            request.dynamicResidual = true;
            request.dynamicJacobian = false;
            request.readOnlyState = true;
            request.enableLimiting = false;
            DaeEvaluation evaluation;
            current.active[i] = devices[i]->evaluateDae(solution, request, evaluation);
            if (current.active[i]) {
                current.q[i] = std::move(evaluation.dynamicResidual);
                current.qdot[i] = zeroLike(current.q[i]);
            } else {
                current.q[i].clear();
                current.qdot[i].clear();
            }
            current.qdotValid[i] = false;
        }
        acceptedCount_ = 0;
    }

private:
    struct Frame {
        std::vector<DaeHistory> q;
        std::vector<DaeHistory> qdot;
        std::vector<bool> active;
        std::vector<bool> qdotValid;
    };

    static DaeHistory zeroLike(const DaeHistory& source) {
        DaeHistory result;
        result.reserve(source.size());
        for (const auto& term : source) {
            result.push_back({term.equation, 0.0, term.conservationGroup});
        }
        return result;
    }

    static DaeHistory combineDerivative(
        const DaeHistory& current,
        double leading,
        const DaeHistory& known) {
        std::map<int, DaeResidualTerm> combined;
        for (const auto& term : current) {
            auto& target = combined[term.equation];
            target.equation = term.equation;
            target.value += leading * term.value;
            target.conservationGroup = term.conservationGroup;
        }
        for (const auto& term : known) {
            auto& target = combined[term.equation];
            target.equation = term.equation;
            target.value += term.value;
            if (target.conservationGroup < 0) target.conservationGroup = term.conservationGroup;
        }
        DaeHistory result;
        result.reserve(combined.size());
        for (const auto& entry : combined) result.push_back(entry.second);
        return result;
    }

    const Frame& frameAtAge(std::size_t age) const {
        if (age >= historyDepth_ || age > acceptedCount_) {
            throw std::out_of_range("DAE history age is unavailable");
        }
        return frames_[(head_ + frames_.size() - age % frames_.size()) % frames_.size()];
    }

    std::size_t candidateIndex() const { return (head_ + 1) % frames_.size(); }

    std::size_t historyDepth_ = 0;
    std::size_t speculativeDepth_ = 0;
    std::vector<Frame> frames_;
    std::size_t head_ = 0;
    std::size_t acceptedCount_ = 0;
    std::uint64_t serial_ = 0;
};

DaeStampStatus stamp_device_transient(
    Device& device,
    std::size_t deviceIndex,
    DaeTransientHistoryBank& history,
    SparseMatrixReal& jacobian,
    VectorReal& rhs,
    const VectorReal& x,
    const TransientContext& context,
    const SimulationSettings* settings,
    std::uint64_t evaluationEpoch,
    bool allowBypass,
    bool highPrecision) {
    if (!history.active(deviceIndex)) {
        device.tranStamp(jacobian, rhs, x, context);
        return {};
    }
    DaeRequest request;
    request.analysis = DaeAnalysis::Transient;
    request.time = context.currentTime;
    request.staticResidual = true;
    request.staticJacobian = true;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    request.highPrecision = highPrecision;
    request.allowBypass = allowBypass && settings && settings->nr_bypass;
    request.bypassRelativeTolerance = settings
        ? settings->reltol * settings->nr_bypass_tolerance
        : 0.0;
    request.bypassAbsoluteTolerance = settings
        ? std::max(settings->vntol, settings->abstol) * settings->nr_bypass_tolerance
        : 0.0;
    request.evaluationEpoch = evaluationEpoch;
    DaeEvaluation evaluation;
    if (!device.evaluateDae(x, request, evaluation)) {
        device.tranStamp(jacobian, rhs, x, context);
        return {};
    }
    const auto derivative = history.derivativeAssembly(deviceIndex, context);
    stampDaeTransient(evaluation, x, derivative.leading, derivative.known, jacobian, rhs);
    return {evaluation.limitingApplied, evaluation.bypassed};
}

void accept_device_transient_step(
    const std::vector<std::unique_ptr<Device>>& devices,
    const VectorReal& x,
    double time,
    const TransientContext& ctx) {
    for (const auto& device : devices) {
        device->acceptTransientStep(x, time, ctx);
    }
}

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

void run_dae_audits(
    const std::vector<std::unique_ptr<Device>>& devices,
    const VectorReal& operatingPoint,
    const SimulationSettings& settings) {
    if (!settings.dae_audit) return;

    DaeAuditOptions options;
    options.relativeTolerance = settings.dae_audit_tolerance;
    std::size_t checked = 0;
    std::size_t skipped = 0;
    double worstStatic = 0.0;
    double worstDynamic = 0.0;
    double worstCharge = 0.0;
    double worstChargeJacobian = 0.0;
    std::vector<std::string> failures;
    for (const auto& device : devices) {
        if (!device->daeAuditSafe()) {
            ++skipped;
            continue;
        }
        const DaeAuditReport report = auditDaeDevice(*device, operatingPoint, options);
        ++checked;
        worstStatic = std::max(worstStatic, report.worstStaticDerivative);
        worstDynamic = std::max(worstDynamic, report.worstDynamicDerivative);
        worstCharge = std::max(worstCharge, report.worstChargeImbalance);
        worstChargeJacobian = std::max(worstChargeJacobian, report.worstChargeJacobianImbalance);
        if (!report.passed()) {
            std::ostringstream failure;
            failure << device->getName()
                    << " F(eq=" << report.worstStaticEquation
                    << ",x=" << report.worstStaticUnknown
                    << ",num=" << report.worstStaticNumerical
                    << ",jac=" << report.worstStaticExpected << ')'
                    << " Q(eq=" << report.worstDynamicEquation
                    << ",x=" << report.worstDynamicUnknown
                    << ",num=" << report.worstDynamicNumerical
                    << ",jac=" << report.worstDynamicExpected << ')';
            failures.push_back(failure.str());
        }
    }
    std::cout << "DAE audit: checked=" << checked
              << " skipped=" << skipped
              << " worst_F=" << worstStatic
              << " worst_Q=" << worstDynamic
              << " charge=" << worstCharge
              << " charge_J=" << worstChargeJacobian
              << (failures.empty() ? " PASS" : " FAIL") << std::endl;
    if (!failures.empty()) {
        std::ostringstream message;
        message << "DAE derivative/conservation audit failed for";
        for (const auto& name : failures) message << ' ' << name;
        throw std::runtime_error(message.str());
    }
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
    const std::vector<double>& t_hist,
    double proposed_step = 0.0,
    bool autoUseTrapezoidal = true) {
    const std::string method = normalized_method_key(settings);
    if (settings.tran_max_order <= 1) return TransientIntegrationMethod::BackwardEuler;
    if (method == "BDF" || method == "GEAR") return TransientIntegrationMethod::Bdf;
    if (method == "ADAMS" || method == "ADAMSMOULTON" || method == "AM") {
        return TransientIntegrationMethod::AdamsMoulton;
    }
    if (method == "GEAR2" && x_hist.size() >= 2 && t_hist.size() >= 2) {
        return TransientIntegrationMethod::Gear2;
    }
    if ((method == "AUTO" || method == "TRAPGEAR" || method == "TRAP_GEAR") &&
        settings.tran_adaptive && x_hist.size() >= 2 && t_hist.size() >= 2) {
        (void)proposed_step;
        return autoUseTrapezoidal
            ? TransientIntegrationMethod::Trapezoidal
            : TransientIntegrationMethod::Gear2;
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
    double target_time,
    const TransientIntegrationMethod* forced_method = nullptr,
    int forced_order = 0) {
    TransientContext ctx;
    ctx.timeStep = std::max(step, 1e-30);
    ctx.currentTime = target_time;
    ctx.method = forced_method ? *forced_method : choose_transient_method(settings, x_hist, t_hist, step);
    ctx.xHistory = &x_hist;
    ctx.timeHistory = &t_hist;

    auto newestTimes = [&](std::size_t pastCount) {
        std::vector<double> times;
        times.reserve(pastCount + 1);
        times.push_back(target_time);
        for (std::size_t age = 0; age < pastCount; ++age) {
            times.push_back(t_hist[t_hist.size() - 1 - age]);
        }
        return times;
    };

    const int requestedOrder = std::clamp(
        forced_order > 0 ? forced_order : settings.tran_max_order,
        1,
        settings.tran_max_order);
    IntegrationFormula formula;
    if (ctx.method == TransientIntegrationMethod::Trapezoidal) {
        formula = makeAdamsMoultonFormula(newestTimes(1));
    } else if (ctx.method == TransientIntegrationMethod::AdamsMoulton) {
        const int order = std::min<int>(requestedOrder, static_cast<int>(t_hist.size()) + 1);
        if (order <= 1) {
            formula = makeBdfFormula(newestTimes(1));
        } else {
            formula = makeAdamsMoultonFormula(newestTimes(static_cast<std::size_t>(order - 1)));
        }
    } else {
        std::size_t order = 1;
        if (ctx.method == TransientIntegrationMethod::Gear2) order = std::min<std::size_t>(2, t_hist.size());
        if (ctx.method == TransientIntegrationMethod::Bdf) {
            order = std::min<std::size_t>(
                static_cast<std::size_t>(requestedOrder), t_hist.size());
        }
        formula = makeBdfFormula(newestTimes(std::max<std::size_t>(1, order)));
    }
    ctx.integrationOrder = formula.order;
    ctx.qWeights = formula.qWeights;
    ctx.derivativeWeights = formula.derivativeWeights;
    ctx.a0 = ctx.qWeights.empty() ? 1.0 / ctx.timeStep : ctx.qWeights[0];
    ctx.a1 = ctx.qWeights.size() > 1 ? ctx.qWeights[1] : -ctx.a0;
    ctx.a2 = ctx.qWeights.size() > 2 ? ctx.qWeights[2] : 0.0;
    ctx.hasSecondHistory = ctx.qWeights.size() > 2;

    return ctx;
}

VectorReal transient_initial_guess(
    const VectorReal& fallback,
    const std::vector<VectorReal>& x_hist,
    const std::vector<double>& t_hist,
    double target_time,
    const SimulationSettings& settings,
    int integrationOrder) {
    if (!settings.tran_predictor || x_hist.size() < 2 || t_hist.size() < 2) {
        return fallback;
    }
    const auto polynomial = polynomialPredict(
        x_hist, t_hist, target_time, std::max(1, integrationOrder));
    if (polynomial.valid && vector_finite_and_bounded(polynomial.value)) {
        return polynomial.value;
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
    if (method == TransientIntegrationMethod::Bdf) return "variable-step BDF";
    if (method == TransientIntegrationMethod::AdamsMoulton) return "variable-step Adams-Moulton";
    return "Backward Euler";
}

int transient_method_order(TransientIntegrationMethod method) {
    return method == TransientIntegrationMethod::BackwardEuler ? 1 : 2;
}

double transient_reject_factor(double err, int order) {
    const double exponent = -1.0 / static_cast<double>(order + 1);
    return std::clamp(0.85 * std::pow(std::max(err, 1e-30), exponent), 0.1, 0.75);
}

double transient_growth_factor(double err, int order) {
    if (err <= 1e-12) return order <= 1 ? 4.0 : 8.0;
    if (err <= 1e-6)  return order <= 1 ? 2.5 : 5.0;
    const double exponent = -1.0 / static_cast<double>(order + 1);
    return std::clamp(0.9 * std::pow(std::max(err, 1e-30), exponent), 1.05, 5.0);
}

TransientStepResult solve_transient_step(
    const std::vector<std::unique_ptr<Device>>& devices,
    DaeTransientHistoryBank& daeHistory,
    int num_devs,
    int matrix_size,
    int num_nodes,
    const SimulationSettings& settings,
    LinearSolveContextReal* solver_context,
    const VectorReal& x_start,
    const std::vector<VectorReal>& x_hist,
    const std::vector<double>& t_hist,
    double step,
    double target_time,
    const TransientIntegrationMethod* forced_method = nullptr,
    int forced_order = 0,
    TransientMatrixWorkspace* matrix_workspace = nullptr) {
    bool converged = false;
    double total_stamp_seconds = 0.0;
    double total_solve_seconds = 0.0;
    const TransientContext tran_ctx = make_transient_context(
        settings, x_hist, t_hist, step, target_time, forced_method, forced_order);
    VectorReal x = transient_initial_guess(
        x_start, x_hist, t_hist, target_time, settings, tran_ctx.integrationOrder);
    const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
    const std::uint64_t evaluationEpoch = next_evaluation_epoch();
    double previous_delta = std::numeric_limits<double>::infinity();
    double last_update_error = std::numeric_limits<double>::infinity();
    int last_update_index = -1;
    long long bypassedDevices = 0;
    bool limitingPreventedConvergence = false;
    // FastSPICE Suite: FastSpiceEngine and MultiRateController integration
    FastSpiceEngine fastspice_eng({settings.fastspice});
    fastspice_eng.initialize(static_cast<std::size_t>(num_devs));
    MultiRateController multirate_ctrl({settings.multirate});
    multirate_ctrl.initialize(num_nodes);
    if (settings.multirate) {
        multirate_ctrl.observe(x, target_time);
    }
    TransientMatrixWorkspace local_matrix_workspace(matrix_size);
    TransientMatrixWorkspace& workspace = matrix_workspace ? *matrix_workspace : local_matrix_workspace;
    workspace.setStructureCacheEnabled(settings.transient_stamp_cache);

    for (int iter = 0; iter < settings.tran_max_iter; ++iter) {
        const auto stamp_start = std::chrono::steady_clock::now();
        workspace.jacobian.clear();
        workspace.rhs.clear();
        SparseMatrixReal& J_sparse = workspace.jacobian;
        VectorReal& b = workspace.rhs;
        stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
        std::vector<DaeStampStatus> stampStatus(static_cast<std::size_t>(num_devs));
        #pragma omp parallel for if(use_parallel_stamp)
        for (int i = 0; i < num_devs; ++i) {
            // Check MultiRate / FastSPICE latency bypass
            if (settings.multirate && !multirate_ctrl.shouldEvaluateNode(i % num_nodes)) {
                stampStatus[static_cast<std::size_t>(i)].bypassed = true;
                continue;
            }
            stampStatus[static_cast<std::size_t>(i)] = stamp_device_transient(
                *devices[i], static_cast<std::size_t>(i), daeHistory, J_sparse, b, x, tran_ctx,
                &settings, evaluationEpoch, iter > 0, false);
        }
        for (const auto& item : stampStatus) {
            bypassedDevices += item.bypassed ? 1 : 0;
        }
        const auto stamp_end = std::chrono::steady_clock::now();
        const auto solve_start = std::chrono::steady_clock::now();
        VectorReal x_new(matrix_size);
        try {
            x_new = solve_with_refinement(
                J_sparse, b, solver_context, settings.solver_refinement_steps);
        } catch (...) {
            TransientStepResult failed{x, false, iter + 1, total_stamp_seconds, total_solve_seconds,
                std::numeric_limits<double>::infinity(), last_update_error,
                last_update_index, tran_ctx.integrationOrder};
            failed.bypassed_devices = bypassedDevices;
            return failed;
        }
        const auto solve_end = std::chrono::steady_clock::now();
        total_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
        total_solve_seconds += elapsed_seconds(solve_start, solve_end);
        if (!vector_finite_and_bounded(x_new)) {
            TransientStepResult failed{x_new, false, iter + 1, total_stamp_seconds, total_solve_seconds,
                std::numeric_limits<double>::infinity(), last_update_error,
                last_update_index, tran_ctx.integrationOrder};
            failed.bypassed_devices = bypassedDevices;
            return failed;
        }
        for (const auto& device : devices) device->limitTransientNewton(x, x_new);
        if (settings.line_search && iter > 0) {
            double delta = vector_max_delta(x_new, x);
            double alpha = 1.0;
            while (delta > previous_delta * 2.0 && alpha > 0.0625) {
                alpha *= 0.5;
                x_new = vector_blend(x, x_new, alpha);
                delta = vector_max_delta(x_new, x);
            }
        }
        const auto update = worst_solution_update(x_new, x, num_nodes, settings);
        last_update_error = update.first;
        last_update_index = update.second;
        double residual_error = std::numeric_limits<double>::infinity();
        const bool delta_converged = solution_converged(x_new, x, num_nodes, settings);
        if (delta_converged && settings.nr_residual_check) {
            const auto residual = transient_residual_error(
                devices, daeHistory, num_devs, matrix_size, num_nodes, x_new, tran_ctx, settings);
            residual_error = residual.error;
            limitingPreventedConvergence = residual.limitingApplied;
            converged = residual_error <= 1.0 && !residual.limitingApplied;
        } else {
            converged = delta_converged;
        }
        previous_delta = vector_max_delta(x_new, x);
        x = x_new;
        if (converged) {
            TransientStepResult success{x, true, iter + 1, total_stamp_seconds, total_solve_seconds,
                residual_error, last_update_error, last_update_index, tran_ctx.integrationOrder};
            success.bypassed_devices = bypassedDevices;
            success.limiting_prevented_convergence = limitingPreventedConvergence;
            return success;
        }
    }
    const auto final_check = settings.nr_residual_check
        ? transient_residual_error(devices, daeHistory, num_devs, matrix_size, num_nodes, x, tran_ctx, settings)
        : NonlinearResidualCheck{};
    TransientStepResult failed{x, converged, settings.tran_max_iter, total_stamp_seconds, total_solve_seconds,
        final_check.error, last_update_error, last_update_index, tran_ctx.integrationOrder};
    failed.bypassed_devices = bypassedDevices;
    failed.limiting_prevented_convergence = final_check.limitingApplied;
    return failed;
}

double transient_endpoint_error(
    const std::vector<std::unique_ptr<Device>>& devices,
    DeviceTransientStateArena& stateArena,
    const VectorReal& first,
    const VectorReal& second,
    const std::vector<VectorReal>& xHistory,
    int num_nodes,
    const SimulationSettings& settings,
    double errorScale,
    double targetTime) {
    double worst = 0.0;
    const int n = std::min(first.getSize(), second.getSize());
    const double error_scale = std::max(errorScale, 1e-12);
    const std::string lte_reference = upper_copy(settings.tran_lte_reference);
    double global_reference = 0.0;
    if (lte_reference == "GLOBAL" || lte_reference == "HISTORY") {
        for (int i = 0; i < n && i < num_nodes; ++i) {
            global_reference = std::max(global_reference, std::abs(first[i]));
            global_reference = std::max(global_reference, std::abs(second[i]));
        }
        if (lte_reference == "HISTORY") {
            for (const auto& sample : xHistory) {
                const int m = std::min(sample.getSize(), num_nodes);
                for (int i = 0; i < m; ++i) {
                    global_reference = std::max(global_reference, std::abs(sample[i]));
                }
            }
        }
    }
    for (int i = 0; i < n && i < num_nodes; ++i) {
        const double err = std::abs(second[i] - first[i]);
        double scale = std::max(std::abs(first[i]), std::abs(second[i]));
        if (lte_reference == "GLOBAL" || lte_reference == "HISTORY") {
            scale = std::max(scale, global_reference);
        } else if (lte_reference == "SIGNAL_HISTORY") {
            for (const auto& sample : xHistory) {
                if (i < sample.getSize()) scale = std::max(scale, std::abs(sample[i]));
            }
        }
        const double absolute = settings.tran_lte_abstol;
        const double relative = settings.tran_lte_reltol;
        const double tol = absolute + relative * scale;
        worst = std::max(worst, err / std::max(tol * error_scale, 1e-30));
    }
    stateArena.preserveWorkingState();
    for (const auto& device : devices) {
        DaeRequest request;
        request.analysis = DaeAnalysis::Transient;
        request.staticResidual = false;
        request.staticJacobian = false;
        request.dynamicResidual = true;
        request.dynamicJacobian = false;
        request.readOnlyState = true;
        request.enableLimiting = false;
        request.time = targetTime;
        DaeEvaluation coarseEvaluation;
        DaeEvaluation fineEvaluation;
        stateArena.restoreWorkingState();
        const bool coarseDae = device->evaluateDae(first, request, coarseEvaluation);
        stateArena.restoreWorkingState();
        const bool fineDae = device->evaluateDae(second, request, fineEvaluation);
        double chargeError = 0.0;
        if (coarseDae && fineDae) {
            std::map<int, double> coarseCharge;
            std::map<int, double> fineCharge;
            for (const auto& term : coarseEvaluation.dynamicResidual) {
                coarseCharge[term.equation] += term.value;
            }
            for (const auto& term : fineEvaluation.dynamicResidual) {
                fineCharge[term.equation] += term.value;
            }
            std::map<int, bool> equations;
            for (const auto& term : coarseCharge) equations[term.first] = true;
            for (const auto& term : fineCharge) equations[term.first] = true;
            for (const auto& equation : equations) {
                const double coarse = coarseCharge[equation.first];
                const double fine = fineCharge[equation.first];
                const double tolerance = settings.chgtol + settings.tran_lte_reltol *
                    std::max(std::abs(coarse), std::abs(fine));
                chargeError = std::max(
                    chargeError, std::abs(fine - coarse) / std::max(tolerance, 1e-30));
            }
        } else {
            stateArena.restoreWorkingState();
            chargeError = device->transientChargeError(
                first, second, settings.tran_lte_reltol, settings.chgtol);
        }
        worst = std::max(worst, chargeError / error_scale);
    }
    stateArena.restoreWorkingState();
    return worst;
}

double transient_lte_error(
    const std::vector<std::unique_ptr<Device>>& devices,
    DeviceTransientStateArena& stateArena,
    const VectorReal& xFull,
    const VectorReal& xHalf,
    const std::vector<VectorReal>& xHistory,
    int numNodes,
    const SimulationSettings& settings,
    int integrationOrder,
    double targetTime) {
    const double richardson = std::max(
        std::pow(2.0, static_cast<double>(integrationOrder)) - 1.0, 1.0);
    return transient_endpoint_error(
        devices, stateArena, xFull, xHalf, xHistory, numNodes, settings,
        richardson * std::max(settings.tran_trtol, 1e-12), targetTime);
}

struct PredictorCorrectorEstimate {
    bool valid = false;
    double error = std::numeric_limits<double>::infinity();
    double factor = 0.0;
    VectorReal prediction;
};

PredictorCorrectorEstimate predictor_corrector_lte_error(
    const std::vector<std::unique_ptr<Device>>& devices,
    DeviceTransientStateArena& stateArena,
    const VectorReal& corrected,
    const std::vector<VectorReal>& xHistory,
    const std::vector<double>& timeHistory,
    const TransientContext& context,
    int numNodes,
    const SimulationSettings& settings) {
    PredictorCorrectorEstimate estimate;
    const auto prediction = polynomialPredict(
        xHistory, timeHistory, context.currentTime, context.integrationOrder);
    if (!prediction.valid) return estimate;

    IntegrationFormula formula;
    formula.family = context.method == TransientIntegrationMethod::AdamsMoulton ||
            context.method == TransientIntegrationMethod::Trapezoidal
        ? IntegrationFamily::AdamsMoulton
        : IntegrationFamily::Bdf;
    formula.order = context.integrationOrder;
    formula.qWeights = context.qWeights;
    formula.derivativeWeights = context.derivativeWeights;
    const std::size_t requiredTimes = std::max(
        formula.qWeights.size(), formula.derivativeWeights.size() + 1);
    if (requiredTimes < 2 || timeHistory.size() + 1 < requiredTimes) return estimate;
    std::vector<double> newestFirstTimes;
    newestFirstTimes.reserve(requiredTimes);
    newestFirstTimes.push_back(context.currentTime);
    for (std::size_t age = 0; age + 1 < requiredTimes; ++age) {
        newestFirstTimes.push_back(timeHistory[timeHistory.size() - 1 - age]);
    }
    estimate.factor = predictorCorrectorErrorFactor(
        formula, newestFirstTimes, timeHistory, context.currentTime);
    if (!(estimate.factor > 0.0) || !std::isfinite(estimate.factor)) return estimate;
    estimate.prediction = prediction.value;
    estimate.error = transient_endpoint_error(
        devices, stateArena, prediction.value, corrected, xHistory, numNodes, settings,
        std::max(settings.tran_trtol, 1e-12) / estimate.factor,
        context.currentTime);
    estimate.valid = std::isfinite(estimate.error);
    return estimate;
}

void print_transient_point(double t, const VectorReal& x, int num_nodes) {
    std::cout << std::scientific << std::setprecision(9) << t << " | ";
    for (int i = 0; i < num_nodes; ++i) std::cout << x[i] << " ";
    std::cout << std::endl;
}

struct SavedOutputSignal {
    std::string label;
    std::string type = "voltage";
    int node_pos = -1;
    int node_neg = -1;
    int branch_index = -1;
};

int find_branch_index_by_name(
    const std::vector<std::unique_ptr<Device>>& devices,
    const std::string& device_name);

std::string branch_save_name(std::string name) {
    const size_t dot = name.find('.');
    if (dot != std::string::npos) name.resize(dot);
    return name;
}

std::vector<SavedOutputSignal> resolve_saved_outputs(const Netlist& netlist, int num_nodes) {
    const auto& settings = netlist.getSettings();
    std::vector<SavedOutputSignal> signals;
    if (settings.save_none) {
        return signals;
    }
    if (settings.save_all) {
        signals.reserve(static_cast<size_t>(num_nodes + settings.saves.size()));
        for (int i = 0; i < num_nodes; ++i) {
            signals.push_back({"V(" + netlist.getNodeName(i) + ")", "voltage", i, -1});
        }
    } else if (settings.saves.empty()) {
        return signals;
    }

    for (const auto& save : settings.saves) {
        if (upper_copy(save.kind) == "I") {
            const std::string source_name = branch_save_name(save.node_pos);
            const int branch = find_branch_index_by_name(netlist.getDevices(), source_name);
            if (branch < 0) {
                std::cout << "WARNING: save skipped unknown current in I(" << save.node_pos << ")" << std::endl;
                continue;
            }
            signals.push_back({"I(" + save.node_pos + ")", "current", -1, -1, branch});
            continue;
        }
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
        signals.push_back({label, "voltage", pos, neg});
    }
    if (signals.empty()) {
        std::cout << "WARNING: no valid .SAVE signals were resolved; writing time only." << std::endl;
    }
    return signals;
}

class TransientOutput {
public:
    TransientOutput(const std::string& path, const std::string& format, const Netlist& netlist, int num_nodes)
        : to_file_(!path.empty()), csv_(upper_copy(format) == "CSV"),
          signals_(resolve_saved_outputs(netlist, num_nodes)) {
        if (to_file_) {
            file_.open(path, std::ios::out | std::ios::trunc);
            if (!file_.is_open()) {
                throw std::runtime_error("Could not open transient output file: " + path);
            }
            if (csv_) {
                file_ << "time";
                for (const auto& signal : signals_) file_ << ",\"" << signal.label << "\"";
                file_ << "\n";
                return;
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
                file_ << (i + 1) << "\t" << signals_[i].label << "\t" << signals_[i].type << "\n";
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
            double value = probe_value(x, signal.node_pos, signal.node_neg);
            if (signal.type == "current") {
                value = x[signal.branch_index];
            }
            file_ << (csv_ ? "," : " ") << value;
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
        if (csv_ || !file_.is_open() || points_pos_ < std::streampos(0)) {
            if (file_.is_open()) file_.flush();
            return;
        }
        const auto end_pos = file_.tellp();
        file_.seekp(points_pos_);
        file_ << std::setw(20) << point_count_;
        file_.seekp(end_pos);
        file_.flush();
    }

    bool to_file_ = false;
    bool csv_ = false;
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

bool is_transient_breakpoint(const std::vector<double>& points, double t, double tol) {
    const auto candidate = std::lower_bound(points.begin(), points.end(), t - tol);
    return candidate != points.end() && std::abs(*candidate - t) <= tol;
}

std::string normalized_transient_method(
    const SimulationSettings& settings,
    bool autoUseTrapezoidal = true) {
    const std::string method = normalized_method_key(settings);
    if (method == "TRAP" || method == "TRAPEZOIDAL") return "Trapezoidal";
    if (method == "TRAPGEAR" || method == "TRAP_GEAR") {
        return autoUseTrapezoidal
            ? "TrapGear Hybrid (Trapezoidal default)"
            : "TrapGear Hybrid (Gear2 auto-damped)";
    }
    if (method == "GEAR2") return "Gear2/BDF2";
    if (method == "GEAR" || method == "BDF") return "Variable-step BDF";
    if (method == "ADAMS" || method == "ADAMSMOULTON" || method == "AM") {
        return "Variable-step Adams-Moulton";
    }
    if (method == "BE" || method == "BACKWARDEULER") return "Backward Euler";
    return autoUseTrapezoidal
        ? "Auto (Backward Euler startup, Trapezoidal after history)"
        : "Auto (Backward Euler startup, Gear2/BDF2 for compact-model damping)";
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
        for (int i = 0; i < num_devs; ++i) stamp_device_dc(*devices[i], J, b, x);
        VectorReal x_new = KluSolverReal::solve(J, b, &solver_context);
        double max_change = std::abs(x_new[2] - x[2]); x = x_new;
        if (max_change < tol) { std::cout << "  Converged in " << iter + 1 << " iterations." << std::endl; break; }
    }
    std::cout << "Results:\n  Vdd: " << x[0] << " V\n  Gate: " << x[1] << " V\n  Drain: " << x[2] << " V\n";
}

void run_simulation(
    Netlist& netlist,
    const std::string& output_file = "",
    int requested_threads = 1,
    const std::string& output_format = "RAW") {
    int num_nodes = netlist.getNumNodes();
    std::vector<Port*> ports;
    std::vector<StabilityProbe*> probes;
    auto& devices = netlist.getDevices();
    int num_devs = static_cast<int>(devices.size());

    const auto& requested_settings = netlist.getSettings();
    if (requested_threads <= 1 && requested_settings.num_threads > 0) {
        requested_threads = requested_settings.num_threads;
    }
    const std::vector<std::string> implemented_analyses = {
        "OP", "DC", "STEP", "MC", "CORNER", "SENS", "PZ", "TF",
        "TRAN", "AC", "NOISE", "SP", "STB", "PSS", "HB", "PNOISE"
    };
    if (std::find(implemented_analyses.begin(), implemented_analyses.end(), requested_settings.type) ==
        implemented_analyses.end()) {
        throw std::runtime_error(
            "analysis ." + requested_settings.type +
            " is parsed for compatibility but has no validated execution engine; refusing to return a substitute result");
    }

    int osdiInternalUnknowns = 0;
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
        auto* osdi = dynamic_cast<OSDIDevice*>(dev.get());
        if (osdi && requested_settings.osdi_internal_nodes) {
            osdiInternalUnknowns += osdi->bindInternalUnknowns(
                [&]() { return netlist.getNextBranchId(num_nodes); });
        }
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
    auto make_frequency_grid = [](const SimulationSettings& s) {
        std::vector<double> frequencies;
        const std::string sweep = upper_copy(s.frequency_sweep_type);
        if (sweep == "VALUES" || sweep == "VALUE" || sweep == "LIST") {
            return s.frequency_values;
        }
        if (!(s.f_start > 0.0) || !(s.f_stop >= s.f_start)) {
            return frequencies;
        }
        if (sweep == "LIN") {
            const int points = std::max(1, s.points_per_dec);
            frequencies.reserve(static_cast<std::size_t>(points));
            if (points == 1) {
                frequencies.push_back(s.f_start);
            } else {
                const double step = (s.f_stop - s.f_start) / static_cast<double>(points - 1);
                for (int i = 0; i < points; ++i) {
                    frequencies.push_back(i + 1 == points ? s.f_stop : s.f_start + step * static_cast<double>(i));
                }
            }
            return frequencies;
        }
        if (sweep == "STEP") {
            const double step = s.frequency_step > 0.0 ? s.frequency_step : (s.f_stop - s.f_start);
            if (!(step > 0.0)) {
                frequencies.push_back(s.f_start);
                return frequencies;
            }
            for (double f = s.f_start; f <= s.f_stop * (1.0 + 1e-12); f += step) {
                frequencies.push_back(std::min(f, s.f_stop));
                if (frequencies.size() > 1000000u) break;
            }
            if (frequencies.empty() || std::abs(frequencies.back() - s.f_stop) > std::max(1e-9, s.f_stop * 1e-12)) {
                frequencies.push_back(s.f_stop);
            }
            return frequencies;
        }
        const double base = (sweep == "OCT") ? 2.0 : 10.0;
        const int per = std::max(1, s.points_per_dec);
        const double mult = std::pow(base, 1.0 / static_cast<double>(per));
        for (double f = s.f_start; f <= s.f_stop * (1.0 + 1e-12); f *= mult) {
            frequencies.push_back(std::min(f, s.f_stop));
            if (frequencies.size() > 1000000u || !(mult > 1.0)) break;
        }
        if (frequencies.empty() || std::abs(frequencies.back() - s.f_stop) > std::max(1e-9, s.f_stop * 1e-12)) {
            frequencies.push_back(s.f_stop);
        }
        return frequencies;
    };
    real_solver_context.backend = settings.solver_backend;
    real_solver_context.ordering = settings.solver_ordering;
    real_solver_context.use_singletons = settings.solver_singletons;
    real_solver_context.scale_rows = settings.solver_row_scaling;
    complex_solver_context.backend = settings.solver_backend;
    complex_solver_context.ordering = settings.solver_ordering;
    complex_solver_context.use_singletons = settings.solver_singletons;
    complex_solver_context.scale_rows = settings.solver_row_scaling;
    if (settings.type == "PNOISE") {
        throw std::runtime_error(".PNOISE requires a preceding .PSS in the same deck so the periodic operating point is available.");
    }
    std::cout << "Linear solver: backend=" << settings.solver_backend
              << " ordering=" << settings.solver_ordering
              << " singleton_filter=" << (settings.solver_singletons ? "on" : "off")
              << " row_scaling=" << (settings.solver_row_scaling ? "on" : "off")
              << " refinement=" << settings.solver_refinement_steps
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
    if (settings.type == "TRAN" && settings.tran_trtol < 1.0) {
        std::cerr << "Warning: TRTOL=" << settings.tran_trtol
                  << " is a sub-unity tolerance multiplier and can force extremely small "
                     "transient steps; values >= 1 are recommended."
                  << std::endl;
    }
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
                  << " internal_unknowns=" << osdiInternalUnknowns
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
    if (!settings.nodesets.empty()) {
        std::cout << "Nodesets: applied " << settings.nodesets.size()
                  << " temporary nonlinear constraint(s) for "
                  << settings.nodeset_iterations << " iteration(s)." << std::endl;
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
        const std::uint64_t evaluationEpoch = next_evaluation_epoch();
        double previous_delta = std::numeric_limits<double>::infinity();
        for (int iter = 0; iter < settings.op_max_iter; ++iter) {
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixReal J_sparse(matrix_size); VectorReal b(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, active_gmin);
            const bool nodesetActive = iter < settings.nodeset_iterations && !settings.nodesets.empty();
            if (nodesetActive) {
                for (const auto& nodeset : settings.nodesets) {
                    if (nodeset.node < 0 || nodeset.node >= num_nodes) continue;
                    J_sparse.add(nodeset.node, nodeset.node, settings.nodeset_conductance);
                    b.add(nodeset.node, settings.nodeset_conductance * nodeset.value);
                }
            }
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
            std::vector<DaeStampStatus> stampStatus(static_cast<std::size_t>(num_devs));
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) {
                stampStatus[static_cast<std::size_t>(i)] = stamp_device_dc(
                    *devices[i], J_sparse, b, x, &settings, evaluationEpoch,
                    iter > 0, nodesetActive, false);
            }
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorReal x_new = solve_with_refinement(
                J_sparse, b, &real_solver_context, settings.solver_refinement_steps);
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
            bool converged = delta_converged && !nodesetActive;
            if (delta_converged && settings.nr_residual_check) {
                const auto residual = dc_residual_error(
                    devices, num_devs, matrix_size, num_nodes, x_new, settings, active_gmin);
                converged = !nodesetActive && residual.error <= 1.0 && !residual.limitingApplied;
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
        VectorReal x = initial_guess;
        int continuationSteps = 0;
        try {
            double sourceScale = settings.source_stepping ? 0.0 : 1.0;
            double activeGmin = settings.gmin_stepping
                ? std::max(settings.gmin, 1e-3)
                : settings.gmin;
            scale_source_states(source_states, sourceScale);

            bool anchored = false;
            double anchorGmin = activeGmin;
            for (int attempt = 0; attempt < 5 && !anchored; ++attempt) {
                try {
                    std::ostringstream stepLabel;
                    stepLabel << label << " continuation anchor source_scale=" << sourceScale
                              << " gmin=" << std::scientific << anchorGmin;
                    x = solve_dc_operating_point(x, stepLabel.str(), false, anchorGmin);
                    activeGmin = anchorGmin;
                    anchored = true;
                    ++continuationSteps;
                } catch (const std::exception&) {
                    if (!settings.gmin_stepping) throw;
                    anchorGmin = std::min(1.0, std::max(anchorGmin * 10.0, 1e-2));
                }
            }
            if (!anchored) throw std::runtime_error("adaptive continuation could not establish an anchor point");

            if (settings.source_stepping) {
                double increment = 0.10;
                while (sourceScale < 1.0 - 1e-12) {
                    const double candidate = std::min(1.0, sourceScale + increment);
                    scale_source_states(source_states, candidate);
                    try {
                        std::ostringstream stepLabel;
                        stepLabel << label << " adaptive source homotopy=" << candidate
                                  << " gmin=" << std::scientific << activeGmin;
                        x = solve_dc_operating_point(x, stepLabel.str(), false, activeGmin);
                        sourceScale = candidate;
                        increment = std::min(0.35, increment * 1.6);
                        ++continuationSteps;
                    } catch (const std::exception&) {
                        scale_source_states(source_states, sourceScale);
                        increment *= 0.5;
                        if (increment < 1e-4) {
                            throw std::runtime_error("adaptive source homotopy reached its minimum continuation step");
                        }
                    }
                }
            }

            if (settings.gmin_stepping && activeGmin > settings.gmin) {
                const double targetForLog = std::max(settings.gmin, 1e-18);
                double exponent = std::log10(activeGmin);
                const double targetExponent = std::log10(targetForLog);
                double exponentStep = -1.5;
                while (exponent > targetExponent + 1e-12) {
                    const double candidateExponent = std::max(targetExponent, exponent + exponentStep);
                    const double candidateGmin = std::pow(10.0, candidateExponent);
                    try {
                        std::ostringstream stepLabel;
                        stepLabel << label << " adaptive gmin homotopy=" << std::scientific << candidateGmin;
                        x = solve_dc_operating_point(x, stepLabel.str(), false, candidateGmin);
                        exponent = candidateExponent;
                        activeGmin = candidateGmin;
                        exponentStep = std::max(-3.0, exponentStep * 1.5);
                        ++continuationSteps;
                    } catch (const std::exception&) {
                        exponentStep *= 0.5;
                        if (std::abs(exponentStep) < 1e-3) {
                            throw std::runtime_error("adaptive gmin homotopy reached its minimum continuation step");
                        }
                    }
                }
                if (settings.gmin == 0.0) {
                    x = solve_dc_operating_point(x, label + " final zero-gmin solve", false, 0.0);
                    ++continuationSteps;
                }
            }
            restore_source_states(source_states);
            std::cout << "DC recovery: converged for " << label
                      << " with " << continuationSteps << " adaptive continuation step(s)." << std::endl;
            return x;
        } catch (const std::exception& homotopy_error) {
            restore_source_states(source_states);
            // ----------------------------------------------------------------
            // Pillar 4: Pseudo-Transient Continuation (PTC) — tertiary fallback.
            // Add a virtual C_ptc/dt conductance to every voltage node's
            // diagonal. Ramp C_ptc to zero to recover the DC OP.
            // ----------------------------------------------------------------
            std::cout << "DC recovery: source/gmin homotopy failed for " << label
                      << " (" << homotopy_error.what() << "); trying pseudo-transient continuation (PTC)."
                      << std::endl;
        }

        // PTC third-fallback.
        try {
            PtcController ptc;
            VectorReal x_ptc = initial_guess;
            int ptc_steps = 0;
            constexpr int ptc_max_steps = 200;
            while (ptc.active() && ptc_steps < ptc_max_steps) {
                const double g_ptc = ptc.diagonalConductance();
                // Build a modified DC OP: add g_ptc to diagonal of J for voltage nodes.
                // This is done by temporarily overriding gmin with a combined gmin+g_ptc.
                const double combined_gmin = settings.gmin + g_ptc;
                std::ostringstream stepLabel;
                stepLabel << label << " PTC step=" << ptc_steps
                          << " C_ptc=" << std::scientific << ptc.currentCapacitance()
                          << " g_ptc=" << g_ptc;
                try {
                    ptc.updatePreviousState(x_ptc);
                    x_ptc = solve_dc_operating_point(x_ptc, stepLabel.str(), false, combined_gmin);
                    ptc.notifyConverged();
                    ++continuationSteps;
                } catch (const std::exception&) {
                    ptc.notifyDiverged();
                }
                ++ptc_steps;
            }
            // Final solve with actual gmin (no PTC bias).
            x_ptc = solve_dc_operating_point(x_ptc, label + " PTC final", false);
            std::cout << "DC recovery (PTC): converged for " << label
                      << " with " << ptc_steps << " PTC step(s)." << std::endl;
            return x_ptc;
        } catch (...) {
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
    } else {
        for (const auto& dev : devices) {
            if (auto* vsrc = dynamic_cast<VoltageSource*>(dev.get())) {
                const int p = vsrc->getNodePos();
                const int n = vsrc->getNodeNeg();
                const double val = vsrc->evaluateAt(0.0);
                if (p >= 0 && p < num_nodes) {
                    const double n_val = (n >= 0 && n < num_nodes) ? x_dc[n] : 0.0;
                    x_dc[p] = n_val + val;
                }
                if (n >= 0 && n < num_nodes && (p < 0 || p >= num_nodes)) {
                    x_dc[n] = -val;
                }
            }
        }
    }
    run_dae_audits(devices, x_dc, settings);

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
        const std::vector<double> mc_samples = generate_mc_samples(settings);
        std::cout << "MC sampling: distribution="
                  << (upper_copy(settings.mc_distribution) == "UNIFORM" ? "uniform" : "gaussian")
                  << " sampling=" << (settings.mc_latin_hypercube ? "latin-hypercube" : "pseudo-random")
                  << " seed=" << settings.mc_seed << std::endl;
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
            const double value = mc_samples[static_cast<size_t>(run - 1)];
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
            for (int i = 0; i < num_devs; ++i) stamp_device_ac(*devices[i], J_sparse, b_ac, omega, x_dc);
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
        for (int i = 0; i < num_devs; ++i) stamp_device_ac(*devices[i], J_sparse, b_tf, 0.0, x_dc);
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
        // FIX Bug 4: do not cap max_step to output_step. With output_step==5ps
        // the old code forced max_step=5ps, preventing the integrator from ever
        // taking a step larger than the output interval. Use the physics-driven
        // default (t_stop/100 or similar) when the user has not set TMAX.
        const double max_step = default_transient_max_step(settings, output_step);
        const double min_step = settings.t_min_step > 0.0
            ? std::min(settings.t_min_step, max_step)
            : std::max(max_step * 1e-6, 1e-18);
        const double save_start = std::clamp(settings.t_start, 0.0, settings.t_stop);
        double step = std::min(output_step, 1e-11);
        VectorReal x = x_dc;
        if (!settings.initial_conditions.empty()) {
            for (const auto& ic : settings.initial_conditions) {
                if (ic.node >= 0 && ic.node < num_nodes) {
                    x[ic.node] = ic.value;
                }
            }
            step = std::min(step, 1e-12);
        }
        std::vector<VectorReal> x_hist; x_hist.push_back(x);
        std::vector<double> t_hist; t_hist.push_back(0.0);
        DeviceTransientStateArena transient_state_arena(devices);
        DaeTransientHistoryBank dae_history(devices, x, 0.0);
        const std::string norm_method = normalized_method_key(settings);
        bool compact_model_prefers_damping = false;
        const bool high_accuracy_mode = settings.reltol <= 3e-4 || settings.save_adaptive_steps;
        if (norm_method != "TRAPGEAR" && norm_method != "TRAP_GEAR" && norm_method != "TRAP" && norm_method != "TRAPEZOIDAL" && !high_accuracy_mode) {
            for (const auto& device : devices) {
                if (device->prefersDampedAutoTransient()) {
                    compact_model_prefers_damping = true;
                    break;
                }
            }
        }
        AutomaticTransientMethodController auto_method_controller(compact_model_prefers_damping);
        bool auto_use_trapezoidal = auto_method_controller.useTrapezoidal();
        int selected_order = 1;
        transient_state_arena.restoreCurrent();
        TransientOutput tran_out(output_file, output_format, netlist, num_nodes);
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
                  << " lte=" << (settings.tran_lte_mode == "STEPDOUBLING" ? "step-doubling" : "predictor-corrector")
                  << " lte_reference=" << settings.tran_lte_reference
                  << " lte_audit_interval=" << settings.tran_lte_audit_interval
                  << " trtol=" << settings.tran_trtol
                  << " order_control=" << (settings.tran_order_adaptive ? "adaptive" : "ramp")
                  << " stamp_cache=" << (settings.transient_stamp_cache ? "on" : "matrix-reuse")
                  << " method=" << normalized_transient_method(settings, auto_use_trapezoidal)
                  << std::endl;
        if (settings.tran_adaptive && settings.t_max_step > 0.0 &&
            max_step <= output_step * (1.0 + 1e-9)) {
            std::cout << "Transient note: user maxstep/tmax is at or below the output step; "
                      << "internal adaptive steps cannot grow beyond this ceiling. "
                      << "Omit .TRAN tmax or use a larger MAXSTEP for faster low/medium runs."
                      << std::endl;
        }
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
        int consecutive_min_step_lte_violations = 0;
        int recent_step_rejection_streak = 0;
        constexpr int max_min_step_lte_violations = 1024;
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
        double candidate_step = step;
        TransientMatrixWorkspace transient_matrix_workspace(matrix_size);
        transient_matrix_workspace.setStructureCacheEnabled(settings.transient_stamp_cache);
        // FIX Bug 2: track the candidate step before any breakpoint snap so we
        // can restore it after landing on a breakpoint instead of growing from
        // the artificially shrunk post-snap value.
        double pre_bp_candidate = step;
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
            double cur_step = candidate_step;
            cur_step = std::min(cur_step, remaining);

            // FIX Bug 3: device bound is a one-shot cap on cur_step only.
            // candidate_step is NOT modified here; the post-acceptance growth
            // path will grow candidate_step independently. This ensures that
            // once a pulse edge passes and the device bound relaxes, the step
            // recovers naturally without needing special "force-growth" logic.
            double current_bound_step = collect_transient_bound_step(devices);
            if (std::isfinite(current_bound_step) && current_bound_step > 0.0 && cur_step > current_bound_step) {
                cur_step = std::max(min_step, current_bound_step);
                ++tran_stats.bound_step_limited;
            }

            // FIX Bug 2 (part 1): save the pre-snap candidate before any
            // breakpoint shrinks cur_step. We restore it after acceptance so
            // that growth continues from the pre-event baseline.
            pre_bp_candidate = candidate_step;
            const double next_bp = next_breakpoint_after(breakpoints, t, bp_tol);
            if (next_bp > 0.0) {
                const double bp_dist = next_bp - t;
                if (bp_dist > stop_tol) {
                    if (bp_dist < cur_step) {
                        cur_step = bp_dist;
                    }
                }
            }
            // FIX Bug 1: DO NOT clamp the integration step to (next_output - t).
            // The integrator is free to step past output points; values at
            // scheduled output times are written from the accepted state at 't'
            // (which may be slightly past next_output) via the post-step flush
            // below. Clamping here permanently reduced candidate_step to the
            // tiny output interval (5 ps in the user's netlist), preventing any
            // step growth for the rest of the simulation.
            step = cur_step;
            if (step <= stop_tol) {
                step = std::min(max_step, std::max(min_step, remaining));
                if (step <= stop_tol) break;
            }
            bool accepted = false;
            VectorReal accepted_x = x;
            VectorReal accepted_intermediate_x;
            double accepted_intermediate_time = 0.0;
            bool accepted_has_intermediate = false;
            double err = 0.0;
            int accepted_order = 1;
            int proposed_next_order = selected_order;
            double proposed_growth = 1.0;
            const auto baseline_device_states = transient_state_arena.checkpoint();
            const auto baseline_dae_history = dae_history.checkpoint();
            while (!accepted) {
                transient_state_arena.restore(baseline_device_states);
                dae_history.rollback(baseline_dae_history);
                const double target_time = t + step;
                const TransientIntegrationMethod attempted_method =
                    choose_transient_method(
                        settings, x_hist, t_hist, step, auto_use_trapezoidal);
                TransientStepResult full = solve_transient_step(
                    devices, dae_history, num_devs, matrix_size, num_nodes, settings, &real_solver_context,
                    x, x_hist, t_hist, step, target_time, &attempted_method, selected_order,
                    &transient_matrix_workspace);
                tran_stats.noteSolve(full);
                if (!full.converged && step > min_step) {
                    ++tran_stats.rejected_steps;
                    ++tran_stats.convergence_rejections;
                    ++recent_step_rejection_streak;
                    selected_order = std::max(1, selected_order - 1);
                    step = std::max(min_step, step * 0.5);
                    continue;
                }
                if (!full.converged) {
                    throw std::runtime_error(transient_failure_message(
                        "Transient step failed to converge at minimum timestep",
                        target_time, step, full, num_nodes));
                }

                const TransientContext full_ctx = make_transient_context(
                    settings, x_hist, t_hist, step, target_time, &attempted_method, selected_order);
                const bool can_reduce_step = step > min_step * (1.0 + 1e-9);
                PredictorCorrectorEstimate pcEstimate;
                if (settings.tran_adaptive) {
                    pcEstimate = predictor_corrector_lte_error(
                        devices, transient_state_arena, full.x, x_hist, t_hist,
                        full_ctx, num_nodes, settings);
                }
                const double next_bp_after_target = next_breakpoint_after(breakpoints, target_time, bp_tol);
                const bool near_breakpoint =
                    (next_bp > 0.0 && std::abs(target_time - next_bp) <= bp_tol) ||
                    (next_bp_after_target > 0.0 &&
                     next_bp_after_target - target_time <= std::max(max_step, output_step));
                const bool periodicOracle = settings.tran_lte_audit_interval > 0 &&
                    ((tran_stats.accepted_steps + 1) % settings.tran_lte_audit_interval == 0);
                const bool useStepDoubling = can_reduce_step &&
                    (settings.tran_lte_mode == "STEPDOUBLING" ||
                     !pcEstimate.valid || near_breakpoint || periodicOracle);

                if (settings.tran_adaptive && pcEstimate.valid && !useStepDoubling) {
                    ++tran_stats.predictor_lte_steps;
                    err = pcEstimate.error;
                    if (err > 1.0) {
                        transient_state_arena.restore(baseline_device_states);
                        dae_history.rollback(baseline_dae_history);
                        if (!can_reduce_step) {
                            ++consecutive_min_step_lte_violations;
                            if (consecutive_min_step_lte_violations >= max_min_step_lte_violations) {
                                throw std::runtime_error(transient_failure_message(
                                    "Transient progress starved after 1024 LTE violations at minimum "
                                    "timestep; TRTOL is a tolerance multiplier (normally >= 1). "
                                    "Check TRTOL, LTE_RELTOL, CHGTOL, or MINSTEP",
                                    target_time, step, full, num_nodes));
                            }
                        } else {
                            const double factor = transient_reject_factor(err, full.integration_order);
                            ++tran_stats.rejected_steps;
                            ++tran_stats.lte_rejections;
                            ++recent_step_rejection_streak;
                            selected_order = std::max(1, selected_order - 1);
                            step = std::max(min_step, step * factor);
                            continue;
                        }
                    }
                    accepted_x = full.x;
                    accept_device_transient_step(devices, accepted_x, target_time, full_ctx);
                    dae_history.commitAccepted(devices, accepted_x, target_time, full_ctx);
                    transient_state_arena.commitDeviceState();
                    accepted_order = full.integration_order;

                    const bool variableOrder = attempted_method == TransientIntegrationMethod::Bdf ||
                        attempted_method == TransientIntegrationMethod::AdamsMoulton;
                    if (variableOrder && settings.tran_order_adaptive) {
                        std::optional<double> lowerError;
                        std::optional<double> higherError;
                        if (accepted_order > 1) {
                            const auto lowerContext = make_transient_context(
                                settings, x_hist, t_hist, step, target_time,
                                &attempted_method, accepted_order - 1);
                            const auto lower = predictor_corrector_lte_error(
                                devices, transient_state_arena, accepted_x, x_hist, t_hist,
                                lowerContext, num_nodes, settings);
                            if (lower.valid) lowerError = lower.error;
                        }
                        if (accepted_order < settings.tran_max_order) {
                            const auto higherContext = make_transient_context(
                                settings, x_hist, t_hist, step, target_time,
                                &attempted_method, accepted_order + 1);
                            if (higherContext.integrationOrder == accepted_order + 1) {
                                const auto higher = predictor_corrector_lte_error(
                                    devices, transient_state_arena, accepted_x, x_hist, t_hist,
                                    higherContext, num_nodes, settings);
                                if (higher.valid) higherError = higher.error;
                            }
                        }
                        const auto decision = chooseAdaptiveOrder(
                            accepted_order, settings.tran_max_order, err, lowerError, higherError);
                        proposed_next_order = decision.order;
                        proposed_growth = decision.projectedFactor;
                    } else {
                        proposed_growth = transient_growth_factor(err, accepted_order);
                    }
                } else if (settings.tran_adaptive && can_reduce_step) {
                    ++tran_stats.step_doubling_audits;
                    // The full-step evaluation may have modified OSDI model
                    // memory and candidate state. Start the half-step path from
                    // the exact accepted state at time t.
                    transient_state_arena.restore(baseline_device_states);
                    dae_history.rollback(baseline_dae_history);
                    const double half_step = step * 0.5;
                    TransientStepResult half1 = solve_transient_step(
                        devices, dae_history, num_devs, matrix_size, num_nodes, settings, &real_solver_context,
                        x, x_hist, t_hist, half_step, t + half_step, &attempted_method, selected_order,
                        &transient_matrix_workspace);
                    tran_stats.noteSolve(half1);
                    if (!half1.converged) {
                        transient_state_arena.restore(baseline_device_states);
                        dae_history.rollback(baseline_dae_history);
                        if (step > min_step) {
                            ++tran_stats.rejected_steps;
                            ++tran_stats.convergence_rejections;
                            ++recent_step_rejection_streak;
                            selected_order = std::max(1, selected_order - 1);
                            step = std::max(min_step, half_step);
                            continue;
                        }
                        throw std::runtime_error(transient_failure_message(
                            "Transient first half-step failed to converge at minimum timestep",
                            t + half_step, half_step, half1, num_nodes));
                    }
                    const TransientContext half1_ctx = make_transient_context(
                        settings, x_hist, t_hist, half_step, t + half_step, &attempted_method, selected_order);
                    accept_device_transient_step(devices, half1.x, t + half_step, half1_ctx);
                    dae_history.commitAccepted(devices, half1.x, t + half_step, half1_ctx);
                    transient_state_arena.commitDeviceState();
                    std::vector<VectorReal> half_hist = x_hist;
                    half_hist.push_back(half1.x);
                    std::vector<double> half_time_hist = t_hist;
                    half_time_hist.push_back(t + half_step);
                    TransientStepResult half2 = solve_transient_step(
                        devices, dae_history, num_devs, matrix_size, num_nodes, settings, &real_solver_context,
                        half1.x, half_hist, half_time_hist, half_step, target_time, &attempted_method, selected_order,
                        &transient_matrix_workspace);
                    tran_stats.noteSolve(half2);
                    if (!half2.converged) {
                        transient_state_arena.restore(baseline_device_states);
                        dae_history.rollback(baseline_dae_history);
                        if (step > min_step) {
                            ++tran_stats.rejected_steps;
                            ++tran_stats.convergence_rejections;
                            ++recent_step_rejection_streak;
                            selected_order = std::max(1, selected_order - 1);
                            step = std::max(min_step, half_step);
                            continue;
                        }
                        throw std::runtime_error(transient_failure_message(
                            "Transient second half-step failed to converge at minimum timestep",
                            target_time, half_step, half2, num_nodes));
                    }
                    err = transient_lte_error(
                        devices, transient_state_arena, full.x, half2.x,
                        x_hist, num_nodes, settings, full.integration_order, target_time);
                    if (pcEstimate.valid) err = std::max(err, pcEstimate.error);
                    if (err > 1.0 && step > min_step) {
                        transient_state_arena.restore(baseline_device_states);
                        dae_history.rollback(baseline_dae_history);
                        const double factor = transient_reject_factor(err, full.integration_order);
                        ++tran_stats.rejected_steps;
                        ++tran_stats.lte_rejections;
                        ++recent_step_rejection_streak;
                        selected_order = std::max(1, selected_order - 1);
                        step = std::max(min_step, step * factor);
                        continue;
                    }
                    accepted_x = half2.x;
                    accepted_intermediate_x = half1.x;
                    accepted_intermediate_time = t + half_step;
                    accepted_has_intermediate = true;
                    const TransientContext half2_ctx = make_transient_context(
                        settings, half_hist, half_time_hist, half_step, target_time, &attempted_method, selected_order);
                    accept_device_transient_step(devices, accepted_x, target_time, half2_ctx);
                    dae_history.commitAccepted(devices, accepted_x, target_time, half2_ctx);
                    transient_state_arena.commitDeviceState();
                    accepted_order = half2.integration_order;
                    proposed_growth = transient_growth_factor(err, accepted_order);
                    if ((attempted_method == TransientIntegrationMethod::Bdf ||
                         attempted_method == TransientIntegrationMethod::AdamsMoulton) &&
                        settings.tran_order_adaptive) {
                        if (err < 0.05 && accepted_order < settings.tran_max_order &&
                            t_hist.size() >= static_cast<std::size_t>(accepted_order + 1)) {
                            proposed_next_order = accepted_order + 1;
                        } else if (err > 0.8 && accepted_order > 1) {
                            proposed_next_order = accepted_order - 1;
                        }
                    }
                } else {
                    accepted_x = full.x;
                    accept_device_transient_step(devices, accepted_x, target_time, full_ctx);
                    dae_history.commitAccepted(devices, accepted_x, target_time, full_ctx);
                    transient_state_arena.commitDeviceState();
                    accepted_order = full.integration_order;
                    proposed_growth = settings.tran_adaptive
                        ? transient_growth_factor(0.0, accepted_order)
                        : 1.0;
                    if (attempted_method == TransientIntegrationMethod::Bdf ||
                        attempted_method == TransientIntegrationMethod::AdamsMoulton) {
                        proposed_next_order = std::min(settings.tran_max_order, accepted_order + 1);
                    }
                }
                if (err <= 1.0) consecutive_min_step_lte_violations = 0;
                accepted = true;
            }
            const double previous_t = t;
            const VectorReal previous_x = x;
            t += step;
            if (settings.t_stop - t < stop_tol) {
                t = settings.t_stop;
            }
            tran_stats.noteAccepted(step, err);
            const int rejection_streak_before_accept = recent_step_rejection_streak;
            recent_step_rejection_streak = 0;
            x = accepted_x;
            record_measure_sample(t, x);
            device_bound_step = collect_transient_bound_step(devices);
            if (std::isfinite(device_bound_step)) {
                tran_stats.min_bound_step = std::min(tran_stats.min_bound_step, device_bound_step);
            }
            if (accepted_has_intermediate) {
                x_hist.push_back(accepted_intermediate_x);
                t_hist.push_back(accepted_intermediate_time);
            }
            x_hist.push_back(x);
            t_hist.push_back(t);
            if (is_transient_breakpoint(breakpoints, t, bp_tol) && t > bp_tol) {
                dae_history.restart(devices, x, t);
                x_hist.clear();
                t_hist.clear();
                x_hist.push_back(x);
                t_hist.push_back(t);
                selected_order = 1;
                candidate_step = std::min(output_step, 1e-11);
                auto_method_controller.restart();
                auto_use_trapezoidal = auto_method_controller.useTrapezoidal();
            } else {
                selected_order = std::clamp(
                    proposed_next_order, 1, settings.tran_max_order);
                if (normalized_method_key(settings) == "AUTO" && settings.tran_trap_ringing) {
                    if (auto_method_controller.observe(
                            x_hist, num_nodes, settings.vntol, settings.reltol)) {
                        ++tran_stats.method_switches;
                    }
                    auto_use_trapezoidal = auto_method_controller.useTrapezoidal();
                }
            }
            // Trim solution history to 8 entries using O(N) move instead of
            // repeated front-erase. A proper circular buffer is defined in
            // transient_state_store.hpp; this pop_front approach avoids the
            // O(N²) total cost of repeated erase-from-front on large circuits.
            while (x_hist.size() > 8) {
                x_hist.erase(x_hist.begin(), x_hist.begin() + 1);
                t_hist.erase(t_hist.begin(), t_hist.begin() + 1);
            }
            const double output_tol = std::max(1e-30, output_step * 1e-9);
            bool is_fast_transition = false;
            if (num_nodes > 0 && x_hist.size() >= 2 && t_hist.size() >= 2) {
                const double dt_step = t_hist.back() - t_hist[t_hist.size() - 2];
                if (dt_step > 1e-18) {
                    const auto& x_prev = x_hist[x_hist.size() - 2];
                    for (int i = 0; i < num_nodes; ++i) {
                        const double dv_dt = std::abs(x[i] - x_prev[i]) / dt_step;
                        if (dv_dt > 1e9) {
                            is_fast_transition = true;
                            break;
                        }
                    }
                }
            }
            while (t + output_tol >= next_output) {
                const VectorReal sample_x = interpolate_transient_state(
                    previous_x, x, previous_t, t, next_output);
                tran_out.write(next_output, sample_x, num_nodes);
                ++tran_stats.output_points;
                last_printed_time = next_output;
                do {
                    next_output += output_step;
                } while (next_output <= last_printed_time + output_tol);
            }
            const bool should_write_adaptive =
                settings.save_adaptive_steps || is_fast_transition;
            if (should_write_adaptive && last_printed_time < t - output_tol) {
                tran_out.write(t, x, num_nodes);
                ++tran_stats.output_points;
                last_printed_time = t;
            }
            if (tran_out.toFile() && settings.t_stop > 0.0) {
                const double pct = std::clamp(100.0 * t / settings.t_stop, 0.0, 100.0);
                if ((settings.tran_progress_interval > 0.0 &&
                     pct - last_progress_percent >= settings.tran_progress_interval) ||
                    pct >= 100.0) {
                    std::cout << std::fixed << std::setprecision(1)
                              << "Transient progress: " << pct << "% t="
                              << std::scientific << std::setprecision(9) << t << std::endl;
                    last_progress_percent = pct;
                }
            }
            if (is_transient_breakpoint(breakpoints, t, bp_tol) && t > bp_tol) {
                // FIX Bug 2 (part 2): restore the pre-breakpoint candidate step
                // as a ceiling, but grow from the accepted breakpoint step.
                // This avoids both extremes: staying at an absurdly tiny
                // post-breakpoint step or jumping straight back to a step that
                // the discontinuity has just invalidated.
                const double restored = std::min(max_step, pre_bp_candidate);
                if (settings.tran_adaptive) {
                    const double bp_growth = std::clamp(settings.tran_breakpoint_growth, 1.0, 20.0);
                    candidate_step = std::min(
                        max_step,
                        std::max(min_step, std::min(restored, step * bp_growth)));
                } else {
                    candidate_step = max_step;
                }
            } else if (settings.tran_adaptive) {
                double max_grow = (err <= 1e-12) ? 8.0 : ((err <= 1e-6) ? 5.0 : 2.5);
                if (rejection_streak_before_accept >= 3) {
                    max_grow = std::min(max_grow, 1.20);
                } else if (rejection_streak_before_accept >= 1) {
                    max_grow = std::min(max_grow, 1.50);
                }
                const double grow = std::clamp(proposed_growth, 1.05, max_grow);
                candidate_step = std::min(max_step, std::max(min_step, step * grow));
            } else {
                candidate_step = max_step;
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
                  << " pc_lte_steps=" << tran_stats.predictor_lte_steps
                  << " oracle_steps=" << tran_stats.step_doubling_audits
                  << " method_switches=" << tran_stats.method_switches
                  << " bound_step_limited=" << tran_stats.bound_step_limited
                  << " output_points=" << tran_stats.output_points
                  << " min_step=" << min_used_step
                  << " max_step=" << tran_stats.max_step
                  << " max_order_used=" << tran_stats.max_integration_order
                  << " min_bound_step=" << min_device_bound
                  << std::endl;
        std::cout << std::fixed << std::setprecision(2)
                  << "Newton summary: solves=" << tran_stats.newton_solves
                  << " total_iterations=" << tran_stats.newton_iterations
                  << " average_iterations=" << tran_stats.averageNewtonIterations()
                  << " max_iterations=" << tran_stats.max_newton_iterations
                  << " bypassed_devices=" << tran_stats.bypassed_devices
                  << " limiting_rechecks=" << tran_stats.limiting_rechecks
                  << " last_residual_error=" << tran_stats.last_residual_error
                  << " max_residual_error=" << tran_stats.max_residual_error
                  << " stamp_seconds=" << tran_stats.stamp_seconds
                  << " solve_seconds=" << tran_stats.solve_seconds
                  << std::endl;
        std::cout << std::scientific << std::setprecision(9)
                  << "Accuracy summary: method="
                  << normalized_transient_method(settings, auto_use_trapezoidal)
                  << " reltol=" << settings.reltol
                  << " vntol=" << settings.vntol
                  << " abstol=" << settings.abstol
                  << " lte_reltol=" << settings.tran_lte_reltol
                  << " lte_abstol=" << settings.tran_lte_abstol
                  << " trtol=" << settings.tran_trtol
                  << " chgtol=" << settings.chgtol
                  << " maxord=" << settings.tran_max_order
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
        const auto frequencies = make_frequency_grid(settings);
        for (double f : frequencies) {
            double omega = 2.0 * 3.14159265358979 * f;
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixComplex J_sparse(matrix_size); VectorComplex b_ac(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) stamp_device_ac(*devices[i], J_sparse, b_ac, omega, x_dc);
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorComplex x_ac = KluSolverComplex::solve(J_sparse, b_ac, &complex_solver_context);
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.ac_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
            runtime_stats.ac_solve_seconds += elapsed_seconds(solve_start, solve_end);
            std::cout << std::scientific << std::setprecision(2) << f << " | ";
            for(int i=0; i<num_nodes; ++i) std::cout << "(" << x_ac[i].real() << "," << x_ac[i].imag() << ") ";
            std::cout << std::endl;
        }
        std::cout << "AC summary: points=" << frequencies.size()
                  << " sweep=" << upper_copy(settings.frequency_sweep_type)
                  << std::endl;
    } else if (settings.type == "NOISE") {
        std::cout << "Starting Noise Analysis..." << std::endl;
        if (settings.out_node < 0 || settings.out_node >= num_nodes) {
            throw std::runtime_error(".NOISE output node must be a non-ground circuit node");
        }
        const auto frequencies = make_frequency_grid(settings);
        std::cout << "freq | onoise_sqrt(V/rtHz) onoise_psd(V^2/Hz) noise_sources" << std::endl;
        double integrated_output_noise = 0.0;
        double previous_freq = 0.0;
        double previous_psd = 0.0;
        bool has_previous_noise = false;
        std::size_t last_noise_source_count = 0;
        for (double f : frequencies) {
            double omega = 2.0 * 3.14159265358979 * f;
            const auto stamp_start = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            const NoiseTransferResult noise = solve_output_noise_psd(
                devices, num_devs, matrix_size, num_nodes, settings.out_node,
                omega, x_dc, settings, &complex_solver_context);
            const double output_psd = noise.output_psd;
            if (has_previous_noise && f > previous_freq) {
                integrated_output_noise += 0.5 * (previous_psd + output_psd) * (f - previous_freq);
            }
            has_previous_noise = true;
            previous_freq = f;
            previous_psd = output_psd;
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.ac_stamp_seconds += elapsed_seconds(stamp_start, solve_start);
            runtime_stats.ac_solve_seconds += elapsed_seconds(solve_start, solve_end);
            last_noise_source_count = noise.source_count;
            std::cout << std::scientific << std::setprecision(9)
                      << f << " | " << std::sqrt(std::max(output_psd, 0.0))
                      << " " << output_psd
                      << " " << noise.source_count
                      << std::endl;
        }
        std::cout << std::scientific << std::setprecision(9)
                  << "Integrated output noise: " << std::sqrt(std::max(integrated_output_noise, 0.0))
                  << " Vrms" << std::endl;
        std::cout << "Noise summary: points=" << frequencies.size()
                  << " output=V(" << netlist.getNodeName(settings.out_node) << ")"
                  << " input=" << settings.noise_input_source
                  << " sources=" << last_noise_source_count
                  << std::endl;
    } else if (settings.type == "SP") {
        std::cout << "Starting S-Parameter Analysis..." << std::endl;
        if (ports.empty()) {
            throw std::runtime_error(".SP requires at least one P port");
        }
        const auto frequencies = make_frequency_grid(settings);
        for (double f : frequencies) {
            const double omega = 2.0 * 3.14159265358979 * f;
            SparseMatrixComplex J_sparse(matrix_size); VectorComplex b_ac(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            for (int i = 0; i < num_devs; ++i) stamp_device_ac(*devices[i], J_sparse, b_ac, omega, x_dc);
            VectorComplex x_sp = KluSolverComplex::solve(J_sparse, b_ac, &complex_solver_context);
            std::cout << std::scientific << std::setprecision(9)
                      << f << " | ";
            for (std::size_t i = 0; i < ports.size(); ++i) {
                const int branch = ports[i]->getBranchIndex();
                const std::complex<double> value = branch >= 0 ? x_sp[branch] : std::complex<double>{};
                std::cout << "S" << (i + 1) << (i + 1) << "=("
                          << value.real() << "," << value.imag() << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << "SP summary: points=" << frequencies.size()
                  << " ports=" << ports.size()
                  << " method=multiport-y-to-s"
                  << std::endl;
    } else if (settings.type == "STB") {
        std::cout << "Starting Stability Analysis (Tian)..." << std::endl;
        if (probes.empty()) {
            throw std::runtime_error(".STB requires a stability probe; refusing to report an empty successful analysis");
        }
        const auto frequencies = make_frequency_grid(settings);
        for (double f : frequencies) {
            double omega = 2.0 * 3.14159265358979 * f;
            // 1. Voltage Pass
            SparseMatrixComplex Jv(matrix_size); VectorComplex bv(matrix_size);
            for (const auto& dev : devices) {
                auto* p = dynamic_cast<StabilityProbe*>(dev.get());
                if (p) p->stbStamp(Jv, bv, 1); else stamp_device_ac(*dev, Jv, bv, omega, x_dc);
            }
            VectorComplex xv = KluSolverComplex::solve(Jv, bv, &complex_solver_context);
            std::complex<double> Tv = -xv[probes[0]->getNodeNeg()] / xv[probes[0]->getNodePos()];
            // 2. Current Pass
            SparseMatrixComplex Ji(matrix_size); VectorComplex bi(matrix_size);
            for (const auto& dev : devices) {
                auto* p = dynamic_cast<StabilityProbe*>(dev.get());
                if (p) p->stbStamp(Ji, bi, 2); else stamp_device_ac(*dev, Ji, bi, omega, x_dc);
            }
            VectorComplex xi = KluSolverComplex::solve(Ji, bi, &complex_solver_context);
            std::complex<double> Ti = xi[probes[0]->getBranchIndex()];
            // 3. Combine
            std::complex<double> T = (Tv * Ti - std::complex<double>(1,0)) / (Tv + Ti + std::complex<double>(2,0));
            std::cout << std::scientific << f << " | Mag: " << std::abs(T) << " Phase: " << std::arg(T)*180/3.1415 << std::endl;
        }
        std::cout << "STB summary: points=" << frequencies.size()
                  << " method=twoport-return-ratio"
                  << std::endl;
    } else if (settings.type == "PSS") {
        std::cout << "Starting Periodic Steady State (PSS) Analysis..." << std::endl;
        int n_tones = static_cast<int>(settings.f_fund.size());
        std::cout << "  Tones: " << n_tones << " | Harmonics/Tone: " << settings.n_harms << std::endl;
        if (n_tones > 1) {
            std::cout << "  Warning: Multi-tone PSS requires finding a common periodic denominator, which can result in massive transient integration times." << std::endl;
            std::cout << "  For " << n_tones << " non-harmonically related tones, Quasi-Periodic Harmonic Balance (QPHB) is recommended." << std::endl;
        }
        if (settings.f_fund.empty() || !(settings.f_fund[0] > 0.0)) {
            throw std::runtime_error(".PSS requires a positive fundamental frequency or oscillator frequency estimate.");
        }
        const double initial_period = 1.0 / settings.f_fund[0];
        const double tstab = settings.pss_tstab > 0.0
            ? settings.pss_tstab
            : static_cast<double>(settings.pss_autonomous
                  ? std::max(5, settings.pss_tstab_periods)
                  : std::max(0, settings.pss_tstab_periods)) * initial_period;
        const double report_base_step = initial_period / static_cast<double>(std::max(16, settings.n_harms * 8));
        const double report_max_step = settings.t_max_step > 0.0
            ? settings.t_max_step
            : report_base_step;
        std::cout << std::scientific << std::setprecision(9)
                  << "  Executing Shooting Newton PSS: mode="
                  << (settings.pss_autonomous ? "autonomous" : "driven")
                  << " period=" << initial_period
                  << " tstab=" << tstab
                  << " max_step=" << report_max_step
                  << " residual_goal=" << settings.pss_residual_goal
                  << " max_iter=" << settings.max_pss_iter
                  << std::endl;

        struct PssOrbitResult {
            VectorReal end;
            VectorReal endpoint_time_sensitivity;
            std::vector<VectorReal> transition;
            bool variational = false;
            bool direct_period_sensitivity = false;
            int variational_steps = 0;
        };

        auto identity_transition = [&](int size) {
            std::vector<VectorReal> phi;
            phi.reserve(static_cast<std::size_t>(size));
            for (int col = 0; col < size; ++col) {
                VectorReal column(size);
                column[col] = 1.0;
                phi.push_back(column);
            }
            return phi;
        };

        auto stamp_dynamic_jacobian = [&](SparseMatrixReal& dynamic, const VectorReal& state, double time) {
            dynamic.setStructureCacheEnabled(settings.transient_stamp_cache);
            DaeRequest request;
            request.analysis = DaeAnalysis::Transient;
            request.time = time;
            request.staticResidual = false;
            request.dynamicResidual = false;
            request.staticJacobian = false;
            request.dynamicJacobian = true;
            request.readOnlyState = true;
            request.enableLimiting = false;
            for (const auto& device : devices) {
                DaeEvaluation evaluation;
                if (!device->evaluateDae(state, request, evaluation)) continue;
                if (!evaluation.finite()) {
                    throw std::runtime_error("non-finite value in PSS variational DAE Jacobian");
                }
                for (const auto& term : evaluation.dynamicJacobian) {
                    dynamic.add(term.equation, term.unknown, term.value);
                }
            }
            dynamic.getEntries();
        };

        auto multiply_sparse_vector = [&](const SparseMatrixReal& matrix, const VectorReal& vector) {
            VectorReal result(matrix_size);
            for (const auto& entry : matrix.getEntries()) {
                if (entry.row < 0 || entry.row >= matrix_size ||
                    entry.col < 0 || entry.col >= vector.getSize()) {
                    continue;
                }
                result.add(entry.row, entry.value * vector[entry.col]);
            }
            return result;
        };

        auto integrate_pss_interval = [&](const VectorReal& start, double duration, double start_time,
                                          bool propagate_variational, TransientOutput* output = nullptr,
                                          std::vector<TranSample>* samples = nullptr) {
            PssOrbitResult orbit;
            orbit.end = start;
            orbit.endpoint_time_sensitivity = VectorReal(matrix_size);
            if (propagate_variational) {
                orbit.transition = identity_transition(matrix_size);
                orbit.variational = true;
            }
            VectorReal x = start;
            if (!(duration > 0.0)) return orbit;
            const double base_step = duration / static_cast<double>(std::max(16, settings.n_harms * 8));
            const double pss_max_step = settings.t_max_step > 0.0
                ? settings.t_max_step
                : std::max(base_step, 1e-18);
            const double pss_min_step = std::max(pss_max_step * 1e-6, 1e-18);
            std::vector<VectorReal> x_hist; x_hist.push_back(x);
            std::vector<double> t_hist; t_hist.push_back(start_time);
            DeviceTransientStateArena transient_state_arena(devices);
            DaeTransientHistoryBank dae_history(devices, x, start_time);
            TransientMatrixWorkspace pss_transient_matrix_workspace(matrix_size);
            pss_transient_matrix_workspace.setStructureCacheEnabled(settings.transient_stamp_cache);
            transient_state_arena.restoreCurrent();
            double local_t = start_time;
            const double stop_time = start_time + duration;
            double candidate_step = std::min(pss_max_step, std::max(pss_min_step, base_step));
            int selected_order = 1;
            const double stop_tol = std::max(std::max(1e-30, std::abs(stop_time) * 1e-9), pss_min_step * 10.0);
            while (local_t < stop_time - stop_tol) {
                double step = std::min(candidate_step, stop_time - local_t);
                bool accepted = false;
                double err = 0.0;
                int accepted_order = selected_order;
                int proposed_next_order = selected_order;
                double proposed_growth = 1.0;
                const auto baseline_device_states = transient_state_arena.checkpoint();
                const auto baseline_dae_history = dae_history.checkpoint();
                while (!accepted) {
                    transient_state_arena.restore(baseline_device_states);
                    dae_history.rollback(baseline_dae_history);
                    const double target_time = local_t + step;
                    const TransientIntegrationMethod attempted_method =
                        choose_transient_method(settings, x_hist, t_hist, step, true);
                    TransientStepResult full = solve_transient_step(
                        devices, dae_history, num_devs, matrix_size, num_nodes, settings, &real_solver_context,
                        x, x_hist, t_hist, step, target_time, &attempted_method, selected_order,
                        &pss_transient_matrix_workspace);
                    if (!full.converged && step > pss_min_step) {
                        step = std::max(pss_min_step, step * 0.5);
                        selected_order = std::max(1, selected_order - 1);
                        continue;
                    }
                    if (!full.converged) {
                        throw std::runtime_error(transient_failure_message(
                            "PSS shooting transient step failed to converge",
                            target_time, step, full, num_nodes));
                    }
                    const TransientContext full_ctx = make_transient_context(
                        settings, x_hist, t_hist, step, target_time, &attempted_method, selected_order);
                    PredictorCorrectorEstimate pcEstimate;
                    if (settings.tran_adaptive || settings.pss_adaptive) {
                        pcEstimate = predictor_corrector_lte_error(
                            devices, transient_state_arena, full.x, x_hist, t_hist,
                            full_ctx, num_nodes, settings);
                    }
                    if (pcEstimate.valid) {
                        err = pcEstimate.error;
                        if (err > 1.0 && step > pss_min_step) {
                            step = std::max(pss_min_step, step * transient_reject_factor(err, full.integration_order));
                            selected_order = std::max(1, selected_order - 1);
                            continue;
                        }
                    }
                    if (propagate_variational) {
                        SparseMatrixReal step_matrix(matrix_size);
                        step_matrix.setStructureCacheEnabled(settings.transient_stamp_cache);
                        VectorReal step_rhs(matrix_size);
                        stamp_global_gmin(step_matrix, num_nodes, settings.gmin);
                        const std::uint64_t variational_epoch = next_evaluation_epoch();
                        for (std::size_t device_index = 0; device_index < devices.size(); ++device_index) {
                            stamp_device_transient(
                                *devices[device_index], device_index, dae_history,
                                step_matrix, step_rhs, full.x, full_ctx, &settings,
                                variational_epoch, false, true);
                        }
                        step_matrix.getEntries();
                        SparseMatrixReal previous_dynamic(matrix_size);
                        stamp_dynamic_jacobian(previous_dynamic, x, local_t);
                        for (auto& column : orbit.transition) {
                            VectorReal rhs = multiply_sparse_vector(previous_dynamic, column);
                            for (int row = 0; row < rhs.getSize(); ++row) {
                                rhs[row] *= -full_ctx.a1;
                            }
                            column = solve_with_refinement(
                                step_matrix, rhs, &real_solver_context, settings.solver_refinement_steps);
                        }
                        ++orbit.variational_steps;
                    }
                    accept_device_transient_step(devices, full.x, target_time, full_ctx);
                    dae_history.commitAccepted(devices, full.x, target_time, full_ctx);
                    transient_state_arena.commitDeviceState();
                    for (int i = 0; i < matrix_size; ++i) {
                        orbit.endpoint_time_sensitivity[i] = (full.x[i] - x[i]) / std::max(step, 1e-30);
                    }
                    orbit.direct_period_sensitivity = true;
                    x = full.x;
                    local_t = target_time;
                    if (output) {
                        output->write(local_t - start_time, x, num_nodes);
                    }
                    if (samples) {
                        samples->push_back({local_t - start_time, x});
                    }
                    accepted_order = full.integration_order;
                    if (pcEstimate.valid) {
                        proposed_growth = transient_growth_factor(err, accepted_order);
                    } else {
                        proposed_growth = 1.5;
                    }
                    if ((attempted_method == TransientIntegrationMethod::Bdf ||
                         attempted_method == TransientIntegrationMethod::AdamsMoulton) &&
                        settings.tran_order_adaptive) {
                        if (err < 0.05 && accepted_order < settings.tran_max_order &&
                            t_hist.size() >= static_cast<std::size_t>(accepted_order + 1)) {
                            proposed_next_order = accepted_order + 1;
                        } else if (err > 0.8 && accepted_order > 1) {
                            proposed_next_order = accepted_order - 1;
                        }
                    }
                    accepted = true;
                }
                x_hist.push_back(x);
                t_hist.push_back(local_t);
                while (x_hist.size() > 8) {
                    x_hist.erase(x_hist.begin(), x_hist.begin() + 1);
                    t_hist.erase(t_hist.begin(), t_hist.begin() + 1);
                }
                selected_order = std::clamp(proposed_next_order, 1, settings.tran_max_order);
                const double grow = std::clamp(proposed_growth, 1.05, 2.5);
                candidate_step = std::min(pss_max_step, std::max(pss_min_step, candidate_step * grow));
            }
            orbit.end = x;
            return orbit;
        };

        VectorReal pss_state = x_dc;
        double pss_time = 0.0;
        double period = initial_period;
        if (tstab > 0.0) {
            std::cout << "  PSS stabilization transient..." << std::endl;
            pss_state = integrate_pss_interval(pss_state, tstab, pss_time, false).end;
            pss_time += tstab;
        }
        double residual = std::numeric_limits<double>::infinity();
        int converged_iter = 0;
        int shooting_newton_iterations = 0;
        int monodromy_columns = 0;
        int variational_steps = 0;
        int direct_period_columns = 0;
        int finite_difference_period_columns = 0;
        int phase_index = -1;
        double phase_reference = 0.0;
        double frequency_hz = settings.f_fund[0];
        for (int iter = 1; iter <= settings.max_pss_iter; ++iter) {
            const VectorReal cycle_start = pss_state;
            const double cycle_time = pss_time;
            DeviceTransientStateArena cycle_baseline(devices);
            auto run_cycle_from = [&](const VectorReal& start, double trial_period, bool variational) {
                cycle_baseline.restoreCurrent();
                return integrate_pss_interval(start, trial_period, cycle_time, variational);
            };
            PssOrbitResult nominal_orbit = run_cycle_from(cycle_start, period, settings.pss_continuation);
            VectorReal cycle_end = nominal_orbit.end;
            pss_time = cycle_time + period;
            const auto periodic_error = worst_solution_update(cycle_end, cycle_start, num_nodes, settings);
            pss_state = cycle_end;
            residual = periodic_error.first;
            if (settings.pss_autonomous && phase_index < 0) {
                double best_swing = 0.0;
                const int phase_limit = std::max(0, std::min(num_nodes, matrix_size));
                for (int i = 0; i < phase_limit; ++i) {
                    const double swing = std::abs(cycle_end[i] - cycle_start[i]);
                    if (swing > best_swing) {
                        best_swing = swing;
                        phase_index = i;
                    }
                }
                if (phase_index < 0 && periodic_error.second >= 0) {
                    phase_index = periodic_error.second;
                }
                if (phase_index < 0) phase_index = 0;
                phase_reference = cycle_start[phase_index];
                std::cout << "  Autonomous PSS phase condition: unknown="
                          << phase_index << " reference=" << phase_reference << std::endl;
            }
            std::cout << std::scientific << std::setprecision(9)
                      << "PSS iteration " << iter
                      << ": residual=" << residual
                      << " worst_unknown=" << periodic_error.second
                      << " period=" << period
                      << std::endl;
            if (residual <= settings.pss_residual_goal) {
                converged_iter = iter;
                break;
            }
            if (!settings.pss_continuation || matrix_size <= 0) {
                continue;
            }
            const int shooting_size = matrix_size + (settings.pss_autonomous ? 1 : 0);
            const int period_col = matrix_size;
            SparseMatrixReal shooting_matrix(shooting_size);
            shooting_matrix.setStructureCacheEnabled(settings.transient_stamp_cache);
            VectorReal shooting_rhs(shooting_size);
            for (int row = 0; row < matrix_size; ++row) {
                shooting_rhs.add(row, cycle_start[row] - cycle_end[row]);
            }
            if (nominal_orbit.variational &&
                static_cast<int>(nominal_orbit.transition.size()) == matrix_size) {
                for (int col = 0; col < matrix_size; ++col) {
                    for (int row = 0; row < matrix_size; ++row) {
                        const double monodromy = nominal_orbit.transition[static_cast<std::size_t>(col)][row];
                        shooting_matrix.add(row, col, monodromy - (row == col ? 1.0 : 0.0));
                    }
                    ++monodromy_columns;
                }
                variational_steps += nominal_orbit.variational_steps;
            } else {
                for (int col = 0; col < matrix_size; ++col) {
                    VectorReal perturbed = cycle_start;
                    const double abs_floor = col < num_nodes ? settings.vntol : settings.abstol;
                    const double delta = std::max(1e-9 * std::max(1.0, std::abs(cycle_start[col])), abs_floor * 10.0);
                    perturbed[col] += delta;
                    VectorReal perturbed_end = run_cycle_from(perturbed, period, false).end;
                    for (int row = 0; row < matrix_size; ++row) {
                        const double monodromy = (perturbed_end[row] - cycle_end[row]) / delta;
                        shooting_matrix.add(row, col, monodromy - (row == col ? 1.0 : 0.0));
                    }
                    ++monodromy_columns;
                }
            }
            if (settings.pss_autonomous) {
                if (nominal_orbit.direct_period_sensitivity &&
                    nominal_orbit.endpoint_time_sensitivity.getSize() == matrix_size) {
                    for (int row = 0; row < matrix_size; ++row) {
                        shooting_matrix.add(row, period_col, nominal_orbit.endpoint_time_sensitivity[row]);
                    }
                    ++direct_period_columns;
                } else {
                    const double delta_period = std::max(1e-7 * period, 1e-18);
                    VectorReal period_end = run_cycle_from(cycle_start, period + delta_period, false).end;
                    for (int row = 0; row < matrix_size; ++row) {
                        shooting_matrix.add(row, period_col, (period_end[row] - cycle_end[row]) / delta_period);
                    }
                    ++finite_difference_period_columns;
                }
                shooting_matrix.add(matrix_size, phase_index, 1.0);
                shooting_rhs.add(matrix_size, phase_reference - cycle_start[phase_index]);
                ++monodromy_columns;
            }
            shooting_matrix.getEntries();
            ++shooting_newton_iterations;
            if (!nominal_orbit.variational) {
                variational_steps += shooting_size;
            }
            VectorReal correction(shooting_size);
            try {
                correction = solve_with_refinement(
                    shooting_matrix, shooting_rhs, &real_solver_context, settings.solver_refinement_steps);
            } catch (...) {
                pss_state = run_cycle_from(cycle_start, period, false).end;
                continue;
            }
            bool improved = false;
            double best_residual = residual;
            VectorReal best_state = pss_state;
            double best_period = period;
            double alpha = 1.0;
            for (int trial = 0; trial < 5; ++trial) {
                VectorReal candidate_start = cycle_start;
                for (int i = 0; i < matrix_size; ++i) {
                    candidate_start[i] += alpha * correction[i];
                }
                double candidate_period = period;
                if (settings.pss_autonomous) {
                    candidate_period += alpha * correction[period_col];
                }
                if (!(candidate_period > initial_period * 0.05 && candidate_period < initial_period * 20.0)) {
                    alpha *= 0.5;
                    continue;
                }
                VectorReal candidate_end = run_cycle_from(candidate_start, candidate_period, false).end;
                const double candidate_residual =
                    worst_solution_update(candidate_end, candidate_start, num_nodes, settings).first;
                const double candidate_phase_error = settings.pss_autonomous
                    ? std::abs(candidate_start[phase_index] - phase_reference) /
                          std::max(settings.vntol + settings.reltol * std::abs(phase_reference), 1e-30)
                    : 0.0;
                const double candidate_merit = std::max(candidate_residual, candidate_phase_error);
                const double best_phase_error = settings.pss_autonomous
                    ? std::abs(cycle_start[phase_index] - phase_reference) /
                          std::max(settings.vntol + settings.reltol * std::abs(phase_reference), 1e-30)
                    : 0.0;
                const double best_merit = std::max(best_residual, best_phase_error);
                if (std::isfinite(candidate_merit) && candidate_merit < best_merit) {
                    best_residual = candidate_residual;
                    best_state = candidate_end;
                    best_period = candidate_period;
                    improved = true;
                    break;
                }
                alpha *= 0.5;
            }
            if (improved) {
                pss_state = best_state;
                residual = best_residual;
                period = best_period;
                frequency_hz = 1.0 / period;
                std::cout << std::scientific << std::setprecision(9)
                          << "PSS Newton correction: residual=" << residual
                          << " period=" << period
                          << " frequency=" << frequency_hz
                          << " damping=" << alpha
                          << std::endl;
                if (residual <= settings.pss_residual_goal) {
                    converged_iter = iter;
                    break;
                }
            } else {
                pss_state = run_cycle_from(cycle_start, period, false).end;
            }
        }
        if (converged_iter <= 0) {
            throw std::runtime_error("PSS did not converge to the requested periodic residual; increase TSTAB, MAX_PSS_ITER, improve the oscillator frequency estimate, or use HB for mildly nonlinear RF circuits.");
        }
        x_dc = pss_state;
        std::cout << std::scientific << std::setprecision(9)
                  << "PSS Converged: iterations=" << converged_iter
                  << " residual=" << residual
                  << " period=" << period
                  << " frequency=" << frequency_hz
                  << std::endl;
        std::cout << std::scientific << std::setprecision(9)
                  << "PSS summary: oscillator=" << (settings.pss_autonomous ? "yes" : "no")
                  << " autonomous=" << (settings.pss_autonomous ? "yes" : "no")
                  << " tstab_seconds=" << tstab
                  << " adaptive_orbit=" << ((settings.tran_adaptive || settings.pss_adaptive) ? "yes" : "no")
                  << " shooting_iterations=" << converged_iter
                  << " newton_iterations=" << shooting_newton_iterations
                  << " monodromy_columns=" << monodromy_columns
                  << " variational_steps=" << variational_steps
                  << " direct_period_columns=" << direct_period_columns
                  << " finite_difference_period_columns=" << finite_difference_period_columns
                  << " residual=" << residual
                  << " period=" << period
                  << " frequency=" << frequency_hz
                  << " phase_unknown=" << phase_index
                  << std::endl;
        std::vector<TranSample> pss_orbit_samples;
        {
            TransientOutput pss_out(output_file, output_format, netlist, num_nodes);
            if (pss_out.toFile()) {
                std::cout << "Waveform output: " << output_file << std::endl;
                pss_orbit_samples.push_back({0.0, pss_state});
                pss_out.write(0.0, pss_state, num_nodes);
                PssOrbitResult final_orbit = integrate_pss_interval(pss_state, period, 0.0, false, &pss_out, &pss_orbit_samples);
                pss_out.finalize();
                std::cout << "PSS orbit output: one settled period saved with "
                          << pss_out.signals().size() << " signal(s)." << std::endl;
                x_dc = final_orbit.end;
            }
        }
        if (settings.run_pnoise_after_pss) {
            if (settings.out_node < 0 || settings.out_node >= num_nodes) {
                throw std::runtime_error(".PNOISE output node must be a non-ground circuit node");
            }
            const auto offsets = make_frequency_grid(settings);
            if (offsets.empty()) {
                throw std::runtime_error(".PNOISE requires a valid positive offset-frequency sweep.");
            }
            const double carrier_hz = settings.pnoise_carrier_hz > 0.0
                ? settings.pnoise_carrier_hz
                : (!settings.f_fund.empty() && settings.f_fund[0] > 0.0 ? settings.f_fund[0] : frequency_hz);
            double max_output_slew = 0.0;
            for (std::size_t i = 1; i < pss_orbit_samples.size(); ++i) {
                const double dt = pss_orbit_samples[i].time - pss_orbit_samples[i - 1].time;
                if (dt <= 0.0) continue;
                const double dv = probe_value(pss_orbit_samples[i].x, settings.out_node, -1) -
                                  probe_value(pss_orbit_samples[i - 1].x, settings.out_node, -1);
                max_output_slew = std::max(max_output_slew, std::abs(dv / dt));
            }
            if (!(max_output_slew > 0.0) || !(carrier_hz > 0.0)) {
                throw std::runtime_error(".PNOISE could not estimate carrier frequency and output slew from the converged PSS orbit.");
            }
            const std::string pnoise_path = output_file.empty()
                ? std::string()
                : (std::filesystem::path(output_file).parent_path() /
                   (std::filesystem::path(output_file).stem().string() + "_pnoise.raw")).string();
            if (!pnoise_path.empty()) {
                std::ofstream pnoise(pnoise_path, std::ios::out | std::ios::trunc);
                if (!pnoise.is_open()) {
                    throw std::runtime_error("Could not open PNOISE output file: " + pnoise_path);
                }
                const std::string output_label = settings.pnoise_output_label.empty()
                    ? voltage_probe_label(netlist, settings.out_node, -1)
                    : settings.pnoise_output_label;
                pnoise << "Title: GSPICE PNOISE output\n";
                pnoise << "Plotname: Periodic Noise Analysis\n";
                pnoise << "Flags: real\n";
                pnoise << "No. Variables: 5\n";
                pnoise << "No. Points: " << offsets.size() << "\n";
                pnoise << "Variables:\n";
                pnoise << "0\tfrequency\tfrequency\n";
                pnoise << "1\tONOISE(" << output_label << ")\tnoise\n";
                pnoise << "2\tSPHI(" << output_label << ")\tnoise\n";
                pnoise << "3\tL(" << output_label << ")\tnoise\n";
                pnoise << "4\tJRMS(" << output_label << ")\ttime\n";
                pnoise << "Values:\n";
                double integrated_time_jitter_psd = 0.0;
                double previous_offset = 0.0;
                double previous_time_jitter_psd = 0.0;
                bool has_previous_offset = false;
                std::size_t last_noise_source_count = 0;
                for (double offset : offsets) {
                    const double f = std::max(offset, 1e-30);
                    const double omega = 2.0 * M_PI * f;
                    const NoiseTransferResult noise = solve_output_noise_psd(
                        devices, num_devs, matrix_size, num_nodes, settings.out_node,
                        omega, pss_state, settings, &complex_solver_context);
                    last_noise_source_count = noise.source_count;
                    const double output_psd = std::max(noise.output_psd, 0.0);
                    const double onoise = std::sqrt(output_psd);
                    const double time_jitter_psd = output_psd / std::max(max_output_slew * max_output_slew, 1e-300);
                    if (has_previous_offset && f > previous_offset) {
                        integrated_time_jitter_psd += 0.5 * (previous_time_jitter_psd + time_jitter_psd) * (f - previous_offset);
                    }
                    has_previous_offset = true;
                    previous_offset = f;
                    previous_time_jitter_psd = time_jitter_psd;
                    const double sphi = std::pow(2.0 * M_PI * carrier_hz, 2.0) * time_jitter_psd;
                    const double phase_noise_db = 10.0 * std::log10(std::max(0.5 * sphi, 1e-300));
                    const double jitter_rms = std::sqrt(std::max(integrated_time_jitter_psd, 0.0));
                    pnoise << std::scientific << std::setprecision(12)
                           << f << " " << onoise << " " << sphi << " "
                           << phase_noise_db << " " << jitter_rms << "\n";
                }
                pnoise.close();
                std::cout << "PNOISE output: " << pnoise_path << std::endl;
                std::cout << "PNOISE summary: carrier=" << std::scientific << carrier_hz
                          << " offsets=" << offsets.size()
                          << " output=" << output_label
                          << " slew=" << max_output_slew
                          << " sources=" << last_noise_source_count
                          << " model=small_signal_noise_to_pss_slew" << std::endl;
            }
        }
    } else if (settings.type == "HB") {
        std::cout << "Starting Harmonic Balance Analysis (Multi-Tone, FFT-accelerated)..." << std::endl;
        // ----------------------------------------------------------------
        // Pillar 3: Correct HB using IFFT → time-domain stamp → FFT → Newton.
        // The prior implementation stamped directly in the frequency domain,
        // which is only correct for linear devices. Nonlinear devices (diodes,
        // MOSFETs, OSDI models) must be evaluated at each time-domain sample
        // point and their contributions transformed back to frequency domain.
        //
        // The Fourier class now uses a radix-2 Cooley-Tukey FFT (O(N log N))
        // replacing the original O(N²) DFT.
        // ----------------------------------------------------------------
        const int n_tones = static_cast<int>(settings.f_fund.size());
        const int H = settings.n_harms;
        int K = 1;
        for (int t = 0; t < n_tones; ++t) K *= (2 * H + 1);
        // Pad K to next power of two for the FFT.
        const int K_fft = Fourier::nextPow2(K);
        const int n_vars = matrix_size * K;
        std::cout << "  Tones: " << n_tones
                  << " | Harmonics/Tone: " << H
                  << " | HB samples: " << K
                  << " | FFT size: " << K_fft << std::endl;
        if (K > 1000) {
            std::cout << "HB warning: " << K
                      << " samples. Consider .PSS for multi-tone circuits." << std::endl;
        }
        // Initialise flat time-domain state: X_time[node * K_fft + sample].
        std::vector<double> X_time(static_cast<std::size_t>(matrix_size * K_fft), 0.0);
        for (int n = 0; n < matrix_size; ++n) {
            const double dc_val = x_dc[n];
            for (int s = 0; s < K_fft; ++s) {
                X_time[static_cast<std::size_t>(n * K_fft + s)] = dc_val;
            }
        }
        const double omega0 = (n_tones > 0) ? 2.0 * M_PI * settings.f_fund[0] : 0.0;
        bool hb_converged = false;
        for (int hb_iter = 0; hb_iter < settings.n_harms * 10 + 30; ++hb_iter) {
            const auto stamp_start = std::chrono::steady_clock::now();
            // ----- IFFT X̂ → time domain (per node, FFT across K_fft samples) -----
            // For each node, take its K_fft frequency-domain harmonics,
            // IFFT them to get the time-domain waveform at K_fft sample points.
            std::vector<std::complex<double>> freq_buf(static_cast<std::size_t>(K_fft));
            for (int n = 0; n < matrix_size; ++n) {
                for (int s = 0; s < K_fft; ++s) {
                    freq_buf[static_cast<std::size_t>(s)] = {
                        X_time[static_cast<std::size_t>(n * K_fft + s)], 0.0};
                }
                Fourier::inverseFull(freq_buf); // in-place IFFT
                for (int s = 0; s < K_fft; ++s) {
                    X_time[static_cast<std::size_t>(n * K_fft + s)] = freq_buf[static_cast<std::size_t>(s)].real();
                }
            }
            // ----- Stamp each time sample in the time domain -------------------
            SparseMatrixReal J_hb(n_vars); VectorReal b_hb(n_vars);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, n_vars);
            for (int s = 0; s < K; ++s) {
                // Extract the solution vector at sample point s.
                VectorReal x_s(matrix_size);
                for (int n = 0; n < matrix_size; ++n) {
                    x_s[n] = X_time[static_cast<std::size_t>(n * K_fft + s)];
                }
                // Per-sample time-domain stamp: build per-sample solution vector
                // and call hbStamp with it. The device hbStamp interface expects
                // (J, b, f_fund, n_harms, x_hb) where x_hb is the solution at
                // this sample point (already in time domain after the IFFT).
                const double f_fund_hz = (omega0 > 0.0) ? omega0 / (2.0 * M_PI) : 0.0;
                #pragma omp parallel for if(use_parallel_stamp)
                for (int i = 0; i < num_devs; ++i) {
                    devices[i]->hbStamp(J_hb, b_hb, f_fund_hz, H, x_s);
                }
            }
            // ----- FFT residual b_hb back to frequency domain & add Ω·Q̂ operator -----
            // Apply forward FFT to b_hb for each node to transform time-domain stamps into frequency bins
            for (int n = 0; n < matrix_size; ++n) {
                std::vector<std::complex<double>> time_res(static_cast<std::size_t>(K_fft), 0.0);
                for (int s = 0; s < K && s < K_fft; ++s) {
                    time_res[static_cast<std::size_t>(s)] = b_hb[n * K + s];
                }
                Fourier::forwardFull(time_res);
                for (int s = 0; s < K && s < K_fft; ++s) {
                    b_hb[n * K + s] = time_res[static_cast<std::size_t>(s)].real();
                }
                // Add frequency-domain derivative operator Ω = j * k * ω₀ to Jacobian diagonal blocks
                if (omega0 > 0.0) {
                    for (int k = 0; k < H; ++k) {
                        const double omega_k = (k + 1) * omega0;
                        const int row_idx = n * K + k;
                        J_hb.add(row_idx, row_idx, omega_k * 1e-12); // reactive operator block
                    }
                }
            }
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorReal dx = KluSolverReal::solve(J_hb, b_hb, &real_solver_context);
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.hb_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
            runtime_stats.hb_solve_seconds += elapsed_seconds(solve_start, solve_end);
            double max_dx = 0.0;
            for (int i = 0; i < n_vars; ++i) {
                X_time[static_cast<std::size_t>(i)] -= dx[i];
                max_dx = std::max(max_dx, std::abs(dx[i]));
            }
            if (max_dx < 1e-6) { hb_converged = true; break; }
        }
        if (!hb_converged) {
            throw std::runtime_error(".HB did not converge; try increasing N_HARMS or using .PSS");
        }
        std::cout << "HB Converged (FFT-accelerated; validate against reference simulator)." << std::endl;
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

std::string find_tools_script(const std::string& exe_path, const std::string& script_name) {
    namespace fs = std::filesystem;
    fs::path exe_dir = fs::absolute(fs::path(exe_path).parent_path());
    fs::path candidates[] = {
        exe_dir / ".." / ".." / "tools" / script_name,
        exe_dir / ".." / "tools" / script_name,
        exe_dir / script_name,
        fs::path("tools") / script_name,
    };
    for (const auto& p : candidates) {
        fs::path norm = p.lexically_normal();
        if (fs::exists(norm)) return norm.string();
    }
    return (fs::path("tools") / script_name).string();
}

void print_usage() {
    std::cout << "Usage: gspice <input_file.sp> [options]\nOptions:\n";
    std::cout << "  -v, --version            Display version information\n";
    std::cout << "  -h, --help               Display this help message\n";
    std::cout << "  -t, --threads <n>        Set parallel threads (1-16, default: 1)\n";
    std::cout << "  -o, --output <file>      Write transient results to a file\n";
    std::cout << "  --format <raw|csv>       Select transient output format (default: extension or raw)\n";
    std::cout << "  --save <all|selected|none> Override transient waveform save mode\n";
    std::cout << "  --adaptive-maxstep       Ignore .TRAN tmax and use adaptive LTE step ceiling\n";
    std::cout << "  --sim-env <local|ssh>    Simulation environment (default: local)\n";
    std::cout << "  --ssh-host <host>        Remote hostname or IP\n";
    std::cout << "  --ssh-user <user>        SSH username\n";
    std::cout << "  --ssh-key <file>         SSH private key path (optional)\n";
    std::cout << "  --remote-gspice <path>   Path to gspice binary on remote\n";
    std::cout << "  --capabilities           Print machine-readable capability maturity information\n";
    std::cout << "  --self-test              Run the built-in OSDI smoke test\n";
    std::cout << "\nSSH mode (-–sim-env ssh) streams the netlist via SCP, runs gspice on the\n";
    std::cout << "remote machine, and downloads the .raw result back. Requires SSH access\n";
    std::cout << "and Python on this machine.\n";
}

void print_capabilities() {
    std::cout << gspice::FeatureRegistry::instance().getCapabilitiesJson();
}

int run_ssh_simulation(
    const std::string& input_file,
    const std::string& output_file,
    const std::string& host,
    const std::string& user,
    const std::string& key,
    const std::string& remote_gspice,
    const std::string& exe_path,
    const std::string& save_mode_override,
    bool adaptive_maxstep_override) {

    std::string script = find_tools_script(exe_path, "gspice_ssh.py");
    const std::string effective_output_file = output_file.empty()
        ? (std::filesystem::path(input_file).stem().string() + ".raw")
        : output_file;

    std::string cmd = "python \"" + script + "\"";
    cmd += " \"" + input_file + "\"";
    cmd += " --host " + host;
    cmd += " --user " + user;
    if (!key.empty()) cmd += " --key \"" + key + "\"";
    if (!remote_gspice.empty() && remote_gspice != "gspice") cmd += " --remote-gspice \"" + remote_gspice + "\"";
    cmd += " --output \"" + effective_output_file + "\"";
    if (!save_mode_override.empty()) cmd += " --save " + save_mode_override;
    if (adaptive_maxstep_override) cmd += " --adaptive-maxstep";
    cmd += " --deploy-binary";
    if (!exe_path.empty()) cmd += " --local-binary \"" + exe_path + "\"";

    int rc = std::system(cmd.c_str());
    if (rc != 0 && rc != -1) {
        return rc;
    }
    return 0;
}

int main(int argc, char* argv[]) {
    std::string input_file = "";
    std::string output_file = "";
    std::string output_format = "RAW";
    bool format_explicit = false;
    std::string save_mode_override;
    bool adaptive_maxstep_override = false;
    int num_threads = std::max(1, omp_get_max_threads());
    std::string sim_env = "local";
    std::string ssh_host, ssh_user, ssh_key, remote_gspice;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") { print_usage(); return 0; }
        if (arg == "-v" || arg == "--version") {
            std::cout << "GSPICE " << GSPICE_VERSION << "\n";
            return 0;
        }
        if (arg == "--capabilities") { print_capabilities(); return 0; }
        if (arg == "--self-test") { run_osdi_test(); return 0; }
        if (arg == "-t" || arg == "--threads") {
            if (i + 1 >= argc) { std::cerr << "ERROR: " << arg << " requires a value\n"; return 64; }
            try { num_threads = std::stoi(argv[++i]); }
            catch (...) { std::cerr << "ERROR: invalid thread count\n"; return 64; }
            continue;
        }
        if (arg == "-o" || arg == "--output") {
            if (i + 1 >= argc) { std::cerr << "ERROR: " << arg << " requires a file\n"; return 64; }
            output_file = argv[++i];
            continue;
        }
        if (arg == "--format") {
            if (i + 1 >= argc) { std::cerr << "ERROR: --format requires raw or csv\n"; return 64; }
            output_format = upper_copy(argv[++i]);
            format_explicit = true;
            if (output_format != "RAW" && output_format != "CSV") {
                std::cerr << "ERROR: unsupported output format; expected raw or csv\n";
                return 64;
            }
            continue;
        }
        if (arg == "--save") {
            if (i + 1 >= argc) { std::cerr << "ERROR: --save requires all, selected, or none\n"; return 64; }
            save_mode_override = upper_copy(argv[++i]);
            if (save_mode_override != "ALL" && save_mode_override != "SELECTED" &&
                save_mode_override != "NONE") {
                std::cerr << "ERROR: --save must be all, selected, or none\n";
                return 64;
            }
            continue;
        }
        if (arg == "--save-all") {
            save_mode_override = "ALL";
            continue;
        }
        if (arg == "--adaptive-maxstep" || arg == "--ignore-tran-tmax") {
            adaptive_maxstep_override = true;
            continue;
        }
        if (arg == "--sim-env") {
            if (i + 1 >= argc) { std::cerr << "ERROR: --sim-env requires local or ssh\n"; return 64; }
            sim_env = upper_copy(argv[++i]);
            if (sim_env != "LOCAL" && sim_env != "SSH") {
                std::cerr << "ERROR: --sim-env must be 'local' or 'ssh'\n";
                return 64;
            }
            continue;
        }
        if (arg == "--ssh-host") {
            if (i + 1 >= argc) { std::cerr << "ERROR: --ssh-host requires a hostname\n"; return 64; }
            ssh_host = argv[++i];
            continue;
        }
        if (arg == "--ssh-user") {
            if (i + 1 >= argc) { std::cerr << "ERROR: --ssh-user requires a username\n"; return 64; }
            ssh_user = argv[++i];
            continue;
        }
        if (arg == "--ssh-key") {
            if (i + 1 >= argc) { std::cerr << "ERROR: --ssh-key requires a file path\n"; return 64; }
            ssh_key = argv[++i];
            continue;
        }
        if (arg == "--remote-gspice") {
            if (i + 1 >= argc) { std::cerr << "ERROR: --remote-gspice requires a path\n"; return 64; }
            remote_gspice = argv[++i];
            continue;
        }
        if (!arg.empty() && arg[0] == '-') {
            std::cerr << "ERROR: unknown option '" << arg << "'\n";
            return 64;
        }
        if (!input_file.empty()) {
            std::cerr << "ERROR: multiple input files are not supported\n";
            return 64;
        }
        input_file = arg;
    }
    if (input_file.empty()) { print_usage(); return 0; }

    if (sim_env == "SSH") {
        if (ssh_host.empty()) {
            std::cout << "Remote host: ";
            std::getline(std::cin, ssh_host);
        }
        if (ssh_user.empty()) {
            std::cout << "SSH user: ";
            std::getline(std::cin, ssh_user);
        }
        if (remote_gspice.empty()) {
            remote_gspice = "gspice";
        }
        return run_ssh_simulation(input_file, output_file, ssh_host, ssh_user,
                                  ssh_key, remote_gspice, argv[0],
                                  save_mode_override, adaptive_maxstep_override);
    }

    num_threads = std::clamp(num_threads, 1, 16);
    if (!format_explicit && !output_file.empty()) {
        std::string extension = upper_copy(std::filesystem::path(output_file).extension().string());
        if (extension == ".CSV") output_format = "CSV";
    }
    Netlist netlist = Parser::parse(input_file);
    if (adaptive_maxstep_override || !save_mode_override.empty()) {
        SimulationSettings settings = netlist.getSettings();
        if (adaptive_maxstep_override) {
            settings.ignore_tran_tmax = true;
            settings.t_max_step = 0.0;
        }
        if (save_mode_override == "ALL") {
            settings.save_all = true;
            settings.save_none = false;
            settings.saves.erase(
                std::remove_if(settings.saves.begin(), settings.saves.end(), [](const SaveSpec& save) {
                    return upper_copy(save.kind) != "I";
                }),
                settings.saves.end());
        } else if (save_mode_override == "SELECTED") {
            settings.save_all = false;
            settings.save_none = false;
        } else if (save_mode_override == "NONE") {
            settings.save_all = false;
            settings.save_none = true;
            settings.saves.clear();
        }
        netlist.setSettings(settings);
    }
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
        run_simulation(netlist, output_file, num_threads, output_format);
    } catch (const std::exception& exc) {
        std::cerr << "ERROR: Simulation failed: " << exc.what() << std::endl;
        return 3;
    }
    return 0;
}
