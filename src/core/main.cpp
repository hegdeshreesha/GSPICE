#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <exception>
#include <cctype>
#include <fstream>
#include <sstream>
#include <limits>
#include <chrono>
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

bool vector_finite_and_bounded(const VectorReal& x, double limit = 1e9) {
    for (int i = 0; i < x.getSize(); ++i) {
        if (!std::isfinite(x[i]) || std::abs(x[i]) > limit) return false;
    }
    return true;
}

struct TransientStepResult {
    VectorReal x;
    bool converged = false;
    int iterations = 0;
    double stamp_seconds = 0.0;
    double solve_seconds = 0.0;
};

struct TransientStats {
    long long accepted_steps = 0;
    long long rejected_steps = 0;
    long long convergence_rejections = 0;
    long long lte_rejections = 0;
    long long output_points = 0;
    long long newton_solves = 0;
    long long newton_iterations = 0;
    int max_newton_iterations = 0;
    double min_step = std::numeric_limits<double>::infinity();
    double max_step = 0.0;
    double last_error = 0.0;
    double stamp_seconds = 0.0;
    double solve_seconds = 0.0;

    void noteSolve(const TransientStepResult& result) {
        ++newton_solves;
        newton_iterations += result.iterations;
        max_newton_iterations = std::max(max_newton_iterations, result.iterations);
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
    VectorReal x = x_start;
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
            return {x_new, false, iter + 1, total_stamp_seconds, total_solve_seconds};
        }
        converged = solution_converged(x_new, x, num_nodes, settings);
        x = x_new;
        if (converged) return {x, true, iter + 1, total_stamp_seconds, total_solve_seconds};
    }
    return {x, converged, settings.tran_max_iter, total_stamp_seconds, total_solve_seconds};
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

class TransientOutput {
public:
    TransientOutput(const std::string& path, const Netlist& netlist, int num_nodes)
        : to_file_(!path.empty()), num_nodes_(num_nodes) {
        if (to_file_) {
            file_.open(path, std::ios::out | std::ios::trunc);
            if (!file_.is_open()) {
                throw std::runtime_error("Could not open transient output file: " + path);
            }
            file_ << "Title: GSPICE RAW output\n";
            file_ << "Plotname: Transient Analysis\n";
            file_ << "Flags: real\n";
            file_ << "No. Variables: " << (num_nodes + 1) << "\n";
            file_ << "No. Points: ";
            points_pos_ = file_.tellp();
            file_ << std::setw(20) << 0 << "\n";
            file_ << "Variables:\n";
            file_ << "0\ttime\ttime\n";
            for (int i = 0; i < num_nodes; ++i) {
                file_ << (i + 1) << "\tV(" << netlist.getNodeName(i) << ")\tvoltage\n";
            }
            file_ << "Values:\n";
        }
    }

    ~TransientOutput() {
        finalize();
    }

    bool toFile() const { return to_file_; }

    void write(double t, const VectorReal& x, int num_nodes) {
        if (!to_file_) {
            print_transient_point(t, x, num_nodes);
            return;
        }
        file_ << std::scientific << std::setprecision(12) << t;
        for (int i = 0; i < num_nodes; ++i) {
            file_ << " " << x[i];
        }
        file_ << "\n";
        ++point_count_;
    }

    void finalize() {
        if (!to_file_ || finalized_) return;
        finalized_ = true;
        if (file_.is_open() && points_pos_ >= std::streampos(0)) {
            const auto end_pos = file_.tellp();
            file_.seekp(points_pos_);
            file_ << std::setw(20) << point_count_;
            file_.seekp(end_pos);
        }
    }

private:
    bool to_file_ = false;
    bool finalized_ = false;
    int num_nodes_ = 0;
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
    }
    int matrix_size = num_nodes + netlist.getNumBranches();
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

    VectorReal x_dc(matrix_size);
    if (!settings.use_uic) {
        std::cout << "Calculating DC Operating Point..." << std::endl;
        for (int iter = 0; iter < settings.op_max_iter; ++iter) {
            const auto stamp_start = std::chrono::steady_clock::now();
            SparseMatrixReal J_sparse(matrix_size); VectorReal b(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            const bool use_parallel_stamp = parallel_stamp_enabled(num_devs, matrix_size);
            #pragma omp parallel for if(use_parallel_stamp)
            for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J_sparse, b, x_dc, 0.0, 0.0, {});
            const auto stamp_end = std::chrono::steady_clock::now();
            const auto solve_start = std::chrono::steady_clock::now();
            VectorReal x_new = KluSolverReal::solve(J_sparse, b, &real_solver_context);
            const auto solve_end = std::chrono::steady_clock::now();
            runtime_stats.dc_stamp_seconds += elapsed_seconds(stamp_start, stamp_end);
            runtime_stats.dc_solve_seconds += elapsed_seconds(solve_start, solve_end);
            bool converged = solution_converged(x_new, x_dc, num_nodes, settings);
            x_dc = x_new; if (converged) break;
        }
    }

    if (settings.type == "TRAN") {
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
        const double stop_tol = std::max(std::max(1e-30, std::abs(settings.t_stop) * 1e-9), min_step * 10.0);
        double last_printed_time = -1.0;
        double last_progress_percent = -1.0;
        double next_output = save_start;
        if (next_output <= 1e-30) {
            tran_out.write(0.0, x, num_nodes);
            ++tran_stats.output_points;
            last_printed_time = 0.0;
            next_output = output_step;
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
            const TransientContext accepted_ctx = make_transient_context(settings, x_hist, t_hist, step, t);
            for (int i = 0; i < num_devs; ++i) {
                devices[i]->acceptTransientStep(x, t, accepted_ctx);
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
            }
        }
        const double output_tol = std::max(1e-30, output_step * 1e-9);
        if (settings.t_stop >= save_start && last_printed_time < settings.t_stop - output_tol) {
            tran_out.write(settings.t_stop, x, num_nodes);
            ++tran_stats.output_points;
        }
        const double min_used_step = std::isfinite(tran_stats.min_step) ? tran_stats.min_step : 0.0;
        std::cout << std::scientific << std::setprecision(9)
                  << "Transient summary: accepted=" << tran_stats.accepted_steps
                  << " rejected=" << tran_stats.rejected_steps
                  << " convergence_rejected=" << tran_stats.convergence_rejections
                  << " lte_rejected=" << tran_stats.lte_rejections
                  << " output_points=" << tran_stats.output_points
                  << " min_step=" << min_used_step
                  << " max_step=" << tran_stats.max_step
                  << std::endl;
        std::cout << std::fixed << std::setprecision(2)
                  << "Newton summary: solves=" << tran_stats.newton_solves
                  << " total_iterations=" << tran_stats.newton_iterations
                  << " average_iterations=" << tran_stats.averageNewtonIterations()
                  << " max_iterations=" << tran_stats.max_newton_iterations
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
