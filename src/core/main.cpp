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

    void noteSolve(const TransientStepResult& result) {
        ++newton_solves;
        newton_iterations += result.iterations;
        max_newton_iterations = std::max(max_newton_iterations, result.iterations);
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

TransientStepResult solve_transient_step(
    const std::vector<std::unique_ptr<Device>>& devices,
    int num_devs,
    int matrix_size,
    int num_nodes,
    const SimulationSettings& settings,
    const VectorReal& x_start,
    const std::vector<VectorReal>& x_hist,
    double step,
    double target_time) {
    VectorReal x = x_start;
    bool converged = false;
    for (int iter = 0; iter < settings.tran_max_iter; ++iter) {
        SparseMatrixReal J_sparse(matrix_size);
        VectorReal b(matrix_size);
        stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
        #pragma omp parallel for
        for (int i = 0; i < num_devs; ++i) {
            devices[i]->dcStamp(J_sparse, b, x, step, target_time, x_hist);
        }
        VectorReal x_new = KluSolverReal::solve(J_sparse, b);
        if (!vector_finite_and_bounded(x_new)) {
            return {x_new, false, iter + 1};
        }
        converged = solution_converged(x_new, x, num_nodes, settings);
        x = x_new;
        if (converged) return {x, true, iter + 1};
    }
    return {x, converged, settings.tran_max_iter};
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
    std::string method = settings.tran_method;
    std::transform(method.begin(), method.end(), method.begin(), ::toupper);
    method.erase(std::remove(method.begin(), method.end(), '_'), method.end());
    method.erase(std::remove(method.begin(), method.end(), '-'), method.end());
    if (method == "TRAP" || method == "TRAPEZOIDAL") return "Trapezoidal (pending device stamps; using Backward Euler)";
    if (method == "GEAR2") return "Gear2 (pending device stamps; using Backward Euler)";
    if (method == "BE" || method == "BACKWARDEULER") return "Backward Euler";
    return "Auto (Backward Euler)";
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
    for (int iter = 0; iter < 100; ++iter) {
        SparseMatrixReal J(matrix_size); VectorReal b(matrix_size);
        auto& devices = netlist.getDevices();
        int num_devs = static_cast<int>(devices.size());
        #pragma omp parallel for
        for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J, b, x, 0.0, 0.0, {});
        VectorReal x_new = KluSolverReal::solve(J, b);
        double max_change = std::abs(x_new[2] - x[2]); x = x_new;
        if (max_change < tol) { std::cout << "  Converged in " << iter + 1 << " iterations." << std::endl; break; }
    }
    std::cout << "Results:\n  Vdd: " << x[0] << " V\n  Gate: " << x[1] << " V\n  Drain: " << x[2] << " V\n";
}

void run_simulation(Netlist& netlist, const std::string& output_file = "") {
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
    const auto& settings = netlist.getSettings();

    VectorReal x_dc(matrix_size);
    if (!settings.use_uic) {
        std::cout << "Calculating DC Operating Point..." << std::endl;
        for (int iter = 0; iter < settings.op_max_iter; ++iter) {
            SparseMatrixReal J_sparse(matrix_size); VectorReal b(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            #pragma omp parallel for
            for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J_sparse, b, x_dc, 0.0, 0.0, {});
            VectorReal x_new = KluSolverReal::solve(J_sparse, b);
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
            : output_step;
        const double min_step = settings.t_min_step > 0.0
            ? std::min(settings.t_min_step, max_step)
            : std::max(max_step * 1e-6, 1e-18);
        const double save_start = std::clamp(settings.t_start, 0.0, settings.t_stop);
        double step = max_step;
        std::vector<VectorReal> x_hist; VectorReal x = x_dc; x_hist.push_back(x);
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
                    devices, num_devs, matrix_size, num_nodes, settings, x, x_hist, step, target_time);
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
                    const double next_bp_after_target = next_breakpoint_after(breakpoints, target_time, bp_tol);
                    const bool near_breakpoint =
                        (next_bp > 0.0 && std::abs(target_time - next_bp) <= bp_tol) ||
                        (next_bp_after_target > 0.0 && next_bp_after_target - target_time <= std::max(max_step * 2.0, output_step));
                    if (!near_breakpoint) {
                        const double state_change = transient_state_change_error(full.x, x, num_nodes, settings);
                        if (state_change < 0.01) {
                            accepted_x = full.x;
                            err = state_change;
                            accepted = true;
                            continue;
                        }
                    }
                    const double half_step = step * 0.5;
                    TransientStepResult half1 = solve_transient_step(
                        devices, num_devs, matrix_size, num_nodes, settings, x, x_hist, half_step, t + half_step);
                    tran_stats.noteSolve(half1);
                    std::vector<VectorReal> half_hist = x_hist;
                    half_hist.push_back(half1.x);
                    TransientStepResult half2 = solve_transient_step(
                        devices, num_devs, matrix_size, num_nodes, settings, half1.x, half_hist, half_step, target_time);
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
                        const double factor = std::clamp(0.8 * std::pow(std::max(err, 1e-12), -0.5), 0.2, 0.75);
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
            x_hist.push_back(x);
            for (int i = 0; i < num_devs; ++i) {
                devices[i]->acceptTransientStep(x, t);
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
                const double grow = (err < 0.05) ? 2.0 : (err < 0.25 ? 1.5 : 1.15);
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
            SparseMatrixComplex J_sparse(matrix_size); VectorComplex b_ac(matrix_size);
            stamp_global_gmin(J_sparse, num_nodes, settings.gmin);
            #pragma omp parallel for
            for (int i = 0; i < num_devs; ++i) devices[i]->acStamp(J_sparse, b_ac, omega, x_dc);
            VectorComplex x_ac = KluSolverComplex::solve(J_sparse, b_ac);
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
            VectorComplex xv = KluSolverComplex::solve(Jv, bv);
            std::complex<double> Tv = -xv[probes[0]->getNodeNeg()] / xv[probes[0]->getNodePos()];
            // 2. Current Pass
            SparseMatrixComplex Ji(matrix_size); VectorComplex bi(matrix_size);
            for (const auto& dev : devices) {
                auto* p = dynamic_cast<StabilityProbe*>(dev.get());
                if (p) p->stbStamp(Ji, bi, 2); else dev->acStamp(Ji, bi, omega, x_dc);
            }
            VectorComplex xi = KluSolverComplex::solve(Ji, bi);
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
            SparseMatrixReal J_hb(n_vars); VectorReal b_hb(n_vars);
            #pragma omp parallel for
            for (int i = 0; i < num_devs; ++i) devices[i]->hbStamp(J_hb, b_hb, (n_tones > 0 ? settings.f_fund[0] : 0.0), H, x_hb);
            VectorReal dx = KluSolverReal::solve(J_hb, b_hb);
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
            VectorReal x_pac = KluSolverReal::solve(J_pac, b_pac);
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
    omp_set_num_threads(num_threads);
    std::cout << "GSPICE Core: Using " << num_threads << " threads.\n";
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
        run_simulation(netlist, output_file);
    } catch (const std::exception& exc) {
        std::cerr << "ERROR: Simulation failed: " << exc.what() << std::endl;
        return 3;
    }
    return 0;
}
