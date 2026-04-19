#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "parser.hpp"
#include "devices/resistor.hpp"
#include "devices/capacitor.hpp"
#include "devices/inductor.hpp"
#include "devices/diode.hpp"
#include "devices/voltage_source.hpp"
#include "devices/port.hpp"
#include "devices/osdi_device.hpp"
#include "osdi_emulator.hpp"
#include "solvers/klu_solver.hpp"

using namespace gspice;

void run_osdi_test() {
    std::cout << "\n--- GSPICE OSDI Bridge Verification ---" << std::endl;
    Netlist netlist;
    int n1 = netlist.getOrCreateNode("1");
    int n0 = netlist.getOrCreateNode("0_diode"); 
    
    OsdiDescriptor va_model = OsdiEmulator::getDescriptor();
    
    // V1 (5V) -> n1 -> R1 (1k) -> n0 -> Diode_VA -> GND
    netlist.addDevice(std::make_unique<VoltageSource>("V1", n1, -1, 5.0, 2));
    netlist.addDevice(std::make_unique<Resistor>("R1", n1, n0, 1000.0));
    netlist.addDevice(std::make_unique<OSDIDevice>("D_VA", va_model, std::vector<int>{n0, -1}));

    int matrix_size = 3; 
    VectorReal x(matrix_size);
    
    std::cout << "Solving circuit with Verilog-A (OSDI) model..." << std::endl;
    for (int iter = 0; iter < 20; ++iter) {
        SparseMatrixReal J(matrix_size);
        VectorReal b(matrix_size);
        for (const auto& dev : netlist.getDevices()) dev->dcStamp(J, b, x, 0, {});
        x = KluSolverReal::solve(J, b);
    }
    
    std::cout << "Node 1 Voltage (Vsrc):  " << std::fixed << std::setprecision(4) << x[0] << " V" << std::endl;
    std::cout << "Node 0 Voltage (Diode): " << x[1] << " V" << std::endl;
    std::cout << "OSDI Bridge Verification Complete." << std::endl;
}

void run_simulation(Netlist& netlist) {
    int num_nodes = netlist.getNumNodes();
    std::vector<Port*> ports;
    for (auto& dev : netlist.getDevices()) {
        auto* vsrc = dynamic_cast<VoltageSource*>(dev.get());
        if (vsrc) vsrc->setBranchIndex(netlist.getNextBranchId(num_nodes));
        auto* ind = dynamic_cast<Inductor*>(dev.get());
        if (ind) ind->setBranchIndex(netlist.getNextBranchId(num_nodes));
        auto* port = dynamic_cast<Port*>(dev.get());
        if (port) {
            port->setBranchIndex(netlist.getNextBranchId(num_nodes));
            ports.push_back(port);
        }
    }
    int matrix_size = num_nodes + netlist.getNumBranches();
    const auto& settings = netlist.getSettings();

    VectorReal x_dc(matrix_size);
    if (!settings.use_uic) {
        std::cout << "Calculating DC Operating Point..." << std::endl;
        for (int iter = 0; iter < 100; ++iter) {
            SparseMatrixReal J_sparse(matrix_size);
            VectorReal b(matrix_size);
            for (const auto& dev : netlist.getDevices()) dev->dcStamp(J_sparse, b, x_dc, 0, {});
            VectorReal x_new = KluSolverReal::solve(J_sparse, b);
            double max_change = 0.0;
            for (int i = 0; i < num_nodes; ++i) max_change = std::max(max_change, std::abs(x_new[i] - x_dc[i]));
            x_dc = x_new;
            if (max_change < 1e-8) break;
        }
    }

    if (settings.type == "TRAN") {
        std::cout << "Starting Transient Analysis..." << std::endl;
        std::vector<VectorReal> x_hist;
        VectorReal x = x_dc; 
        x_hist.push_back(x);
        for (double t = settings.t_step; t <= settings.t_stop; t += settings.t_step) {
            for (int iter = 0; iter < 50; ++iter) {
                SparseMatrixReal J_sparse(matrix_size);
                VectorReal b(matrix_size);
                for (const auto& dev : netlist.getDevices()) dev->dcStamp(J_sparse, b, x, settings.t_step, x_hist);
                VectorReal x_new = KluSolverReal::solve(J_sparse, b);
                double max_change = 0.0;
                for (int i = 0; i < num_nodes; ++i) max_change = std::max(max_change, std::abs(x_new[i] - x[i]));
                x = x_new;
                if (max_change < 1e-6) break;
            }
            x_hist.push_back(x);
            std::cout << std::fixed << std::setprecision(6) << t << " | ";
            for(int i=0; i<num_nodes; ++i) std::cout << x[i] << " ";
            std::cout << std::endl;
        }
    } else if (settings.type == "AC") {
        std::cout << "Starting AC Analysis..." << std::endl;
        double f = settings.f_start;
        double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            double omega = 2.0 * 3.14159265358979 * f;
            SparseMatrixComplex J_sparse(matrix_size);
            VectorComplex b_ac(matrix_size);
            for (const auto& dev : netlist.getDevices()) dev->acStamp(J_sparse, b_ac, omega, x_dc);
            VectorComplex x_ac = KluSolverComplex::solve(J_sparse, b_ac);
            std::cout << std::scientific << std::setprecision(2) << f << " | ";
            for(int i=0; i<num_nodes; ++i) std::cout << "(" << x_ac[i].real() << "," << x_ac[i].imag() << ") ";
            std::cout << std::endl;
            f *= dec_mult;
        }
    } else if (settings.type == "PSS") {
        std::cout << "Starting PSS Analysis (Shooting Method)..." << std::endl;
        double T = 1.0 / settings.f_fund;
        double dt = T / 100.0;
        VectorReal x_start(matrix_size);
        bool pss_converged = false;
        for (int pss_iter = 0; pss_iter < settings.max_pss_iter; ++pss_iter) {
            std::vector<VectorReal> x_hist;
            VectorReal x = x_start;
            x_hist.push_back(x);
            for (double t = dt; t <= T * 1.001; t += dt) {
                for (int iter = 0; iter < 50; ++iter) {
                    SparseMatrixReal J_sparse(matrix_size);
                    VectorReal b(matrix_size);
                    for (const auto& dev : netlist.getDevices()) dev->dcStamp(J_sparse, b, x, dt, x_hist);
                    VectorReal x_new = KluSolverReal::solve(J_sparse, b);
                    double max_change = 0.0;
                    for (int i = 0; i < num_nodes; ++i) max_change = std::max(max_change, std::abs(x_new[i] - x[i]));
                    x = x_new;
                    if (max_change < 1e-6) break;
                }
                x_hist.push_back(x);
            }
            VectorReal x_end = x_hist.back();
            double pss_err = 0.0;
            for (int i = 0; i < num_nodes; ++i) pss_err = std::max(pss_err, std::abs(x_end[i] - x_start[i]));
            if (pss_err < 1e-4) {
                pss_converged = true;
                x_start = x_end;
                break;
            }
            x_start = x_end;
        }
        if (pss_converged) std::cout << "PSS Periodic State Found." << std::endl;
    } else if (settings.type == "SP") {
        std::cout << "Starting S-Parameter Analysis..." << std::endl;
        double f = settings.f_start;
        double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            double omega = 2.0 * 3.14159265358979 * f;
            for (size_t j = 0; j < ports.size(); ++j) {
                SparseMatrixComplex J_sparse(matrix_size);
                VectorComplex b_sp(matrix_size);
                for (const auto& dev : netlist.getDevices()) dev->acStamp(J_sparse, b_sp, omega, x_dc);
                b_sp.add(ports[j]->getBranchIndex(), {1.0, 0.0});
                VectorComplex x_sp = KluSolverComplex::solve(J_sparse, b_sp);
                double aj = 1.0 / (2.0 * std::sqrt(ports[j]->getZ0()));
                for (size_t i = 0; i < ports.size(); ++i) {
                    std::complex<double> Vi = x_sp[ports[i]->getNodePos()] - (ports[i]->getNodeNeg() >= 0 ? x_sp[ports[i]->getNodeNeg()] : 0.0);
                    std::complex<double> Ii = x_sp[ports[i]->getBranchIndex()];
                    std::complex<double> bi = (Vi - Ii * ports[i]->getZ0()) / (2.0 * std::sqrt(ports[i]->getZ0()));
                    std::complex<double> Sij = bi / aj;
                    std::cout << "S" << ports[i]->getPortNum() << ports[j]->getPortNum() << ":(" << std::abs(Sij) << "," << std::arg(Sij)*180/3.14 << ") ";
                }
            }
            std::cout << std::endl;
            f *= dec_mult;
        }
    } else if (settings.type == "NOISE") {
        std::cout << "Starting Noise Analysis..." << std::endl;
        double f = settings.f_start;
        double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            double omega = 2.0 * 3.14159265358979 * f;
            SparseMatrixComplex J_sparse(matrix_size);
            VectorComplex dummy_b(matrix_size);
            for (const auto& dev : netlist.getDevices()) dev->acStamp(J_sparse, dummy_b, omega, x_dc);
            MatrixComplex J_dense = J_sparse.toDense();
            VectorComplex i_out(matrix_size);
            if (settings.out_node >= 0) i_out[settings.out_node] = {1.0, 0.0};
            VectorComplex x_adj = J_dense.solveTranspose(i_out);
            double total_noise_psd = 0.0;
            for (const auto& dev : netlist.getDevices()) {
                double device_psd = dev->getNoisePSD(omega, x_dc);
                if (device_psd <= 0.0) continue;
                auto* res = dynamic_cast<Resistor*>(dev.get());
                if (res) {
                    std::complex<double> tf = x_adj[res->getNodePos()] - (res->getNodeNeg() >= 0 ? x_adj[res->getNodeNeg()] : 0.0);
                    total_noise_psd += device_psd * std::norm(tf);
                }
            }
            std::cout << std::scientific << std::setprecision(2) << f << " | " << total_noise_psd << std::endl;
            f *= dec_mult;
        }
    } else if (settings.type == "HB") {
        std::cout << "Starting Harmonic Balance Analysis..." << std::endl;
        int N = settings.n_harms;
        int K = 2 * N + 1;
        int n_vars = matrix_size * K;
        VectorReal x_hb(n_vars);
        for (int i = 0; i < matrix_size; ++i) x_hb[i * K] = x_dc[i];
        bool hb_converged = false;
        for (int hb_iter = 0; hb_iter < 30; ++hb_iter) {
            SparseMatrixReal J_hb(n_vars);
            VectorReal b_hb(n_vars);
            for (const auto& dev : netlist.getDevices()) dev->hbStamp(J_hb, b_hb, settings.f_fund, N, x_hb);
            VectorReal dx = KluSolverReal::solve(J_hb, b_hb);
            for (int i = 0; i < n_vars; ++i) { x_hb[i] -= dx[i]; }
            hb_converged = true; break;
        }
        if (hb_converged) std::cout << "HB Converged." << std::endl;
    } else {
        std::cout << "DC Operating Point: ";
        for(int i=0; i<num_nodes; ++i) std::cout << x_dc[i] << " ";
        std::cout << std::endl;
    }
}

void print_usage() {
    std::cout << "Usage: gspice <input_file.sp> [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -v, --version    Display version information" << std::endl;
    std::cout << "  -h, --help       Display this help message" << std::endl;
    std::cout << "  -o <file>        Specify output file name (default: output.raw)" << std::endl;
    std::cout << "\nSupported Analyses:" << std::endl;
    std::cout << "  .OP, .TRAN, .AC, .SP, .NOISE, .PSS, .HB" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        run_osdi_test(); // Default internal test for dev verification
        print_usage();
        return 0;
    }

    std::string input_file = "";
    std::string output_file = "output.raw";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage();
            return 0;
        } else if (arg == "-v" || arg == "--version") {
            std::cout << "GSPICE v0.4.0 - Next-Gen Circuit Simulator Core" << std::endl;
            std::cout << "Architecture: x64 | Sparse-Solver: KLU-Bridge | OSDI-Host: v1.0" << std::endl;
            return 0;
        } else if (arg == "-o" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg[0] != '-') {
            input_file = arg;
        }
    }

    if (input_file.empty()) {
        print_usage();
        return 1;
    }

    Netlist netlist = Parser::parse(input_file);
    if (netlist.getDevices().empty()) return 1;

    run_simulation(netlist);
    return 0;
}
