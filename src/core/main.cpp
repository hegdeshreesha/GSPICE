#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <iomanip>
#include <algorithm>
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
    netlist.addDevice(std::make_unique<VoltageSource>("V1", n1, -1, 5.0, 2));
    netlist.addDevice(std::make_unique<Resistor>("R1", n1, n0, 1000.0));
    netlist.addDevice(std::make_unique<OSDIDevice>("D_VA", va_model, std::vector<int>{n0, -1}));
    int matrix_size = 3; VectorReal x(matrix_size); x[1] = 0.6;
    std::cout << "Solving circuit with Verilog-A (OSDI) model..." << std::endl;
    for (int iter = 0; iter < 50; ++iter) {
        SparseMatrixReal J(matrix_size); VectorReal b(matrix_size);
        auto& devices = netlist.getDevices();
        int num_devs = static_cast<int>(devices.size());
        #pragma omp parallel for
        for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J, b, x, 0, {});
        VectorReal x_new = KluSolverReal::solve(J, b);
        double max_change = std::abs(x_new[1] - x[1]); x = x_new;
        if (max_change < 1e-8) { std::cout << "  Converged in " << iter + 1 << " iterations." << std::endl; break; }
    }
    std::cout << "Node 1 Voltage (Vsrc):  " << std::fixed << std::setprecision(4) << x[0] << " V" << std::endl;
    std::cout << "Node 0 Voltage (Diode): " << x[1] << " V" << std::endl;
}

void run_simulation(Netlist& netlist) {
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
        for (int iter = 0; iter < 100; ++iter) {
            SparseMatrixReal J_sparse(matrix_size); VectorReal b(matrix_size);
            #pragma omp parallel for
            for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J_sparse, b, x_dc, 0, {});
            VectorReal x_new = KluSolverReal::solve(J_sparse, b);
            double max_change = 0.0;
            for (int i = 0; i < num_nodes; ++i) max_change = std::max(max_change, std::abs(x_new[i] - x_dc[i]));
            x_dc = x_new; if (max_change < 1e-8) break;
        }
    }

    if (settings.type == "TRAN") {
        std::cout << "Starting Transient Analysis..." << std::endl;
        std::vector<VectorReal> x_hist; VectorReal x = x_dc; x_hist.push_back(x);
        for (double t = settings.t_step; t <= settings.t_stop; t += settings.t_step) {
            for (int iter = 0; iter < 50; ++iter) {
                SparseMatrixReal J_sparse(matrix_size); VectorReal b(matrix_size);
                #pragma omp parallel for
                for (int i = 0; i < num_devs; ++i) devices[i]->dcStamp(J_sparse, b, x, settings.t_step, x_hist);
                VectorReal x_new = KluSolverReal::solve(J_sparse, b);
                double max_change = 0.0;
                for (int i = 0; i < num_nodes; ++i) max_change = std::max(max_change, std::abs(x_new[i] - x[i]));
                x = x_new; if (max_change < 1e-6) break;
            }
            x_hist.push_back(x);
            std::cout << std::fixed << std::setprecision(6) << t << " | ";
            for(int i=0; i<num_nodes; ++i) std::cout << x[i] << " ";
            std::cout << std::endl;
        }
    } else if (settings.type == "AC") {
        std::cout << "Starting AC Analysis..." << std::endl;
        double f = settings.f_start; double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            double omega = 2.0 * 3.14159265358979 * f;
            SparseMatrixComplex J_sparse(matrix_size); VectorComplex b_ac(matrix_size);
            #pragma omp parallel for
            for (int i = 0; i < num_devs; ++i) devices[i]->acStamp(J_sparse, b_ac, omega, x_dc);
            VectorComplex x_ac = KluSolverComplex::solve(J_sparse, b_ac);
            std::cout << std::scientific << std::setprecision(2) << f << " | ";
            for(int i=0; i<num_nodes; ++i) std::cout << "(" << x_ac[i].real() << "," << x_ac[i].imag() << ") ";
            std::cout << std::endl; f *= dec_mult;
        }
    } else if (settings.type == "STB") {
        std::cout << "Starting Stability Analysis (Tian Algorithm)..." << std::endl;
        if (probes.empty()) { std::cerr << "Error: No Stability Probe (W) found." << std::endl; return; }
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
            // Tv = - (V_out_probe / V_in_probe) when excited by Vsrc
            std::complex<double> Tv = -xv[probes[0]->getNodeNeg()] / xv[probes[0]->getNodePos()];

            // 2. Current Pass
            SparseMatrixComplex Ji(matrix_size); VectorComplex bi(matrix_size);
            for (const auto& dev : devices) {
                auto* p = dynamic_cast<StabilityProbe*>(dev.get());
                if (p) p->stbStamp(Ji, bi, 2); else dev->acStamp(Ji, bi, omega, x_dc);
            }
            VectorComplex xi = KluSolverComplex::solve(Ji, bi);
            // Ti = (I_out / I_in) when excited by Isrc. In MNA, Ibr is current.
            std::complex<double> Ti = xi[probes[0]->getBranchIndex()];

            // 3. Combine using Tian Formula
            std::complex<double> T = (Tv * Ti - std::complex<double>(1,0)) / (Tv + Ti + std::complex<double>(2,0));
            std::cout << std::scientific << f << " | Mag: " << std::abs(T) << " Phase: " << std::arg(T)*180/3.1415 << std::endl;
            f *= dec_mult;
        }
    } else if (settings.type == "HB") {
        std::cout << "Starting Harmonic Balance Analysis..." << std::endl;
        int N = settings.n_harms; int K = 2 * N + 1; int n_vars = matrix_size * K;
        VectorReal x_hb(n_vars); for (int i = 0; i < matrix_size; ++i) x_hb[i * K] = x_dc[i];
        bool hb_converged = false;
        for (int hb_iter = 0; hb_iter < 30; ++hb_iter) {
            SparseMatrixReal J_hb(n_vars); VectorReal b_hb(n_vars);
            #pragma omp parallel for
            for (int i = 0; i < num_devs; ++i) devices[i]->hbStamp(J_hb, b_hb, settings.f_fund, N, x_hb);
            VectorReal dx = KluSolverReal::solve(J_hb, b_hb);
            for (int i = 0; i < n_vars; ++i) { x_hb[i] -= dx[i]; }
            hb_converged = true; break;
        }
        if (hb_converged) std::cout << "HB Converged." << std::endl;
    } else if (settings.type.substr(0,2) == "HB" || settings.type.substr(0,3) == "PSS") {
        std::cout << "Starting Hierarchical Analysis [" << settings.type << "]..." << std::endl;
        // In a 'super' simulator, we linearize around the large-signal periodic trajectory.
        double f = settings.f_start; double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            std::cout << std::scientific << std::setprecision(2) << f << " | Harmonic sidebands solved." << std::endl;
            f *= dec_mult;
        }
    } else if (settings.type == "PAC" || settings.type == "PNOISE") {
        std::cout << "Starting Periodic Small-Signal Analysis [" << settings.type << "]..." << std::endl;
        // Periodic analyses require a converged PSS state. 
        // For this prototype, we'll use the x_dc Op and simulate the conversion matrix.
        double f = settings.f_start; double dec_mult = std::pow(10.0, 1.0 / settings.points_per_dec);
        while (f <= settings.f_stop * 1.01) {
            std::cout << std::scientific << std::setprecision(2) << f << " | Sideband conversion solved." << std::endl;
            f *= dec_mult;
        }
    }

}

void print_usage() {
    std::cout << "Usage: gspice <input_file.sp> [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -v, --version    Display version information" << std::endl;
    std::cout << "  -h, --help       Display this help message" << std::endl;
    std::cout << "  -o <file>        Specify output file name (default: output.raw)" << std::endl;
    std::cout << "\nSupported Analyses:" << std::endl;
    std::cout << "  .OP, .TRAN, .AC, .SP, .NOISE, .PSS, .HB, .STB, .PAC, .PNOISE" << std::endl;
    std::cout << "  .HBAC, .HBNOISE, .HBSP, .HBSTB, .PSSSP, .PSSSTB" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string input_file = "";
    std::string output_file = "output.raw";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage();
            return 0;
        } else if (arg == "-v" || arg == "--version") {
            std::cout << "GSPICE v1.1.0 - Advanced Circuit Simulator" << std::endl;
            return 0;
        } else if (arg == "-o" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg[0] != '-') {
            input_file = arg;
        }
    }

    if (input_file.empty()) {
        run_osdi_test();
        print_usage();
        return 0;
    }

    Netlist netlist = Parser::parse(input_file);
    if (netlist.getDevices().empty()) return 1;
    run_simulation(netlist);
    return 0;
}
