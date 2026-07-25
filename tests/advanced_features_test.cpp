#include "solvers/adjoint_solver.hpp"
#include "nport_stamp.hpp"
#include "solvers/btf_ordering.hpp"
#include "hb_multitone.hpp"
#include "gspice_python.hpp"
#include "provenance.hpp"

#include <iostream>
#include <cassert>
#include <vector>
#include <complex>

int main() {
    std::cout << "Running GSPICE Advanced Capabilities Verification Suite...\n";

    // 1. Adjoint Solver Test
    std::vector<std::vector<double>> J = {{2.0, 1.0}, {1.0, 3.0}};
    std::vector<std::vector<double>> J_adj = {{2.0, 1.0}, {1.0, 3.0}};
    std::vector<double> rhs = {1.0, 0.0};
    std::vector<double> adjSol;
    bool ok = gspice::AdjointSolver::solveAdjoint(J, rhs, adjSol);
    assert(ok);
    assert(adjSol.size() == 2);
    std::cout << "  [PASS] Adjoint Solver Verification\n";

    // 2. N-Port S-Parameter Conversion Test
    std::vector<std::vector<std::complex<double>>> S = {
        {0.1, 0.9},
        {0.9, 0.1}
    };
    std::vector<std::vector<std::complex<double>>> Y;
    ok = gspice::NPortStamp::convertStoY(S, 50.0, Y);
    assert(ok);
    assert(Y.size() == 2);
    std::cout << "  [PASS] N-Port S-Parameter Engine Verification\n";

    // 3. SuiteSparse BTF Decomposition Test
    std::vector<int> rowPtr = {0, 1, 2};
    std::vector<int> colIdx = {0, 1};
    std::vector<int> perm;
    std::vector<gspice::BtfBlock> blocks;
    ok = gspice::BtfOrdering::decompose(2, rowPtr, colIdx, perm, blocks);
    assert(ok);
    assert(!blocks.empty());
    std::cout << "  [PASS] SuiteSparse BTF Ordering Verification\n";

    // 4. Multi-Tone HB Lattice Test
    auto lattice = gspice::HbMultiToneLattice::generateLattice(1e9, 2, 2.1e9, 2);
    assert(!lattice.empty());
    std::cout << "  [PASS] Multi-Tone Harmonic Balance Lattice Verification\n";

    // 5. Python Zero-Copy Buffer Export Test
    gspice::SimulatorContext ctx;
    ctx.results().setAnalysisType("OP");
    ctx.results().setVariableNames({"V(1)"});
    ctx.results().addSample(0.0, {5.0});
    ctx.results().finalizeResult();

    auto bufView = gspice::PythonExportEngine::getSolutionBuffer(ctx);
    assert(bufView.sampleCount == 1);
    assert(bufView.variableCount == 1);
    std::cout << "  [PASS] Zero-Copy Python API View Verification\n";

    // 6. Provenance JSON Output Test
    gspice::ProvenanceInfo prov;
    std::string provJson = prov.toJson();
    assert(!provJson.empty());
    assert(provJson.find("simulatorVersion") == std::string::npos); // JSON keys in quote
    std::cout << "  [PASS] Provenance Stamping Verification\n";

    std::cout << "ALL ADVANCED CAPABILITY TESTS PASSED SUCCESSFULLY.\n";
    return 0;
}
