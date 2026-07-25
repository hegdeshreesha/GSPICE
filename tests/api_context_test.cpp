#include "gspice_api.h"
#include <iostream>
#include <cassert>
#include <cmath>

int main() {
    std::cout << "Running GSPICE C API SimulatorContext Integration Test...\n";

    gspice_context_t ctx = gspice_create_context();
    assert(ctx != nullptr);

    const char* netlist = 
        "Simple RC Circuit Test\n"
        "V1 1 0 5.0\n"
        "R1 1 2 1k\n"
        "C1 2 0 1u\n"
        ".op\n"
        ".end\n";

    int ret = gspice_load_netlist_string(ctx, netlist);
    assert(ret == 0);

    int node_count = gspice_get_node_count(ctx);
    assert(node_count > 0);

    ret = gspice_run_simulation(ctx);
    assert(ret == 0);

    int node1_idx = gspice_get_node_index(ctx, "1");
    assert(node1_idx >= 0);

    gspice_destroy_context(ctx);
    std::cout << "C API Integration Test PASSED successfully.\n";
    return 0;
}
