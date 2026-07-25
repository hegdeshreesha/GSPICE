#include "gspice_api.h"
#include "simulator_context.hpp"
#include "circuit_elaborator.hpp"
#include "analysis_manager.hpp"

#include <iostream>
#include <string>

struct gspice_context_impl {
    gspice::SimulatorContext context;
    std::string last_error;
};

extern "C" {

gspice_context_t gspice_create_context(void) {
    return new gspice_context_impl();
}

void gspice_destroy_context(gspice_context_t ctx) {
    if (ctx) {
        delete static_cast<gspice_context_impl*>(ctx);
    }
}

int gspice_load_netlist_string(gspice_context_t ctx, const char* netlist_text) {
    if (!ctx || !netlist_text) return -1;
    auto* impl = static_cast<gspice_context_impl*>(ctx);
    if (!gspice::CircuitElaborator::elaborate(impl->context, std::string(netlist_text))) {
        impl->last_error = impl->context.getLastError();
        return -1;
    }
    return 0;
}

int gspice_load_netlist_file(gspice_context_t ctx, const char* filepath) {
    if (!ctx || !filepath) return -1;
    auto* impl = static_cast<gspice_context_impl*>(ctx);
    std::ifstream file(filepath);
    if (!file.is_open()) {
        impl->last_error = "Could not open netlist file: " + std::string(filepath);
        return -1;
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return gspice_load_netlist_string(ctx, buffer.str().c_str());
}

int gspice_run_simulation(gspice_context_t ctx) {
    if (!ctx) return -1;
    auto* impl = static_cast<gspice_context_impl*>(ctx);
    if (impl->context.settings().type == "TRAN") {
        if (!gspice::AnalysisManager::runTransient(impl->context)) return -1;
    } else {
        if (!gspice::AnalysisManager::runOperatingPoint(impl->context)) return -1;
    }
    return 0;
}

int gspice_get_node_count(gspice_context_t ctx) {
    if (!ctx) return 0;
    auto* impl = static_cast<gspice_context_impl*>(ctx);
    return impl->context.nodeCount();
}

double gspice_get_node_voltage(gspice_context_t ctx, int node_index) {
    if (!ctx) return 0.0;
    auto* impl = static_cast<gspice_context_impl*>(ctx);
    if (node_index < 0 || node_index >= impl->context.nodeCount()) return 0.0;
    if (impl->context.solution().getSize() > node_index) {
        return impl->context.solution()[node_index];
    }
    return 0.0;
}

int gspice_get_node_index(gspice_context_t ctx, const char* node_name) {
    if (!ctx || !node_name) return -1;
    auto* impl = static_cast<gspice_context_impl*>(ctx);
    return impl->context.netlist().findNode(std::string(node_name));
}

const char* gspice_get_last_error(gspice_context_t ctx) {
    if (!ctx) return "Invalid context pointer.";
    auto* impl = static_cast<gspice_context_impl*>(ctx);
    return impl->last_error.c_str();
}

} // extern "C"
