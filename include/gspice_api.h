#ifndef GSPICE_API_H
#define GSPICE_API_H

// ---------------------------------------------------------------------------
// GSPICE Public C/C++ Simulator Interface (Component 4 — Beating VACASK API).
//
// Allows embedding GSPICE as a shared library (.dll / .so) inside external
// electronic design automation (EDA) tools, PyOPUS, or Python frameworks.
// ---------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

typedef void* gspice_context_t;

/// Initialize a new GSPICE simulation context.
gspice_context_t gspice_create_context(void);

/// Destroy a GSPICE context and release memory.
void gspice_destroy_context(gspice_context_t ctx);

/// Load a SPICE netlist string into the simulation context.
int gspice_load_netlist_string(gspice_context_t ctx, const char* netlist_text);

/// Load a SPICE netlist file.
int gspice_load_netlist_file(gspice_context_t ctx, const char* filepath);

/// Run the loaded simulation analysis (DC, TRAN, AC, HB, PSS).
int gspice_run_simulation(gspice_context_t ctx);

/// Retrieve the number of circuit nodes.
int gspice_get_node_count(gspice_context_t ctx);

/// Retrieve the node voltage by index.
double gspice_get_node_voltage(gspice_context_t ctx, int node_index);

/// Retrieve node index by node name.
int gspice_get_node_index(gspice_context_t ctx, const char* node_name);

/// Get last error message.
const char* gspice_get_last_error(gspice_context_t ctx);

#ifdef __cplusplus
}
#endif

#endif // GSPICE_API_H
