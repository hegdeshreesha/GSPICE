#ifndef GSPICE_NETLIST_HPP
#define GSPICE_NETLIST_HPP

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <utility>
#include "device.hpp"
#include "compact_model.hpp"

namespace gspice {

struct SweepSpec {
    std::string source;
    double start = 0.0;
    double stop = 0.0;
    double step = 0.0;
};

struct MeasureSpec {
    std::string analysis = "TRAN";
    std::string name;
    std::string op;
    int node_pos = -1;
    int node_neg = -1;
    bool has_at = false;
    bool has_from = false;
    bool has_to = false;
    double at = 0.0;
    double from = 0.0;
    double to = 0.0;
};

struct CornerSpec {
    std::string name;
    std::vector<std::pair<std::string, double>> source_values;
};

struct OutputSpec {
    std::string name;
    int node_pos = -1;
    int node_neg = -1;
    bool has_min = false;
    bool has_max = false;
    double min_value = 0.0;
    double max_value = 0.0;
};

struct SaveSpec {
    std::string kind = "V";
    std::string node_pos;
    std::string node_neg = "0";
};

struct InitialConditionSpec {
    int node = -1;
    double value = 0.0;
};

struct SimulationSettings {
    std::string type = "OP"; // Default to DC Operating Point
    double t_stop = 0.0;
    double t_step = 0.0;
    double t_start = 0.0;
    double t_max_step = 0.0;
    double t_min_step = 0.0;
    bool ignore_tran_tmax = false;
    bool tran_adaptive = true;
    bool tran_predictor = true;
    std::string tran_method = "AUTO";
    std::string tran_lte_mode = "PREDICTOR";
    int tran_lte_audit_interval = 0;
    int tran_lte_charge_interval = 16;
    std::string tran_lte_reference = "HISTORY";
    double tran_breakpoint_growth = 4.0;
    bool tran_order_adaptive = true;
    bool tran_trap_ringing = true;
    double tran_lte_reltol = 5e-3;
    double tran_lte_abstol = 1e-6;
    double tran_trtol = 1.0;
    double chgtol = 1e-14;
    double fluxtol = 1e-12;
    int tran_max_order = 2;
    bool save_adaptive_steps = false;
    double tran_progress_interval = 1.0;
    double f_start = 0.0;
    double f_stop = 0.0;
    int points_per_dec = 0;
    std::string frequency_sweep_type = "DEC";
    double frequency_step = 0.0;
    std::vector<double> frequency_values;
    bool use_uic = false;
    double temperature_c = 27.0;
    double nominal_temperature_c = 27.0;
    double geometry_scale = 1.0;

    // DC sweep parameters
    std::string dc_sweep_source;
    double dc_start = 0.0;
    double dc_stop = 0.0;
    double dc_step = 0.0;
    std::vector<SweepSpec> dc_sweeps;
    std::vector<SweepSpec> step_sweeps;

    // Monte Carlo source-variation parameters
    int mc_runs = 0;
    unsigned int mc_seed = 1;
    std::string mc_source;
    std::string mc_distribution = "GAUSSIAN";
    double mc_mean = 0.0;
    double mc_sigma = 0.0;
    double mc_lower = 0.0;
    double mc_upper = 0.0;
    bool mc_latin_hypercube = false;
    std::vector<CornerSpec> corners;
    std::vector<OutputSpec> output_specs;

    // Transfer function parameters
    std::string tf_input_source;
    int tf_out_pos = -1;
    int tf_out_neg = -1;

    // Sensitivity parameters
    std::string sens_source;
    int sens_out_pos = -1;
    int sens_out_neg = -1;

    // Solver / numerical controls, SPICE-like defaults
    double reltol = 1e-3;
    double vntol = 1e-6;
    double abstol = 1e-12;
    double gmin = 1e-12;
    double minr = 1e-12;
    int op_max_iter = 100;
    int tran_max_iter = 50;
    std::string solver_backend = "AUTO";
    std::string solver_ordering = "AUTO";
    bool solver_singletons = true;
    bool solver_row_scaling = true;
    int solver_refinement_steps = 1;
    bool source_stepping = true;
    bool gmin_stepping = true;
    bool line_search = true;
    bool nr_residual_check = true;
    bool nr_bypass = true;
    double nr_bypass_tolerance = 0.1;
    int nodeset_iterations = 2;
    double nodeset_conductance = 1e6;
    bool dae_audit = false;
    double dae_audit_tolerance = 2e-4;
    bool fastspice = false;
    bool multirate = false;
    bool transient_stamp_cache = true;
    bool parallel_solve = false;
    bool ticer = false;
    double ticer_fmax = 1e9;
    int num_threads = 0;

    // PSS / HB Parameters
    std::vector<double> f_fund; // List of fundamental frequencies (e.g., f1, f2, f3, f4)
    int n_harms = 0;             // Number of harmonics per tone
    int max_pss_iter = 10;       // Max shooting iterations
    bool pss_autonomous = false; // Solve period as an unknown with a phase condition
    double pss_tstab = 0.0;      // Stabilization/runup time before shooting
    int pss_tstab_periods = 0;   // Stabilization/runup periods before shooting
    bool pss_adaptive = false;   // Use adaptive transient substeps inside each PSS sample interval
    bool pss_continuation = true;
    int pss_continuation_steps = 3;
    double pss_residual_goal = 1.0;
    bool pnoise_phase_noise = false;
    bool pnoise_jitter = false;
    bool run_pnoise_after_pss = false;
    double pnoise_carrier_hz = 0.0;
    std::string pnoise_output_label;

    // Noise Parameters
    int out_node = -1;
    std::string noise_input_source;
    bool transient_noise = false;
    double transient_noise_fmin = 0.0;
    double transient_noise_fmax = 0.0;
    double transient_noise_scale = 1.0;
    unsigned int transient_noise_seed = 1;
    std::string transient_noise_mode = "ZOH";
    int transient_noise_oversample = 6;
    int transient_noise_colored_tones_per_dec = 4;

    // Measurements
    std::vector<MeasureSpec> measures;
    std::vector<InitialConditionSpec> initial_conditions;
    std::vector<InitialConditionSpec> nodesets;
    bool save_all = true;
    bool save_none = false;
    std::vector<SaveSpec> saves;

    // Hierarchical / Periodic Small-Signal Flags
    bool is_periodic = false; // Set if PAC, HBAC, etc.
};

struct ModelCard {
    std::string name;
    std::string type;
    std::unordered_map<std::string, std::string> params;
};

class Netlist {
public:
    Netlist() {
        // Node "0" is always Ground (-1)
        node_map_["0"] = -1;
    }

    void addDevice(std::unique_ptr<Device> dev) {
        devices_.push_back(std::move(dev));
    }

    int getOrCreateNode(const std::string& name) {
        if (node_map_.count(name)) {
            return node_map_[name];
        }
        int new_id = next_node_id_++;
        node_map_[name] = new_id;
        node_names_[new_id] = name;
        return new_id;
    }

    int getNumNodes() const { return next_node_id_; }

    std::string getNodeName(int index) const {
        auto it = node_names_.find(index);
        if (it != node_names_.end()) {
            return it->second;
        }
        return "node" + std::to_string(index);
    }

    int findNode(const std::string& name) const {
        auto it = node_map_.find(name);
        if (it != node_map_.end()) return it->second;
        std::string lower = normalizeKey(name);
        for (const auto& [key, value] : node_map_) {
            if (normalizeKey(key) == lower) return value;
        }
        return -2;
    }
    
    const std::vector<std::unique_ptr<Device>>& getDevices() const {
        return devices_;
    }

    void setSettings(const SimulationSettings& settings) {
        settings_ = settings;
    }

    const SimulationSettings& getSettings() const {
        return settings_;
    }

    void addModelCard(const ModelCard& model) {
        model_cards_[model.name] = model;
        model_cards_[normalizeKey(model.name)] = model;
        std::string upper = model.name;
        std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        model_cards_[upper] = model;
    }

    const ModelCard* findModelCard(const std::string& name) const {
        auto it = model_cards_.find(name);
        if (it != model_cards_.end()) return &it->second;
        it = model_cards_.find(normalizeKey(name));
        if (it != model_cards_.end()) return &it->second;
        std::string upper = name;
        std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        it = model_cards_.find(upper);
        if (it != model_cards_.end()) return &it->second;
        return nullptr;
    }

    void addWarning(const std::string& message) {
        warnings_.push_back(message);
    }

    void addError(const std::string& message) {
        errors_.push_back(message);
    }

    void addModelStatus(const std::string& message) {
        model_status_.push_back(message);
    }

    const std::vector<std::string>& getWarnings() const {
        return warnings_;
    }

    const std::vector<std::string>& getErrors() const {
        return errors_;
    }

    const std::vector<std::string>& getModelStatus() const {
        return model_status_;
    }

private:
    std::vector<std::unique_ptr<Device>> devices_;
    std::unordered_map<std::string, int> node_map_;
    std::map<int, std::string> node_names_;
    std::unordered_map<std::string, ModelCard> model_cards_;
    std::vector<std::string> model_status_;
    std::vector<std::string> warnings_;
    std::vector<std::string> errors_;
    int next_node_id_ = 0;
    int next_branch_id_ = 0; // We'll need this for voltage sources
    SimulationSettings settings_;

    static std::string normalizeKey(std::string key) {
        std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return key;
    }

public:
    // Helper to get branch indices for voltage sources
    int getNextBranchId(int total_nodes) {
        return total_nodes + next_branch_id_++;
    }
    int getNumBranches() const { return next_branch_id_; }
};

} // namespace gspice

#endif // GSPICE_NETLIST_HPP
