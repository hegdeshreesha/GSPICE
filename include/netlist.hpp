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
    std::string kind = "V";
    std::string device_name;
    int node_pos = -1;
    int node_neg = -1;
    bool has_at = false;
    bool has_from = false;
    bool has_to = false;
    bool has_when_value = false;
    double at = 0.0;
    double from = 0.0;
    double to = 0.0;
    double when_value = 0.0;
    std::string crossing = "ANY";
    int crossing_count = 1;
    bool has_target = false;
    std::string target_kind = "V";
    std::string target_device_name;
    int target_node_pos = -1;
    int target_node_neg = -1;
    bool target_has_when_value = false;
    double target_when_value = 0.0;
    std::string target_crossing = "ANY";
    int target_crossing_count = 1;
};

struct FourSpec {
    double frequency = 0.0;
    int harmonics = 9;
    int node_pos = -1;
    int node_neg = -1;
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
    bool tran_max_step_auto = false;
    double t_min_step = 0.0;
    bool tran_adaptive = true;
    bool tran_predictor = true;
    std::string tran_method = "AUTO";
    std::string tran_lte_mode = "PREDICTOR";
    int tran_lte_audit_interval = 0;
    bool tran_order_adaptive = true;
    bool tran_trap_ringing = true;
    double tran_lte_reltol = 5e-3;
    double tran_lte_abstol = 1e-6;
    double tran_trtol = 1.0;
    double chgtol = 1e-14;
    double cshunt = -1.0;  // <0 means AUTO; 0 explicitly disables the shunt floor.
    int tran_max_order = 2;
    bool save_adaptive_steps = false;
    bool transient_noise = false;
    unsigned int transient_noise_seed = 1;
    double transient_noise_scale = 0.0;
    double transient_noise_fmax = 0.0;
    std::string transient_noise_mode = "WHITE";
    double f_start = 0.0;
    double f_stop = 0.0;
    int points_per_dec = 0;
    std::string f_sweep_type = "DEC";
    std::vector<double> f_values;
    bool use_uic = false;
    double temperature_c = 27.0;

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
    bool tran_verbose_debug = false;
    bool fastspice = false;
    bool multirate = false;
    bool parallel_solve = false;
    bool ticer = false;
    double ticer_fmax = 1e9;
    int num_threads = 0;

    // PSS / HB Parameters
    std::vector<double> f_fund; // List of fundamental frequencies (e.g., f1, f2, f3, f4)
    int n_harms = 0;             // Number of harmonics per tone
    int max_pss_iter = 10;       // Max shooting iterations
    bool pss_requested = false;
    double pss_tstab = 0.0;
    int pss_tstab_periods = 0;
    double pss_residual_goal = 1.0;
    bool hb_native_required = false;

    // Noise Parameters
    int out_node = -1;
    std::string pnoise_input_source;
    bool pnoise_phase_noise = false;
    bool pnoise_jitter = false;
    double pnoise_carrier = 0.0;

    // Small-signal transfer output for EXTERNAL_ORACLE-style ACXF/DCXF aliases.
    int xf_out_pos = -1;
    int xf_out_neg = -1;

    // Measurements
    std::vector<MeasureSpec> measures;
    std::vector<FourSpec> fours;
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
    CompactModelInfo compact;
};

struct GsdiArtifactInfo {
    std::string model_type;
    std::string path;
    int terminal_count = 0;
    std::vector<std::string> parameter_names;
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

    // Allocates a hidden model-internal node (e.g. series-resistance node of a
    // compact model). The column participates in the MNA matrix like any other
    // node; the generated name is namespaced so it cannot collide with user
    // node names and stays probeable by name for debugging.
    int createInternalNode(const std::string& device_name, std::size_t index) {
        return getOrCreateNode("%internal." + device_name + "." + std::to_string(index));
    }

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
        ModelCard cached = model;
        cached.compact = CompactModelRegistry::instance().classify(model.type);
        model_cards_[normalizeKey(model.name)] = cached;
    }

    const ModelCard* findModelCard(const std::string& name) const {
        auto it = model_cards_.find(normalizeKey(name));
        if (it != model_cards_.end()) return &it->second;
        return nullptr;
    }

    void addGsdiArtifact(const GsdiArtifactInfo& artifact) {
        gsdi_artifacts_[normalizeKey(artifact.model_type)] = artifact;
    }

    const GsdiArtifactInfo* findGsdiArtifact(const std::string& modelType) const {
        auto it = gsdi_artifacts_.find(normalizeKey(modelType));
        if (it != gsdi_artifacts_.end()) return &it->second;
        return nullptr;
    }

    void addWarning(const std::string& message) {
        if (std::find(warnings_.begin(), warnings_.end(), message) != warnings_.end()) return;
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
    std::unordered_map<std::string, GsdiArtifactInfo> gsdi_artifacts_;
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
