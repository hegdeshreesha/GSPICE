#ifndef GSPICE_NETLIST_HPP
#define GSPICE_NETLIST_HPP

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <map>
#include "device.hpp"

namespace gspice {

struct SimulationSettings {
    std::string type = "OP"; // Default to DC Operating Point
    double t_stop = 0.0;
    double t_step = 0.0;
    double f_start = 0.0;
    double f_stop = 0.0;
    int points_per_dec = 0;
    bool use_uic = false;

    // PSS / HB Parameters
    std::vector<double> f_fund; // List of fundamental frequencies (e.g., f1, f2, f3, f4)
    int n_harms = 0;             // Number of harmonics per tone
    int max_pss_iter = 10;       // Max shooting iterations

    // Noise Parameters
    int out_node = -1;

    // Hierarchical / Periodic Small-Signal Flags
    bool is_periodic = false; // Set if PAC, HBAC, etc.
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
    
    const std::vector<std::unique_ptr<Device>>& getDevices() const {
        return devices_;
    }

    void setSettings(const SimulationSettings& settings) {
        settings_ = settings;
    }

    const SimulationSettings& getSettings() const {
        return settings_;
    }

    void addOsdiModel(const std::string& name, const std::string& path) {
        model_lib_map_[name] = path;
    }

    void addWarning(const std::string& message) {
        warnings_.push_back(message);
    }

    void addError(const std::string& message) {
        errors_.push_back(message);
    }

    const std::vector<std::string>& getWarnings() const {
        return warnings_;
    }

    const std::vector<std::string>& getErrors() const {
        return errors_;
    }

private:
    std::vector<std::unique_ptr<Device>> devices_;
    std::unordered_map<std::string, int> node_map_;
    std::map<int, std::string> node_names_;
    std::unordered_map<std::string, std::string> model_lib_map_;
    std::vector<std::string> warnings_;
    std::vector<std::string> errors_;
    int next_node_id_ = 0;
    int next_branch_id_ = 0; // We'll need this for voltage sources
    SimulationSettings settings_;

public:
    // Helper to get branch indices for voltage sources
    int getNextBranchId(int total_nodes) {
        return total_nodes + next_branch_id_++;
    }
    int getNumBranches() const { return next_branch_id_; }
};

} // namespace gspice

#endif // GSPICE_NETLIST_HPP
