#ifndef GSPICE_NETLIST_HPP
#define GSPICE_NETLIST_HPP

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
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
        return new_id;
    }

    int getNumNodes() const { return next_node_id_; }
    
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

private:
    std::vector<std::unique_ptr<Device>> devices_;
    std::unordered_map<std::string, int> node_map_;
    std::unordered_map<std::string, std::string> model_lib_map_;
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
