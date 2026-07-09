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
#include "device.hpp"
#include "osdi_loader.hpp"

namespace gspice {

struct SimulationSettings {
    std::string type = "OP"; // Default to DC Operating Point
    double t_stop = 0.0;
    double t_step = 0.0;
    double t_start = 0.0;
    double t_max_step = 0.0;
    double t_min_step = 0.0;
    bool tran_adaptive = true;
    std::string tran_method = "AUTO";
    double tran_lte_reltol = 5e-3;
    double tran_lte_abstol = 1e-6;
    double f_start = 0.0;
    double f_stop = 0.0;
    int points_per_dec = 0;
    bool use_uic = false;

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

    // PSS / HB Parameters
    std::vector<double> f_fund; // List of fundamental frequencies (e.g., f1, f2, f3, f4)
    int n_harms = 0;             // Number of harmonics per tone
    int max_pss_iter = 10;       // Max shooting iterations

    // Noise Parameters
    int out_node = -1;

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
        legacy_osdi_model_paths_[name] = path;
    }

    bool loadOsdiLibrary(const std::string& path, std::string& errorOut) {
        auto loader = std::make_unique<OSDILoader>(path);
        if (!loader->isLoaded()) {
            errorOut = loader->getError();
            return false;
        }
        for (const auto& desc : loader->getAvailableModels()) {
            const char* modelName = desc.name ? desc.name : desc.model_name;
            if (modelName) {
                osdi_descriptors_[modelName] = desc;
                osdi_descriptors_[normalizeKey(modelName)] = desc;
            }
        }
        osdi_loaders_.push_back(std::move(loader));
        return true;
    }

    void addModelCard(const ModelCard& model) {
        model_cards_[model.name] = model;
    }

    const ModelCard* findModelCard(const std::string& name) const {
        auto it = model_cards_.find(name);
        if (it != model_cards_.end()) return &it->second;
        return nullptr;
    }

    const OsdiDescriptor* findOsdiDescriptor(const std::string& modelType) const {
        auto it = osdi_descriptors_.find(modelType);
        if (it != osdi_descriptors_.end()) return &it->second;
        it = osdi_descriptors_.find(normalizeKey(modelType));
        if (it != osdi_descriptors_.end()) return &it->second;
        return nullptr;
    }

    bool tryAutoLoadOsdiForModelType(
        const std::string& modelType,
        const std::vector<std::filesystem::path>& searchRoots,
        std::string& errorOut) {
        if (findOsdiDescriptor(modelType)) return true;

        std::vector<std::string> names = osdiCandidateNames(modelType);
        std::vector<std::string> attempted;
        for (const auto& root : searchRoots) {
            if (root.empty()) continue;
            for (const auto& name : names) {
                std::filesystem::path candidate = root / name;
                attempted.push_back(candidate.string());
                if (!std::filesystem::exists(candidate)) continue;
                std::string loadError;
                if (loadOsdiLibrary(candidate.string(), loadError) && findOsdiDescriptor(modelType)) {
                    return true;
                }
                if (!loadError.empty()) {
                    errorOut = loadError;
                    return false;
                }
            }
        }
        errorOut = "No OSDI library found for model type '" + modelType + "'. Tried: ";
        for (size_t i = 0; i < attempted.size(); ++i) {
            if (i) errorOut += "; ";
            errorOut += attempted[i];
        }
        return false;
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
    std::vector<std::unique_ptr<OSDILoader>> osdi_loaders_;
    std::unordered_map<std::string, int> node_map_;
    std::map<int, std::string> node_names_;
    std::unordered_map<std::string, std::string> legacy_osdi_model_paths_;
    std::unordered_map<std::string, ModelCard> model_cards_;
    std::unordered_map<std::string, OsdiDescriptor> osdi_descriptors_;
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

    static std::vector<std::string> osdiCandidateNames(const std::string& modelType) {
        std::vector<std::string> names;
        auto addName = [&names](const std::string& stem) {
            if (stem.empty()) return;
            const std::string osdi = stem + ".osdi";
            if (std::find(names.begin(), names.end(), osdi) == names.end()) {
                names.push_back(osdi);
            }
        };
        const std::string lower = normalizeKey(modelType);
        addName(modelType);
        addName(lower);
        std::string upper = modelType;
        std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        addName(upper);
        if (lower == "psp103va") {
            addName("psp103");
            addName("PSP103");
        } else if (lower == "pspnqs103va") {
            addName("psp103_nqs");
            addName("PSP103_NQS");
        } else if (lower.size() > 2 && lower.substr(lower.size() - 2) == "va") {
            addName(lower.substr(0, lower.size() - 2));
        }
        return names;
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
