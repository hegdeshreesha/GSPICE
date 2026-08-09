#ifndef GSPICE_GMC_HPP
#define GSPICE_GMC_HPP

#include "device.hpp"

#include <cctype>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace gspice {

struct GmcModelDefinition {
    std::string name;
    std::string type;
    std::unordered_map<std::string, std::string> model_params;
    std::unordered_map<std::string, std::string> instance_params;
    std::vector<int> nodes;
    double temperature_c = 27.0;
    // Global columns for non-collapsible hidden unknowns (own MNA columns the
    // elaborator allocates; trailing so old aggregate initializers stay valid).
    std::vector<int> internal_nodes;
};

class GmcRegistry {
public:
    using Factory = std::function<std::unique_ptr<Device>(const GmcModelDefinition&)>;

    static GmcRegistry& instance() {
        static GmcRegistry registry;
        return registry;
    }

    void registerModel(const std::string& type, Factory factory) {
        registerModel(type, std::move(factory), 0);
    }

    // Models with hidden unknowns (non-collapsible internal nodes) declare how
    // many extra global columns the elaborator must allocate for each instance.
    void registerModel(const std::string& type, Factory factory, std::size_t internal_node_count) {
        const std::string key = normalize(type);
        factories_[key] = std::move(factory);
        internal_counts_[key] = internal_node_count;
    }

    // Number of hidden columns the elaborator must allocate for an instance of
    // this model type; 0 when the type is not registered.
    std::size_t internalNodeCount(const std::string& type) const {
        const auto it = internal_counts_.find(normalize(type));
        return it == internal_counts_.end() ? 0 : it->second;
    }

    bool hasModel(const std::string& type) const {
        return factories_.find(normalize(type)) != factories_.end();
    }

    std::unique_ptr<Device> create(const GmcModelDefinition& definition) const {
        auto it = factories_.find(normalize(definition.type));
        if (it == factories_.end()) return nullptr;
        return it->second(definition);
    }

private:
    static std::string normalize(std::string value) {
        for (char& c : value) {
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
        return value;
    }

    std::unordered_map<std::string, Factory> factories_;
    std::unordered_map<std::string, std::size_t> internal_counts_;
};

} // namespace gspice

#endif // GSPICE_GMC_HPP
