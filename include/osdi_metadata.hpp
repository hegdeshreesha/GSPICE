#ifndef GSPICE_OSDI_METADATA_HPP
#define GSPICE_OSDI_METADATA_HPP

#include "osdi.h"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace gspice {

enum class OsdiParameterKind {
    Model,
    Instance,
    Opvar
};

enum class OsdiValueType {
    Real,
    Integer,
    String
};

struct OsdiParameterInfo {
    uint32_t id = 0;
    std::string canonical_name;
    std::vector<std::string> aliases;
    OsdiParameterKind kind = OsdiParameterKind::Model;
    OsdiValueType type = OsdiValueType::Real;
    uint32_t length = 1;
    std::string units;
    std::string description;
};

struct OsdiNodeInfo {
    uint32_t index = 0;
    std::string name;
    std::string units;
    std::string residual_units;
    bool is_flow = false;
    bool has_resistive_residual = false;
    bool has_reactive_residual = false;
};

struct OsdiNoiseInfo {
    uint32_t index = 0;
    std::string name;
    OsdiNodePair nodes{UINT32_MAX, UINT32_MAX};
    uint32_t type = UINT32_MAX;
};

class OsdiDescriptorMetadata {
public:
    explicit OsdiDescriptorMetadata(const OsdiDescriptor& desc) {
        descriptor_name = desc.name ? desc.name : (desc.model_name ? desc.model_name : "<unnamed>");
        has_bound_step = desc.bound_step_offset != UINT32_MAX;
        has_collapsible_nodes = desc.num_collapsible != 0 && desc.collapsible != nullptr;
        has_noise = desc.num_noise_src != 0 && desc.noise_sources != nullptr;
        has_opvars = desc.num_opvars != 0;
        has_limiting_rhs = desc.load_limit_rhs_resist != nullptr || desc.load_limit_rhs_react != nullptr;
        has_spice_rhs = desc.load_spice_rhs_dc != nullptr || desc.load_spice_rhs_tran != nullptr;
        has_reactive = desc.load_residual_react != nullptr ||
            desc.load_jacobian_react != nullptr ||
            desc.write_jacobian_array_react != nullptr ||
            desc.load_jacobian_tran != nullptr;
        uses_abstime = (desc.module_flags & 1u) != 0;

        buildNodes(desc);
        buildParameters(desc);
        buildJacobianLists(desc);
        buildNoise(desc);
    }

    std::string descriptor_name;
    std::vector<OsdiParameterInfo> parameters;
    std::vector<uint32_t> model_parameter_ids;
    std::vector<uint32_t> instance_parameter_ids;
    std::vector<uint32_t> opvar_ids;
    std::vector<OsdiNodeInfo> nodes;
    std::vector<uint32_t> nonzero_resistive_residual_nodes;
    std::vector<uint32_t> nonzero_reactive_residual_nodes;
    std::vector<uint32_t> nonzero_resistive_jacobian_entries;
    std::vector<uint32_t> nonzero_reactive_jacobian_entries;
    std::vector<OsdiNoiseInfo> noise_sources;
    std::unordered_map<std::string, uint32_t> node_name_to_index;
    std::unordered_map<std::string, uint32_t> parameter_name_to_id;
    std::unordered_map<std::string, std::vector<uint32_t>> noise_name_to_indices;
    bool has_bound_step = false;
    bool uses_abstime = false;
    bool has_collapsible_nodes = false;
    bool has_noise = false;
    bool has_opvars = false;
    bool has_limiting_rhs = false;
    bool has_spice_rhs = false;
    bool has_reactive = false;

    const OsdiParameterInfo* findParameter(const std::string& name) const {
        const auto it = parameter_name_to_id.find(normalize(name));
        if (it == parameter_name_to_id.end() || it->second >= parameters.size()) return nullptr;
        return &parameters[it->second];
    }

    const OsdiNodeInfo* findNode(const std::string& name) const {
        const auto it = node_name_to_index.find(normalize(name));
        if (it == node_name_to_index.end() || it->second >= nodes.size()) return nullptr;
        return &nodes[it->second];
    }

    std::string summary() const {
        std::ostringstream out;
        out << descriptor_name
            << "{nodes=" << nodes.size()
            << ",terms=" << terminal_count_
            << ",params=" << parameters.size()
            << ",inst_params=" << instance_parameter_ids.size()
            << ",opvars=" << opvar_ids.size()
            << ",jac=" << (nonzero_resistive_jacobian_entries.size() + nonzero_reactive_jacobian_entries.size());
        if (has_reactive) out << ",reactive";
        if (has_bound_step) out << ",bound_step";
        if (has_limiting_rhs) out << ",limiting_rhs";
        if (has_spice_rhs) out << ",spice_rhs";
        if (has_noise) out << ",noise=" << noise_sources.size();
        if (has_collapsible_nodes) out << ",collapsible";
        out << "}";
        return out.str();
    }

private:
    uint32_t terminal_count_ = 0;

    static std::string normalize(std::string value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return value;
    }

    static std::string cstr(const char* value) {
        return value ? std::string(value) : std::string();
    }

    void buildNodes(const OsdiDescriptor& desc) {
        terminal_count_ = desc.num_terminals;
        if (!desc.nodes) return;
        nodes.reserve(desc.num_nodes);
        for (uint32_t i = 0; i < desc.num_nodes; ++i) {
            const OsdiNode& node = desc.nodes[i];
            OsdiNodeInfo info;
            info.index = i;
            info.name = node.name ? node.name : ("node" + std::to_string(i));
            info.units = cstr(node.units);
            info.residual_units = cstr(node.residual_units);
            info.is_flow = node.is_flow;
            info.has_resistive_residual = desc.load_residual_resist != nullptr && node.resist_residual_off != UINT32_MAX;
            info.has_reactive_residual = desc.load_residual_react != nullptr && node.react_residual_off != UINT32_MAX;
            if (info.has_resistive_residual) nonzero_resistive_residual_nodes.push_back(i);
            if (info.has_reactive_residual) nonzero_reactive_residual_nodes.push_back(i);
            node_name_to_index[normalize(info.name)] = i;
            nodes.push_back(std::move(info));
        }
    }

    void buildParameters(const OsdiDescriptor& desc) {
        if (!desc.param_opvar) return;
        parameters.reserve(desc.num_params);
        for (uint32_t id = 0; id < desc.num_params; ++id) {
            const OsdiParamOpvar& raw = desc.param_opvar[id];
            OsdiParameterInfo info;
            info.id = id;
            info.kind = parameterKind(raw.flags);
            info.type = valueType(raw.flags);
            info.length = raw.len == 0 ? 1 : raw.len;
            info.units = cstr(raw.units);
            info.description = cstr(raw.description);
            const uint32_t name_count = raw.name ? (1 + raw.num_alias) : 0;
            for (uint32_t n = 0; n < name_count; ++n) {
                if (!raw.name[n]) continue;
                std::string name(raw.name[n]);
                if (info.canonical_name.empty()) info.canonical_name = name;
                info.aliases.push_back(name);
                parameter_name_to_id[normalize(name)] = id;
            }
            if (info.canonical_name.empty()) {
                info.canonical_name = "param" + std::to_string(id);
                parameter_name_to_id[normalize(info.canonical_name)] = id;
            }

            if (info.kind == OsdiParameterKind::Instance) {
                instance_parameter_ids.push_back(id);
            } else if (info.kind == OsdiParameterKind::Opvar) {
                opvar_ids.push_back(id);
            } else {
                model_parameter_ids.push_back(id);
            }
            parameters.push_back(std::move(info));
        }
    }

    static OsdiParameterKind parameterKind(uint32_t flags) {
        const uint32_t kind = flags & PARA_KIND_MASK;
        if (kind == PARA_KIND_INST) return OsdiParameterKind::Instance;
        if (kind == PARA_KIND_OPVAR) return OsdiParameterKind::Opvar;
        return OsdiParameterKind::Model;
    }

    static OsdiValueType valueType(uint32_t flags) {
        const uint32_t type = flags & PARA_TY_MASK;
        if (type == PARA_TY_INT) return OsdiValueType::Integer;
        if (type == PARA_TY_STR) return OsdiValueType::String;
        return OsdiValueType::Real;
    }

    void buildJacobianLists(const OsdiDescriptor& desc) {
        if (!desc.jacobian_entries) return;
        for (uint32_t i = 0; i < desc.num_jacobian_entries; ++i) {
            const uint32_t flags = desc.jacobian_entries[i].flags;
            if ((flags & (JACOBIAN_ENTRY_RESIST | JACOBIAN_ENTRY_RESIST_CONST)) != 0) {
                nonzero_resistive_jacobian_entries.push_back(i);
            }
            if ((flags & (JACOBIAN_ENTRY_REACT | JACOBIAN_ENTRY_REACT_CONST)) != 0) {
                nonzero_reactive_jacobian_entries.push_back(i);
            }
        }
    }

    void buildNoise(const OsdiDescriptor& desc) {
        if (!desc.noise_sources) return;
        noise_sources.reserve(desc.num_noise_src);
        for (uint32_t i = 0; i < desc.num_noise_src; ++i) {
            OsdiNoiseInfo info;
            info.index = i;
            info.name = desc.noise_sources[i].name ? desc.noise_sources[i].name : ("noise" + std::to_string(i));
            info.nodes = desc.noise_sources[i].nodes;
            if (desc.noise_source_type) info.type = desc.noise_source_type[i];
            noise_name_to_indices[normalize(info.name)].push_back(i);
            noise_sources.push_back(std::move(info));
        }
    }
};

} // namespace gspice

#endif // GSPICE_OSDI_METADATA_HPP
