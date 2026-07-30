#ifndef GSPICE_OSDI_DEVICE_HPP
#define GSPICE_OSDI_DEVICE_HPP

#include "device.hpp"
#include "circuit_topology.hpp"
#include "osdi.h"
#include "osdi_metadata.hpp"
#include "osdi_model_state.hpp"
#include "utils.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace gspice {

class OSDIDevice : public Device {
public:
    using ParamMap = std::unordered_map<std::string, std::string>;

    OSDIDevice(
        const std::string& name,
        const OsdiDescriptor& desc,
        const std::vector<int>& nodes,
        const ParamMap& modelParams = {},
        const ParamMap& instanceParams = {},
        double temperatureC = 27.0,
        bool useLimitingRhs = false,
        bool useTranJacobian = false,
        bool bindFullModelParams = false,
        bool useSpiceRhs = false,
        const OsdiDescriptorMetadata* loadedMetadata = nullptr,
        std::shared_ptr<OSDISharedModelState> sharedModel = nullptr,
        double simulationGmin = 1e-12,
        double simulationMinResistance = 1e-12,
        double nominalTemperatureC = 27.0,
        double geometryScale = 1.0,
        double relativeTolerance = 1e-3,
        double voltageTolerance = 1e-6,
        double absoluteTolerance = 1e-12,
        double chargeTolerance = 1e-14,
        double fluxTolerance = 1e-12)
        : Device(name),
          desc_(desc),
          metadata_(loadedMetadata ? *loadedMetadata : OsdiDescriptorMetadata(desc)),
          nodes_(nodes),
          terminal_nodes_(nodes),
          use_limiting_rhs_(useLimitingRhs),
          use_tran_jacobian_(useTranJacobian),
          bind_full_model_params_override_(bindFullModelParams),
          use_spice_rhs_(useSpiceRhs),
          simulation_gmin_(simulationGmin),
          simulation_minr_(simulationMinResistance),
          nominal_temperature_c_(nominalTemperatureC),
          geometry_scale_(geometryScale),
          relative_tolerance_(relativeTolerance),
          voltage_tolerance_(voltageTolerance),
          absolute_tolerance_(absoluteTolerance),
          charge_tolerance_(chargeTolerance),
          flux_tolerance_(fluxTolerance),
          temperature_k_(temperatureC + 273.15) {
        if (desc_.legacy_evaluate) {
            instance_data_ = desc_.legacy_create_instance ? desc_.legacy_create_instance(nullptr) : nullptr;
            return;
        }
        validateRealOsdiDescriptor();
        bind_full_model_params_ = shouldBindFullModelParams(modelParams);
        instance_storage_.resize(wordsFor(desc_.instance_size));
        instance_data_ = instance_storage_.empty() ? nullptr : instance_storage_.data();

        const char* descriptorName = desc_.name ? desc_.name : desc_.model_name;
        if (sharedModel) {
            if (!sharedModel->initialized ||
                sharedModel->model_size != desc_.model_size ||
                sharedModel->descriptor_name != std::string(descriptorName ? descriptorName : "") ||
                sharedModel->simulation_gmin != simulation_gmin_ ||
                sharedModel->simulation_minr != simulation_minr_ ||
                sharedModel->nominal_temperature_c != nominal_temperature_c_ ||
                sharedModel->geometry_scale != geometry_scale_) {
                throw std::runtime_error("OSDI shared model state does not match the descriptor");
            }
            shared_model_ = std::move(sharedModel);
        } else {
            shared_model_ = std::make_shared<OSDISharedModelState>();
            shared_model_->storage.resize(wordsFor(desc_.model_size));
            shared_model_->descriptor_name = descriptorName ? descriptorName : "";
            shared_model_->model_size = desc_.model_size;
            shared_model_->simulation_gmin = simulation_gmin_;
            shared_model_->simulation_minr = simulation_minr_;
            shared_model_->nominal_temperature_c = nominal_temperature_c_;
            shared_model_->geometry_scale = geometry_scale_;
            model_data_ = shared_model_->data();
            applyParams(model_data_, modelParams, false, shared_model_->string_params);
            setupModel();
            shared_model_->initialized = true;
        }
        model_data_ = shared_model_->data();
        applyParams(instance_data_, instanceParams, true, instance_string_params_);
        setupInstance();
        initializeStateIndexTable();
        buildNodeMapping();
        prev_state_.assign(desc_.num_states, 0.0);
        next_state_.assign(desc_.num_states, 0.0);
        prev_react_.assign(local_node_count_, 0.0);
        prev2_react_.assign(local_node_count_, 0.0);
        prev_react_derivative_.assign(local_node_count_, 0.0);
        read_only_snapshot_.resize(transientStateBytes());
    }

    // Reserve one global unknown for every internal descriptor node. Collapsed
    // aliases are anchored while inactive. Keeping the reservation stable lets
    // setup_instance change an internal collapse pattern without changing the
    // matrix dimension or invalidating every other device's unknown indices.
    int bindInternalUnknowns(const std::function<int()>& allocateUnknown) {
        if (desc_.legacy_evaluate || internal_unknowns_bound_) return 0;
        if (terminal_nodes_.size() < desc_.num_terminals) {
            throw std::logic_error("OSDI terminals are not bound");
        }
        reserved_internal_unknowns_.assign(desc_.num_nodes, -2);
        int allocated = 0;
        for (uint32_t node = desc_.num_terminals; node < desc_.num_nodes; ++node) {
            reserved_internal_unknowns_[node] = allocateUnknown();
            ++allocated;
        }
        internal_unknowns_bound_ = true;
        applyCollapseTopology();
        resetTopologyDependentState();
        return allocated;
    }

    // Re-run setup_instance after an environment or parameter update. Internal
    // collapse changes reuse the unknown reservations above. External-terminal
    // changes are published through collectNodeCollapses() for circuit-wide
    // equation re-elaboration at the surrounding analysis boundary.
    bool refreshInstanceSetupAndTopology() {
        if (desc_.legacy_evaluate) return false;
        const auto before = collapsePattern();
        setupInstance();
        initializeStateIndexTable();
        const auto after = collapsePattern();
        if (before == after) {
            invalidateDaeCache();
            return false;
        }
        if (internal_unknowns_bound_) {
            applyCollapseTopology();
            resetTopologyDependentState();
        } else {
            buildNodeMapping();
            invalidateDaeCache();
        }
        ++topology_revision_;
        return true;
    }

    std::uint64_t topologyRevision() const {
        return topology_revision_;
    }

    void collectNodeCollapses(std::vector<NodeCollapse>& collapses) const override {
        if (desc_.legacy_evaluate || desc_.num_terminals == 0 ||
            desc_.num_collapsible == 0 || !desc_.collapsible ||
            desc_.collapsed_offset == UINT32_MAX || !instance_data_) {
            return;
        }
        std::vector<uint32_t> alias(desc_.num_nodes);
        for (uint32_t node = 0; node < desc_.num_nodes; ++node) alias[node] = node;
        const auto* active = reinterpret_cast<const bool*>(
            reinterpret_cast<const unsigned char*>(instance_data_) + desc_.collapsed_offset);
        for (uint32_t i = 0; i < desc_.num_collapsible; ++i) {
            if (!active[i]) continue;
            const uint32_t from = desc_.collapsible[i].node_1;
            const uint32_t to = desc_.collapsible[i].node_2;
            if (from < alias.size()) alias[from] = to;
        }
        const auto resolve = [&](uint32_t start) {
            uint32_t current = start;
            for (uint32_t guard = 0; guard <= desc_.num_nodes; ++guard) {
                if (current == UINT32_MAX || current >= alias.size()) return UINT32_MAX;
                if (alias[current] == current) return current;
                current = alias[current];
            }
            throw std::runtime_error("OSDI collapsible-node cycle detected");
        };
        std::unordered_map<uint32_t, std::vector<int>> terminalGroups;
        for (uint32_t terminal = 0;
             terminal < desc_.num_terminals && terminal < terminal_nodes_.size();
             ++terminal) {
            terminalGroups[resolve(terminal)].push_back(terminal_nodes_[terminal]);
        }
        for (const auto& [root, terminals] : terminalGroups) {
            const bool grounded = root == UINT32_MAX ||
                std::find(terminals.begin(), terminals.end(), -1) != terminals.end();
            int representative = -1;
            for (int terminal : terminals) {
                if (terminal < 0) continue;
                if (grounded) {
                    collapses.push_back({terminal, -1});
                } else if (representative < 0) {
                    representative = terminal;
                } else if (terminal != representative) {
                    collapses.push_back({representative, terminal});
                }
            }
        }
    }

    bool reconfigureInstanceTemperature(double temperatureC) {
        const double updated = temperatureC + 273.15;
        if (!std::isfinite(updated) || updated <= 0.0) {
            throw std::invalid_argument("OSDI instance temperature must be above absolute zero");
        }
        if (updated == temperature_k_) return false;
        temperature_k_ = updated;
        return refreshInstanceSetupAndTopology();
    }

private:
    void applyCollapseTopology() {
        nodes_.resize(desc_.num_nodes, -2);
        for (uint32_t terminal = 0;
             terminal < desc_.num_terminals && terminal < terminal_nodes_.size();
             ++terminal) {
            nodes_[terminal] = terminal_nodes_[terminal];
        }
        std::vector<uint32_t> alias(desc_.num_nodes, UINT32_MAX);
        for (uint32_t i = 0; i < desc_.num_nodes; ++i) alias[i] = i;
        if (desc_.num_collapsible && desc_.collapsible && desc_.collapsed_offset != UINT32_MAX) {
            const bool* collapsed = reinterpret_cast<const bool*>(
                reinterpret_cast<const unsigned char*>(instance_data_) + desc_.collapsed_offset);
            for (uint32_t i = 0; i < desc_.num_collapsible; ++i) {
                if (!collapsed[i]) continue;
                const uint32_t from = desc_.collapsible[i].node_1;
                const uint32_t to = desc_.collapsible[i].node_2;
                if (from < alias.size()) alias[from] = to;
            }
        }
        const auto resolve = [&](uint32_t start) {
            uint32_t current = start;
            for (uint32_t guard = 0; guard <= desc_.num_nodes; ++guard) {
                if (current == UINT32_MAX || current >= alias.size()) return UINT32_MAX;
                const uint32_t next = alias[current];
                if (next == current) return current;
                current = next;
            }
            throw std::runtime_error("OSDI collapsible-node cycle detected");
        };

        std::unordered_map<uint32_t, int> rootUnknown;
        for (uint32_t terminal = 0; terminal < desc_.num_terminals; ++terminal) {
            const uint32_t root = resolve(terminal);
            if (root == UINT32_MAX) continue;
            const int terminalUnknown = nodes_[terminal] < 0 ? -1 : nodes_[terminal];
            const auto found = rootUnknown.find(root);
            // Circuit-wide external aliases are projected after assembly. Use
            // one representative locally so this OSDI instance also sees the
            // collapsed topology before the global projection is applied.
            if (found == rootUnknown.end()) rootUnknown[root] = terminalUnknown;
        }
        for (uint32_t node = 0; node < desc_.num_nodes; ++node) {
            const uint32_t root = resolve(node);
            if (root == UINT32_MAX) {
                nodes_[node] = -1;
                continue;
            }
            auto found = rootUnknown.find(root);
            if (found == rootUnknown.end()) {
                if (root < desc_.num_terminals ||
                    root >= reserved_internal_unknowns_.size() ||
                    reserved_internal_unknowns_[root] < 0) {
                    throw std::runtime_error(
                        "OSDI instance '" + name_ +
                        "' requested an invalid internal collapse root");
                }
                const int unknown = reserved_internal_unknowns_[root];
                rootUnknown[root] = unknown;
                nodes_[node] = unknown;
            } else {
                nodes_[node] = found->second;
            }
        }
        buildNodeMapping();
    }

    void resetTopologyDependentState() {
        prev_react_.assign(local_node_count_, 0.0);
        prev2_react_.assign(local_node_count_, 0.0);
        prev_react_derivative_.assign(local_node_count_, 0.0);
        prev_react_derivative_valid_ = false;
        read_only_snapshot_.resize(transientStateBytes());
        invalidateDaeCache();
    }

public:

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        if (desc_.legacy_evaluate || !desc_.eval || use_spice_rhs_) return false;

        if (canBypassDae(x, request)) {
            evaluation = cached_dae_evaluation_;
            evaluation.bypassed = true;
            return true;
        }

        struct ReadOnlyStateRestore {
            OSDIDevice& device;
            bool enabled = false;

            ReadOnlyStateRestore(OSDIDevice& target, bool shouldRestore)
                : device(target), enabled(shouldRestore) {
                if (enabled && !device.read_only_snapshot_.empty()) {
                    device.saveTransientStateBytes(
                        device.read_only_snapshot_.data(), device.read_only_snapshot_.size());
                }
            }

            ~ReadOnlyStateRestore() noexcept {
                if (!enabled || device.read_only_snapshot_.empty()) return;
                try {
                    device.restoreTransientStateBytes(
                        device.read_only_snapshot_.data(), device.read_only_snapshot_.size());
                } catch (...) {
                }
            }
        } readOnlyState(*this, request.readOnlyState);

        uint32_t flags = 0;
        if (request.staticResidual) flags |= CALC_RESIST_RESIDUAL;
        if (request.dynamicResidual) flags |= CALC_REACT_RESIDUAL;
        if (request.staticJacobian) flags |= CALC_RESIST_JACOBIAN;
        if (request.dynamicJacobian) flags |= CALC_REACT_JACOBIAN;
        const bool requestLimiting = request.enableLimiting && osdiLimitingRhsEnabled() &&
            request.analysis != DaeAnalysis::SmallSignal;
        if (requestLimiting && request.staticResidual && desc_.load_limit_rhs_resist) {
            flags |= CALC_RESIST_LIM_RHS;
        }
        if (requestLimiting && request.dynamicResidual && desc_.load_limit_rhs_react) {
            flags |= CALC_REACT_LIM_RHS;
        }
        switch (request.analysis) {
        case DaeAnalysis::OperatingPoint:
            flags |= CALC_OP | ANALYSIS_DC | ANALYSIS_STATIC;
            break;
        case DaeAnalysis::SmallSignal:
            flags |= ANALYSIS_AC | ANALYSIS_STATIC;
            break;
        case DaeAnalysis::Transient:
            flags |= ANALYSIS_TRAN;
            break;
        case DaeAnalysis::HarmonicBalance:
            return false;
        }
        if (request.nodeset) flags |= ANALYSIS_NODESET;
        if (requestLimiting) flags = limitingEvalFlags(flags);

        std::vector<double> solve = localSolveFromGlobal(x);
        scratch_state_.clear();
        double* nextState = next_state_.empty() ? nullptr : next_state_.data();
        if (request.readOnlyState && desc_.num_states != 0) {
            scratch_state_.assign(desc_.num_states, 0.0);
            nextState = scratch_state_.data();
        }
        const char* analysisName = "op";
        const char* analysisType = "static";
        if (request.analysis == DaeAnalysis::Transient) {
            analysisName = "tran";
            analysisType = "tran";
        } else if (request.analysis == DaeAnalysis::SmallSignal) {
            analysisName = "ac";
            analysisType = "ac";
        }
        OsdiSimParas paras = simulationParameters(
            request.simulationGmin,
            request.newtonIteration,
            requestLimiting && !limiting_initialized_,
            request.sourceScaleFactor,
            analysisName,
            analysisType);
        OsdiSimInfo info{
            paras,
            request.time,
            solve.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            nextState,
            flags
        };
        resetBoundStepRequest();
        const uint32_t result = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        if (!request.readOnlyState && requestLimiting) noteLimitingEval();
        if ((result & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }
        evaluation.limitingApplied = requestLimiting && (result & EVAL_RET_FLAG_LIM) != 0;
        evaluation.finishRequested = (result & EVAL_RET_FLAG_FINISH) != 0;
        evaluation.stopRequested = (result & EVAL_RET_FLAG_STOP) != 0;
        evaluation.maximumTimeStep = transientBoundStep();

        scratch_static_residual_.assign(local_node_count_, 0.0);
        scratch_dynamic_residual_.assign(local_node_count_, 0.0);
        if (request.staticResidual && desc_.load_residual_resist) {
            desc_.load_residual_resist(instance_data_, model_data_, scratch_static_residual_.data());
        }
        if (request.dynamicResidual && desc_.load_residual_react) {
            desc_.load_residual_react(instance_data_, model_data_, scratch_dynamic_residual_.data());
        }
        if (evaluation.limitingApplied && request.staticResidual && desc_.load_limit_rhs_resist) {
            scratch_limit_correction_.assign(local_node_count_, 0.0);
            desc_.load_limit_rhs_resist(instance_data_, model_data_, scratch_limit_correction_.data());
            for (size_t i = 0; i < scratch_static_residual_.size(); ++i) {
                scratch_static_residual_[i] -= scratch_limit_correction_[i];
            }
        }
        if (evaluation.limitingApplied && request.dynamicResidual && desc_.load_limit_rhs_react) {
            scratch_limit_correction_.assign(local_node_count_, 0.0);
            desc_.load_limit_rhs_react(instance_data_, model_data_, scratch_limit_correction_.data());
            for (size_t i = 0; i < scratch_dynamic_residual_.size(); ++i) {
                scratch_dynamic_residual_[i] -= scratch_limit_correction_[i];
            }
        }
        double staticGroundResidual = 0.0;
        double dynamicGroundResidual = 0.0;
        if (request.staticResidual) {
            evaluation.staticResidual.reserve(local_node_count_ + inactive_global_nodes_.size() + 1);
        }
        if (request.dynamicResidual) {
            evaluation.dynamicResidual.reserve(local_node_count_ + 1);
        }
        for (uint32_t local = 1; local < local_node_count_; ++local) {
            const int equation = globalNodeForLocal(local);
            if (equation < 0) continue;
            const bool flowEquation = localIsFlow(local);
            if (request.staticResidual) {
                evaluation.staticResidual.push_back({equation, scratch_static_residual_[local]});
                if (!flowEquation) staticGroundResidual -= scratch_static_residual_[local];
            }
            if (request.dynamicResidual) {
                const int conservationGroup = flowEquation ? -1 : 0;
                evaluation.dynamicResidual.push_back(
                    {equation, scratch_dynamic_residual_[local], conservationGroup});
                if (!flowEquation) dynamicGroundResidual -= scratch_dynamic_residual_[local];
            }
        }
        // OSDI reserves local node zero for ground, but residual loaders do not
        // materialize its discarded MNA row. Reconstruct it even when ground
        // is only used by an internal branch and is not an explicit terminal.
        // Normal MNA stamping ignores equation -1.
        if (request.staticResidual) {
            evaluation.staticResidual.push_back({-1, staticGroundResidual});
        }
        if (request.dynamicResidual) {
            evaluation.dynamicResidual.push_back({-1, dynamicGroundResidual, 0});
        }

        scratch_static_jacobian_.assign(
            std::max<uint32_t>(desc_.num_resistive_jacobian_entries, 1), 0.0);
        scratch_dynamic_jacobian_.assign(
            std::max<uint32_t>(desc_.num_reactive_jacobian_entries, 1), 0.0);
        bool haveStaticJacobian = false;
        bool haveDynamicJacobian = false;
        if (request.staticJacobian && desc_.write_jacobian_array_resist) {
            desc_.write_jacobian_array_resist(instance_data_, model_data_, scratch_static_jacobian_.data());
            haveStaticJacobian = true;
        } else if (request.staticJacobian && desc_.load_jacobian_resist &&
                   desc_.jacobian_ptr_resist_offset != UINT32_MAX) {
            populateJacobianPointers(scratch_static_jacobian_);
            desc_.load_jacobian_resist(instance_data_, model_data_);
            haveStaticJacobian = true;
        }
        if (request.dynamicJacobian && desc_.write_jacobian_array_react) {
            desc_.write_jacobian_array_react(instance_data_, model_data_, scratch_dynamic_jacobian_.data());
            haveDynamicJacobian = true;
        } else if (request.dynamicJacobian && desc_.load_jacobian_react) {
            populateReactiveJacobianPointers(scratch_dynamic_jacobian_);
            desc_.load_jacobian_react(instance_data_, model_data_, 1.0);
            haveDynamicJacobian = true;
        }

        size_t staticIndex = 0;
        size_t dynamicIndex = 0;
        if (request.staticJacobian) {
            evaluation.staticJacobian.reserve(
                desc_.num_jacobian_entries + inactive_global_nodes_.size());
        }
        if (request.dynamicJacobian) {
            evaluation.dynamicJacobian.reserve(desc_.num_jacobian_entries);
        }
        for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
            const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
            const bool hasStatic = (entry.flags & JACOBIAN_ENTRY_RESIST) != 0;
            const bool hasDynamic = (entry.flags & JACOBIAN_ENTRY_REACT) != 0;
            double staticValue = 0.0;
            double dynamicValue = 0.0;
            if (hasStatic) {
                if (haveStaticJacobian && staticIndex < scratch_static_jacobian_.size()) {
                    staticValue = scratch_static_jacobian_[staticIndex];
                }
                ++staticIndex;
            }
            if (hasDynamic) {
                if (haveDynamicJacobian && dynamicIndex < scratch_dynamic_jacobian_.size()) {
                    dynamicValue = scratch_dynamic_jacobian_[dynamicIndex];
                }
                ++dynamicIndex;
            }
            if (entry.nodes.node_2 >= node_mapping_.size() ||
                entry.nodes.node_1 >= node_mapping_.size()) continue;
            const uint32_t rowLocal = node_mapping_[entry.nodes.node_1];
            const uint32_t columnLocal = node_mapping_[entry.nodes.node_2];
            const int equation = globalNodeForLocal(rowLocal);
            const int unknown = globalNodeForLocal(columnLocal);
            if (rowLocal == 0 || equation < 0 || unknown < 0) continue;
            const bool flowEquation = localIsFlow(rowLocal);
            if (request.staticJacobian && hasStatic && staticValue != 0.0) {
                evaluation.staticJacobian.push_back({equation, unknown, staticValue});
                if (!flowEquation) evaluation.staticJacobian.push_back({-1, unknown, -staticValue});
            }
            if (request.dynamicJacobian && hasDynamic && dynamicValue != 0.0) {
                const int conservationGroup = flowEquation ? -1 : 0;
                evaluation.dynamicJacobian.push_back(
                    {equation, unknown, dynamicValue, conservationGroup});
                if (!flowEquation) {
                    evaluation.dynamicJacobian.push_back({-1, unknown, -dynamicValue, 0});
                }
            }
        }
        // Collapsed descriptor nodes may already have been allocated as hidden
        // global unknowns while the instance topology was being elaborated.
        // They no longer participate in the compact model. Anchor those unused
        // hidden unknowns independently so exact node collapse cannot leave a
        // singular global matrix or perturb the active terminal equations.
        for (int orphan : inactive_global_nodes_) {
            if (orphan < 0 || orphan >= x.getSize()) continue;
            if (request.staticResidual) evaluation.staticResidual.push_back({orphan, x[orphan]});
            if (request.staticJacobian) evaluation.staticJacobian.push_back({orphan, orphan, 1.0});
        }
        storeDaeCache(x, request, evaluation);
        return true;
    }

    bool daeAuditSafe() const override {
        return !desc_.legacy_evaluate && desc_.eval && !use_spice_rhs_;
    }

    bool daeAuditUnknown(int unknown) const override {
        return !globalIsFlow(unknown);
    }

    void annotateDaeUnknowns(std::vector<DaeUnknownKind>& kinds) const override {
        for (uint32_t i = 0; i < node_mapping_.size() && i < nodes_.size(); ++i) {
            const int global = nodes_[i];
            if (global < 0 || global >= static_cast<int>(kinds.size()) || !desc_.nodes) continue;
            const DaeUnknownKind candidate = descriptorNodeIsFlow(i)
                ? DaeUnknownKind::Flow
                : DaeUnknownKind::Potential;
            DaeUnknownKind& current = kinds[static_cast<std::size_t>(global)];
            if (current == DaeUnknownKind::Unspecified || current == candidate) {
                current = candidate;
            } else if (candidate == DaeUnknownKind::Potential) {
                // A collapsed alias shared with an electrical terminal remains
                // a potential unknown even if an auxiliary flow descriptor
                // points at the same global index.
                current = DaeUnknownKind::Potential;
            }
        }
    }

    void annotateDaeTolerances(
        std::vector<double>& unknownAbsolute,
        std::vector<double>& residualAbsolute) const override {
        for (uint32_t i = 0; i < node_mapping_.size() &&
             i < nodes_.size() && i < metadata_.nodes.size(); ++i) {
            const int global = nodes_[i];
            if (global < 0 || global >= static_cast<int>(unknownAbsolute.size()) ||
                global >= static_cast<int>(residualAbsolute.size())) {
                continue;
            }
            const OsdiNodeInfo& node = metadata_.nodes[i];
            if (std::isfinite(node.unknown_abstol) && node.unknown_abstol > 0.0) {
                unknownAbsolute[static_cast<std::size_t>(global)] = std::min(
                    unknownAbsolute[static_cast<std::size_t>(global)],
                    node.unknown_abstol);
            }
            if (std::isfinite(node.residual_abstol) && node.residual_abstol > 0.0) {
                residualAbsolute[static_cast<std::size_t>(global)] = std::min(
                    residualAbsolute[static_cast<std::size_t>(global)],
                    node.residual_abstol);
            }
        }
    }

    bool prefersDampedAutoTransient() const override { return true; }

    bool canCacheTransientStamp() const override {
        return !metadata_.has_bound_step && !metadata_.uses_abstime;
    }

    void collectTransientFootprint(std::vector<int>& unknowns) const override {
        unknowns.reserve(unknowns.size() + nodes_.size());
        for (int node : nodes_) {
            if (node >= 0 && std::find(unknowns.begin(), unknowns.end(), node) == unknowns.end()) {
                unknowns.push_back(node);
            }
        }
    }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        if (desc_.legacy_evaluate) {
            legacyDcStamp(J, b, x);
            return;
        }
        if (
            timeStep > 0.0 &&
            osdiTransientEnabled() &&
            desc_.load_residual_react &&
            (desc_.load_jacobian_react || desc_.write_jacobian_array_react) &&
            !x_hist.empty()
        ) {
            TransientContext ctx;
            ctx.timeStep = timeStep;
            ctx.currentTime = currentTime;
            ctx.method = TransientIntegrationMethod::BackwardEuler;
            ctx.a0 = 1.0 / timeStep;
            ctx.a1 = -ctx.a0;
            ctx.xHistory = &x_hist;
            transientStamp(J, b, x, ctx);
            return;
        }

        std::vector<double> solve(local_node_count_, 0.0);
        for (size_t i = 0; i < nodes_.size() && i < node_mapping_.size(); ++i) {
            const uint32_t local = node_mapping_[i];
            if (local == UINT32_MAX || local >= solve.size()) continue;
            solve[local] = nodes_[i] >= 0 ? x[nodes_[i]] : 0.0;
        }

        OsdiSimParas paras = simulationParameters(simulation_gmin_);
        OsdiSimInfo info{
            paras,
            currentTime,
            solve.data(),
            nullptr,
            nullptr,
            limitingEvalFlags(CALC_RESIST_RESIDUAL | CALC_RESIST_JACOBIAN |
                (osdiLimitingRhsEnabled() ? CALC_RESIST_LIM_RHS : 0u) |
                CALC_OP | ANALYSIS_DC | ANALYSIS_STATIC)
        };
        resetBoundStepRequest();
        uint32_t ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        noteLimitingEval();
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }

        std::vector<double> residual(local_node_count_, 0.0);
        if (desc_.load_residual_resist) {
            desc_.load_residual_resist(instance_data_, model_data_, residual.data());
        }
        if (limitingApplied(ret) && desc_.load_limit_rhs_resist) {
            std::vector<double> limited(local_node_count_, 0.0);
            desc_.load_limit_rhs_resist(instance_data_, model_data_, limited.data());
            for (size_t i = 0; i < residual.size() && i < limited.size(); ++i) {
                residual[i] -= limited[i];
            }
        }
        std::vector<double> spiceRhs(local_node_count_, 0.0);
        const bool useSpiceRhs = use_spice_rhs_ && desc_.load_spice_rhs_dc;
        if (useSpiceRhs) {
            desc_.load_spice_rhs_dc(instance_data_, model_data_, spiceRhs.data(), solve.data());
        }

        std::vector<double> jacobian(std::max<uint32_t>(desc_.num_resistive_jacobian_entries, desc_.num_jacobian_entries), 0.0);
        bool wroteArray = false;
        if (desc_.load_jacobian_resist && desc_.jacobian_ptr_resist_offset != UINT32_MAX) {
            populateJacobianPointers(jacobian);
            desc_.load_jacobian_resist(instance_data_, model_data_);
            wroteArray = true;
        }
        if (desc_.write_jacobian_array_resist) {
            desc_.write_jacobian_array_resist(instance_data_, model_data_, jacobian.data());
            wroteArray = true;
        }

        for (uint32_t localRow = 1; localRow < local_node_count_; ++localRow) {
            const int globalRow = globalNodeForLocal(localRow);
            if (globalRow < 0 || localRow >= residual.size()) continue;

            double rhs = useSpiceRhs ? spiceRhs[localRow] : -residual[localRow];
            size_t jacIndex = 0;
            for (uint32_t e = 0; wroteArray && e < desc_.num_jacobian_entries; ++e) {
                const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
                const bool isResistive = (entry.flags & JACOBIAN_ENTRY_RESIST) != 0;
                if (!isResistive) continue;
                const double g = jacIndex < jacobian.size() ? jacobian[jacIndex] : 0.0;
                ++jacIndex;

                if (entry.nodes.node_2 >= node_mapping_.size() || entry.nodes.node_1 >= node_mapping_.size()) continue;
                const uint32_t rowLocal = node_mapping_[entry.nodes.node_1];
                const uint32_t colLocal = node_mapping_[entry.nodes.node_2];
                if (rowLocal != localRow || colLocal == UINT32_MAX) continue;

                const int colGlobal = globalNodeForLocal(colLocal);
                const double vcol = colLocal < solve.size() ? solve[colLocal] : 0.0;
                if (!useSpiceRhs) {
                    rhs += g * vcol;
                }
                if (colGlobal >= 0) {
                    J.add(globalRow, colGlobal, g);
                }
            }
            b.add(globalRow, rhs);
        }
    }

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, const TransientContext& ctx) override {
        if (desc_.legacy_evaluate) {
            legacyDcStamp(J, b, x);
            return;
        }
        if (
            ctx.timeStep > 0.0 &&
            ctx.xHistory &&
            !ctx.xHistory->empty() &&
            osdiTransientEnabled() &&
            desc_.load_residual_react &&
            (desc_.load_jacobian_react || desc_.write_jacobian_array_react)
        ) {
            transientStamp(J, b, x, ctx);
            return;
        }

        static const std::vector<VectorReal> empty_history;
        dcStamp(J, b, x, 0.0, ctx.currentTime, empty_history);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        if (desc_.legacy_evaluate) {
            const int n = static_cast<int>(nodes_.size());
            std::vector<double> voltages(n, 0.0);
            std::vector<double> currents(n, 0.0);
            std::vector<double> charges(n, 0.0);
            std::vector<double> jacobian(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
            for (int i = 0; i < n; ++i) {
                voltages[i] = nodes_[i] >= 0 ? x_dc[nodes_[i]] : 0.0;
            }
            desc_.legacy_evaluate(instance_data_, voltages.data(), currents.data(), charges.data(), jacobian.data());
            for (int row = 0; row < n; ++row) {
                if (nodes_[row] < 0) continue;
                for (int col = 0; col < n; ++col) {
                    if (nodes_[col] < 0) continue;
                    J.add(nodes_[row], nodes_[col], {jacobian[static_cast<size_t>(row) * n + col], 0.0});
                }
            }
            return;
        }
        if (!desc_.eval) return;

        std::vector<double> solve = localSolveFromGlobal(x_dc);
        OsdiSimParas paras = simulationParameters(simulation_gmin_);
        OsdiSimInfo info{
            paras,
            0.0,
            solve.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            next_state_.empty() ? nullptr : next_state_.data(),
            CALC_RESIST_JACOBIAN | CALC_REACT_JACOBIAN | ANALYSIS_AC | ANALYSIS_STATIC
        };
        const uint32_t ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }

        std::vector<double> jacResist(std::max<uint32_t>(desc_.num_resistive_jacobian_entries, 1), 0.0);
        bool wroteResistArray = false;
        if (desc_.write_jacobian_array_resist) {
            desc_.write_jacobian_array_resist(instance_data_, model_data_, jacResist.data());
            wroteResistArray = true;
        } else if (desc_.load_jacobian_resist && desc_.jacobian_ptr_resist_offset != UINT32_MAX) {
            populateJacobianPointers(jacResist);
            desc_.load_jacobian_resist(instance_data_, model_data_);
            wroteResistArray = true;
        }

        std::vector<double> jacReact(std::max<uint32_t>(desc_.num_reactive_jacobian_entries, 1), 0.0);
        bool wroteReactArray = false;
        if (desc_.write_jacobian_array_react) {
            desc_.write_jacobian_array_react(instance_data_, model_data_, jacReact.data());
            wroteReactArray = true;
        } else if (desc_.load_jacobian_react) {
            populateReactiveJacobianPointers(jacReact);
            desc_.load_jacobian_react(instance_data_, model_data_, 1.0);
            wroteReactArray = true;
        }

        size_t resistIndex = 0;
        size_t reactIndex = 0;
        for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
            const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
            double g = 0.0;
            double c = 0.0;
            if ((entry.flags & JACOBIAN_ENTRY_RESIST) != 0) {
                if (wroteResistArray && resistIndex < jacResist.size()) g = jacResist[resistIndex];
                ++resistIndex;
            }
            if ((entry.flags & JACOBIAN_ENTRY_REACT) != 0) {
                if (wroteReactArray && reactIndex < jacReact.size()) c = jacReact[reactIndex];
                ++reactIndex;
            }
            if (g == 0.0 && c == 0.0) continue;
            if (entry.nodes.node_2 >= node_mapping_.size() || entry.nodes.node_1 >= node_mapping_.size()) continue;
            const uint32_t rowLocal = node_mapping_[entry.nodes.node_1];
            const uint32_t colLocal = node_mapping_[entry.nodes.node_2];
            const int row = globalNodeForLocal(rowLocal);
            const int col = globalNodeForLocal(colLocal);
            if (row < 0 || col < 0) continue;
            J.add(row, col, {g, omega * c});
        }
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc, std::vector<NoiseSource>& sources) const override {
        if (desc_.legacy_evaluate || !desc_.load_noise || desc_.num_noise_src == 0 || !desc_.noise_sources) return;
        std::vector<double> solve = localSolveFromGlobal(x_dc);
        OsdiSimParas paras = simulationParameters(simulation_gmin_);
        std::vector<double> prevState = prev_state_;
        std::vector<double> scratchState(desc_.num_states, 0.0);
        OsdiSimInfo info{
            paras,
            0.0,
            solve.data(),
            prevState.empty() ? nullptr : prevState.data(),
            scratchState.empty() ? nullptr : scratchState.data(),
            CALC_NOISE | ANALYSIS_NOISE
        };
        const uint32_t ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }
        std::vector<double> densities(desc_.num_noise_src, 0.0);
        std::vector<double> powers(desc_.num_noise_src, 0.0);
        std::vector<double> exponents(desc_.num_noise_src, 0.0);
        const double freq = std::max(omega / (2.0 * 3.14159265358979323846), 0.0);
        desc_.load_noise(instance_data_, model_data_, freq, densities.data());
        if (desc_.load_noise_params) {
            desc_.load_noise_params(
                instance_data_, model_data_, powers.data(), exponents.data());
        }
        for (uint32_t i = 0; i < desc_.num_noise_src; ++i) {
            if (densities[i] <= 0.0) continue;
            const int nodePos = globalNodeForDescriptorNode(desc_.noise_sources[i].nodes.node_1);
            const int nodeNeg = globalNodeForDescriptorNode(desc_.noise_sources[i].nodes.node_2);
            if (nodePos < 0 && nodeNeg < 0) continue;
            const char* noiseName = desc_.noise_sources[i].name ? desc_.noise_sources[i].name : "noise";
            const uint32_t noiseType = desc_.noise_source_type
                ? desc_.noise_source_type[i]
                : UINT32_MAX;
            sources.push_back({
                name_ + "." + noiseName,
                nodePos,
                nodeNeg,
                densities[i],
                noiseType,
                powers[i],
                exponents[i]});
        }
    }

    void acceptTransientStep(const VectorReal& x, double currentTime) override {
        (void)x;
        (void)currentTime;
        if (!next_state_.empty()) {
            prev_state_ = next_state_;
        }
        invalidateDaeCache();
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        if (!next_state_.empty()) {
            prev_state_ = next_state_;
        }
        if (!desc_.legacy_evaluate && desc_.load_residual_react && ctx.timeStep > 0.0) {
            accept_react_now_.assign(local_node_count_, 0.0);
            if (loadReactiveResidualAt(x, currentTime, accept_react_now_)) {
                accept_react_derivative_.assign(local_node_count_, 0.0);
                const bool useTrap =
                    ctx.method == TransientIntegrationMethod::Trapezoidal &&
                    prev_react_derivative_valid_ &&
                    prev_react_.size() == accept_react_now_.size();
                const bool useSecond =
                    ctx.hasSecondHistory &&
                    prev_react_.size() == accept_react_now_.size() &&
                    prev2_react_.size() == accept_react_now_.size();
                if (useTrap) {
                    for (size_t i = 0; i < accept_react_now_.size(); ++i) {
                        accept_react_derivative_[i] =
                            ctx.a0 * accept_react_now_[i] +
                            ctx.a1 * prev_react_[i] -
                            prev_react_derivative_[i];
                    }
                } else if (useSecond) {
                    for (size_t i = 0; i < accept_react_now_.size(); ++i) {
                        accept_react_derivative_[i] =
                            ctx.a0 * accept_react_now_[i] +
                            ctx.a1 * prev_react_[i] +
                            ctx.a2 * prev2_react_[i];
                    }
                } else {
                    const double alpha = 1.0 / ctx.timeStep;
                    for (size_t i = 0; i < accept_react_now_.size(); ++i) {
                        const double prev = prev_react_.size() == accept_react_now_.size()
                            ? prev_react_[i]
                            : accept_react_now_[i];
                        accept_react_derivative_[i] = alpha * (accept_react_now_[i] - prev);
                    }
                }
                prev2_react_ = prev_react_;
                prev_react_.swap(accept_react_now_);
                prev_react_derivative_.swap(accept_react_derivative_);
                prev_react_derivative_valid_ = true;
            }
        }
        invalidateDaeCache();
    }

    bool acceptedDaeDynamicResidual(DaeHistory& residual) const override {
        residual.clear();
        if (desc_.legacy_evaluate || !desc_.load_residual_react ||
            prev_react_.size() != local_node_count_) {
            return false;
        }
        residual.reserve(local_node_count_ + 1);
        double dynamicGroundResidual = 0.0;
        for (uint32_t local = 1; local < local_node_count_; ++local) {
            const int equation = globalNodeForLocal(local);
            if (equation < 0) continue;
            const bool flowEquation = localIsFlow(local);
            const int conservationGroup = flowEquation ? -1 : 0;
            residual.push_back({equation, prev_react_[local], conservationGroup});
            if (!flowEquation) dynamicGroundResidual -= prev_react_[local];
        }
        residual.push_back({-1, dynamicGroundResidual, 0});
        return true;
    }

    void limitTransientNewton(const VectorReal& previous, VectorReal& candidate) const override {
        if (nodes_.empty() || candidate.getSize() <= 0) return;
        const int n = static_cast<int>(nodes_.size());
        for (int i = 0; i < n; ++i) {
            const int node = nodes_[i];
            if (node < 0 || node >= candidate.getSize() || node >= previous.getSize()) continue;
            const double v_prev = previous[node];
            const double v_cand = candidate[node];
            const double delta = v_cand - v_prev;
            if (std::abs(delta) > 0.5) {
                const double sign = (delta > 0.0) ? 1.0 : -1.0;
                candidate[node] = v_prev + sign * (0.5 + 0.1 * std::log(1.0 + (std::abs(delta) - 0.5)));
            }
        }
    }

    double transientChargeError(
        const VectorReal& coarse,
        const VectorReal& fine,
        double reltol,
        double chgtol) override {
        if (!desc_.load_residual_react || local_node_count_ <= 1) return 0.0;
        std::vector<double> q_coarse(local_node_count_, 0.0);
        std::vector<double> q_fine(local_node_count_, 0.0);
        if (!const_cast<OSDIDevice*>(this)->loadReactiveResidualAt(coarse, 0.0, q_coarse) ||
            !const_cast<OSDIDevice*>(this)->loadReactiveResidualAt(fine, 0.0, q_fine)) {
            return 0.0;
        }
        double max_err = 0.0;
        for (std::size_t i = 1; i < local_node_count_; ++i) {
            const double dq = std::abs(q_fine[i] - q_coarse[i]);
            const double tol = reltol * std::abs(q_fine[i]) + chgtol;
            if (tol > 0.0) {
                max_err = std::max(max_err, dq / tol);
            }
        }
        return max_err;
    }

    std::size_t transientStateBytes() const override {
        return instance_storage_.size() * sizeof(std::max_align_t) +
            (prev_state_.size() + next_state_.size() + prev_react_.size() +
             prev2_react_.size() + prev_react_derivative_.size()) * sizeof(double) +
            2 * sizeof(std::uint8_t);
    }

    void saveTransientStateBytes(std::byte* destination, std::size_t size) const override {
        if (size != transientStateBytes()) throw std::invalid_argument("OSDI transient state size");
        std::size_t offset = 0;
        const auto write = [&](const void* source, std::size_t bytes) {
            if (bytes != 0) std::memcpy(destination + offset, source, bytes);
            offset += bytes;
        };
        write(instance_storage_.data(), instance_storage_.size() * sizeof(std::max_align_t));
        write(prev_state_.data(), prev_state_.size() * sizeof(double));
        write(next_state_.data(), next_state_.size() * sizeof(double));
        write(prev_react_.data(), prev_react_.size() * sizeof(double));
        write(prev2_react_.data(), prev2_react_.size() * sizeof(double));
        write(prev_react_derivative_.data(), prev_react_derivative_.size() * sizeof(double));
        const std::uint8_t derivativeValid = prev_react_derivative_valid_ ? 1u : 0u;
        const std::uint8_t limitingInitialized = limiting_initialized_ ? 1u : 0u;
        write(&derivativeValid, sizeof(derivativeValid));
        write(&limitingInitialized, sizeof(limitingInitialized));
    }

    void restoreTransientStateBytes(const std::byte* source, std::size_t size) override {
        if (size != transientStateBytes()) throw std::invalid_argument("OSDI transient state size");
        std::size_t offset = 0;
        const auto read = [&](void* destination, std::size_t bytes) {
            if (bytes != 0) std::memcpy(destination, source + offset, bytes);
            offset += bytes;
        };
        read(instance_storage_.data(), instance_storage_.size() * sizeof(std::max_align_t));
        read(prev_state_.data(), prev_state_.size() * sizeof(double));
        read(next_state_.data(), next_state_.size() * sizeof(double));
        read(prev_react_.data(), prev_react_.size() * sizeof(double));
        read(prev2_react_.data(), prev2_react_.size() * sizeof(double));
        read(prev_react_derivative_.data(), prev_react_derivative_.size() * sizeof(double));
        std::uint8_t derivativeValid = 0;
        std::uint8_t limitingInitialized = 0;
        read(&derivativeValid, sizeof(derivativeValid));
        read(&limitingInitialized, sizeof(limitingInitialized));
        prev_react_derivative_valid_ = derivativeValid != 0;
        limiting_initialized_ = limitingInitialized != 0;
        invalidateDaeCache();
    }

    const OsdiDescriptorMetadata& metadata() const {
        return metadata_;
    }

    std::shared_ptr<OSDISharedModelState> sharedModelState() const {
        return shared_model_;
    }

    std::vector<std::string> getOpvarNames() const {
        std::vector<std::string> names;
        names.reserve(metadata_.opvar_ids.size());
        for (uint32_t id : metadata_.opvar_ids) {
            if (id < metadata_.parameters.size()) {
                names.push_back(metadata_.parameters[id].canonical_name);
            }
        }
        return names;
    }

    bool readOpvar(const std::string& name, double& value) const {
        if (desc_.legacy_evaluate || !desc_.access || !instance_data_) return false;
        const OsdiParameterInfo* info = metadata_.findParameter(name);
        if (!info || info->kind != OsdiParameterKind::Opvar) return false;
        const uint32_t flags = info->kind == OsdiParameterKind::Instance
            ? (ACCESS_FLAG_READ | ACCESS_FLAG_INSTANCE)
            : (ACCESS_FLAG_READ | ACCESS_FLAG_INSTANCE);
        void* raw = desc_.access(instance_data_, model_data_, info->id, flags);
        if (!raw) return false;
        if (info->type == OsdiValueType::Real) {
            value = *reinterpret_cast<double*>(raw);
            return true;
        }
        if (info->type == OsdiValueType::Integer) {
            value = static_cast<double>(*reinterpret_cast<int32_t*>(raw));
            return true;
        }
        return false;
    }

    void collectOperatingPointVariables(
        const VectorReal& x,
        std::vector<OperatingPointVariable>& variables) override {
        if (desc_.legacy_evaluate || !desc_.eval || !desc_.access ||
            metadata_.opvar_ids.empty()) {
            return;
        }
        std::vector<double> solve = localSolveFromGlobal(x);
        std::vector<double> scratchState(desc_.num_states, 0.0);
        OsdiSimParas paras = simulationParameters(
            simulation_gmin_, 0, false, 1.0, "op", "static");
        OsdiSimInfo info{
            paras,
            0.0,
            solve.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            scratchState.empty() ? nullptr : scratchState.data(),
            CALC_OP | CALC_RESIST_RESIDUAL | ANALYSIS_DC | ANALYSIS_STATIC};
        const uint32_t result =
            desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        if ((result & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error(
                "OSDI device '" + name_ +
                "' aborted while evaluating operating-point variables");
        }

        for (uint32_t id : metadata_.opvar_ids) {
            if (id >= metadata_.parameters.size()) continue;
            const OsdiParameterInfo& parameter = metadata_.parameters[id];
            void* raw = desc_.access(
                instance_data_, model_data_, id,
                ACCESS_FLAG_READ | ACCESS_FLAG_INSTANCE);
            if (!raw) continue;
            OperatingPointVariable variable;
            variable.name = name_ + "." + parameter.canonical_name;
            variable.units = parameter.units;
            const std::size_t length = std::max<std::size_t>(parameter.length, 1);
            if (parameter.type == OsdiValueType::Real) {
                const auto* values = reinterpret_cast<const double*>(raw);
                variable.numericValues.assign(values, values + length);
            } else if (parameter.type == OsdiValueType::Integer) {
                const auto* values = reinterpret_cast<const int32_t*>(raw);
                variable.numericValues.reserve(length);
                for (std::size_t i = 0; i < length; ++i) {
                    variable.numericValues.push_back(
                        static_cast<double>(values[i]));
                }
            } else {
                const char* const* value =
                    reinterpret_cast<const char* const*>(raw);
                if (value && *value) variable.stringValue = *value;
            }
            variables.push_back(std::move(variable));
        }
        invalidateDaeCache();
    }

    double boundStep() const {
        if (!instance_data_ || desc_.bound_step_offset == UINT32_MAX) return 0.0;
        const auto* raw = reinterpret_cast<const double*>(
            reinterpret_cast<const unsigned char*>(instance_data_) + desc_.bound_step_offset);
        if (!raw || !std::isfinite(*raw) || *raw <= 0.0) return 0.0;
        return *raw;
    }

    double transientBoundStep() const override {
        return boundStep();
    }

private:
    OsdiDescriptor desc_{};
    OsdiDescriptorMetadata metadata_;
    std::vector<int> nodes_;
    std::vector<int> terminal_nodes_;
    std::vector<int> reserved_internal_unknowns_;
    std::shared_ptr<OSDISharedModelState> shared_model_;
    std::vector<std::max_align_t> instance_storage_;
    void* model_data_ = nullptr;
    void* instance_data_ = nullptr;
    std::vector<uint32_t> node_mapping_;
    std::vector<int> inactive_global_nodes_;
    size_t local_node_count_ = 0;
    std::vector<std::unique_ptr<char[]>> instance_string_params_;
    std::vector<double> prev_state_;
    std::vector<double> next_state_;
    std::vector<double> prev_react_;
    std::vector<double> prev2_react_;
    std::vector<double> prev_react_derivative_;
    std::vector<double> scratch_state_;
    std::vector<double> scratch_static_residual_;
    std::vector<double> scratch_dynamic_residual_;
    std::vector<double> scratch_limit_correction_;
    std::vector<double> scratch_static_jacobian_;
    std::vector<double> scratch_dynamic_jacobian_;
    std::vector<double> accept_react_now_;
    std::vector<double> accept_react_derivative_;
    std::vector<double> scratch_reactive_solve_;
    std::vector<double> scratch_reactive_state_;
    std::vector<std::byte> read_only_snapshot_;
    bool prev_react_derivative_valid_ = false;
    bool bind_full_model_params_ = false;
    bool use_limiting_rhs_ = false;
    bool use_tran_jacobian_ = false;
    bool bind_full_model_params_override_ = false;
    bool use_spice_rhs_ = false;
    bool limiting_initialized_ = false;
    bool internal_unknowns_bound_ = false;
    std::uint64_t topology_revision_ = 0;
    double simulation_gmin_ = 1e-12;
    double simulation_minr_ = 1e-12;
    double nominal_temperature_c_ = 27.0;
    double geometry_scale_ = 1.0;
    double relative_tolerance_ = 1e-3;
    double voltage_tolerance_ = 1e-6;
    double absolute_tolerance_ = 1e-12;
    double charge_tolerance_ = 1e-14;
    double flux_tolerance_ = 1e-12;
    double temperature_k_ = 300.15;
    mutable std::array<char*, 16> simparam_names_{
        const_cast<char*>("iniLim"),
        const_cast<char*>("gmin"),
        const_cast<char*>("gdev"),
        const_cast<char*>("tnom"),
        const_cast<char*>("minr"),
        const_cast<char*>("scale"),
        const_cast<char*>("iteration"),
        const_cast<char*>("simulatorVersion"),
        const_cast<char*>("simulatorSubversion"),
        const_cast<char*>("sourceScaleFactor"),
        const_cast<char*>("reltol"),
        const_cast<char*>("vntol"),
        const_cast<char*>("abstol"),
        const_cast<char*>("chgtol"),
        const_cast<char*>("fluxtol"),
        nullptr};
    mutable std::array<double, 15> simparam_values_{};
    mutable std::array<char*, 4> simparam_string_names_{
        const_cast<char*>("analysis_name"),
        const_cast<char*>("analysis_type"),
        const_cast<char*>("cwd"),
        nullptr};
    mutable std::array<char*, 4> simparam_string_values_{nullptr, nullptr, nullptr, nullptr};
    mutable std::string simparam_analysis_name_ = "setup";
    mutable std::string simparam_analysis_type_ = "static";
    mutable std::string simparam_cwd_ = std::filesystem::current_path().string();
    bool dae_cache_valid_ = false;
    DaeRequest cached_dae_request_{};
    DaeEvaluation cached_dae_evaluation_{};
    std::vector<double> cached_dae_inputs_;

    static size_t wordsFor(uint32_t bytes) {
        if (bytes == 0) return 0;
        return (static_cast<size_t>(bytes) + sizeof(std::max_align_t) - 1) / sizeof(std::max_align_t);
    }

    static void* osdiHandle() {
        static char handle[] = "gspice";
        return handle;
    }

    static bool osdiTransientEnabled() {
        const char* value = std::getenv("GSPICE_ENABLE_OSDI_TRAN");
        if (!value) return true;
        std::string text(value);
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return !(text == "0" || text == "false" || text == "no" || text == "off");
    }

    bool fullOsdiModelParamBindingEnabled() const {
        if (bind_full_model_params_override_) return true;
        const char* value = std::getenv("GSPICE_BIND_FULL_OSDI_MODEL_PARAMS");
        if (!value) return false;
        std::string text(value);
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return text == "1" || text == "true" || text == "yes" || text == "on";
    }

    bool standardOsdiTranJacobianEnabled() const {
        if (use_tran_jacobian_) return true;
        const char* value = std::getenv("GSPICE_USE_OSDI_TRAN_JACOBIAN");
        if (!value) return false;
        std::string text(value);
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return text == "1" || text == "true" || text == "yes" || text == "on";
    }

    bool osdiLimitingRhsEnabled() const {
        if (use_limiting_rhs_) return true;
        const char* value = std::getenv("GSPICE_USE_OSDI_LIMITING_RHS");
        if (!value) return false;
        std::string text(value);
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return text == "1" || text == "true" || text == "yes" || text == "on";
    }

    static std::string lower(std::string value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return value;
    }

    void validateRealOsdiDescriptor() const {
        if (!desc_.name || !desc_.setup_model || !desc_.setup_instance || !desc_.eval ||
            !desc_.load_residual_resist || !desc_.jacobian_entries || desc_.num_terminals == 0) {
            throw std::runtime_error("OSDI descriptor is incomplete");
        }
        if (nodes_.size() < desc_.num_terminals) {
            throw std::runtime_error("OSDI device has fewer circuit nodes than the model requires");
        }
    }

    void setupModel() {
        OsdiSimParas paras = simulationParameters(
            simulation_gmin_, 0, true, 1.0, "setup_model", "static");
        OsdiInitInfo res{0, 0, nullptr};
        desc_.setup_model(osdiHandle(), model_data_, &paras, &res);
        checkInitResult(res, "model");
    }

    void setupInstance() {
        if (instance_data_ && desc_.num_collapsible != 0 &&
            desc_.collapsed_offset != UINT32_MAX) {
            auto* collapsed = reinterpret_cast<bool*>(
                reinterpret_cast<unsigned char*>(instance_data_) +
                desc_.collapsed_offset);
            std::fill(
                collapsed,
                collapsed + desc_.num_collapsible,
                false);
        }
        OsdiSimParas paras = simulationParameters(
            simulation_gmin_, 0, true, 1.0, "setup_instance", "static");
        OsdiInitInfo res{0, 0, nullptr};
        desc_.setup_instance(osdiHandle(), instance_data_, model_data_, temperature_k_, desc_.num_terminals, &paras, &res);
        checkInitResult(res, "instance");
    }

    std::vector<bool> collapsePattern() const {
        std::vector<bool> pattern(desc_.num_collapsible, false);
        if (!instance_data_ || desc_.num_collapsible == 0 ||
            desc_.collapsed_offset == UINT32_MAX) {
            return pattern;
        }
        const auto* collapsed = reinterpret_cast<const bool*>(
            reinterpret_cast<const unsigned char*>(instance_data_) +
            desc_.collapsed_offset);
        for (uint32_t i = 0; i < desc_.num_collapsible; ++i) {
            pattern[i] = collapsed[i];
        }
        return pattern;
    }

    void initializeStateIndexTable() {
        if (!instance_data_ || desc_.state_idx_off == UINT32_MAX || desc_.num_states == 0) return;
        auto* stateIndices = reinterpret_cast<uint32_t*>(
            reinterpret_cast<unsigned char*>(instance_data_) + desc_.state_idx_off);
        for (uint32_t i = 0; i < desc_.num_states; ++i) {
            stateIndices[i] = i;
        }
    }

    void resetBoundStepRequest() {
        if (!instance_data_ || desc_.bound_step_offset == UINT32_MAX) return;
        auto* raw = reinterpret_cast<double*>(
            reinterpret_cast<unsigned char*>(instance_data_) + desc_.bound_step_offset);
        *raw = std::numeric_limits<double>::infinity();
    }

    bool limitingApplied(uint32_t evalResult) const {
        return osdiLimitingRhsEnabled() && (evalResult & EVAL_RET_FLAG_LIM) != 0;
    }

    uint32_t limitingEvalFlags(uint32_t baseFlags) const {
        if (!osdiLimitingRhsEnabled()) return baseFlags;
        uint32_t flags = baseFlags | ENABLE_LIM;
        if (!limiting_initialized_) flags |= INIT_LIM;
        return flags;
    }

    void noteLimitingEval() {
        if (osdiLimitingRhsEnabled()) limiting_initialized_ = true;
    }

    void checkInitResult(const OsdiInitInfo& res, const std::string& stage) const {
        if ((res.flags & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI " + stage + " setup aborted with $fatal");
        }
        if (res.num_errors != 0) {
            std::string message =
                "OSDI " + stage + " setup reported " +
                std::to_string(res.num_errors) + " initialization error(s)";
            if (res.errors) {
                for (uint32_t i = 0; i < res.num_errors; ++i) {
                    const OsdiInitError& error = res.errors[i];
                    message += i == 0 ? ": " : ", ";
                    if (error.code == INIT_ERR_OUT_OF_BOUNDS &&
                        desc_.param_opvar &&
                        error.payload.parameter_id <
                            desc_.num_params + desc_.num_opvars) {
                        const OsdiParamOpvar& parameter =
                            desc_.param_opvar[error.payload.parameter_id];
                        const char* name =
                            parameter.name && parameter.name[0]
                            ? parameter.name[0]
                            : "<unnamed>";
                        message += "parameter '" + std::string(name) +
                            "' is out of bounds";
                    } else {
                        message += "error code " +
                            std::to_string(error.code);
                    }
                }
            }
            throw std::runtime_error(message);
        }
    }

    void applyParams(
        void* data,
        const ParamMap& params,
        bool instanceParams,
        std::vector<std::unique_ptr<char[]>>& stringParams) {
        if (!data || !desc_.access || !desc_.param_opvar) return;
        for (uint32_t id = 0; id < desc_.num_params; ++id) {
            const OsdiParamOpvar& param = desc_.param_opvar[id];
            const bool isInstance = (param.flags & PARA_KIND_INST) != 0;
            if (isInstance != instanceParams) {
                continue;
            }

            auto found = findParamValue(param, params);
            if (found == params.end()) {
                continue;
            }
            if (!instanceParams && !bind_full_model_params_ && !isTypeParam(param)) {
                continue;
            }
            const uint32_t flags = instanceParams ? (ACCESS_FLAG_SET | ACCESS_FLAG_INSTANCE) : ACCESS_FLAG_SET;
            void* raw = instanceParams
                ? desc_.access(data, nullptr, id, flags)
                : desc_.access(nullptr, data, id, flags);
            if (!raw) continue;
            const uint32_t ty = param.flags & PARA_TY_MASK;
            try {
                std::string normalizedValue = found->second;
                normalizedValue.erase(
                    std::remove_if(normalizedValue.begin(), normalizedValue.end(), [](char c) {
                        return c == '\'' || c == '"' || c == '{' || c == '}';
                    }),
                    normalizedValue.end());
                if (ty == PARA_TY_REAL) {
                    *reinterpret_cast<double*>(raw) = Utils::parseValue(normalizedValue);
                } else if (ty == PARA_TY_INT) {
                    *reinterpret_cast<int32_t*>(raw) = static_cast<int32_t>(Utils::parseValue(normalizedValue));
                } else if (ty == PARA_TY_STR) {
                    auto text = std::make_unique<char[]>(found->second.size() + 1);
                    std::memcpy(text.get(), found->second.c_str(), found->second.size() + 1);
                    *reinterpret_cast<char**>(raw) = text.get();
                    stringParams.push_back(std::move(text));
                }
            } catch (const std::exception& ex) {
                const char* parameterName =
                    param.name && param.name[0] ? param.name[0] : "<unnamed>";
                throw std::runtime_error(
                    "OSDI parameter '" + std::string(parameterName) +
                    "' has an unresolved or invalid value '" + found->second +
                    "': " + ex.what());
            }
        }
    }

    ParamMap::const_iterator findParamValue(const OsdiParamOpvar& param, const ParamMap& params) const {
        if (!param.name) return params.end();
        const uint32_t nameCount = 1 + param.num_alias;
        for (uint32_t i = 0; i < nameCount; ++i) {
            if (!param.name[i]) continue;
            auto it = params.find(lower(param.name[i]));
            if (it != params.end()) return it;
        }
        return params.end();
    }

    bool isTypeParam(const OsdiParamOpvar& param) const {
        if (!param.name) return false;
        const uint32_t nameCount = 1 + param.num_alias;
        for (uint32_t i = 0; i < nameCount; ++i) {
            if (param.name[i] && lower(param.name[i]) == "type") return true;
        }
        return false;
    }

    bool shouldBindFullModelParams(const ParamMap& params) const {
        (void)params;
        // Binding is an explicit all-or-selector policy. With full binding
        // enabled every recognized value must parse or construction fails;
        // disabling it binds only topology selectors such as `type`.
        return fullOsdiModelParamBindingEnabled();
    }

    void buildNodeMapping() {
        node_mapping_.assign(desc_.num_nodes, UINT32_MAX);
        inactive_global_nodes_.clear();
        uint32_t nextLocal = 1; // OSDI convention reserves node 0 for ground.
        for (uint32_t i = 0; i < desc_.num_nodes; ++i) {
            if (i < nodes_.size() && nodes_[i] < 0) {
                node_mapping_[i] = 0;
            } else {
                node_mapping_[i] = nextLocal++;
            }
        }

        if (desc_.num_collapsible && desc_.collapsible && desc_.collapsed_offset != UINT32_MAX) {
            bool* collapsed = reinterpret_cast<bool*>(reinterpret_cast<unsigned char*>(instance_data_) + desc_.collapsed_offset);
            for (uint32_t i = 0; i < desc_.num_collapsible; ++i) {
                if (!collapsed[i]) continue;
                const uint32_t from = desc_.collapsible[i].node_1;
                const uint32_t to = desc_.collapsible[i].node_2;
                if (from >= node_mapping_.size()) continue;
                const uint32_t mappedTo = (to == UINT32_MAX || to >= node_mapping_.size()) ? 0 : node_mapping_[to];
                const uint32_t mappedFrom = node_mapping_[from];
                for (auto& mapped : node_mapping_) {
                    if (mapped == mappedFrom) mapped = mappedTo;
                }
            }
        }

        local_node_count_ = 0;
        for (uint32_t mapped : node_mapping_) {
            if (mapped != UINT32_MAX) local_node_count_ = std::max(local_node_count_, static_cast<size_t>(mapped) + 1);
        }
        local_node_count_ = std::max(local_node_count_, static_cast<size_t>(desc_.num_terminals) + 1);

        if (desc_.node_mapping_offset != UINT32_MAX) {
            auto* osdiNodeMap = reinterpret_cast<uint32_t*>(
                reinterpret_cast<unsigned char*>(instance_data_) + desc_.node_mapping_offset);
            for (uint32_t i = 0; i < desc_.num_nodes; ++i) {
                osdiNodeMap[i] = node_mapping_[i];
            }
        }

        for (uint32_t i = desc_.num_terminals; i < node_mapping_.size() && i < nodes_.size(); ++i) {
            if (nodes_[i] < 0) continue;
            if (globalNodeForLocal(node_mapping_[i]) != nodes_[i]) {
                inactive_global_nodes_.push_back(nodes_[i]);
            }
        }
        // Reservations for currently collapsed internal nodes no longer
        // appear in nodes_. Anchor them explicitly until a later setup refresh
        // activates them.
        for (uint32_t i = desc_.num_terminals;
             i < reserved_internal_unknowns_.size();
             ++i) {
            const int reserved = reserved_internal_unknowns_[i];
            if (reserved < 0) continue;
            if (std::find(nodes_.begin(), nodes_.end(), reserved) == nodes_.end() &&
                std::find(
                    inactive_global_nodes_.begin(),
                    inactive_global_nodes_.end(),
                    reserved) == inactive_global_nodes_.end()) {
                inactive_global_nodes_.push_back(reserved);
            }
        }
    }

    std::vector<double> daeBypassInputs(const VectorReal& x) const {
        std::vector<double> inputs;
        if (desc_.inputs && desc_.num_inputs != 0) {
            inputs.reserve(desc_.num_inputs);
            for (uint32_t i = 0; i < desc_.num_inputs; ++i) {
                const int first = globalNodeForDescriptorNode(desc_.inputs[i].node_1);
                const int second = globalNodeForDescriptorNode(desc_.inputs[i].node_2);
                const double firstValue = first >= 0 && first < x.getSize() ? x[first] : 0.0;
                const double secondValue = second >= 0 && second < x.getSize() ? x[second] : 0.0;
                inputs.push_back(firstValue - secondValue);
            }
            return inputs;
        }
        inputs.reserve(node_mapping_.size());
        for (uint32_t descriptorNode = 0; descriptorNode < node_mapping_.size(); ++descriptorNode) {
            const int global = globalNodeForDescriptorNode(descriptorNode);
            inputs.push_back(global >= 0 && global < x.getSize() ? x[global] : 0.0);
        }
        return inputs;
    }

    static bool sameDaeRequest(const DaeRequest& lhs, const DaeRequest& rhs) {
        return lhs.analysis == rhs.analysis && lhs.time == rhs.time &&
            lhs.staticResidual == rhs.staticResidual && lhs.dynamicResidual == rhs.dynamicResidual &&
            lhs.staticJacobian == rhs.staticJacobian && lhs.dynamicJacobian == rhs.dynamicJacobian &&
            lhs.enableLimiting == rhs.enableLimiting && lhs.nodeset == rhs.nodeset &&
            lhs.simulationGmin == rhs.simulationGmin &&
            lhs.sourceScaleFactor == rhs.sourceScaleFactor;
    }

    bool canBypassDae(const VectorReal& x, const DaeRequest& request) const {
        if (!request.allowBypass || request.highPrecision || request.readOnlyState ||
            metadata_.has_bound_step || metadata_.uses_abstime ||
            !dae_cache_valid_ || cached_dae_evaluation_.limitingApplied ||
            cached_dae_evaluation_.finishRequested ||
            cached_dae_evaluation_.stopRequested ||
            !sameDaeRequest(request, cached_dae_request_)) {
            return false;
        }
        const auto inputs = daeBypassInputs(x);
        if (inputs.size() != cached_dae_inputs_.size()) return false;
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            const double scale = std::max(std::abs(inputs[i]), std::abs(cached_dae_inputs_[i]));
            const double tolerance = request.bypassAbsoluteTolerance +
                request.bypassRelativeTolerance * scale;
            if (std::abs(inputs[i] - cached_dae_inputs_[i]) > tolerance) return false;
        }
        return true;
    }

    void storeDaeCache(
        const VectorReal& x,
        const DaeRequest& request,
        const DaeEvaluation& evaluation) {
        if (!request.allowBypass || request.highPrecision || request.readOnlyState ||
            evaluation.limitingApplied ||
            evaluation.finishRequested ||
            evaluation.stopRequested) {
            return;
        }
        cached_dae_request_ = request;
        cached_dae_evaluation_ = evaluation;
        cached_dae_evaluation_.bypassed = false;
        cached_dae_inputs_ = daeBypassInputs(x);
        dae_cache_valid_ = true;
    }

    void invalidateDaeCache() {
        dae_cache_valid_ = false;
        cached_dae_inputs_.clear();
        cached_dae_evaluation_.clear();
    }

    int terminalForLocal(uint32_t local) const {
        for (uint32_t term = 0; term < desc_.num_terminals && term < node_mapping_.size(); ++term) {
            if (node_mapping_[term] == local) return static_cast<int>(term);
        }
        return -1;
    }

    int globalNodeForLocal(uint32_t local) const {
        if (local == UINT32_MAX || local == 0) return -1;
        for (uint32_t i = 0; i < node_mapping_.size() && i < nodes_.size(); ++i) {
            if (node_mapping_[i] == local) return nodes_[static_cast<size_t>(i)];
        }
        return -1;
    }

    int globalNodeForDescriptorNode(uint32_t descriptorNode) const {
        if (descriptorNode >= node_mapping_.size()) return -1;
        const uint32_t local = node_mapping_[descriptorNode];
        return globalNodeForLocal(local);
    }

    OsdiSimParas simulationParameters(
        double requestedGmin,
        int iteration = 0,
        bool initializeLimiting = false,
        double sourceScaleFactor = 1.0,
        const char* analysisName = "op",
        const char* analysisType = "static") const {
        simparam_values_[0] = initializeLimiting ? 1.0 : 0.0;
        simparam_values_[1] =
            std::isfinite(requestedGmin) && requestedGmin >= 0.0
            ? requestedGmin
            : simulation_gmin_;
        simparam_values_[2] = 0.0;
        simparam_values_[3] = nominal_temperature_c_;
        simparam_values_[4] = simulation_minr_;
        simparam_values_[5] = geometry_scale_;
        simparam_values_[6] = static_cast<double>(std::max(iteration, 0));
        simparam_values_[7] = 0.0;
        simparam_values_[8] = 1.0;
        simparam_values_[9] = std::clamp(sourceScaleFactor, 0.0, 1.0);
        simparam_values_[10] = relative_tolerance_;
        simparam_values_[11] = voltage_tolerance_;
        simparam_values_[12] = absolute_tolerance_;
        simparam_values_[13] = charge_tolerance_;
        simparam_values_[14] = flux_tolerance_;
        simparam_analysis_name_ = analysisName ? analysisName : "";
        simparam_analysis_type_ = analysisType ? analysisType : "";
        simparam_string_values_[0] = simparam_analysis_name_.data();
        simparam_string_values_[1] = simparam_analysis_type_.data();
        simparam_string_values_[2] = simparam_cwd_.data();
        return OsdiSimParas{
            simparam_names_.data(),
            simparam_values_.data(),
            simparam_string_names_.data(),
            simparam_string_values_.data()};
    }

    bool localIsFlow(uint32_t local) const {
        bool foundFlow = false;
        for (uint32_t i = 0; i < node_mapping_.size(); ++i) {
            if (node_mapping_[i] != local || !desc_.nodes) continue;
            if (!descriptorNodeIsFlow(i)) return false;
            foundFlow = true;
        }
        return foundFlow;
    }

    bool globalIsFlow(int global) const {
        bool foundFlow = false;
        for (uint32_t i = 0; i < node_mapping_.size() && i < nodes_.size(); ++i) {
            if (nodes_[i] != global || !desc_.nodes) continue;
            if (!descriptorNodeIsFlow(i)) return false;
            foundFlow = true;
        }
        return foundFlow;
    }

    bool descriptorNodeIsFlow(uint32_t descriptorNode) const {
        if (descriptorNode >= desc_.num_nodes || !desc_.nodes) return false;
        if (desc_.unknown_nature) {
            const OsdiNatureRef ref = desc_.unknown_nature[descriptorNode];
            if (ref.ref_type == NATREF_DISCIPLINE_FLOW) return true;
            if (ref.ref_type == NATREF_DISCIPLINE_POTENTIAL) return false;
        }
        return desc_.nodes[descriptorNode].is_flow;
    }

    void populateJacobianPointers(std::vector<double>& jacobian) {
        auto** ptrs = reinterpret_cast<double**>(reinterpret_cast<unsigned char*>(instance_data_) + desc_.jacobian_ptr_resist_offset);
        size_t jacIndex = 0;
        for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
            const bool isResistive = (desc_.jacobian_entries[e].flags & JACOBIAN_ENTRY_RESIST) != 0;
            ptrs[e] = isResistive && jacIndex < jacobian.size() ? &jacobian[jacIndex++] : nullptr;
        }
    }

    void populateTranJacobianPointers(std::vector<double>& jacobian) {
        auto** ptrs = reinterpret_cast<double**>(reinterpret_cast<unsigned char*>(instance_data_) + desc_.jacobian_ptr_resist_offset);
        size_t jacIndex = 0;
        for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
            ptrs[e] = jacIndex < jacobian.size() ? &jacobian[jacIndex++] : nullptr;
            if (desc_.jacobian_entries[e].react_ptr_off != UINT32_MAX) {
                auto** reactPtr = reinterpret_cast<double**>(
                    reinterpret_cast<unsigned char*>(instance_data_) + desc_.jacobian_entries[e].react_ptr_off);
                *reactPtr = ptrs[e];
            }
        }
    }

    void populateReactiveJacobianPointers(std::vector<double>& jacobian) {
        size_t jacIndex = 0;
        for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
            const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
            const bool isReactive = (entry.flags & JACOBIAN_ENTRY_REACT) != 0;
            if (!isReactive) continue;
            if (entry.react_ptr_off != UINT32_MAX) {
                auto** reactPtr = reinterpret_cast<double**>(
                    reinterpret_cast<unsigned char*>(instance_data_) + entry.react_ptr_off);
                *reactPtr = jacIndex < jacobian.size() ? &jacobian[jacIndex] : nullptr;
            }
            ++jacIndex;
        }
    }

    std::vector<double> localSolveFromGlobal(const VectorReal& x) const {
        std::vector<double> solve(local_node_count_, 0.0);
        fillLocalSolveFromGlobal(x, solve);
        return solve;
    }

    void fillLocalSolveFromGlobal(const VectorReal& x, std::vector<double>& solve) const {
        solve.assign(local_node_count_, 0.0);
        for (size_t i = 0; i < nodes_.size() && i < node_mapping_.size(); ++i) {
            const uint32_t local = node_mapping_[i];
            if (local == UINT32_MAX || local >= solve.size()) continue;
            solve[local] = nodes_[i] >= 0 ? x[nodes_[i]] : 0.0;
        }
    }

    bool loadReactiveResidualAt(
        const VectorReal& state,
        double time,
        std::vector<double>& reactive) {
        if (!desc_.load_residual_react) return false;
        fillLocalSolveFromGlobal(state, scratch_reactive_solve_);
        scratch_reactive_state_.assign(desc_.num_states, 0.0);
        OsdiSimParas paras = simulationParameters(simulation_gmin_);
        OsdiSimInfo info{
            paras,
            time,
            scratch_reactive_solve_.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            scratch_reactive_state_.empty() ? nullptr : scratch_reactive_state_.data(),
            CALC_REACT_RESIDUAL | ANALYSIS_TRAN
        };
        resetBoundStepRequest();
        const uint32_t ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }
        desc_.load_residual_react(instance_data_, model_data_, reactive.data());
        return true;
    }

    void transientStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        const TransientContext& ctx) {
        const double timeStep = ctx.timeStep;
        const double currentTime = ctx.currentTime;
        std::vector<double> solve = localSolveFromGlobal(x);

        const auto& history = *ctx.xHistory;
        const VectorReal& xPrev = history.back();
        const VectorReal& xPrev2 = (ctx.hasSecondHistory && history.size() >= 2)
            ? history[history.size() - 2]
            : history.back();
        double prevTime = currentTime - timeStep;
        double prev2Time = prevTime - timeStep;
        if (ctx.timeHistory && !ctx.timeHistory->empty()) {
            prevTime = ctx.timeHistory->back();
            if (ctx.timeHistory->size() >= 2) {
                prev2Time = (*ctx.timeHistory)[ctx.timeHistory->size() - 2];
            }
        }

        std::vector<double> reactPrev(local_node_count_, 0.0);
        std::vector<double> reactPrev2(local_node_count_, 0.0);
        loadReactiveResidualAt(xPrev, prevTime, reactPrev);
        if (ctx.hasSecondHistory && history.size() >= 2) {
            loadReactiveResidualAt(xPrev2, prev2Time, reactPrev2);
        }

        OsdiSimParas paras = simulationParameters(simulation_gmin_);

        OsdiSimInfo info{
            paras,
            currentTime,
            solve.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            next_state_.empty() ? nullptr : next_state_.data(),
            limitingEvalFlags(CALC_RESIST_RESIDUAL | CALC_RESIST_JACOBIAN |
                CALC_REACT_RESIDUAL | CALC_REACT_JACOBIAN |
                (osdiLimitingRhsEnabled() ? (CALC_RESIST_LIM_RHS | CALC_REACT_LIM_RHS) : 0u) |
                ANALYSIS_TRAN)
        };
        resetBoundStepRequest();
        uint32_t ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        noteLimitingEval();
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }

        double alpha = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && history.size() >= 2;
        const bool useTrap =
            ctx.method == TransientIntegrationMethod::Trapezoidal &&
            prev_react_derivative_valid_ &&
            prev_react_derivative_.size() == local_node_count_;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            alpha = 1.0 / timeStep;
            a1 = -alpha;
            a2 = 0.0;
            useSecond = false;
        }

        std::vector<double> residualResist(local_node_count_, 0.0);
        std::vector<double> residualReact(local_node_count_, 0.0);
        if (desc_.load_residual_resist) {
            desc_.load_residual_resist(instance_data_, model_data_, residualResist.data());
        }
        if (desc_.load_residual_react) {
            desc_.load_residual_react(instance_data_, model_data_, residualReact.data());
        }
        std::vector<double> limitResist(local_node_count_, 0.0);
        std::vector<double> limitReact(local_node_count_, 0.0);
        if (limitingApplied(ret)) {
            if (desc_.load_limit_rhs_resist) {
                desc_.load_limit_rhs_resist(instance_data_, model_data_, limitResist.data());
            }
            if (desc_.load_limit_rhs_react) {
                desc_.load_limit_rhs_react(instance_data_, model_data_, limitReact.data());
            }
        }
        std::vector<double> spiceRhs(local_node_count_, 0.0);
        const bool useSpiceRhs = use_spice_rhs_ && desc_.load_spice_rhs_tran;
        if (useSpiceRhs) {
            desc_.load_spice_rhs_tran(instance_data_, model_data_, spiceRhs.data(), solve.data(), alpha);
        }

        std::vector<double> jacTran(std::max<uint32_t>(desc_.num_jacobian_entries, 1), 0.0);
        bool wroteTranArray = false;
        if (
            standardOsdiTranJacobianEnabled() &&
            desc_.load_jacobian_tran &&
            desc_.jacobian_ptr_resist_offset != UINT32_MAX
        ) {
            populateTranJacobianPointers(jacTran);
            desc_.load_jacobian_tran(instance_data_, model_data_, alpha);
            wroteTranArray = true;
        }

        std::vector<double> jacResist(std::max<uint32_t>(desc_.num_resistive_jacobian_entries, 1), 0.0);
        bool wroteResistArray = false;
        if (!wroteTranArray && desc_.write_jacobian_array_resist) {
            desc_.write_jacobian_array_resist(instance_data_, model_data_, jacResist.data());
            wroteResistArray = true;
        } else if (!wroteTranArray && desc_.load_jacobian_resist && desc_.jacobian_ptr_resist_offset != UINT32_MAX) {
            populateJacobianPointers(jacResist);
            desc_.load_jacobian_resist(instance_data_, model_data_);
            wroteResistArray = true;
        }

        std::vector<double> jacReact(std::max<uint32_t>(desc_.num_reactive_jacobian_entries, 1), 0.0);
        bool wroteReactArray = false;
        if (!wroteTranArray && desc_.write_jacobian_array_react) {
            desc_.write_jacobian_array_react(instance_data_, model_data_, jacReact.data());
            for (double& value : jacReact) value *= alpha;
            wroteReactArray = true;
        } else if (!wroteTranArray && desc_.load_jacobian_react) {
            populateReactiveJacobianPointers(jacReact);
            desc_.load_jacobian_react(instance_data_, model_data_, alpha);
            wroteReactArray = true;
        }

        for (uint32_t localRow = 1; localRow < local_node_count_; ++localRow) {
            const int globalRow = globalNodeForLocal(localRow);
            if (globalRow < 0 || localRow >= residualResist.size()) continue;

            double reactiveResidual = alpha * residualReact[localRow] + a1 * reactPrev[localRow];
            if (useSecond && localRow < reactPrev2.size()) {
                reactiveResidual += a2 * reactPrev2[localRow];
            }
            if (useTrap) {
                reactiveResidual -= prev_react_derivative_[localRow];
            }
            const double limitingCorrection = limitResist[localRow] + alpha * limitReact[localRow];
            const double residual = residualResist[localRow] + reactiveResidual - limitingCorrection;
            double rhs = useSpiceRhs ? spiceRhs[localRow] : -residual;
            size_t resistIndex = 0;
            size_t reactIndex = 0;
            for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
                const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
                double g = 0.0;
                const bool hasResistive = (entry.flags & JACOBIAN_ENTRY_RESIST) != 0;
                const bool hasReactive = (entry.flags & JACOBIAN_ENTRY_REACT) != 0;
                if (wroteTranArray) {
                    if ((hasResistive || hasReactive) && e < jacTran.size()) {
                        g = jacTran[e];
                    }
                } else if (hasResistive) {
                    if (wroteResistArray && resistIndex < jacResist.size()) {
                        g += jacResist[resistIndex];
                    }
                    ++resistIndex;
                }
                if (!wroteTranArray && hasReactive) {
                    if (wroteReactArray && reactIndex < jacReact.size()) {
                        g += jacReact[reactIndex];
                    }
                    ++reactIndex;
                }
                if (entry.nodes.node_2 >= node_mapping_.size() || entry.nodes.node_1 >= node_mapping_.size()) continue;
                const uint32_t rowLocal = node_mapping_[entry.nodes.node_1];
                const uint32_t colLocal = node_mapping_[entry.nodes.node_2];
                if (rowLocal != localRow) continue;
                const int colGlobal = globalNodeForLocal(colLocal);
                const double vcol = colLocal < solve.size() ? solve[colLocal] : 0.0;
                if (!useSpiceRhs) {
                    rhs += g * vcol;
                }
                if (colGlobal >= 0) {
                    J.add(globalRow, colGlobal, g);
                }
            }
            b.add(globalRow, rhs);
        }
    }

    void legacyDcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x) {
        int n = static_cast<int>(nodes_.size());
        std::vector<double> voltages(n);
        for (int i = 0; i < n; ++i) {
            voltages[i] = (nodes_[i] >= 0) ? x[nodes_[i]] : 0.0;
        }

        std::vector<double> currents(n, 0.0);
        std::vector<double> charges(n, 0.0);
        std::vector<double> jacobian(n * n, 0.0);
        desc_.legacy_evaluate(instance_data_, voltages.data(), currents.data(), charges.data(), jacobian.data());

        for (int i = 0; i < n; ++i) {
            if (nodes_[i] < 0) continue;
            for (int j = 0; j < n; ++j) {
                if (nodes_[j] < 0) continue;
                J.add(nodes_[i], nodes_[j], jacobian[i * n + j]);
            }
            double rhs = -currents[i];
            for (int j = 0; j < n; ++j) {
                rhs += jacobian[i * n + j] * voltages[j];
            }
            b.add(nodes_[i], rhs);
        }
    }
};

} // namespace gspice

#endif
