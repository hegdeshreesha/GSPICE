#ifndef GSPICE_OSDI_DEVICE_HPP
#define GSPICE_OSDI_DEVICE_HPP

#include "device.hpp"
#include "osdi.h"
#include "osdi_metadata.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
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
        bool useSpiceRhs = false)
        : Device(name),
          desc_(desc),
          metadata_(desc),
          nodes_(nodes),
          use_limiting_rhs_(useLimitingRhs),
          use_tran_jacobian_(useTranJacobian),
          bind_full_model_params_override_(bindFullModelParams),
          use_spice_rhs_(useSpiceRhs),
          temperature_k_(temperatureC + 273.15) {
        if (desc_.legacy_evaluate) {
            instance_data_ = desc_.legacy_create_instance ? desc_.legacy_create_instance(nullptr) : nullptr;
            return;
        }
        validateRealOsdiDescriptor();
        bind_full_model_params_ = shouldBindFullModelParams(modelParams);
        model_storage_.resize(wordsFor(desc_.model_size));
        instance_storage_.resize(wordsFor(desc_.instance_size));
        model_data_ = model_storage_.empty() ? nullptr : model_storage_.data();
        instance_data_ = instance_storage_.empty() ? nullptr : instance_storage_.data();

        applyParams(model_data_, modelParams, false);
        setupModel();
        applyParams(instance_data_, instanceParams, true);
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

    // Allocate only the uncollapsed internal OSDI unknowns after ordinary
    // circuit nodes are known. This is the same MNA namespace used by branch
    // currents, so hidden compact-model unknowns do not leak into waveform
    // output and collapsed nodes never leave unused global rows behind.
    int bindInternalUnknowns(const std::function<int()>& allocateUnknown) {
        if (desc_.legacy_evaluate || internal_unknowns_bound_) return 0;
        if (nodes_.size() < desc_.num_terminals) {
            throw std::logic_error("OSDI terminals are not bound");
        }
        nodes_.resize(desc_.num_nodes, -2);
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
            if (nodes_[terminal] < 0) {
                rootUnknown[root] = -1;
            } else if (rootUnknown.find(root) == rootUnknown.end()) {
                rootUnknown[root] = nodes_[terminal];
            }
        }
        int allocated = 0;
        for (uint32_t node = 0; node < desc_.num_nodes; ++node) {
            const uint32_t root = resolve(node);
            if (root == UINT32_MAX) {
                nodes_[node] = -1;
                continue;
            }
            auto found = rootUnknown.find(root);
            if (found == rootUnknown.end()) {
                const int unknown = allocateUnknown();
                rootUnknown[root] = unknown;
                nodes_[node] = unknown;
                ++allocated;
            } else {
                nodes_[node] = found->second;
            }
        }
        internal_unknowns_bound_ = true;
        buildNodeMapping();
        prev_react_.assign(local_node_count_, 0.0);
        prev2_react_.assign(local_node_count_, 0.0);
        prev_react_derivative_.assign(local_node_count_, 0.0);
        read_only_snapshot_.resize(transientStateBytes());
        invalidateDaeCache();
        return allocated;
    }

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
        std::vector<double> scratchState;
        double* nextState = next_state_.empty() ? nullptr : next_state_.data();
        if (request.readOnlyState && desc_.num_states != 0) {
            scratchState.assign(desc_.num_states, 0.0);
            nextState = scratchState.data();
        }
        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};
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
        if (!request.readOnlyState) noteLimitingEval();
        if ((result & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }
        evaluation.limitingApplied = requestLimiting && (result & EVAL_RET_FLAG_LIM) != 0;
        evaluation.maximumTimeStep = transientBoundStep();

        std::vector<double> staticResidual(local_node_count_, 0.0);
        std::vector<double> dynamicResidual(local_node_count_, 0.0);
        if (request.staticResidual && desc_.load_residual_resist) {
            desc_.load_residual_resist(instance_data_, model_data_, staticResidual.data());
        }
        if (request.dynamicResidual && desc_.load_residual_react) {
            desc_.load_residual_react(instance_data_, model_data_, dynamicResidual.data());
        }
        if (evaluation.limitingApplied && request.staticResidual && desc_.load_limit_rhs_resist) {
            std::vector<double> correction(local_node_count_, 0.0);
            desc_.load_limit_rhs_resist(instance_data_, model_data_, correction.data());
            for (size_t i = 0; i < staticResidual.size(); ++i) staticResidual[i] -= correction[i];
        }
        if (evaluation.limitingApplied && request.dynamicResidual && desc_.load_limit_rhs_react) {
            std::vector<double> correction(local_node_count_, 0.0);
            desc_.load_limit_rhs_react(instance_data_, model_data_, correction.data());
            for (size_t i = 0; i < dynamicResidual.size(); ++i) dynamicResidual[i] -= correction[i];
        }
        double staticGroundResidual = 0.0;
        double dynamicGroundResidual = 0.0;
        bool hasGroundTerminal = false;
        for (uint32_t terminal = 0;
             terminal < desc_.num_terminals && terminal < node_mapping_.size();
             ++terminal) {
            hasGroundTerminal = hasGroundTerminal || node_mapping_[terminal] == 0;
        }
        for (uint32_t local = 1; local < local_node_count_; ++local) {
            const int equation = globalNodeForLocal(local);
            if (equation < 0) continue;
            if (request.staticResidual) {
                evaluation.staticResidual.push_back({equation, staticResidual[local]});
                staticGroundResidual -= staticResidual[local];
            }
            if (request.dynamicResidual) {
                evaluation.dynamicResidual.push_back({equation, dynamicResidual[local], 0});
                dynamicGroundResidual -= dynamicResidual[local];
            }
        }
        // OSDI reserves local node zero for ground, but residual loaders do not
        // need to materialize its discarded MNA row. Reconstruct that row for
        // conservation/derivative audits; normal stamping ignores equation -1.
        if (request.staticResidual && hasGroundTerminal) {
            evaluation.staticResidual.push_back({-1, staticGroundResidual});
        }
        if (request.dynamicResidual && hasGroundTerminal) {
            evaluation.dynamicResidual.push_back({-1, dynamicGroundResidual, 0});
        }

        std::vector<double> staticJacobian(
            std::max<uint32_t>(desc_.num_resistive_jacobian_entries, 1), 0.0);
        std::vector<double> dynamicJacobian(
            std::max<uint32_t>(desc_.num_reactive_jacobian_entries, 1), 0.0);
        bool haveStaticJacobian = false;
        bool haveDynamicJacobian = false;
        if (request.staticJacobian && desc_.write_jacobian_array_resist) {
            desc_.write_jacobian_array_resist(instance_data_, model_data_, staticJacobian.data());
            haveStaticJacobian = true;
        } else if (request.staticJacobian && desc_.load_jacobian_resist &&
                   desc_.jacobian_ptr_resist_offset != UINT32_MAX) {
            populateJacobianPointers(staticJacobian);
            desc_.load_jacobian_resist(instance_data_, model_data_);
            haveStaticJacobian = true;
        }
        if (request.dynamicJacobian && desc_.write_jacobian_array_react) {
            desc_.write_jacobian_array_react(instance_data_, model_data_, dynamicJacobian.data());
            haveDynamicJacobian = true;
        } else if (request.dynamicJacobian && desc_.load_jacobian_react) {
            populateReactiveJacobianPointers(dynamicJacobian);
            desc_.load_jacobian_react(instance_data_, model_data_, 1.0);
            haveDynamicJacobian = true;
        }

        size_t staticIndex = 0;
        size_t dynamicIndex = 0;
        std::unordered_map<int, double> staticGroundJacobian;
        std::unordered_map<int, double> dynamicGroundJacobian;
        for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
            const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
            const bool hasStatic = (entry.flags & JACOBIAN_ENTRY_RESIST) != 0;
            const bool hasDynamic = (entry.flags & JACOBIAN_ENTRY_REACT) != 0;
            double staticValue = 0.0;
            double dynamicValue = 0.0;
            if (hasStatic) {
                if (haveStaticJacobian && staticIndex < staticJacobian.size()) {
                    staticValue = staticJacobian[staticIndex];
                }
                ++staticIndex;
            }
            if (hasDynamic) {
                if (haveDynamicJacobian && dynamicIndex < dynamicJacobian.size()) {
                    dynamicValue = dynamicJacobian[dynamicIndex];
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
            if (request.staticJacobian && hasStatic && staticValue != 0.0) {
                evaluation.staticJacobian.push_back({equation, unknown, staticValue});
                staticGroundJacobian[unknown] -= staticValue;
            }
            if (request.dynamicJacobian && hasDynamic && dynamicValue != 0.0) {
                evaluation.dynamicJacobian.push_back({equation, unknown, dynamicValue, 0});
                dynamicGroundJacobian[unknown] -= dynamicValue;
            }
        }
        if (request.staticJacobian && hasGroundTerminal) {
            for (const auto& entry : staticGroundJacobian) {
                evaluation.staticJacobian.push_back({-1, entry.first, entry.second});
            }
        }
        if (request.dynamicJacobian && hasGroundTerminal) {
            for (const auto& entry : dynamicGroundJacobian) {
                evaluation.dynamicJacobian.push_back({-1, entry.first, entry.second, 0});
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

    bool prefersDampedAutoTransient() const override { return true; }

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

        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};
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
        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};
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
        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};
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
        const double freq = std::max(omega / (2.0 * 3.14159265358979323846), 0.0);
        desc_.load_noise(instance_data_, model_data_, freq, densities.data());
        for (uint32_t i = 0; i < desc_.num_noise_src; ++i) {
            if (densities[i] <= 0.0) continue;
            const int nodePos = globalNodeForDescriptorNode(desc_.noise_sources[i].nodes.node_1);
            const int nodeNeg = globalNodeForDescriptorNode(desc_.noise_sources[i].nodes.node_2);
            if (nodePos < 0 && nodeNeg < 0) continue;
            const char* noiseName = desc_.noise_sources[i].name ? desc_.noise_sources[i].name : "noise";
            sources.push_back({name_ + "." + noiseName, nodePos, nodeNeg, densities[i]});
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
            std::vector<double> reactNow(local_node_count_, 0.0);
            if (loadReactiveResidualAt(x, currentTime, reactNow)) {
                std::vector<double> deriv(local_node_count_, 0.0);
                const bool useTrap =
                    ctx.method == TransientIntegrationMethod::Trapezoidal &&
                    prev_react_derivative_valid_ &&
                    prev_react_.size() == reactNow.size();
                const bool useSecond =
                    ctx.hasSecondHistory &&
                    prev_react_.size() == reactNow.size() &&
                    prev2_react_.size() == reactNow.size();
                if (useTrap) {
                    for (size_t i = 0; i < reactNow.size(); ++i) {
                        deriv[i] = ctx.a0 * reactNow[i] + ctx.a1 * prev_react_[i] - prev_react_derivative_[i];
                    }
                } else if (useSecond) {
                    for (size_t i = 0; i < reactNow.size(); ++i) {
                        deriv[i] = ctx.a0 * reactNow[i] + ctx.a1 * prev_react_[i] + ctx.a2 * prev2_react_[i];
                    }
                } else {
                    const double alpha = 1.0 / ctx.timeStep;
                    for (size_t i = 0; i < reactNow.size(); ++i) {
                        const double prev = prev_react_.size() == reactNow.size() ? prev_react_[i] : reactNow[i];
                        deriv[i] = alpha * (reactNow[i] - prev);
                    }
                }
                prev2_react_ = prev_react_;
                prev_react_ = reactNow;
                prev_react_derivative_ = deriv;
                prev_react_derivative_valid_ = true;
            }
        }
        invalidateDaeCache();
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
        return model_storage_.size() * sizeof(std::max_align_t) +
            instance_storage_.size() * sizeof(std::max_align_t) +
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
        write(model_storage_.data(), model_storage_.size() * sizeof(std::max_align_t));
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
        read(model_storage_.data(), model_storage_.size() * sizeof(std::max_align_t));
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
    std::vector<std::max_align_t> model_storage_;
    std::vector<std::max_align_t> instance_storage_;
    void* model_data_ = nullptr;
    void* instance_data_ = nullptr;
    std::vector<uint32_t> node_mapping_;
    std::vector<int> inactive_global_nodes_;
    size_t local_node_count_ = 0;
    std::vector<std::unique_ptr<char[]>> string_params_;
    std::vector<double> prev_state_;
    std::vector<double> next_state_;
    std::vector<double> prev_react_;
    std::vector<double> prev2_react_;
    std::vector<double> prev_react_derivative_;
    std::vector<std::byte> read_only_snapshot_;
    bool prev_react_derivative_valid_ = false;
    bool bind_full_model_params_ = false;
    bool use_limiting_rhs_ = false;
    bool use_tran_jacobian_ = false;
    bool bind_full_model_params_override_ = false;
    bool use_spice_rhs_ = false;
    bool limiting_initialized_ = false;
    bool internal_unknowns_bound_ = false;
    double temperature_k_ = 300.15;
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
        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};
        OsdiInitInfo res{0, 0, nullptr};
        desc_.setup_model(osdiHandle(), model_data_, &paras, &res);
        checkInitResult(res, "model");
    }

    void setupInstance() {
        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};
        OsdiInitInfo res{0, 0, nullptr};
        desc_.setup_instance(osdiHandle(), instance_data_, model_data_, temperature_k_, desc_.num_terminals, &paras, &res);
        checkInitResult(res, "instance");
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
            throw std::runtime_error("OSDI " + stage + " setup reported " + std::to_string(res.num_errors) + " initialization error(s)");
        }
    }

    void applyParams(void* data, const ParamMap& params, bool instanceParams) {
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
                    string_params_.push_back(std::move(text));
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
            lhs.evaluationEpoch == rhs.evaluationEpoch;
    }

    bool canBypassDae(const VectorReal& x, const DaeRequest& request) const {
        if (!request.allowBypass || request.highPrecision || request.readOnlyState ||
            !dae_cache_valid_ || cached_dae_evaluation_.limitingApplied ||
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
            evaluation.limitingApplied) {
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
        for (size_t i = 0; i < nodes_.size() && i < node_mapping_.size(); ++i) {
            const uint32_t local = node_mapping_[i];
            if (local == UINT32_MAX || local >= solve.size()) continue;
            solve[local] = nodes_[i] >= 0 ? x[nodes_[i]] : 0.0;
        }
        return solve;
    }

    bool loadReactiveResidualAt(
        const VectorReal& state,
        double time,
        std::vector<double>& reactive) {
        if (!desc_.load_residual_react) return false;
        std::vector<double> solve = localSolveFromGlobal(state);
        std::vector<double> scratchState(desc_.num_states, 0.0);
        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};
        OsdiSimInfo info{
            paras,
            time,
            solve.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            scratchState.empty() ? nullptr : scratchState.data(),
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

        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};

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
