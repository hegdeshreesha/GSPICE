#ifndef GSPICE_OSDI_DEVICE_HPP
#define GSPICE_OSDI_DEVICE_HPP

#include "device.hpp"
#include "osdi.h"
#include "utils.hpp"
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <cstring>
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
        const ParamMap& instanceParams = {})
        : Device(name), desc_(desc), nodes_(nodes) {
        if (desc_.legacy_evaluate) {
            instance_data_ = desc_.legacy_create_instance ? desc_.legacy_create_instance(nullptr) : nullptr;
            return;
        }
        validateRealOsdiDescriptor();
        model_storage_.resize(wordsFor(desc_.model_size));
        instance_storage_.resize(wordsFor(desc_.instance_size));
        model_data_ = model_storage_.empty() ? nullptr : model_storage_.data();
        instance_data_ = instance_storage_.empty() ? nullptr : instance_storage_.data();

        applyParams(model_data_, modelParams, false);
        setupModel();
        applyParams(instance_data_, instanceParams, true);
        setupInstance();
        buildNodeMapping();
        prev_state_.assign(desc_.num_states, 0.0);
        next_state_.assign(desc_.num_states, 0.0);
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
            transientStamp(J, b, x, timeStep, currentTime, x_hist.back());
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
            CALC_RESIST_RESIDUAL | CALC_RESIST_JACOBIAN | CALC_OP | ANALYSIS_DC | ANALYSIS_STATIC
        };
        uint32_t ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }

        std::vector<double> residual(local_node_count_, 0.0);
        if (desc_.load_residual_resist) {
            desc_.load_residual_resist(instance_data_, model_data_, residual.data());
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

        for (size_t term = 0; term < nodes_.size() && term < node_mapping_.size(); ++term) {
            const int globalRow = nodes_[term];
            const uint32_t localRow = node_mapping_[term];
            if (globalRow < 0 || localRow == UINT32_MAX || localRow >= residual.size()) continue;

            double rhs = -residual[localRow];
            size_t jacIndex = 0;
            for (uint32_t e = 0; wroteArray && e < desc_.num_jacobian_entries; ++e) {
                const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
                const bool isResistive = (entry.flags & (JACOBIAN_ENTRY_RESIST | JACOBIAN_ENTRY_RESIST_CONST)) != 0;
                if (!isResistive) continue;
                const double g = jacIndex < jacobian.size() ? jacobian[jacIndex] : 0.0;
                ++jacIndex;

                if (entry.nodes.node_2 >= node_mapping_.size() || entry.nodes.node_1 >= node_mapping_.size()) continue;
                const uint32_t rowLocal = node_mapping_[entry.nodes.node_2];
                const uint32_t colLocal = node_mapping_[entry.nodes.node_1];
                if (rowLocal != localRow || colLocal == UINT32_MAX) continue;

                const int colTerm = terminalForLocal(colLocal);
                const double vcol = colLocal < solve.size() ? solve[colLocal] : 0.0;
                rhs -= g * vcol;
                if (colTerm >= 0 && nodes_[static_cast<size_t>(colTerm)] >= 0) {
                    J.add(globalRow, nodes_[static_cast<size_t>(colTerm)], -g);
                }
            }
            b.add(globalRow, rhs);
        }
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)J;
        (void)b;
        (void)omega;
        (void)x_dc;
    }

    void acceptTransientStep(const VectorReal& x, double currentTime) override {
        (void)x;
        (void)currentTime;
        if (!next_state_.empty()) {
            prev_state_ = next_state_;
        }
    }

private:
    OsdiDescriptor desc_{};
    std::vector<int> nodes_;
    std::vector<std::max_align_t> model_storage_;
    std::vector<std::max_align_t> instance_storage_;
    void* model_data_ = nullptr;
    void* instance_data_ = nullptr;
    std::vector<uint32_t> node_mapping_;
    size_t local_node_count_ = 0;
    std::vector<std::unique_ptr<char[]>> string_params_;
    std::vector<double> prev_state_;
    std::vector<double> next_state_;

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
        desc_.setup_instance(osdiHandle(), instance_data_, model_data_, 300.15, desc_.num_terminals, &paras, &res);
        checkInitResult(res, "instance");
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
            if (!instanceParams && !isTypeParam(param)) {
                continue;
            }

            const uint32_t flags = instanceParams ? (ACCESS_FLAG_SET | ACCESS_FLAG_INSTANCE) : ACCESS_FLAG_SET;
            void* raw = instanceParams
                ? desc_.access(data, nullptr, id, flags)
                : desc_.access(nullptr, data, id, flags);
            if (!raw) continue;
            const uint32_t ty = param.flags & PARA_TY_MASK;
            try {
                if (ty == PARA_TY_REAL) {
                    *reinterpret_cast<double*>(raw) = Utils::parseValue(found->second);
                } else if (ty == PARA_TY_INT) {
                    *reinterpret_cast<int32_t*>(raw) = static_cast<int32_t>(Utils::parseValue(found->second));
                } else if (ty == PARA_TY_STR) {
                    auto text = std::make_unique<char[]>(found->second.size() + 1);
                    std::memcpy(text.get(), found->second.c_str(), found->second.size() + 1);
                    *reinterpret_cast<char**>(raw) = text.get();
                    string_params_.push_back(std::move(text));
                }
            } catch (const std::exception&) {
                // PDK cards often contain expressions or unresolved .param symbols.
                // Leave those at the OpenVAF model default until GSPICE has a full expression evaluator.
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

    void buildNodeMapping() {
        node_mapping_.assign(desc_.num_nodes, UINT32_MAX);
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
    }

    int terminalForLocal(uint32_t local) const {
        for (uint32_t term = 0; term < desc_.num_terminals && term < node_mapping_.size(); ++term) {
            if (node_mapping_[term] == local) return static_cast<int>(term);
        }
        return -1;
    }

    void populateJacobianPointers(std::vector<double>& jacobian) {
        auto** ptrs = reinterpret_cast<double**>(reinterpret_cast<unsigned char*>(instance_data_) + desc_.jacobian_ptr_resist_offset);
        size_t jacIndex = 0;
        for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
            const bool isResistive = (desc_.jacobian_entries[e].flags & (JACOBIAN_ENTRY_RESIST | JACOBIAN_ENTRY_RESIST_CONST)) != 0;
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
            const bool isReactive = (entry.flags & (JACOBIAN_ENTRY_REACT | JACOBIAN_ENTRY_REACT_CONST)) != 0;
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

    void transientStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const VectorReal& x_prev) {
        std::vector<double> solve = localSolveFromGlobal(x);
        std::vector<double> solvePrev = localSolveFromGlobal(x_prev);

        char* nullName = nullptr;
        char* nullStrName = nullptr;
        OsdiSimParas paras{&nullName, nullptr, &nullStrName, nullptr};

        std::vector<double> reactPrev(local_node_count_, 0.0);
        std::vector<double> scratchState(desc_.num_states, 0.0);
        OsdiSimInfo prevInfo{
            paras,
            currentTime - timeStep,
            solvePrev.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            scratchState.empty() ? nullptr : scratchState.data(),
            CALC_REACT_RESIDUAL | ANALYSIS_TRAN
        };
        uint32_t ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &prevInfo);
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }
        if (desc_.load_residual_react) {
            desc_.load_residual_react(instance_data_, model_data_, reactPrev.data());
        }

        OsdiSimInfo info{
            paras,
            currentTime,
            solve.data(),
            prev_state_.empty() ? nullptr : prev_state_.data(),
            next_state_.empty() ? nullptr : next_state_.data(),
            CALC_RESIST_RESIDUAL | CALC_RESIST_JACOBIAN |
                CALC_REACT_RESIDUAL | CALC_REACT_JACOBIAN |
                ANALYSIS_TRAN
        };
        ret = desc_.eval(osdiHandle(), instance_data_, model_data_, &info);
        if ((ret & EVAL_RET_FLAG_FATAL) != 0) {
            throw std::runtime_error("OSDI device '" + name_ + "' aborted with $fatal");
        }

        const double alpha = 1.0 / timeStep;
        std::vector<double> residualResist(local_node_count_, 0.0);
        std::vector<double> residualReact(local_node_count_, 0.0);
        if (desc_.load_residual_resist) {
            desc_.load_residual_resist(instance_data_, model_data_, residualResist.data());
        }
        if (desc_.load_residual_react) {
            desc_.load_residual_react(instance_data_, model_data_, residualReact.data());
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
            for (double& value : jacReact) value *= alpha;
            wroteReactArray = true;
        } else if (desc_.load_jacobian_react) {
            populateReactiveJacobianPointers(jacReact);
            desc_.load_jacobian_react(instance_data_, model_data_, alpha);
            wroteReactArray = true;
        }

        for (size_t term = 0; term < nodes_.size() && term < node_mapping_.size(); ++term) {
            const int globalRow = nodes_[term];
            const uint32_t localRow = node_mapping_[term];
            if (globalRow < 0 || localRow == UINT32_MAX || localRow >= residualResist.size()) continue;

            const double residual =
                residualResist[localRow] + alpha * (residualReact[localRow] - reactPrev[localRow]);
            double rhs = -residual;
            size_t resistIndex = 0;
            size_t reactIndex = 0;
            for (uint32_t e = 0; e < desc_.num_jacobian_entries; ++e) {
                const OsdiJacobianEntry& entry = desc_.jacobian_entries[e];
                double g = 0.0;
                if ((entry.flags & (JACOBIAN_ENTRY_RESIST | JACOBIAN_ENTRY_RESIST_CONST)) != 0) {
                    if (wroteResistArray && resistIndex < jacResist.size()) {
                        g += jacResist[resistIndex];
                    }
                    ++resistIndex;
                }
                if ((entry.flags & (JACOBIAN_ENTRY_REACT | JACOBIAN_ENTRY_REACT_CONST)) != 0) {
                    if (wroteReactArray && reactIndex < jacReact.size()) {
                        g += jacReact[reactIndex];
                    }
                    ++reactIndex;
                }
                if (entry.nodes.node_2 >= node_mapping_.size() || entry.nodes.node_1 >= node_mapping_.size()) continue;
                const uint32_t rowLocal = node_mapping_[entry.nodes.node_2];
                const uint32_t colLocal = node_mapping_[entry.nodes.node_1];
                if (rowLocal != localRow) continue;
                const int colTerm = terminalForLocal(colLocal);
                const double vcol = colLocal < solve.size() ? solve[colLocal] : 0.0;
                rhs -= g * vcol;
                if (colTerm >= 0 && nodes_[static_cast<size_t>(colTerm)] >= 0) {
                    J.add(globalRow, nodes_[static_cast<size_t>(colTerm)], -g);
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
