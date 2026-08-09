#ifndef GSPICE_GSDI_HPP
#define GSPICE_GSDI_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace gspice {

using GsdiParamMap = std::unordered_map<std::string, double>;

// Collapse decision sentinels for GsdiCollapseMap. A collapsible node pair is
// a pair of nodes whose difference is forced to zero (Verilog-A `V(x,y) <+ 0`),
// so the simulator merges the two columns into one unknown.
inline constexpr int GsdiNoCollapse = -1;
inline constexpr int GsdiCollapseToGround = -2;

enum class GsdiNodeRole {
    Terminal,    // external device terminal, participates in KCL
    Internal,    // hidden node created by the model (e.g. series resistance)
    Collapsible  // internal node declared collapsible with another node
};

struct GsdiNodeInfo {
    std::string name;
    int index = -1;          // local node index (position inside descriptor.nodes)
    GsdiNodeRole role = GsdiNodeRole::Terminal;
    int collapse_partner = GsdiNoCollapse;  // local index merged into, or GsdiCollapseToGround
};

struct GsdiParameterInfo {
    std::string name;
    double default_value = 0.0;
    std::string units;
    std::string description;
    bool is_model_param = true;  // true: .MODEL card, false: instance
    bool has_min = false;
    double min_value = 0.0;
    bool has_max = false;
    double max_value = 0.0;
};

struct GsdiOpvarInfo {
    std::string name;
    std::string description;
};

enum class GsdiAnalysis {
    OperatingPoint,
    Dc,
    Transient,
    Ac,
    Noise,
    Pss,
    Pnoise
};

struct GsdiModelDescriptor {
    std::string model_type;
    std::string version = "1.0";
    int terminal_count = 0;
    std::vector<GsdiNodeInfo> nodes;                 // terminal nodes first, then internal/collapsible
    std::vector<GsdiParameterInfo> parameters;       // model and instance parameters
    std::vector<GsdiOpvarInfo> opvars;               // computed operating-point variables
    std::vector<std::pair<int, int>> jacobian_pattern;  // static sparsity: (equation, unknown) local indices
    bool supports_op = true;
    bool supports_transient = false;
    bool supports_ac = false;
    bool supports_noise = false;
    bool supports_pss = false;
    bool supports_pnoise = false;
};

struct GsdiModelCard {
    std::string name;
    std::string type;
    GsdiParamMap parameters;
};

struct GsdiEvalRequest {
    GsdiAnalysis analysis = GsdiAnalysis::OperatingPoint;
    const double* solution = nullptr;
    std::size_t solution_size = 0;
    double time = 0.0;
    double time_step = 0.0;
    double omega = 0.0;
    bool residual = true;
    bool jacobian = true;
    bool dynamic_residual = false;
    bool dynamic_jacobian = false;
    bool noise = false;
    bool allow_bypass = true;
    bool read_only_state = false;
};

struct GsdiResidual {
    int equation = -1;
    double value = 0.0;
    int conservation_group = -1;
};

struct GsdiJacobian {
    int equation = -1;
    int unknown = -1;
    double value = 0.0;
    int conservation_group = -1;
};

struct GsdiNoiseSource {
    int node_pos = -1;
    int node_neg = -1;
    double spectral_density = 0.0;
    std::string name;
};

struct GsdiEvalResult {
    std::vector<GsdiResidual> static_residual;
    std::vector<GsdiJacobian> static_jacobian;
    std::vector<GsdiResidual> dynamic_residual;
    std::vector<GsdiJacobian> dynamic_jacobian;
    std::vector<GsdiNoiseSource> noise;
    std::vector<double> opvars;  // computed operating-point variables, descriptor.opvars order
    bool limiting_applied = false;
    bool bypassed = false;

    void clear() {
        static_residual.clear();
        static_jacobian.clear();
        dynamic_residual.clear();
        dynamic_jacobian.clear();
        noise.clear();
        opvars.clear();
        limiting_applied = false;
        bypassed = false;
    }

    bool finite() const {
        const auto finite_residuals = [](const auto& terms) {
            for (const auto& term : terms) {
                if (!std::isfinite(term.value)) return false;
            }
            return true;
        };
        for (const auto& source : noise) {
            if (!std::isfinite(source.spectral_density)) return false;
        }
        return finite_residuals(static_residual) && finite_residuals(static_jacobian) &&
               finite_residuals(dynamic_residual) && finite_residuals(dynamic_jacobian);
    }
};

// Maps a model's local node space (terminals + internal + collapsible nodes)
// onto the collapsed unknown vector the simulator stamps. Collapsed nodes are
// removed; their partner keeps an alias. The simulator passes the collapsed
// solution vector to evaluate(); the map tells the DAE adapter which local
// node index each collapsed-column value belongs to.
class GsdiCollapseMap {
public:
    GsdiCollapseMap() = default;

    // collapse_decisions: one entry per local node. GsdiNoCollapse keeps the
    // node, GsdiCollapseToGround merges it into the reference node, and any
    // other value merges into that partner's local index.
    GsdiCollapseMap(const GsdiModelDescriptor& descriptor,
                    const std::vector<int>& collapse_decisions) {
        const std::size_t local_count = descriptor.nodes.size();
        collapse_.assign(local_count, GsdiNoCollapse);
        value_source_.assign(local_count, GsdiNoCollapse);
        decisions_.assign(local_count, GsdiNoCollapse);
        expand_.clear();

        for (std::size_t local = 0; local < local_count; ++local) {
            const int decision =
                local < collapse_decisions.size() ? collapse_decisions[local] : GsdiNoCollapse;
            decisions_[local] = decision;
        }
        // Resolve keep/merge decisions transitively. A collapsed node's value
        // and its KCL equation are carried by the kept node its merge chain
        // terminates at (or by the reference node when the chain ends at
        // GsdiCollapseToGround). Cycles and out-of-range partners fail closed
        // by treating the node as merged away.
        for (std::size_t local = 0; local < local_count; ++local) {
            if (decisions_[local] == GsdiNoCollapse) {
                collapse_[local] = static_cast<int>(expand_.size());
                expand_.push_back(static_cast<int>(local));
                value_source_[local] = static_cast<int>(local);
            } else {
                collapse_[local] = GsdiNoCollapse;
                value_source_[local] = resolveSource(static_cast<int>(local));
            }
        }
    }

    // Builds the standard collapse plan from the model descriptor: every node
    // declared Collapsible merges into its declared collapse_partner, or into
    // the reference node when the partner is GsdiCollapseToGround. Nodes
    // declared Terminal/Internal stay as their own columns.
    static GsdiCollapseMap standard(const GsdiModelDescriptor& descriptor) {
        std::vector<int> decisions(descriptor.nodes.size(), GsdiNoCollapse);
        for (std::size_t local = 0; local < descriptor.nodes.size(); ++local) {
            const GsdiNodeInfo& node = descriptor.nodes[local];
            if (node.role != GsdiNodeRole::Collapsible) continue;
            const int partner = node.collapse_partner;
            const bool valid_partner =
                partner == GsdiCollapseToGround ||
                (partner >= 0 && static_cast<std::size_t>(partner) < descriptor.nodes.size() &&
                 partner != static_cast<int>(local));
            decisions[local] = valid_partner ? partner : GsdiNoCollapse;
        }
        return GsdiCollapseMap(descriptor, decisions);
    }

    // Collapsed unknown index for a local node, or GsdiNoCollapse(-1) when the
    // node was merged away (into a partner or into the reference node).
    int collapsedIndex(int local_node) const {
        return (local_node >= 0 && static_cast<std::size_t>(local_node) < collapse_.size())
                   ? collapse_[local_node] : GsdiNoCollapse;
    }

    // Reverse: local node that owns a collapsed unknown column.
    int localNodeOf(int collapsed_index) const {
        return (collapsed_index >= 0 && static_cast<std::size_t>(collapsed_index) < expand_.size())
                   ? expand_[collapsed_index] : GsdiNoCollapse;
    }

    // Kept local node whose column carries this node's value, or
    // GsdiNoCollapse(-1) when the node's value is the reference (ground).
    // The model reads its local voltage from nodes_[valueSource(local)].
    int valueSource(int local_node) const {
        return (local_node >= 0 && static_cast<std::size_t>(local_node) < value_source_.size())
                   ? value_source_[local_node] : GsdiNoCollapse;
    }

    // Number of local nodes in the model's node space (terminals + internal +
    // collapsible), i.e. the size of the voltage vector evaluate() receives.
    int localCount() const { return static_cast<int>(collapse_.size()); }

    int unknownCount() const { return static_cast<int>(expand_.size()); }

private:
    // Follows a node's merge chain to the kept node that carries its value,
    // or returns GsdiNoCollapse(-1) when the chain terminates at the
    // reference node / is invalid / cycles.
    int resolveSource(int local) const {
        int hops = 0;
        while (local >= 0 && static_cast<std::size_t>(local) < decisions_.size() &&
               decisions_[local] != GsdiNoCollapse) {
            if (++hops > static_cast<int>(decisions_.size())) return GsdiNoCollapse;
            if (decisions_[local] == GsdiCollapseToGround) return GsdiNoCollapse;
            local = decisions_[local];
        }
        if (local >= 0 && static_cast<std::size_t>(local) < decisions_.size() &&
            decisions_[local] == GsdiNoCollapse) {
            return local;
        }
        return GsdiNoCollapse;
    }

    std::vector<int> collapse_;      // local node -> collapsed unknown index
    std::vector<int> expand_;        // collapsed unknown index -> local node
    std::vector<int> value_source_;  // local node -> kept local node carrying its value
    std::vector<int> decisions_;     // raw merge decisions (kept for chain walking)
};

class GsdiInstance {
public:
    virtual ~GsdiInstance() = default;

    virtual bool evaluate(
        const GsdiEvalRequest& request,
        GsdiEvalResult& result) = 0;

    virtual std::size_t terminalCount() const = 0;
    virtual std::size_t internalNodeCount() const { return 0; }
    virtual std::size_t stateBytes() const { return 0; }
    virtual void saveState(std::byte*, std::size_t) const {}
    virtual void restoreState(const std::byte*, std::size_t) {}
};

class GsdiModel {
public:
    virtual ~GsdiModel() = default;
    virtual const GsdiModelDescriptor& descriptor() const = 0;
    virtual std::unique_ptr<GsdiInstance> createInstance(
        const GsdiModelCard& card,
        const std::vector<int>& terminal_nodes) const = 0;
};

using GsdiModelFactory = std::function<std::unique_ptr<GsdiModel>()>;

// Compatibility names retained for the initial GSDI prototype.
using GsdiEvaluationRequest = GsdiEvalRequest;
using GsdiEvaluationResult = GsdiEvalResult;
using GsdiConductance = GsdiJacobian;
using GsdiCurrent = GsdiResidual;
using GsdiNoise = GsdiNoiseSource;

} // namespace gspice

#endif // GSPICE_GSDI_HPP
