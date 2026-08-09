#ifndef GSPICE_GDI_HPP
#define GSPICE_GDI_HPP

/**
 * GSPICE Device Interface (GDI)
 * Lightweight, permissively licensed (Apache 2.0 / MIT) C++ interface for 
 * GMC-compiled compact models and native GSPICE device plugins.
 */

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <cstdint>

#include "gsdi.hpp"

namespace gspice {

struct GdiParameter {
    std::string name;
    double default_value = 0.0;
    std::string units;
    std::string description;
    bool is_model_param = true; // true = .MODEL, false = instance
};

struct GdiNodeInfo {
    std::string name;
    int index = -1;
    bool is_internal = false;
};

struct GdiEvaluationRequest {
    double time = 0.0;
    double time_step = 0.0;
    bool static_residual = true;
    bool static_jacobian = true;
    bool dynamic_residual = false;
    bool dynamic_jacobian = false;
    bool noise = false;
    bool allow_bypass = true;
};

struct GdiNoiseSource {
    int node_pos = -1;
    int node_neg = -1;
    double spectral_density = 0.0;
    std::string name;
};

struct GdiElementContrib {
    int row = -1;
    int col = -1;
    double value = 0.0;
};

struct GdiResidualContrib {
    int node = -1;
    double value = 0.0;
};

struct GdiEvaluationResult {
    std::vector<GdiResidualContrib> static_residual;
    std::vector<GdiElementContrib>  static_jacobian;
    std::vector<GdiResidualContrib> dynamic_residual;
    std::vector<GdiElementContrib>  dynamic_jacobian;
    std::vector<GdiNoiseSource> noise;
    bool bypassed = false;
    bool limiting_applied = false;

    void clear() {
        static_residual.clear();
        static_jacobian.clear();
        dynamic_residual.clear();
        dynamic_jacobian.clear();
        noise.clear();
        bypassed = false;
        limiting_applied = false;
    }
};

/**
 * Abstract interface for GMC-generated model instances.
 */
class GdiInstanceBase {
public:
    virtual ~GdiInstanceBase() = default;

    virtual bool evaluate(
        const double* voltages,
        const GdiEvaluationRequest& request,
        GdiEvaluationResult& result) = 0;

    virtual std::size_t numTerminals() const = 0;
    virtual std::size_t numInternalNodes() const = 0;
    virtual std::size_t numStates() const = 0;
};

class GsdiGdiAdapter final : public GsdiInstance {
public:
    explicit GsdiGdiAdapter(std::unique_ptr<GdiInstanceBase> instance)
        : instance_(std::move(instance)) {}

    bool evaluate(const GsdiEvalRequest& request, GsdiEvalResult& result) override {
        if (!instance_ || !request.solution) return false;
        GdiEvaluationRequest gdi_request;
        gdi_request.time = request.time;
        gdi_request.time_step = request.time_step;
        gdi_request.static_residual = request.residual;
        gdi_request.static_jacobian = request.jacobian;
        gdi_request.dynamic_residual = request.dynamic_residual;
        gdi_request.dynamic_jacobian = request.dynamic_jacobian;
        gdi_request.noise = request.noise;
        gdi_request.allow_bypass = request.allow_bypass;

        GdiEvaluationResult gdi_result;
        if (!instance_->evaluate(request.solution, gdi_request, gdi_result)) return false;
        result.clear();
        for (const auto& term : gdi_result.static_residual) {
            result.static_residual.push_back({term.node, term.value, -1});
        }
        for (const auto& term : gdi_result.static_jacobian) {
            result.static_jacobian.push_back({term.row, term.col, term.value, -1});
        }
        for (const auto& term : gdi_result.dynamic_residual) {
            result.dynamic_residual.push_back({term.node, term.value, -1});
        }
        for (const auto& term : gdi_result.dynamic_jacobian) {
            result.dynamic_jacobian.push_back({term.row, term.col, term.value, -1});
        }
        for (const auto& source : gdi_result.noise) {
            result.noise.push_back({source.node_pos, source.node_neg, source.spectral_density, source.name});
        }
        result.bypassed = gdi_result.bypassed;
        result.limiting_applied = gdi_result.limiting_applied;
        return true;
    }

    std::size_t terminalCount() const override {
        return instance_ ? instance_->numTerminals() : 0;
    }

    std::size_t internalNodeCount() const override {
        return instance_ ? instance_->numInternalNodes() : 0;
    }

    std::size_t stateBytes() const override {
        return instance_ ? instance_->numStates() * sizeof(double) : 0;
    }

private:
    std::unique_ptr<GdiInstanceBase> instance_;
};

/**
 * Abstract factory interface for GMC-generated model cards.
 */
class GdiModuleBase {
public:
    virtual ~GdiModuleBase() = default;

    virtual std::string modelName() const = 0;
    virtual std::vector<GdiParameter> listParameters() const = 0;
    virtual std::unique_ptr<GdiInstanceBase> createInstance(
        const std::unordered_map<std::string, double>& modelParams,
        const std::unordered_map<std::string, double>& instanceParams) = 0;
};

} // namespace gspice

#endif // GSPICE_GDI_HPP
