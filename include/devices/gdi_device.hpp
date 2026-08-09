#ifndef GSPICE_GDI_DEVICE_HPP
#define GSPICE_GDI_DEVICE_HPP

#include "device.hpp"
#include "gdi.hpp"
#include "gsdi.hpp"
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>
#include <string>

namespace gspice {

/**
 * Adapter class wrapping GMC-compiled GDI model instances into GSPICE's 
 * high-performance DAE device execution engine.
 */
class GdiDevice : public Device {
public:
    GdiDevice(const std::string& name,
              std::unique_ptr<GdiInstanceBase> gdiInstance,
              const std::vector<int>& nodes,
              std::unique_ptr<Device> legacy = nullptr)
        : Device(name),
          gsdi_(std::make_unique<GsdiGdiAdapter>(std::move(gdiInstance))),
          nodes_(nodes),
          legacy_(std::move(legacy)) {}

    GdiDevice(const std::string& name,
              std::unique_ptr<GsdiInstance> gsdiInstance,
              const std::vector<int>& nodes,
              std::unique_ptr<Device> legacy = nullptr)
        : Device(name), gsdi_(std::move(gsdiInstance)), nodes_(nodes), legacy_(std::move(legacy)) {}

    // Constructor with a collapse plan. nodes_ must list only the columns the
    // simulator allocated (terminals + non-collapsible internal nodes). Local
    // node values and stamps for collapsed nodes are routed through the kept
    // node that carries their column, so a model with internal/collapsible
    // nodes can run against the terminal-only MNA space.
    GdiDevice(const std::string& name,
              std::unique_ptr<GsdiInstance> gsdiInstance,
              const std::vector<int>& nodes,
              GsdiCollapseMap collapse,
              std::unique_ptr<Device> legacy = nullptr)
        : Device(name), gsdi_(std::move(gsdiInstance)), nodes_(nodes),
          collapse_(std::move(collapse)), legacy_(std::move(legacy)) {}

    bool daeAuditSafe() const override { return true; }

    bool prefersDampedAutoTransient() const override { return true; }

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        if (!gsdi_) return false;
        evaluation.clear();

        auto localVoltages = localVoltagesOf(x);

        GsdiEvalRequest gsdiReq;
        gsdiReq.analysis = request.analysis == DaeAnalysis::Transient
            ? GsdiAnalysis::Transient : GsdiAnalysis::OperatingPoint;
        gsdiReq.solution = localVoltages.data();
        gsdiReq.solution_size = localVoltages.size();
        gsdiReq.time = request.time;
        gsdiReq.residual = request.staticResidual;
        gsdiReq.jacobian = request.staticJacobian;
        gsdiReq.dynamic_residual = request.dynamicResidual;
        gsdiReq.dynamic_jacobian = request.dynamicJacobian;
        gsdiReq.noise = false;
        gsdiReq.allow_bypass = request.allowBypass;

        GsdiEvalResult gsdiRes;
        if (!gsdi_->evaluate(gsdiReq, gsdiRes)) {
            return false;
        }

        // Map local contributions to global GSPICE DAE matrix/vectors. Local
        // indices of collapsed nodes resolve to the kept node's column; nodes
        // merged into the reference node map to the reference index (-1).
        // Reference-index terms are kept in the evaluation (matching the
        // emission contract of Capacitor/Diode/BJT/BSIM, whose grounded rows
        // carry the opposite charge for a closed conservation group); the DAE
        // stampers already ignore negative indices.
        if (request.staticResidual) {
            for (const auto& res : gsdiRes.static_residual) {
                evaluation.staticResidual.push_back({globalIndex(res.equation), res.value});
            }
        }
        if (request.staticJacobian) {
            for (const auto& elem : gsdiRes.static_jacobian) {
                evaluation.staticJacobian.push_back(
                    {globalIndex(elem.equation), globalIndex(elem.unknown), elem.value});
            }
        }
        if (request.dynamicResidual) {
            for (const auto& res : gsdiRes.dynamic_residual) {
                evaluation.dynamicResidual.push_back({globalIndex(res.equation), res.value, 0});
            }
        }
        if (request.dynamicJacobian) {
            for (const auto& elem : gsdiRes.dynamic_jacobian) {
                evaluation.dynamicJacobian.push_back(
                    {globalIndex(elem.equation), globalIndex(elem.unknown), elem.value, 0});
            }
        }

        evaluation.bypassed = gsdiRes.bypassed;
        evaluation.limitingApplied = gsdiRes.limiting_applied;
        return true;
    }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x,
                 double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        (void)timeStep; (void)currentTime; (void)x_hist;
        DaeRequest req;
        req.staticResidual = true;
        req.staticJacobian = true;
        DaeEvaluation eval;
        if (evaluateDae(x, req, eval)) {
            for (const auto& jac : eval.staticJacobian) {
                J.add(jac.equation, jac.unknown, jac.value);
            }
            for (const auto& res : eval.staticResidual) {
                b.add(res.equation, -res.value);
            }
        }
    }

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, const TransientContext& ctx) override {
        if (legacy_) {
            legacy_->tranStamp(J, b, x, ctx);
            return;
        }
        DaeRequest request;
        request.analysis = DaeAnalysis::Transient;
        request.time = ctx.currentTime;
        request.staticResidual = true;
        request.staticJacobian = true;
        request.dynamicResidual = ctx.timeStep > 0.0;
        request.dynamicJacobian = ctx.timeStep > 0.0;

        DaeEvaluation current;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) {
            request.dynamicResidual = false;
            request.dynamicJacobian = false;
            if (evaluateDae(x, request, current)) stampDaeStatic(current, x, J, b);
            return;
        }
        if (!evaluateDae(x, request, current)) return;

        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        if (ctx.method == TransientIntegrationMethod::Trapezoidal) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }
        DaeEvaluation previous;
        if (!evaluateDae((*ctx.xHistory)[ctx.xHistory->size() - 1], request, previous)) return;
        DaeHistory history;
        appendScaledDaeResidual(history, previous.dynamicResidual, a1);
        if (useSecond) {
            DaeEvaluation previous2;
            if (!evaluateDae((*ctx.xHistory)[ctx.xHistory->size() - 2], request, previous2)) return;
            appendScaledDaeResidual(history, previous2.dynamicResidual, a2);
        }
        stampDaeTransient(current, x, a0, history, J, b);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        if (legacy_) {
            legacy_->acStamp(J, b, omega, x_dc);
            return;
        }
        (void)b;
        auto voltages = localVoltagesOf(x_dc);
        GsdiEvalRequest request;
        request.analysis = GsdiAnalysis::Ac;
        request.solution = voltages.data();
        request.solution_size = voltages.size();
        request.omega = omega;
        request.residual = false;
        request.jacobian = true;
        request.dynamic_residual = false;
        request.dynamic_jacobian = true;
        GsdiEvalResult result;
        if (!gsdi_ || !gsdi_->evaluate(request, result)) return;
        for (const auto& elem : result.static_jacobian) {
            const int r = globalIndex(elem.equation);
            const int c = globalIndex(elem.unknown);
            if (r >= 0 && c >= 0) J.add(r, c, {elem.value, 0.0});
        }
        for (const auto& elem : result.dynamic_jacobian) {
            const int r = globalIndex(elem.equation);
            const int c = globalIndex(elem.unknown);
            if (r >= 0 && c >= 0) J.add(r, c, {0.0, omega * elem.value});
        }
    }

    double getNoisePSD(double omega, const VectorReal& x_dc) override {
        if (legacy_) return legacy_->getNoisePSD(omega, x_dc);
        auto voltages = localVoltagesOf(x_dc);
        GsdiEvalRequest request;
        request.solution = voltages.data();
        request.solution_size = voltages.size();
        request.omega = omega;
        request.residual = false;
        request.jacobian = false;
        request.noise = true;
        GsdiEvalResult result;
        return gsdi_ && gsdi_->evaluate(request, result)
            ? [&] { double total = 0.0; for (const auto& source : result.noise) total += source.spectral_density; return total; }()
            : 0.0;
    }

    void collectNoiseSources(
        double omega,
        const VectorReal& x_dc,
        std::vector<NoiseSource>& sources) const override {
        if (legacy_) legacy_->collectNoiseSources(omega, x_dc, sources);
        if (!gsdi_) return;
        auto voltages = localVoltagesOf(x_dc);
        GsdiEvalRequest request;
        request.solution = voltages.data();
        request.solution_size = voltages.size();
        request.omega = omega;
        request.residual = false;
        request.jacobian = false;
        request.noise = true;
        GsdiEvalResult result;
        if (!gsdi_->evaluate(request, result)) return;
        for (const auto& source : result.noise) {
            const int pos = globalIndex(source.node_pos);
            const int neg = globalIndex(source.node_neg);
            if (pos >= 0 && neg >= 0) {
                sources.push_back({source.name, pos, neg, source.spectral_density,
                                   UINT32_MAX, 0.0, 0.0});
            }
        }
    }

private:
    // Global MNA column for a local node index. With a collapse plan the
    // column is the kept node carrying the local node's value (partner of the
    // merge chain), or -1 when the node merged into the reference node.
    int globalIndex(int local) const {
        const int source = collapse_ ? collapse_->valueSource(local) : local;
        return (source >= 0 && source < static_cast<int>(nodes_.size())) ? nodes_[source] : -1;
    }

    std::vector<double> localVoltagesOf(const VectorReal& x) const {
        const int count = collapse_ ? collapse_->localCount() : static_cast<int>(nodes_.size());
        std::vector<double> voltages(static_cast<std::size_t>(count), 0.0);
        for (int i = 0; i < count; ++i) {
            const int node = globalIndex(i);
            voltages[i] = (node >= 0 && node < x.getSize()) ? x[node] : 0.0;
        }
        return voltages;
    }

    std::unique_ptr<GsdiInstance> gsdi_;
    std::vector<int> nodes_;
    std::optional<GsdiCollapseMap> collapse_;
    std::unique_ptr<Device> legacy_;
};

} // namespace gspice

#endif // GSPICE_GDI_DEVICE_HPP
