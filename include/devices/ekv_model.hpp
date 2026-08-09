#ifndef GSPICE_EKV_MODEL_HPP
#define GSPICE_EKV_MODEL_HPP

#include "device.hpp"
#include <cmath>
#include <string>
#include <algorithm>
#include <array>

namespace gspice {

/**
 * Native C++ EKV 2.6 Compact MOSFET Model
 * Continuous formulation across subthreshold, moderate, and strong inversion.
 */
class EkvMosfet : public Device {
public:
    EkvMosfet(const std::string& name, int nodeD, int nodeG, int nodeS, int nodeB, int type = 1,
              double W = 1e-6, double L = 1e-6, double Vto = 0.5, double Kp = 50e-6,
              double gamma = 0.7, double phi = 0.6, double lambda = 0.02, double theta = 0.0)
        : Device(name), nodeD_(nodeD), nodeG_(nodeG), nodeS_(nodeS), nodeB_(nodeB),
          type_(type), W_(W), L_(L), Vto_(Vto), Kp_(Kp), gamma_(gamma), phi_(phi),
          lambda_(lambda), theta_(theta) {}

    bool daeAuditSafe() const override { return true; }

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        const std::array<int, 4> nodes = {nodeD_, nodeG_, nodeS_, nodeB_};
        if (request.staticResidual) {
            const auto currents = terminalCurrents(x);
            for (int row = 0; row < 4; ++row) {
                evaluation.staticResidual.push_back({nodes[row], currents[row]});
            }
        }
        if (request.staticJacobian) {
            const auto jacobian = numericalCurrentJacobian(x);
            for (int row = 0; row < 4; ++row) {
                for (int column = 0; column < 4; ++column) {
                    evaluation.staticJacobian.push_back(
                        {nodes[row], nodes[column], jacobian[row][column]});
                }
            }
        }
        return true;
    }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x,
                 double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        (void)timeStep; (void)currentTime; (void)x_hist;
        double Vd = (nodeD_ >= 0) ? x[nodeD_] : 0.0;
        double Vg = (nodeG_ >= 0) ? x[nodeG_] : 0.0;
        double Vs = (nodeS_ >= 0) ? x[nodeS_] : 0.0;
        double Vb = (nodeB_ >= 0) ? x[nodeB_] : 0.0;

        double Vgb = type_ * (Vg - Vb);
        double Vdb = type_ * (Vd - Vb);
        double Vsb = type_ * (Vs - Vb);
        double Vt = 0.02585; // thermal voltage at 300K

        // EKV pinchoff voltage VP calculation
        double vgb_vto = Vgb - Vto_;
        double VP = 0.0;
        if (vgb_vto > 0.0) {
            double gamma_term = gamma_ / (2.0 * std::sqrt(std::max(phi_, 1e-6)));
            double n = 1.0 + gamma_term;
            VP = vgb_vto / n;
        }

        // Normalized forward & reverse currents
        double exp_f = std::exp(std::clamp((VP - Vsb) / (2.0 * Vt), -40.0, 40.0));
        double exp_r = std::exp(std::clamp((VP - Vdb) / (2.0 * Vt), -40.0, 40.0));

        double log_f = std::log(1.0 + exp_f);
        double log_r = std::log(1.0 + exp_r);

        double if_val = log_f * log_f;
        double ir_val = log_r * log_r;

        double Is = 2.0 * 1.3 * Kp_ * (W_ / std::max(L_, 1e-12)) * Vt * Vt;
        double Vds = Vd - Vs;
        double Ids_raw = Is * (if_val - ir_val) * (1.0 + lambda_ * std::abs(Vds));
        double Ids = type_ * Ids_raw;

        // Analytical / numerical small-signal conductances
        double gm = Is * (2.0 * log_f / (1.0 + exp_f) * exp_f / (2.0 * Vt)) * type_;
        double gds = std::max(Is * (1.0 + lambda_ * std::abs(Vds)) * 1e-3, 1e-12);

        double Ieq = Ids - gm * (Vg - Vs) - gds * (Vd - Vs);

        J.add(nodeD_, nodeD_, gds);
        J.add(nodeD_, nodeG_, gm);
        J.add(nodeD_, nodeS_, -(gm + gds));

        J.add(nodeS_, nodeD_, -gds);
        J.add(nodeS_, nodeG_, -gm);
        J.add(nodeS_, nodeS_, (gm + gds));

        b.add(nodeD_, -Ieq);
        b.add(nodeS_, Ieq);
    }

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, const TransientContext& ctx) override {
        static const std::vector<VectorReal> empty_history;
        dcStamp(J, b, x, 0.0, ctx.currentTime, empty_history);
    }

private:
    std::array<double, 4> terminalCurrents(const VectorReal& x) const {
        double Vd = (nodeD_ >= 0) ? x[nodeD_] : 0.0;
        double Vg = (nodeG_ >= 0) ? x[nodeG_] : 0.0;
        double Vs = (nodeS_ >= 0) ? x[nodeS_] : 0.0;
        double Vb = (nodeB_ >= 0) ? x[nodeB_] : 0.0;

        double Vgb = type_ * (Vg - Vb);
        double Vdb = type_ * (Vd - Vb);
        double Vsb = type_ * (Vs - Vb);
        double Vt = 0.02585;

        double vgb_vto = Vgb - Vto_;
        double VP = (vgb_vto > 0.0) ? vgb_vto / (1.0 + gamma_ / (2.0 * std::sqrt(std::max(phi_, 1e-6)))) : 0.0;

        double exp_f = std::exp(std::clamp((VP - Vsb) / (2.0 * Vt), -40.0, 40.0));
        double exp_r = std::exp(std::clamp((VP - Vdb) / (2.0 * Vt), -40.0, 40.0));

        double log_f = std::log(1.0 + exp_f);
        double log_r = std::log(1.0 + exp_r);

        double Is = 2.0 * 1.3 * Kp_ * (W_ / std::max(L_, 1e-12)) * Vt * Vt;
        double Ids = type_ * Is * (log_f * log_f - log_r * log_r) * (1.0 + lambda_ * std::abs(Vd - Vs));

        return {Ids, 0.0, -Ids, 0.0};
    }

    std::array<std::array<double, 4>, 4> numericalCurrentJacobian(const VectorReal& x) const {
        std::array<std::array<double, 4>, 4> J{};
        const double eps = 1e-7;
        VectorReal x_pert = x;
        const std::array<int, 4> nodes = {nodeD_, nodeG_, nodeS_, nodeB_};
        const auto baseCurrents = terminalCurrents(x);

        for (int j = 0; j < 4; ++j) {
            int nj = nodes[j];
            if (nj < 0) continue;
            double orig = x_pert[nj];
            x_pert[nj] = orig + eps;
            const auto pertCurrents = terminalCurrents(x_pert);
            x_pert[nj] = orig;
            for (int i = 0; i < 4; ++i) {
                J[i][j] = (pertCurrents[i] - baseCurrents[i]) / eps;
            }
        }
        return J;
    }

    int nodeD_, nodeG_, nodeS_, nodeB_;
    int type_; // 1 for NMOS, -1 for PMOS
    double W_, L_;
    double Vto_, Kp_, gamma_, phi_, lambda_, theta_;
};

} // namespace gspice

#endif // GSPICE_EKV_MODEL_HPP
