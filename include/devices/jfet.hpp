#ifndef GSPICE_JFET_HPP
#define GSPICE_JFET_HPP

#include "device.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace gspice {

class Jfet : public Device {
public:
    Jfet(
        const std::string& name,
        int nodeD,
        int nodeG,
        int nodeS,
        int type,
        double beta,
        double vto,
        double lambda,
        double is,
        double n)
        : Device(name),
          nodeD_(nodeD),
          nodeG_(nodeG),
          nodeS_(nodeS),
          type_(type >= 0 ? 1 : -1),
          beta_(std::max(beta, 0.0)),
          vto_(vto),
          lambda_(std::max(lambda, 0.0)),
          is_(std::max(is, 0.0)),
          n_(std::max(n, 1e-6)) {}

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        const std::array<int, 3> nodes = {nodeD_, nodeG_, nodeS_};
        if (request.staticResidual) {
            const auto currents = terminalCurrents(x);
            for (int row = 0; row < 3; ++row) {
                evaluation.staticResidual.push_back({nodes[row], currents[row]});
            }
        }
        if (request.staticJacobian) {
            const auto jacobian = numericalCurrentJacobian(x);
            for (int row = 0; row < 3; ++row) {
                for (int column = 0; column < 3; ++column) {
                    evaluation.staticJacobian.push_back({nodes[row], nodes[column], jacobian[row][column]});
                }
            }
        }
        return true;
    }

    bool daeAuditSafe() const override { return true; }

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
        DaeRequest request;
        DaeEvaluation evaluation;
        evaluateDae(x, request, evaluation);
        stampDaeStatic(evaluation, x, J, b);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        DaeRequest request;
        request.analysis = DaeAnalysis::SmallSignal;
        request.staticResidual = false;
        DaeEvaluation evaluation;
        evaluateDae(x_dc, request, evaluation);
        stampDaeSmallSignal(evaluation, omega, J);
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc, std::vector<NoiseSource>& sources) const override {
        (void)omega;
        const auto jacobian = numericalCurrentJacobian(x_dc);
        const double k = 1.380649e-23;
        const double t = 298.15;
        const double channelConductance = std::max(jacobian[0][0], 0.0);
        const double psd = 4.0 * k * t * channelConductance;
        if (psd > 0.0) sources.push_back({name_ + ".channel", nodeD_, nodeS_, psd});
    }

    bool probeCurrent(const VectorReal& x, double& current, double time = 0.0) const override {
        (void)time;
        current = terminalCurrents(x)[0];
        return true;
    }

private:
    static double nodeVoltage(const VectorReal& x, int node) {
        return node >= 0 ? x[node] : 0.0;
    }

    static double limitedExp(double arg) {
        return std::exp(std::clamp(arg, -80.0, 40.0));
    }

    double gateDiodeCurrent(double voltage) const {
        const double vt = n_ * 0.02585;
        const double limited = std::clamp(voltage, -2.0, 0.85);
        const double ev = limitedExp(limited / vt);
        const double gd = is_ / vt * ev + 1e-12;
        return is_ * (ev - 1.0) + 1e-12 * limited + gd * (voltage - limited);
    }

    double nChannelCurrent(double vgs, double vds) const {
        if (vds < 0.0) return -nChannelCurrent(vgs - vds, -vds);
        const double vgt = vgs - vto_;
        if (vgt <= 0.0) return 1e-12 * vds;
        if (vds < vgt) {
            return beta_ * (2.0 * vgt * vds - vds * vds) * (1.0 + lambda_ * vds);
        }
        return beta_ * vgt * vgt * (1.0 + lambda_ * vds);
    }

    std::array<double, 3> terminalVoltages(const VectorReal& x) const {
        return {nodeVoltage(x, nodeD_), nodeVoltage(x, nodeG_), nodeVoltage(x, nodeS_)};
    }

    std::array<double, 3> terminalCurrentsFromVoltages(const std::array<double, 3>& voltage) const {
        const double vd = type_ * voltage[0];
        const double vg = type_ * voltage[1];
        const double vs = type_ * voltage[2];
        double ids = type_ * nChannelCurrent(vg - vs, vd - vs);

        const double igs = type_ * gateDiodeCurrent(type_ * (voltage[1] - voltage[2]));
        const double igd = type_ * gateDiodeCurrent(type_ * (voltage[1] - voltage[0]));
        return {ids - igd, igs + igd, -ids - igs};
    }

    std::array<double, 3> terminalCurrents(const VectorReal& x) const {
        return terminalCurrentsFromVoltages(terminalVoltages(x));
    }

    std::array<std::array<double, 3>, 3> numericalCurrentJacobian(const VectorReal& x) const {
        const auto base = terminalVoltages(x);
        std::array<std::array<double, 3>, 3> jacobian{};
        for (int column = 0; column < 3; ++column) {
            auto plus = base;
            auto minus = base;
            const double delta = std::max(1e-7, std::abs(base[column]) * 1e-7);
            plus[column] += delta;
            minus[column] -= delta;
            const auto plusCurrent = terminalCurrentsFromVoltages(plus);
            const auto minusCurrent = terminalCurrentsFromVoltages(minus);
            for (int row = 0; row < 3; ++row) {
                jacobian[row][column] = (plusCurrent[row] - minusCurrent[row]) / (2.0 * delta);
            }
        }
        return jacobian;
    }

    int nodeD_;
    int nodeG_;
    int nodeS_;
    int type_;
    double beta_;
    double vto_;
    double lambda_;
    double is_;
    double n_;
};

} // namespace gspice

#endif // GSPICE_JFET_HPP
