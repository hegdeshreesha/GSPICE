#ifndef GSPICE_RESISTOR_HPP
#define GSPICE_RESISTOR_HPP

#include "device.hpp"
#include <string>

namespace gspice {

class Resistor : public Device {
public:
    Resistor(const std::string& name, int nodePos, int nodeNeg, double value)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), value_(value) {}

    double getConductance() const { return 1.0 / value_; }
    int getNodePos() const { return nodePos_; }
    int getNodeNeg() const { return nodeNeg_; }

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        const double conductance = getConductance();
        const double voltage =
            (nodePos_ >= 0 ? x[nodePos_] : 0.0) -
            (nodeNeg_ >= 0 ? x[nodeNeg_] : 0.0);
        if (request.staticResidual) {
            const double current = conductance * voltage;
            evaluation.staticResidual.push_back({nodePos_, current});
            evaluation.staticResidual.push_back({nodeNeg_, -current});
        }
        if (request.staticJacobian) {
            evaluation.staticJacobian.push_back({nodePos_, nodePos_, conductance});
            evaluation.staticJacobian.push_back({nodePos_, nodeNeg_, -conductance});
            evaluation.staticJacobian.push_back({nodeNeg_, nodePos_, -conductance});
            evaluation.staticJacobian.push_back({nodeNeg_, nodeNeg_, conductance});
        }
        return true;
    }

    bool daeAuditSafe() const override { return true; }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
        DaeEvaluation evaluation;
        DaeRequest request;
        evaluateDae(x, request, evaluation);
        stampDaeStatic(evaluation, x, J, b);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        DaeEvaluation evaluation;
        DaeRequest request;
        request.analysis = DaeAnalysis::SmallSignal;
        request.staticResidual = false;
        evaluateDae(x_dc, request, evaluation);
        stampDaeSmallSignal(evaluation, omega, J);
    }

    double getNoisePSD(double omega, const VectorReal& x_dc) override {
        double k = 1.380649e-23;
        double T = 298.15;
        double G = getConductance();
        return 4.0 * k * T * G;
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc, std::vector<NoiseSource>& sources) const override {
        (void)omega;
        (void)x_dc;
        double k = 1.380649e-23;
        double T = 298.15;
        double G = getConductance();
        if (G <= 0.0) return;
        sources.push_back({name_ + ".thermal", nodePos_, nodeNeg_, 4.0 * k * T * G});
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        double G = getConductance();
        int K = 2 * n_harms + 1;
        for (int k = 0; k < K; ++k) {
            J.add(nodePos_ * K + k, nodePos_ * K + k, G);
            J.add(nodeNeg_ * K + k, nodeNeg_ * K + k, G);
            J.add(nodePos_ * K + k, nodeNeg_ * K + k, -G);
            J.add(nodeNeg_ * K + k, nodePos_ * K + k, -G);
        }
    }

    void pacStamp(SparseMatrixReal& J, VectorReal& b, double f_in, double f_fund, int n_harms, const VectorReal& x_periodic) override {
        hbStamp(J, b, f_fund, n_harms, x_periodic); // Resistor is freq-independent, same as HB
    }

private:
    int nodePos_;
    int nodeNeg_;
    double value_;
};

} // namespace gspice

#endif // GSPICE_RESISTOR_HPP
