#ifndef GSPICE_BJT_HPP
#define GSPICE_BJT_HPP

#include "device.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <string>

namespace gspice {

class Bjt : public Device {
public:
    Bjt(
        const std::string& name,
        int nodeC,
        int nodeB,
        int nodeE,
        int type = 1,
        double is = 1e-16,
        double bf = 100.0,
        double br = 1.0,
        double nf = 1.0,
        double nr = 1.0,
        double area = 1.0)
        : Device(name),
          nodeC_(nodeC),
          nodeB_(nodeB),
          nodeE_(nodeE),
          type_(type >= 0 ? 1 : -1),
          is_(std::max(is * std::max(area, 1e-30), 1e-30)),
          bf_(std::max(bf, 1e-9)),
          br_(std::max(br, 1e-9)),
          nf_(std::max(nf, 1e-9)),
          nr_(std::max(nr, 1e-9)) {}

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
        stampLinearized(J, b, x);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        const auto jac = numericalJacobian(x_dc);
        const std::array<int, 3> nodes = {nodeC_, nodeB_, nodeE_};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                J.add(nodes[r], nodes[c], std::complex<double>{jac[r][c], 0.0});
            }
        }
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc, std::vector<NoiseSource>& sources) const override {
        (void)omega;
        const double q = 1.602176634e-19;
        const auto currents = terminalCurrents(x_dc);
        const double icPsd = 2.0 * q * std::abs(currents[0]);
        const double ibPsd = 2.0 * q * std::abs(currents[1]);
        if (icPsd > 0.0) sources.push_back({name_ + ".collector_shot", nodeC_, nodeE_, icPsd});
        if (ibPsd > 0.0) sources.push_back({name_ + ".base_shot", nodeB_, nodeE_, ibPsd});
    }

private:
    static double nodeVoltage(const VectorReal& x, int node) {
        return node >= 0 ? x[node] : 0.0;
    }

    static double limitedExp(double arg) {
        return std::exp(std::clamp(arg, -80.0, 40.0));
    }

    std::array<double, 3> terminalVoltages(const VectorReal& x) const {
        return {nodeVoltage(x, nodeC_), nodeVoltage(x, nodeB_), nodeVoltage(x, nodeE_)};
    }

    std::array<double, 3> terminalCurrentsFromVoltages(const std::array<double, 3>& v) const {
        const double vc = v[0];
        const double vb = v[1];
        const double ve = v[2];
        const double vbe = type_ * (vb - ve);
        const double vbc = type_ * (vb - vc);
        const double ibe = is_ / bf_ * (limitedExp(vbe / (nf_ * vt_)) - 1.0);
        const double ibc = is_ / br_ * (limitedExp(vbc / (nr_ * vt_)) - 1.0);
        const double itr = is_ * (limitedExp(vbe / (nf_ * vt_)) - limitedExp(vbc / (nr_ * vt_)));

        const double ic = itr - ibc;
        const double ib = ibe + ibc;
        const double ie = -(ic + ib);
        return {type_ * ic, type_ * ib, type_ * ie};
    }

    std::array<double, 3> terminalCurrents(const VectorReal& x) const {
        return terminalCurrentsFromVoltages(terminalVoltages(x));
    }

    std::array<std::array<double, 3>, 3> numericalJacobian(const VectorReal& x) const {
        const auto baseV = terminalVoltages(x);
        std::array<std::array<double, 3>, 3> jac{};
        for (int c = 0; c < 3; ++c) {
            auto plusV = baseV;
            auto minusV = baseV;
            const double dv = std::max(1e-7, std::abs(baseV[c]) * 1e-7);
            plusV[c] += dv;
            minusV[c] -= dv;
            const auto plusI = terminalCurrentsFromVoltages(plusV);
            const auto minusI = terminalCurrentsFromVoltages(minusV);
            for (int r = 0; r < 3; ++r) {
                jac[r][c] = (plusI[r] - minusI[r]) / (2.0 * dv);
            }
        }
        return jac;
    }

    void stampLinearized(SparseMatrixReal& J, VectorReal& b, const VectorReal& x) const {
        const auto currents = terminalCurrents(x);
        const auto jac = numericalJacobian(x);
        const std::array<int, 3> nodes = {nodeC_, nodeB_, nodeE_};
        const auto volts = terminalVoltages(x);
        for (int r = 0; r < 3; ++r) {
            double rhs = -currents[r];
            for (int c = 0; c < 3; ++c) {
                J.add(nodes[r], nodes[c], jac[r][c]);
                rhs += jac[r][c] * volts[c];
            }
            b.add(nodes[r], rhs);
        }
    }

    static constexpr double vt_ = 0.025852;
    int nodeC_;
    int nodeB_;
    int nodeE_;
    int type_;
    double is_;
    double bf_;
    double br_;
    double nf_;
    double nr_;
};

} // namespace gspice

#endif // GSPICE_BJT_HPP
