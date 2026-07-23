#ifndef GSPICE_BJT_HPP
#define GSPICE_BJT_HPP

#include "device.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <string>
#include <cstring>
#include <cstdint>

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
        double area = 1.0,
        double cje = 0.0,
        double cjc = 0.0,
        double tf = 0.0)
        : Device(name),
          nodeC_(nodeC),
          nodeB_(nodeB),
          nodeE_(nodeE),
          type_(type >= 0 ? 1 : -1),
          is_(std::max(is * std::max(area, 1e-30), 1e-30)),
          bf_(std::max(bf, 1e-9)),
          br_(std::max(br, 1e-9)),
          nf_(std::max(nf, 1e-9)),
          nr_(std::max(nr, 1e-9)),
          cje_(std::max(cje * std::max(area, 1e-30), 0.0)),
          cjc_(std::max(cjc * std::max(area, 1e-30), 0.0)),
          tf_(std::max(tf, 0.0)) {}

    bool evaluateDae(
        const VectorReal& x,
        const DaeRequest& request,
        DaeEvaluation& evaluation) override {
        evaluation.clear();
        const std::array<int, 3> nodes = {nodeC_, nodeB_, nodeE_};
        if (request.staticResidual) {
            const auto currents = terminalCurrents(x);
            for (int row = 0; row < 3; ++row) {
                evaluation.staticResidual.push_back({nodes[row], currents[row]});
            }
        }
        if (request.staticJacobian) {
            const auto jacobian = numericalJacobian(x);
            for (int row = 0; row < 3; ++row) {
                for (int column = 0; column < 3; ++column) {
                    evaluation.staticJacobian.push_back(
                        {nodes[row], nodes[column], jacobian[row][column]});
                }
            }
        }
        if (request.dynamicResidual) {
            const auto charges = terminalCharges(x);
            for (int row = 0; row < 3; ++row) {
                evaluation.dynamicResidual.push_back({nodes[row], charges[row], 0});
            }
        }
        if (request.dynamicJacobian) {
            const auto jacobian = numericalChargeJacobian(x);
            for (int row = 0; row < 3; ++row) {
                for (int column = 0; column < 3; ++column) {
                    evaluation.dynamicJacobian.push_back(
                        {nodes[row], nodes[column], jacobian[row][column], 0});
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
        const auto caps = junctionCaps(x_dc);
        for (const auto& cap : caps) {
            const std::complex<double> y{0.0, omega * cap.value};
            J.add(cap.pos, cap.pos, y);
            J.add(cap.neg, cap.neg, y);
            J.add(cap.pos, cap.neg, -y);
            J.add(cap.neg, cap.pos, -y);
        }
    }

    void tranStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, const TransientContext& ctx) override {
        stampLinearized(J, b, x);
        const auto caps = junctionCaps(x);
        for (size_t i = 0; i < caps.size(); ++i) stampCap(J, b, caps[i], ctx, i);
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        (void)currentTime;
        if (ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;
        const auto caps = junctionCaps(x);
        for (size_t i = 0; i < caps.size(); ++i) {
            const auto& cap = caps[i];
            double a0 = ctx.a0;
            double a1 = ctx.a1;
            double a2 = ctx.a2;
            bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
            const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && capCurrentValid_[i];
            if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
                a0 = 1.0 / ctx.timeStep;
                a1 = -a0;
                a2 = 0.0;
                useSecond = false;
            }
            const double vNow = capVoltage(x, cap);
            const double vPrev = capVoltage(ctx.xHistory->back(), cap);
            const double vPrev2 = useSecond
                ? capVoltage((*ctx.xHistory)[ctx.xHistory->size() - 2], cap) : vPrev;
            if (useTrap) {
                capCurrents_[i] = cap.value * (a0 * vNow + a1 * vPrev) - capCurrents_[i];
            } else {
                capCurrents_[i] = cap.value *
                    (a0 * vNow + a1 * vPrev + (useSecond ? a2 * vPrev2 : 0.0));
            }
            capCurrentValid_[i] = true;
        }
    }

    std::size_t transientStateBytes() const override {
        return sizeof(double) * capCurrents_.size() + sizeof(std::uint8_t) * capCurrentValid_.size();
    }

    void saveTransientStateBytes(std::byte* destination, std::size_t size) const override {
        if (size != transientStateBytes()) throw std::invalid_argument("BJT transient state size");
        std::size_t offset = 0;
        std::memcpy(destination, capCurrents_.data(), sizeof(double) * capCurrents_.size());
        offset += sizeof(double) * capCurrents_.size();
        for (bool valid : capCurrentValid_) {
            const std::uint8_t byte = valid ? 1u : 0u;
            std::memcpy(destination + offset++, &byte, sizeof(byte));
        }
    }

    void restoreTransientStateBytes(const std::byte* source, std::size_t size) override {
        if (size != transientStateBytes()) throw std::invalid_argument("BJT transient state size");
        std::size_t offset = 0;
        std::memcpy(capCurrents_.data(), source, sizeof(double) * capCurrents_.size());
        offset += sizeof(double) * capCurrents_.size();
        for (std::size_t i = 0; i < capCurrentValid_.size(); ++i) {
            std::uint8_t byte = 0;
            std::memcpy(&byte, source + offset++, sizeof(byte));
            capCurrentValid_[i] = byte != 0;
        }
    }

    double transientChargeError(
        const VectorReal& coarse,
        const VectorReal& fine,
        double reltol,
        double chgtol) override {
        double worst = 0.0;
        const auto coarseCaps = junctionCaps(coarse);
        const auto fineCaps = junctionCaps(fine);
        for (size_t i = 0; i < coarseCaps.size(); ++i) {
            const double qCoarse = coarseCaps[i].value * capVoltage(coarse, coarseCaps[i]);
            const double qFine = fineCaps[i].value * capVoltage(fine, fineCaps[i]);
            const double tol = chgtol + reltol * std::max(std::abs(qCoarse), std::abs(qFine));
            worst = std::max(worst, std::abs(qFine - qCoarse) / std::max(tol, 1e-30));
        }
        return worst;
    }

    void limitTransientNewton(const VectorReal& previous, VectorReal& candidate) const override {
        limitJunction(previous, candidate, nodeB_, nodeE_, nf_ * vt_);
        limitJunction(previous, candidate, nodeB_, nodeC_, nr_ * vt_);
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
    struct JunctionCap {
        int pos = -1;
        int neg = -1;
        double value = 0.0;
    };

    static double nodeVoltage(const VectorReal& x, int node) {
        return node >= 0 ? x[node] : 0.0;
    }

    static double limitedExp(double arg) {
        return std::exp(std::clamp(arg, -80.0, 40.0));
    }

    void limitJunction(
        const VectorReal& previous,
        VectorReal& candidate,
        int pos,
        int neg,
        double thermal) const {
        const double oldVoltage = type_ * (nodeVoltage(previous, pos) - nodeVoltage(previous, neg));
        double newVoltage = type_ * (nodeVoltage(candidate, pos) - nodeVoltage(candidate, neg));
        const double vt = std::max(thermal, 1e-12);
        const double vcrit = vt * std::log(vt / (std::sqrt(2.0) * std::max(is_, 1e-30)));
        if (newVoltage > vcrit && std::abs(newVoltage - oldVoltage) > 2.0 * vt) {
            if (oldVoltage > 0.0) {
                const double arg = 1.0 + (newVoltage - oldVoltage) / vt;
                newVoltage = arg > 0.0 ? oldVoltage + vt * std::log(arg) : vcrit;
            } else {
                newVoltage = vt * std::log(std::max(newVoltage, vt) / vt);
            }
        }
        const double wanted = type_ * newVoltage;
        const double actual = nodeVoltage(candidate, pos) - nodeVoltage(candidate, neg);
        const double correction = wanted - actual;
        if (pos >= 0 && neg >= 0) {
            candidate[pos] += 0.5 * correction;
            candidate[neg] -= 0.5 * correction;
        } else if (pos >= 0) {
            candidate[pos] += correction;
        } else if (neg >= 0) {
            candidate[neg] -= correction;
        }
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

    std::array<double, 3> terminalChargesFromVoltages(const std::array<double, 3>& v) const {
        const double vc = v[0];
        const double vb = v[1];
        const double ve = v[2];
        const double physicalVbe = vb - ve;
        const double physicalVbc = vb - vc;
        const double normalizedVbe = type_ * physicalVbe;
        const double forwardTransport = type_ * tf_ * is_ *
            (limitedExp(normalizedVbe / (nf_ * vt_)) - 1.0);
        const double qbe = cje_ * physicalVbe + forwardTransport;
        const double qbc = cjc_ * physicalVbc;
        // C, B, E ordering. Pairwise construction guarantees exact terminal
        // charge conservation even when diffusion charge is nonlinear.
        return {-qbc, qbe + qbc, -qbe};
    }

    std::array<double, 3> terminalCharges(const VectorReal& x) const {
        return terminalChargesFromVoltages(terminalVoltages(x));
    }

    std::array<JunctionCap, 2> junctionCaps(const VectorReal& x) const {
        const auto jac = numericalJacobian(x);
        const double gmForward = std::max(std::abs(jac[0][1]), 0.0);
        const double diffusion = tf_ * gmForward;
        return {{{nodeB_, nodeE_, cje_ + diffusion}, {nodeB_, nodeC_, cjc_}}};
    }

    double capVoltage(const VectorReal& x, const JunctionCap& cap) const {
        return nodeVoltage(x, cap.pos) - nodeVoltage(x, cap.neg);
    }

    void stampCap(
        SparseMatrixReal& J,
        VectorReal& b,
        const JunctionCap& cap,
        const TransientContext& ctx,
        size_t index) const {
        if (cap.value <= 0.0 || ctx.timeStep <= 0.0 || !ctx.xHistory || ctx.xHistory->empty()) return;
        double a0 = ctx.a0;
        double a1 = ctx.a1;
        double a2 = ctx.a2;
        bool useSecond = ctx.hasSecondHistory && ctx.xHistory->size() >= 2;
        const bool useTrap = ctx.method == TransientIntegrationMethod::Trapezoidal && capCurrentValid_[index];
        if (ctx.method == TransientIntegrationMethod::Trapezoidal && !useTrap) {
            a0 = 1.0 / ctx.timeStep;
            a1 = -a0;
            a2 = 0.0;
            useSecond = false;
        }
        const double vPrev = capVoltage(ctx.xHistory->back(), cap);
        const double vPrev2 = useSecond
            ? capVoltage((*ctx.xHistory)[ctx.xHistory->size() - 2], cap) : vPrev;
        const double geq = cap.value * a0;
        double ieq = -cap.value * (a1 * vPrev + (useSecond ? a2 * vPrev2 : 0.0));
        if (useTrap) ieq += capCurrents_[index];
        J.add(cap.pos, cap.pos, geq);
        J.add(cap.neg, cap.neg, geq);
        J.add(cap.pos, cap.neg, -geq);
        J.add(cap.neg, cap.pos, -geq);
        b.add(cap.pos, ieq);
        b.add(cap.neg, -ieq);
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

    std::array<std::array<double, 3>, 3> numericalChargeJacobian(const VectorReal& x) const {
        const auto baseV = terminalVoltages(x);
        std::array<std::array<double, 3>, 3> jac{};
        for (int column = 0; column < 3; ++column) {
            auto plusV = baseV;
            auto minusV = baseV;
            const double dv = std::max(1e-7, std::abs(baseV[column]) * 1e-7);
            plusV[column] += dv;
            minusV[column] -= dv;
            const auto plusQ = terminalChargesFromVoltages(plusV);
            const auto minusQ = terminalChargesFromVoltages(minusV);
            for (int row = 0; row < 3; ++row) {
                jac[row][column] = (plusQ[row] - minusQ[row]) / (2.0 * dv);
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
    double cje_;
    double cjc_;
    double tf_;
    std::array<double, 2> capCurrents_{};
    std::array<bool, 2> capCurrentValid_{};
};

} // namespace gspice

#endif // GSPICE_BJT_HPP
