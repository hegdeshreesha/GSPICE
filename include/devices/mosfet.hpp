#ifndef GSPICE_MOSFET_HPP
#define GSPICE_MOSFET_HPP

#include "device.hpp"
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cctype>
#include <array>
#include <cstring>
#include <cstdint>

namespace gspice {

class Mosfet : public Device {
public:
    /**
     * MOSFET Level 1 (Shichman-Hodges)
     * @param nodeD Drain
     * @param nodeG Gate
     * @param nodeS Source
     * @param nodeB Bulk
     * @param type 1 for NMOS, -1 for PMOS
     */
    Mosfet(const std::string& name, int nodeD, int nodeG, int nodeS, int nodeB, int type = 1,
           double W = 1e-6, double L = 1e-6, double Vth = 0.5, double Kp = 100e-6,
           double lambda = 0.05, double gamma = 0.4, double phi = 0.7)
        : Device(name), nodeD_(nodeD), nodeG_(nodeG), nodeS_(nodeS), nodeB_(nodeB),
          type_(type), W_(W), L_(L), Vth_(Vth), Kp_(Kp),
          lambda_(lambda), gamma_(gamma), phi_(phi) {}

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
        if (primitiveMosTransientCapsEnabled() && request.dynamicResidual) {
            const auto charges = terminalCharges(x);
            for (int row = 0; row < 4; ++row) {
                evaluation.dynamicResidual.push_back({nodes[row], charges[row], 0});
            }
        }
        if (primitiveMosTransientCapsEnabled() && request.dynamicJacobian) {
            const auto jacobian = numericalChargeJacobian(x);
            for (int row = 0; row < 4; ++row) {
                for (int column = 0; column < 4; ++column) {
                    evaluation.dynamicJacobian.push_back(
                        {nodes[row], nodes[column], jacobian[row][column], 0});
                }
            }
        }
        return true;
    }

    bool daeAuditSafe() const override { return true; }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        double Vd = (nodeD_ >= 0) ? x[nodeD_] : 0.0;
        double Vg = (nodeG_ >= 0) ? x[nodeG_] : 0.0;
        double Vs = (nodeS_ >= 0) ? x[nodeS_] : 0.0;
        double Vb = (nodeB_ >= 0) ? x[nodeB_] : 0.0;
        
        stampJunctionDiodes(J, b, x);

        // Transform PMOS into NMOS-like equations.
        double Vgs = type_ * (Vg - Vs);
        double Vds = type_ * (Vd - Vs);
        double Vsb = std::max(0.0, type_ * (Vs - Vb));

        const double sqrtPhi = std::sqrt(std::max(phi_, 1e-12));
        const double VthEff = Vth_ + gamma_ * (std::sqrt(std::max(phi_ + Vsb, 1e-12)) - sqrtPhi);
        
        double Ids = 0.0;
        double gm = 0.0;
        double gds = 0.0;
        double beta = Kp_ * (W_ / std::max(L_, 1e-12));
        double Vov = Vgs - VthEff;

        if (Vds <= 0.0) {
            Ids = 0.0;
            gm = 0.0;
            gds = 1e-12;
        } else if (Vov <= 0.0) {
            const double thermal = 0.02585;
            const double n = 1.5;
            const double Id0 = 1e-12 * (W_ / std::max(L_, 1e-12));
            const double expArg = std::clamp(Vov / (n * thermal), -80.0, 40.0);
            const double expTerm = std::exp(expArg);
            const double vdsTerm = 1.0 - std::exp(-Vds / thermal);
            Ids = Id0 * expTerm * vdsTerm;
            gm = Ids / (n * thermal);
            gds = Id0 * expTerm * std::exp(-Vds / thermal) / thermal;
            gds = std::max(gds, 1e-12);
        } else if (Vds < Vov) {
            // Linear
            Ids = beta * (Vov * Vds - 0.5 * Vds * Vds);
            gm = beta * Vds;
            gds = beta * (Vov - Vds);
        } else {
            // Saturation with channel-length modulation.
            const double Idsat = 0.5 * beta * Vov * Vov;
            Ids = Idsat * (1.0 + lambda_ * Vds);
            gm = beta * Vov * (1.0 + lambda_ * Vds);
            gds = std::max(Idsat * lambda_, 1e-12);
        }

        // Apply device type polarity to current only. The transformed
        // derivatives with respect to original terminal voltages remain
        // positive for both NMOS and PMOS in this drain-to-source stamp.
        Ids *= type_;
        // gds is always positive conductant

        // Stamping
        // Drain Row:  [  gds   gm  - (gm+gds) ] [ Vd Vg Vs ] = [ -Ieq ]
        // Source Row: [ -gds  -gm    (gm+gds) ] [ Vd Vg Vs ] = [  Ieq ]
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
        if (!primitiveMosTransientCapsEnabled()) return;
        const auto caps = parasiticCaps();
        for (size_t i = 0; i < caps.size(); ++i) {
            stampCap(J, b, caps[i], ctx, i);
        }
    }

    void acceptTransientStep(const VectorReal& x, double currentTime, const TransientContext& ctx) override {
        (void)currentTime;
        if (!primitiveMosTransientCapsEnabled() || ctx.timeStep <= 0.0 ||
            !ctx.xHistory || ctx.xHistory->empty()) return;
        const auto caps = parasiticCaps();
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
        if (size != transientStateBytes()) throw std::invalid_argument("MOS transient state size");
        std::size_t offset = 0;
        std::memcpy(destination, capCurrents_.data(), sizeof(double) * capCurrents_.size());
        offset += sizeof(double) * capCurrents_.size();
        for (bool valid : capCurrentValid_) {
            const std::uint8_t byte = valid ? 1u : 0u;
            std::memcpy(destination + offset++, &byte, sizeof(byte));
        }
    }

    void restoreTransientStateBytes(const std::byte* source, std::size_t size) override {
        if (size != transientStateBytes()) throw std::invalid_argument("MOS transient state size");
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
        if (!primitiveMosTransientCapsEnabled()) return 0.0;
        double worst = 0.0;
        for (const auto& cap : parasiticCaps()) {
            const double qCoarse = cap.value * capVoltage(coarse, cap);
            const double qFine = cap.value * capVoltage(fine, cap);
            const double tol = chgtol + reltol * std::max(std::abs(qCoarse), std::abs(qFine));
            worst = std::max(worst, std::abs(qFine - qCoarse) / std::max(tol, 1e-30));
        }
        return worst;
    }

    void limitTransientNewton(const VectorReal& previous, VectorReal& candidate) const override {
        // Limit both bulk junctions. Channel-voltage limiting is deliberately
        // left to compact models because Level-1 has no continuous charge
        // formulation from which to derive a reliable fetlim update.
        limitBodyJunction(previous, candidate, nodeB_, nodeD_);
        limitBodyJunction(previous, candidate, nodeB_, nodeS_);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        const auto ss = evaluateSmallSignal(x_dc);
        const std::complex<double> gds = {ss.gds, 0.0};
        const std::complex<double> gm = {ss.gm, 0.0};
        J.add(nodeD_, nodeD_, gds);
        J.add(nodeD_, nodeG_, gm);
        J.add(nodeD_, nodeS_, -(gm + gds));
        J.add(nodeS_, nodeD_, -gds);
        J.add(nodeS_, nodeG_, -gm);
        J.add(nodeS_, nodeS_, gm + gds);
    }

    void collectNoiseSources(double omega, const VectorReal& x_dc, std::vector<NoiseSource>& sources) const override {
        (void)omega;
        const auto ss = evaluateSmallSignal(x_dc);
        const double k = 1.380649e-23;
        const double T = 298.15;
        const double channelConductance = std::max(ss.gds + (2.0 / 3.0) * std::abs(ss.gm), 0.0);
        const double psd = 4.0 * k * T * channelConductance;
        if (psd > 0.0) {
            sources.push_back({name_ + ".channel", nodeD_, nodeS_, psd});
        }
    }

private:
    struct ParasiticCap {
        int pos = -1;
        int neg = -1;
        double value = 0.0;
    };

    struct SmallSignal {
        double gm = 0.0;
        double gds = 0.0;
    };

    int nodeD_, nodeG_, nodeS_, nodeB_;
    int type_;
    double W_, L_, Vth_, Kp_;
    double lambda_, gamma_, phi_;
    std::array<double, 4> capCurrents_{};
    std::array<bool, 4> capCurrentValid_{};

    static bool primitiveMosTransientCapsEnabled() {
        const char* value = std::getenv("GSPICE_ENABLE_PRIMITIVE_MOS_CAPS");
        if (!value) return true;
        std::string setting(value);
        std::transform(setting.begin(), setting.end(), setting.begin(), [](unsigned char ch) {
            return static_cast<char>(std::tolower(ch));
        });
        return !(setting == "0" || setting == "false" || setting == "no" || setting == "off");
    }

    static double nodeVoltage(const VectorReal& x, int node) {
        return node >= 0 ? x[node] : 0.0;
    }

    std::array<double, 4> terminalVoltages(const VectorReal& x) const {
        return {nodeVoltage(x, nodeD_), nodeVoltage(x, nodeG_),
                nodeVoltage(x, nodeS_), nodeVoltage(x, nodeB_)};
    }

    std::array<double, 4> terminalCurrentsFromVoltages(
        const std::array<double, 4>& voltage) const {
        const double vd = voltage[0];
        const double vg = voltage[1];
        const double vs = voltage[2];
        const double vb = voltage[3];
        const double vgs = type_ * (vg - vs);
        const double vds = type_ * (vd - vs);
        const double vsb = std::max(0.0, type_ * (vs - vb));
        const double sqrtPhi = std::sqrt(std::max(phi_, 1e-12));
        const double threshold = Vth_ + gamma_ *
            (std::sqrt(std::max(phi_ + vsb, 1e-12)) - sqrtPhi);
        const double beta = Kp_ * (W_ / std::max(L_, 1e-12));
        const double overdrive = vgs - threshold;
        double ids = 0.0;
        if (vds <= 0.0) {
            ids = 1e-12 * vds;
        } else if (overdrive <= 0.0) {
            constexpr double thermal = 0.02585;
            constexpr double slope = 1.5;
            const double id0 = 1e-12 * (W_ / std::max(L_, 1e-12));
            ids = id0 * std::exp(std::clamp(overdrive / (slope * thermal), -80.0, 40.0)) *
                (1.0 - std::exp(-vds / thermal));
        } else if (vds < overdrive) {
            ids = beta * (overdrive * vds - 0.5 * vds * vds);
        } else {
            ids = 0.5 * beta * overdrive * overdrive * (1.0 + lambda_ * vds);
        }
        ids *= type_;

        std::array<double, 4> current{ids, 0.0, -ids, 0.0};
        const double areaScale = std::max(W_ * L_ / 1e-12, 1e-3);
        const double saturation = 1e-15 * areaScale;
        const auto addDiode = [&](int anode, int cathode) {
            constexpr double thermal = 0.02585;
            const double raw = voltage[anode] - voltage[cathode];
            const double limited = std::clamp(raw, -2.0, 0.85);
            const double exponential = std::exp(limited / thermal);
            const double conductance = saturation / thermal * exponential + 1e-12;
            double diodeCurrent = saturation * (exponential - 1.0) + 1e-12 * limited;
            diodeCurrent += conductance * (raw - limited);
            current[anode] += diodeCurrent;
            current[cathode] -= diodeCurrent;
        };
        if (type_ > 0) {
            addDiode(3, 0);
            addDiode(3, 2);
        } else {
            addDiode(0, 3);
            addDiode(2, 3);
        }
        return current;
    }

    std::array<double, 4> terminalCurrents(const VectorReal& x) const {
        return terminalCurrentsFromVoltages(terminalVoltages(x));
    }

    std::array<double, 4> terminalChargesFromVoltages(
        const std::array<double, 4>& voltage) const {
        const double vd = voltage[0];
        const double vg = voltage[1];
        const double vs = voltage[2];
        const double vb = voltage[3];
        const double vgs = type_ * (vg - vs);
        const double vds = type_ * (vd - vs);
        const double vsb = std::max(0.0, type_ * (vs - vb));
        const double sqrtPhi = std::sqrt(std::max(phi_, 1e-12));
        const double threshold = Vth_ + gamma_ *
            (std::sqrt(std::max(phi_ + vsb, 1e-12)) - sqrtPhi);
        const double overdrive = vgs - threshold;
        constexpr double smoothing = 1e-3;
        const double positiveOverdrive = 0.5 *
            (overdrive + std::sqrt(overdrive * overdrive + smoothing * smoothing));
        const double coxArea = 5.0e-3 * std::max(W_ * L_, 0.0);
        const double gateChannelCharge = type_ * coxArea * positiveOverdrive;
        const double saturationBlend = std::clamp(
            vds / std::max(positiveOverdrive + smoothing, smoothing), 0.0, 1.0);
        const double drainFraction = 0.5 - saturationBlend / 6.0;

        std::array<double, 4> charge{
            -drainFraction * gateChannelCharge,
            gateChannelCharge,
            -(1.0 - drainFraction) * gateChannelCharge,
            0.0
        };
        const double overlap = 3.0e-10 * std::max(W_, 0.0);
        const double junction = std::max(1.0e-18,
            0.5 * 2.0e-3 * std::max(W_ * L_, 0.0));
        const auto addPairCharge = [&](int positive, int negative, double capacitance) {
            const double q = capacitance * (voltage[positive] - voltage[negative]);
            charge[positive] += q;
            charge[negative] -= q;
        };
        addPairCharge(1, 0, overlap);
        addPairCharge(1, 2, overlap);
        addPairCharge(0, 3, junction);
        addPairCharge(2, 3, junction);
        return charge;
    }

    std::array<double, 4> terminalCharges(const VectorReal& x) const {
        return terminalChargesFromVoltages(terminalVoltages(x));
    }

    std::array<std::array<double, 4>, 4> numericalCurrentJacobian(const VectorReal& x) const {
        const auto base = terminalVoltages(x);
        std::array<std::array<double, 4>, 4> jacobian{};
        for (int column = 0; column < 4; ++column) {
            auto plus = base;
            auto minus = base;
            const double delta = std::max(1e-7, std::abs(base[column]) * 1e-7);
            plus[column] += delta;
            minus[column] -= delta;
            const auto plusCurrent = terminalCurrentsFromVoltages(plus);
            const auto minusCurrent = terminalCurrentsFromVoltages(minus);
            for (int row = 0; row < 4; ++row) {
                jacobian[row][column] =
                    (plusCurrent[row] - minusCurrent[row]) / (2.0 * delta);
            }
        }
        return jacobian;
    }

    std::array<std::array<double, 4>, 4> numericalChargeJacobian(const VectorReal& x) const {
        const auto base = terminalVoltages(x);
        std::array<std::array<double, 4>, 4> jacobian{};
        for (int column = 0; column < 4; ++column) {
            auto plus = base;
            auto minus = base;
            const double delta = std::max(1e-7, std::abs(base[column]) * 1e-7);
            plus[column] += delta;
            minus[column] -= delta;
            const auto plusCharge = terminalChargesFromVoltages(plus);
            const auto minusCharge = terminalChargesFromVoltages(minus);
            for (int row = 0; row < 4; ++row) {
                jacobian[row][column] =
                    (plusCharge[row] - minusCharge[row]) / (2.0 * delta);
            }
        }
        return jacobian;
    }

    void limitBodyJunction(
        const VectorReal& previous,
        VectorReal& candidate,
        int anode,
        int cathode) const {
        constexpr double vt = 0.02585;
        constexpr double saturation = 1e-14;
        const double oldVoltage = nodeVoltage(previous, anode) - nodeVoltage(previous, cathode);
        double newVoltage = nodeVoltage(candidate, anode) - nodeVoltage(candidate, cathode);
        const double vcrit = vt * std::log(vt / (std::sqrt(2.0) * saturation));
        if (newVoltage > vcrit && std::abs(newVoltage - oldVoltage) > 2.0 * vt) {
            if (oldVoltage > 0.0) {
                const double arg = 1.0 + (newVoltage - oldVoltage) / vt;
                newVoltage = arg > 0.0 ? oldVoltage + vt * std::log(arg) : vcrit;
            } else {
                newVoltage = vt * std::log(std::max(newVoltage, vt) / vt);
            }
        }
        const double correction = newVoltage -
            (nodeVoltage(candidate, anode) - nodeVoltage(candidate, cathode));
        if (anode >= 0 && cathode >= 0) {
            candidate[anode] += 0.5 * correction;
            candidate[cathode] -= 0.5 * correction;
        } else if (anode >= 0) {
            candidate[anode] += correction;
        } else if (cathode >= 0) {
            candidate[cathode] -= correction;
        }
    }

    SmallSignal evaluateSmallSignal(const VectorReal& x) const {
        const double Vd = nodeVoltage(x, nodeD_);
        const double Vg = nodeVoltage(x, nodeG_);
        const double Vs = nodeVoltage(x, nodeS_);
        const double Vb = nodeVoltage(x, nodeB_);
        const double Vgs = type_ * (Vg - Vs);
        const double Vds = type_ * (Vd - Vs);
        const double Vsb = std::max(0.0, type_ * (Vs - Vb));
        const double sqrtPhi = std::sqrt(std::max(phi_, 1e-12));
        const double VthEff = Vth_ + gamma_ * (std::sqrt(std::max(phi_ + Vsb, 1e-12)) - sqrtPhi);
        const double beta = Kp_ * (W_ / std::max(L_, 1e-12));
        const double Vov = Vgs - VthEff;

        SmallSignal ss;
        if (Vds <= 0.0) {
            ss.gds = 1e-12;
        } else if (Vov <= 0.0) {
            const double thermal = 0.02585;
            const double n = 1.5;
            const double Id0 = 1e-12 * (W_ / std::max(L_, 1e-12));
            const double expArg = std::clamp(Vov / (n * thermal), -80.0, 40.0);
            const double expTerm = std::exp(expArg);
            const double vdsTerm = 1.0 - std::exp(-Vds / thermal);
            const double ids = Id0 * expTerm * vdsTerm;
            ss.gm = ids / (n * thermal);
            ss.gds = Id0 * expTerm * std::exp(-Vds / thermal) / thermal;
            ss.gds = std::max(ss.gds, 1e-12);
        } else if (Vds < Vov) {
            ss.gm = beta * Vds;
            ss.gds = beta * (Vov - Vds);
        } else {
            const double Idsat = 0.5 * beta * Vov * Vov;
            ss.gm = beta * Vov * (1.0 + lambda_ * Vds);
            ss.gds = std::max(Idsat * lambda_, 1e-12);
        }
        return ss;
    }

    static double capVoltage(const VectorReal& x, const ParasiticCap& cap) {
        return nodeVoltage(x, cap.pos) - nodeVoltage(x, cap.neg);
    }

    void stampCap(
        SparseMatrixReal& J,
        VectorReal& b,
        const ParasiticCap& cap,
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

    std::array<ParasiticCap, 4> parasiticCaps() const {
        const double Cox = 5.0e-3;       // F/m^2, rough thin-oxide order of magnitude
        const double Cov = 3.0e-10;      // F/m, overlap order of magnitude
        const double CjArea = 2.0e-3;    // F/m^2, rough junction area capacitance
        const double Cmin = 1e-18;
        const double gateArea = std::max(W_ * L_, 0.0);
        const double Cchannel = Cox * gateArea;
        const double Cgs = std::max(0.5 * Cchannel + Cov * W_, Cmin);
        const double Cgd = std::max(Cov * W_, Cmin);
        const double Cdb = std::max(0.5 * CjArea * gateArea, Cmin);
        const double Csb = std::max(0.5 * CjArea * gateArea, Cmin);

        return {{{nodeG_, nodeS_, Cgs}, {nodeG_, nodeD_, Cgd},
                 {nodeD_, nodeB_, Cdb}, {nodeS_, nodeB_, Csb}}};
    }

    void stampDiode(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        int anode,
        int cathode,
        double saturationCurrent) const {
        const double Vt = 0.02585;
        double Vd = nodeVoltage(x, anode) - nodeVoltage(x, cathode);
        Vd = std::clamp(Vd, -2.0, 0.85);
        const double expV = std::exp(Vd / Vt);
        double Id = saturationCurrent * (expV - 1.0);
        double gd = (saturationCurrent / Vt) * expV + 1e-12;
        Id += 1e-12 * Vd;
        const double Ieq = Id - gd * Vd;

        J.add(anode, anode, gd);
        J.add(cathode, cathode, gd);
        J.add(anode, cathode, -gd);
        J.add(cathode, anode, -gd);
        b.add(anode, -Ieq);
        b.add(cathode, Ieq);
    }

    void stampJunctionDiodes(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x) const {
        const double areaScale = std::max(W_ * L_ / 1e-12, 1e-3);
        const double Is = 1e-15 * areaScale;
        if (type_ > 0) {
            stampDiode(J, b, x, nodeB_, nodeD_, Is);
            stampDiode(J, b, x, nodeB_, nodeS_, Is);
        } else {
            stampDiode(J, b, x, nodeD_, nodeB_, Is);
            stampDiode(J, b, x, nodeS_, nodeB_, Is);
        }
    }
};

} // namespace gspice

#endif // GSPICE_MOSFET_HPP
