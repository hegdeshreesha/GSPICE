#ifndef GSPICE_MOSFET_HPP
#define GSPICE_MOSFET_HPP

#include "device.hpp"
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cctype>

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

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        double Vd = (nodeD_ >= 0) ? x[nodeD_] : 0.0;
        double Vg = (nodeG_ >= 0) ? x[nodeG_] : 0.0;
        double Vs = (nodeS_ >= 0) ? x[nodeS_] : 0.0;
        double Vb = (nodeB_ >= 0) ? x[nodeB_] : 0.0;
        
        if (primitiveMosTransientCapsEnabled()) {
            stampParasiticCaps(J, b, timeStep, x_hist);
        }
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
    struct SmallSignal {
        double gm = 0.0;
        double gds = 0.0;
    };

    int nodeD_, nodeG_, nodeS_, nodeB_;
    int type_;
    double W_, L_, Vth_, Kp_;
    double lambda_, gamma_, phi_;

    static bool primitiveMosTransientCapsEnabled() {
        const char* value = std::getenv("GSPICE_ENABLE_PRIMITIVE_MOS_CAPS");
        if (!value) return false;
        std::string setting(value);
        std::transform(setting.begin(), setting.end(), setting.begin(), [](unsigned char ch) {
            return static_cast<char>(std::tolower(ch));
        });
        return setting == "1" || setting == "true" || setting == "yes" || setting == "on";
    }

    static double nodeVoltage(const VectorReal& x, int node) {
        return node >= 0 ? x[node] : 0.0;
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

    static void stampCap(
        SparseMatrixReal& J,
        VectorReal& b,
        int nodePos,
        int nodeNeg,
        double value,
        double timeStep,
        const std::vector<VectorReal>& x_hist) {
        if (value <= 0.0 || timeStep <= 0.0 || x_hist.empty()) return;
        const double Geq = value / timeStep;
        const VectorReal& x_prev = x_hist.back();
        const double Vprev = nodeVoltage(x_prev, nodePos) - nodeVoltage(x_prev, nodeNeg);
        const double Ieq = Geq * Vprev;
        J.add(nodePos, nodePos, Geq);
        J.add(nodeNeg, nodeNeg, Geq);
        J.add(nodePos, nodeNeg, -Geq);
        J.add(nodeNeg, nodePos, -Geq);
        b.add(nodePos, Ieq);
        b.add(nodeNeg, -Ieq);
    }

    void stampParasiticCaps(
        SparseMatrixReal& J,
        VectorReal& b,
        double timeStep,
        const std::vector<VectorReal>& x_hist) const {
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

        stampCap(J, b, nodeG_, nodeS_, Cgs, timeStep, x_hist);
        stampCap(J, b, nodeG_, nodeD_, Cgd, timeStep, x_hist);
        stampCap(J, b, nodeD_, nodeB_, Cdb, timeStep, x_hist);
        stampCap(J, b, nodeS_, nodeB_, Csb, timeStep, x_hist);
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
