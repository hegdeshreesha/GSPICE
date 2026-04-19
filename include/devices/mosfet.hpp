#ifndef GSPICE_MOSFET_HPP
#define GSPICE_MOSFET_HPP

#include "device.hpp"
#include <cmath>
#include <string>
#include <algorithm>

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
           double W = 1e-6, double L = 1e-6, double Vth = 0.5, double Kp = 100e-6)
        : Device(name), nodeD_(nodeD), nodeG_(nodeG), nodeS_(nodeS), nodeB_(nodeB),
          type_(type), W_(W), L_(L), Vth_(Vth), Kp_(Kp) {}

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        double Vd = (nodeD_ >= 0) ? x[nodeD_] : 0.0;
        double Vg = (nodeG_ >= 0) ? x[nodeG_] : 0.0;
        double Vs = (nodeS_ >= 0) ? x[nodeS_] : 0.0;
        
        // NMOS orientation
        double Vgs = type_ * (Vg - Vs);
        double Vds = type_ * (Vd - Vs);
        
        double Ids = 0.0;
        double gm = 0.0;
        double gds = 0.0;
        double beta = Kp_ * (W_ / L_);

        if (Vgs < Vth_) {
            // Cutoff
            Ids = 0.0; gm = 0.0; gds = 1e-12; // Gmin
        } else if (Vds < (Vgs - Vth_)) {
            // Linear
            Ids = beta * ((Vgs - Vth_) * Vds - 0.5 * Vds * Vds);
            gm = beta * Vds;
            gds = beta * (Vgs - Vth_ - Vds);
        } else {
            // Saturation
            Ids = 0.5 * beta * (Vgs - Vth_) * (Vgs - Vth_);
            gm = beta * (Vgs - Vth_);
            gds = 1e-12; // Basic Level 1 has no lambda
        }

        // Apply device type polarity
        Ids *= type_;
        gm *= type_;
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
        // Reuse DC logic for small-signal gm/gds
        // (Implementation omitted for brevity, similar to dcStamp with complex types)
    }

private:
    int nodeD_, nodeG_, nodeS_, nodeB_;
    int type_;
    double W_, L_, Vth_, Kp_;
};

} // namespace gspice

#endif // GSPICE_MOSFET_HPP
