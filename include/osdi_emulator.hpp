#ifndef GSPICE_OSDI_EMULATOR_HPP
#define GSPICE_OSDI_EMULATOR_HPP

#include "osdi.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace gspice {

class OsdiEmulator {
public:
    /**
     * Emulates a Level 50 industrial MOSFET (BSIM-like complexity).
     */
    static void evaluate_mos50(void* instance_data, OsdiEvaluationData* data) {
        double Vd = data->voltages[0];
        double Vg = data->voltages[1];
        double Vs = data->voltages[2];
        double Vb = data->voltages[3];

        double Vgs = Vg - Vs;
        double Vds = Vd - Vs;
        double Vbs = Vb - Vs;

        // Complex threshold dependency (Level 50 feature)
        double Vth = 0.5 - 0.1 * Vbs + 0.05 * std::sqrt(std::abs(Vbs));
        
        double Ids = 0.0;
        double gm = 0.0;
        double gds = 0.0;
        double beta = 1e-3;

        // Exponential sub-threshold + Quadratic saturation
        if (Vgs < Vth) {
            Ids = 1e-9 * std::exp((Vgs - Vth) / 0.05); // Subthreshold
            gm = Ids / 0.05;
            gds = 1e-12;
        } else {
            double Vov = Vgs - Vth;
            if (Vds < Vov) {
                Ids = beta * (Vov * Vds - 0.5 * Vds * Vds);
                gm = beta * Vds;
                gds = beta * (Vov - Vds);
            } else {
                Ids = 0.5 * beta * Vov * Vov * (1.0 + 0.02 * Vds); // Velocity saturation
                gm = beta * Vov;
                gds = 0.5 * beta * Vov * Vov * 0.02;
            }
        }

        // OSDI Terminal Assignment: 0=D, 1=G, 2=S, 3=B
        data->currents[0] = Ids;  // Drain current IN
        data->currents[1] = 0.0;  // Gate current
        data->currents[2] = -Ids; // Source current OUT
        data->currents[3] = 0.0;  // Bulk current

        // Jacobian Stamping (High-level derivatives)
        data->jacobian[0*4 + 0] = gds;   // dId/dVd
        data->jacobian[0*4 + 1] = gm;    // dId/dVg
        data->jacobian[0*4 + 2] = -(gm + gds); // dId/dVs
        data->jacobian[2*4 + 0] = -gds;
        data->jacobian[2*4 + 1] = -gm;
        data->jacobian[2*4 + 2] = (gm + gds);
    }

    static void* create_instance(void* model_data) { return nullptr; }

    static OsdiDescriptor getDescriptor() {
        static const char* term_names[] = {"d", "g", "s", "b"};
        static OsdiDescriptor desc;
        desc.model_name = "mos_level_50";
        desc.metadata.num_terminals = 4;
        desc.metadata.terminal_names = term_names;
        desc.evaluate = evaluate_mos50;
        desc.create_instance = create_instance;
        return desc;
    }
};

} // namespace gspice

#endif // GSPICE_OSDI_EMULATOR_HPP
