#ifndef GSPICE_OSDI_EMULATOR_HPP
#define GSPICE_OSDI_EMULATOR_HPP

#include "osdi.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace gspice {

class OsdiEmulator {
public:
    static void evaluate_diode(void* instance_data, OsdiEvaluationData* data) {
        double Vd = data->voltages[0] - data->voltages[1];
        
        // Stabilizing limit
        if (Vd > 0.8) Vd = 0.8;
        if (Vd < -2.0) Vd = -2.0;

        double Is = 1e-14;
        double Vt = 0.026;

        double expV = std::exp(Vd / Vt);
        double Id = Is * (expV - 1.0);
        double gd = (Is / Vt) * expV;

        // Standard SPICE Convention: Current flowing INTO the terminal is positive.
        // Terminal 0 (Anode): Current enters => +Id
        // Terminal 1 (Cathode): Current enters => -Id
        data->currents[0] = Id;
        data->currents[1] = -Id;

        // Jacobian: dI_in / dV
        // d(Id)/dV_a = gd,  d(Id)/dV_c = -gd
        // d(-Id)/dV_a = -gd, d(-Id)/dV_c = gd
        data->jacobian[0*2 + 0] = gd;
        data->jacobian[0*2 + 1] = -gd;
        data->jacobian[1*2 + 0] = -gd;
        data->jacobian[1*2 + 1] = gd;
    }

    static void* create_instance(void* model_data) { return nullptr; }

    static OsdiDescriptor getDescriptor() {
        static const char* term_names[] = {"a", "c"};
        static OsdiDescriptor desc;
        desc.model_name = "lumen_diode_va";
        desc.version_major = 1;
        desc.version_minor = 0;
        desc.metadata.name = "lumen_diode";
        desc.metadata.num_terminals = 2;
        desc.metadata.terminal_names = term_names;
        desc.evaluate = evaluate_diode;
        desc.create_instance = create_instance;
        return desc;
    }
};

} // namespace gspice

#endif // GSPICE_OSDI_EMULATOR_HPP
