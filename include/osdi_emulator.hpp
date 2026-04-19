#ifndef GSPICE_OSDI_EMULATOR_HPP
#define GSPICE_OSDI_EMULATOR_HPP

#include "osdi.h"
#include <cmath>
#include <iostream>

namespace gspice {

/**
 * Emulates a compiled OSDI shared library in memory.
 * Used to verify the OSDI bridge when OpenVAF is not present.
 */
class OsdiEmulator {
public:
    static void evaluate_diode(void* instance_data, OsdiEvaluationData* data) {
        double Vd = data->voltages[0] - data->voltages[1];
        double Is = 1e-14;
        double Vt = 0.026;

        double expV = std::exp(Vd / Vt);
        double Id = Is * (expV - 1.0);
        double gd = (Is / Vt) * expV;

        // OSDI Terminal 0 (Anode), Terminal 1 (Cathode)
        data->currents[0] = Id;
        data->currents[1] = -Id;

        // Jacobian: dI/dV
        data->jacobian[0*2 + 0] = gd;  // dI0/dV0
        data->jacobian[0*2 + 1] = -gd; // dI0/dV1
        data->jacobian[1*2 + 0] = -gd; // dI1/dV0
        data->jacobian[1*2 + 1] = gd;  // dI1/dV1
    }

    static void* create_instance(void* model_data) {
        return nullptr; // No instance state needed for this simple diode
    }

    static OsdiDescriptor getDescriptor() {
        static const char* term_names[] = {"a", "c"};
        OsdiDescriptor desc;
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
