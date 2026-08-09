#ifndef GSPICE_PSP103_GSDI_HPP
#define GSPICE_PSP103_GSDI_HPP

#include "gsdi.hpp"

#include <vector>

namespace gspice {

// GSDI model descriptor for the native PSP103.4 evaluator. The device is a
// pure 4-terminal model (drain, gate, source, bulk) with no hidden nodes, so
// the local node space coincides with the terminal space and GsdiCollapseMap::
// standard() collapses nothing (identity plan). Routing PSP103 instances
// through GsdiDaeDeviceAdapter + GdiDevice makes the native evaluator run on
// the same descriptor + collapse pipeline as the GMC-generated models, with
// the DAE stamp layer identical to the direct native path.
inline const GsdiModelDescriptor& psp103GsdiDescriptor() {
    static const GsdiModelDescriptor d = [] {
        GsdiModelDescriptor d;
        d.model_type = "PSP103VA";
        d.version = "1.0";
        d.terminal_count = 4;
        d.nodes = {
            {"d", 0, GsdiNodeRole::Terminal, GsdiNoCollapse},
            {"g", 1, GsdiNodeRole::Terminal, GsdiNoCollapse},
            {"s", 2, GsdiNodeRole::Terminal, GsdiNoCollapse},
            {"b", 3, GsdiNodeRole::Terminal, GsdiNoCollapse},
        };
        d.parameters = {
            {"TYPE", 1.0, "", "channel type (+1 NMOS, -1 PMOS)", true},
            {"VFB", 0.0, "V", "flat-band voltage", true},
            {"VFB0", 0.0, "V", "flat-band voltage alias", true},
            {"PHIB", 0.7, "V", "bulk Fermi potential", true},
            {"PHIBO", 0.7, "V", "bulk Fermi potential alias", true},
            {"TOX", 0.0, "m", "gate oxide thickness", true},
            {"TOXO", 0.0, "m", "gate oxide thickness alias", true},
            {"EPSROX", 3.9, "", "relative gate oxide permittivity", true},
            {"EPSROXO", 3.9, "", "relative gate oxide permittivity alias", true},
            {"BET", 0.0, "A/V^2", "transconductance", true},
            {"BETO", 0.0, "A/V^2", "transconductance alias", true},
            {"U0", 0.0, "cm^2/Vs", "low-field mobility", true},
            {"MOBILITY", 0.0, "cm^2/Vs", "low-field mobility alias", true},
            {"VTO", 0.5, "V", "threshold voltage alias", true},
            {"VT0", 0.5, "V", "threshold voltage alias", true},
            {"W", 1e-6, "m", "instance width", false},
            {"L", 1e-6, "m", "instance length", false},
            {"M", 1.0, "", "instance multiplier", false},
            {"NF", 1.0, "", "instance finger count", false},
        };
        for (int row = 0; row < 4; ++row) {
            for (int column = 0; column < 4; ++column) {
                d.jacobian_pattern.emplace_back(row, column);
            }
        }
        d.supports_op = true;
        d.supports_transient = true;
        d.supports_ac = true;
        d.supports_noise = true;
        return d;
    }();
    return d;
}

} // namespace gspice

#endif // GSPICE_PSP103_GSDI_HPP