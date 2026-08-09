#include "devices/psp103_gsdi.hpp"
#include "devices/psp103_model.hpp"
#include "device.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace {

bool finite(double v) { return std::isfinite(v); }

}  // namespace

int main() {
    const gspice::GsdiParamMap instance = {{"W", 1e-6}, {"L", 1e-6}};
    const auto params = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VFB", -0.85}, {"TOX", 2e-9}, {"PHIB", 0.7}, {"BET", 1e-3}
    });

    gspice::Psp103Mosfet device("M1", 0, 1, 2, 3, params, instance, 27.0);

    // Saturated NMOS bias: strong-inversion channel noise must be nonzero.
    gspice::VectorReal x(4);
    x[0] = 0.5;
    x[1] = 1.2;
    x[2] = 0.0;
    x[3] = 0.0;

    const double omega = 2.0 * 3.141592653589793 * 1.0e3;
    const double psd = device.getNoisePSD(omega, x);
    assert(finite(psd));
    assert(psd > 0.0);

    std::vector<gspice::NoiseSource> sources;
    device.collectNoiseSources(omega, x, sources);
    assert(!sources.empty());
    bool hasChannel = false;
    for (const auto& source : sources) {
        assert(finite(source.currentPsd));
        assert(source.currentPsd >= 0.0);
        if (source.name.find(".drain") != std::string::npos) {
            hasChannel = true;
            assert(source.currentPsd > 0.0);
            assert(source.nodePos >= 0 && source.nodeNeg >= 0 && source.nodePos != source.nodeNeg);
        }
    }
    assert(hasChannel);

    // Off bias: no channel current, so channel noise collapses to shot-only
    // (gate leakage is exponentially suppressed, but the PSD must stay finite).
    gspice::VectorReal off(4);
    off[0] = 0.0;
    off[1] = -0.5;
    off[2] = 0.0;
    off[3] = 0.0;
    assert(finite(device.getNoisePSD(omega, off)));

    // SFL flicker term: doubling the flicker coefficient must raise the
    // low-frequency drain PSD (channel thermal is unchanged between the two
    // devices at identical bias).
    const auto paramsFl = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VFB", -0.85}, {"TOX", 2e-9}, {"PHIB", 0.7}, {"BET", 1e-3},
        {"SFL", 2.0e-6}
    });
    gspice::Psp103Mosfet deviceFl("M2", 0, 1, 2, 3, paramsFl, instance, 27.0);
    const double slow = 2.0 * 3.141592653589793 * 10.0;
    const double psdBase = device.getNoisePSD(slow, x);
    const double psdFl = deviceFl.getNoisePSD(slow, x);
    assert(finite(psdBase));
    assert(finite(psdFl));
    assert(psdFl > psdBase);

    // Foundry PSP103 aliases used by IHP cards must feed the native noise
    // equations rather than only being accepted by the parser.
    const auto paramsAliasNoise = gspice::Psp103ParameterSet::from({
        {"TYPE", 1.0}, {"VFB", -0.85}, {"TOX", 2e-9}, {"PHIB", 0.7}, {"BET", 1e-3},
        {"ALPNOI", 2.0e-6}, {"EFO", 1.0}, {"FNTO", 0.5}
    });
    gspice::Psp103Mosfet aliasNoise("M3", 0, 1, 2, 3, paramsAliasNoise, instance, 27.0);
    const double psdAlias = aliasNoise.getNoisePSD(slow, x);
    assert(finite(psdAlias));
    assert(psdAlias > psdBase);

    return 0;
}
