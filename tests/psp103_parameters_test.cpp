#include "devices/psp103_parameters.hpp"

#include <cassert>
#include <cmath>

int main() {
    using namespace gspice;

    const Psp103ParameterSet model = Psp103ParameterSet::from({
        {"type", -1.0},
        {"TNOM", 25.0}
    });
    const Psp103PreparedModel prepared = model.prepare({
        {"width", 2e-6},
        {"length", 0.5e-6},
        {"NF", 4.0},
        {"M", 2.0}
    }, 85.0);

    assert(prepared.valid());
    assert(prepared.polarity == -1);
    assert(std::abs(prepared.effective_width_m - 16e-6) < 1e-18);
    assert(std::abs(prepared.effective_length_m - 0.5e-6) < 1e-18);
    assert(std::abs(prepared.width_length_ratio - 32.0) < 1e-12);
    assert(std::abs(prepared.temperature.nominal_c - 25.0) < 1e-12);
    assert(std::abs(prepared.temperature.kelvin - 358.15) < 1e-12);
    assert(prepared.temperature.thermal_voltage > 0.030 &&
           prepared.temperature.thermal_voltage < 0.031);

    const auto valid_model = Psp103ParameterSet::from({
        {"vfbo", 0.2}, {"toxo", 2e-9}, {"epsroxo", 3.9},
        {"ndepo", 1e23}, {"beto", 1e-3}, {"phibo", 0.7}});
    assert(valid_model.validateIntrinsicCore());
    assert(!valid_model.featureRequested({"SWIGATE", "SWGATE"}));

    const auto invalid_model = Psp103ParameterSet::from({{"VFB", 0.2}});
    assert(!invalid_model.validateIntrinsicCore());
    assert(invalid_model.validateIntrinsicCore().missing.size() == 5);
    assert(Psp103ParameterSet::handlesModelParameter("vto"));
    assert(Psp103ParameterSet::handlesModelParameter("beta"));
    assert(Psp103ParameterSet::handlesModelParameter("lamda"));
    assert(Psp103ParameterSet::handlesModelParameter("vfb0"));
    assert(Psp103ParameterSet::handlesModelParameter("tox"));
    assert(Psp103ParameterSet::handlesModelParameter("epsrox"));
    assert(Psp103ParameterSet::handlesModelParameter("u0"));
    assert(!Psp103ParameterSet::handlesModelParameter("xj"));
    return 0;
}
