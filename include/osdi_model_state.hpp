#ifndef GSPICE_OSDI_MODEL_STATE_HPP
#define GSPICE_OSDI_MODEL_STATE_HPP

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace gspice {

// OSDI model data belongs to a .MODEL card and is shared by all instances
// referencing that card. Instance data and transient history remain private to
// each OSDIDevice.
struct OSDISharedModelState {
    std::vector<std::max_align_t> storage;
    std::vector<std::unique_ptr<char[]>> string_params;
    std::string descriptor_name;
    std::size_t model_size = 0;
    double simulation_gmin = 1e-12;
    double simulation_minr = 1e-12;
    double nominal_temperature_c = 27.0;
    double geometry_scale = 1.0;
    bool initialized = false;

    void* data() {
        return storage.empty() ? nullptr : storage.data();
    }

    const void* data() const {
        return storage.empty() ? nullptr : storage.data();
    }
};

} // namespace gspice

#endif
