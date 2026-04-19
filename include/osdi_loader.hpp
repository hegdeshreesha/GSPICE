#ifndef GSPICE_OSDI_LOADER_HPP
#define GSPICE_OSDI_LOADER_HPP

#include <windows.h>
#include <string>
#include <vector>
#include <iostream>
#include "osdi.h"

namespace gspice {

class OSDILoader {
public:
    OSDILoader(const std::string& libraryPath) {
        hModule_ = LoadLibraryA(libraryPath.c_str());
        if (!hModule_) {
            std::cerr << "Failed to load OSDI library: " << libraryPath << " (Error " << GetLastError() << ")" << std::endl;
            return;
        }

        // Standard OSDI exports
        auto* num_desc_ptr = (uint32_t*)GetProcAddress(hModule_, "OSDI_NUM_DESCRIPTORS");
        auto* descriptors_ptr = (OsdiDescriptor**)GetProcAddress(hModule_, "OSDI_DESCRIPTORS");

        if (num_desc_ptr && descriptors_ptr) {
            uint32_t num = *num_desc_ptr;
            OsdiDescriptor* descs = *descriptors_ptr;
            for (uint32_t i = 0; i < num; ++i) {
                available_models_.push_back(descs[i]);
            }
        }
    }

    ~OSDILoader() {
        if (hModule_) FreeLibrary(hModule_);
    }

    const std::vector<OsdiDescriptor>& getAvailableModels() const {
        return available_models_;
    }

private:
    HMODULE hModule_ = nullptr;
    std::vector<OsdiDescriptor> available_models_;
};

} // namespace gspice

#endif // GSPICE_OSDI_LOADER_HPP
