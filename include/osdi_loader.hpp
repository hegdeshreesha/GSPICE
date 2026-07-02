#ifndef GSPICE_OSDI_LOADER_HPP
#define GSPICE_OSDI_LOADER_HPP

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include "osdi.h"
#include "osdi_emulator.hpp"

namespace gspice {

class OSDILoader {
public:
    OSDILoader(const std::string& libraryPath) {
        library_path_ = libraryPath;
        if (libraryPath == "builtin:mos_level_50") {
            available_models_.push_back(OsdiEmulator::getDescriptor());
            loaded_ = true;
            return;
        }

        hModule_ = LoadLibraryA(libraryPath.c_str());
        if (!hModule_) {
            error_ = "Failed to load OSDI library: " + libraryPath + " (Windows error " + std::to_string(GetLastError()) + ")";
            std::cerr << error_ << std::endl;
            return;
        }

        // Standard OSDI exports
        auto* num_desc_ptr = reinterpret_cast<uint32_t*>(GetProcAddress(hModule_, "OSDI_NUM_DESCRIPTORS"));
        auto* descriptor_size_ptr = reinterpret_cast<uint32_t*>(GetProcAddress(hModule_, "OSDI_DESCRIPTOR_SIZE"));
        auto* descriptors_base = reinterpret_cast<const char*>(GetProcAddress(hModule_, "OSDI_DESCRIPTORS"));

        if (num_desc_ptr && descriptor_size_ptr && descriptors_base) {
            uint32_t num = *num_desc_ptr;
            uint32_t descriptor_size = *descriptor_size_ptr;
            if (descriptor_size == 0 || descriptor_size > sizeof(OsdiDescriptor)) {
                error_ = "OSDI descriptor size is unsupported in " + libraryPath +
                    " (module size " + std::to_string(descriptor_size) +
                    ", GSPICE size " + std::to_string(sizeof(OsdiDescriptor)) + ")";
                std::cerr << error_ << std::endl;
                return;
            }
            for (uint32_t i = 0; i < num; ++i) {
                OsdiDescriptor desc{};
                std::memcpy(&desc, descriptors_base + static_cast<size_t>(i) * descriptor_size, descriptor_size);
                available_models_.push_back(desc);
            }
            loaded_ = true;
        } else {
            error_ = "OSDI library does not export OSDI_NUM_DESCRIPTORS/OSDI_DESCRIPTOR_SIZE/OSDI_DESCRIPTORS symbols: " + libraryPath;
            std::cerr << error_ << std::endl;
        }
    }

    ~OSDILoader() {
        if (hModule_) FreeLibrary(hModule_);
    }

    const std::vector<OsdiDescriptor>& getAvailableModels() const {
        return available_models_;
    }

    bool isLoaded() const { return loaded_; }
    const std::string& getError() const { return error_; }
    const std::string& getPath() const { return library_path_; }

private:
    HMODULE hModule_ = nullptr;
    bool loaded_ = false;
    std::string library_path_;
    std::string error_;
    std::vector<OsdiDescriptor> available_models_;
};

} // namespace gspice

#endif // GSPICE_OSDI_LOADER_HPP
