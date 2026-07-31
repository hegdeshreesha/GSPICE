#ifndef GSPICE_OSDI_LOADER_HPP
#define GSPICE_OSDI_LOADER_HPP

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include "osdi.h"
#include "osdi_emulator.hpp"
#include "osdi_metadata.hpp"

namespace gspice {

inline void osdiLogMessage(void*, char* message, uint32_t level) {
    std::ostream& stream = (level & LOG_LVL_MASK) >= LOG_LVL_WARN
        ? std::cerr
        : std::cout;
    stream << "OSDI";
    switch (level & LOG_LVL_MASK) {
    case LOG_LVL_DEBUG: stream << "(debug)"; break;
    case LOG_LVL_INFO: stream << "(info)"; break;
    case LOG_LVL_WARN: stream << "(warn)"; break;
    case LOG_LVL_ERR: stream << "(error)"; break;
    case LOG_LVL_FATAL: stream << "(fatal)"; break;
    default: break;
    }
    if ((level & LOG_FMT_ERR) != 0) {
        stream << ": failed to format \"" << (message ? message : "") << "\"\n";
        return;
    }
    stream << ": " << (message ? message : "");
    std::free(message);
}

class OSDILoader {
public:
    OSDILoader(const std::string& libraryPath) {
        library_path_ = libraryPath;
        if (libraryPath == "builtin:mos_level_50") {
            available_models_.push_back(OsdiEmulator::getDescriptor());
            metadata_.emplace_back(available_models_.back());
            loaded_ = true;
            return;
        }

        hModule_ = openLibrary(libraryPath);
        if (!hModule_) {
            error_ = "Failed to load OSDI library: " + libraryPath + " (" + lastLoaderError() + ")";
            std::cerr << error_ << std::endl;
            return;
        }

        auto* major_ptr = reinterpret_cast<uint32_t*>(
            getSymbol(hModule_, "OSDI_VERSION_MAJOR"));
        auto* minor_ptr = reinterpret_cast<uint32_t*>(
            getSymbol(hModule_, "OSDI_VERSION_MINOR"));
        if (!major_ptr || !minor_ptr ||
            (*major_ptr == 0u && *minor_ptr < OSDI_VERSION_MINOR_CURR)) {
            error_ = "Unsupported OSDI interface version in " + libraryPath +
                "; GSPICE requires OSDI 0.4 or newer";
            std::cerr << error_ << std::endl;
            return;
        }
        version_major_ = *major_ptr;
        version_minor_ = *minor_ptr;

        auto* nature_count = reinterpret_cast<uint32_t*>(
            getSymbol(hModule_, "OSDI_NATURES_LEN"));
        auto* discipline_count = reinterpret_cast<uint32_t*>(
            getSymbol(hModule_, "OSDI_DISCIPLINES_LEN"));
        auto* attribute_count = reinterpret_cast<uint32_t*>(
            getSymbol(hModule_, "OSDI_ATTRIBUTES_LEN"));
        const auto* natures = reinterpret_cast<const OsdiNature*>(
            getSymbol(hModule_, "OSDI_NATURES"));
        const auto* disciplines = reinterpret_cast<const OsdiDiscipline*>(
            getSymbol(hModule_, "OSDI_DISCIPLINES"));
        const auto* attributes = reinterpret_cast<const OsdiAttribute*>(
            getSymbol(hModule_, "OSDI_ATTRIBUTES"));
        nature_count_ = nature_count ? *nature_count : 0u;
        discipline_count_ = discipline_count ? *discipline_count : 0u;
        attribute_count_ = attribute_count ? *attribute_count : 0u;

        // OSDI exports a writable function-pointer slot used by generated
        // $display/$warning/$error code.
        if (auto** logger = reinterpret_cast<void**>(
                getSymbol(hModule_, "osdi_log"))) {
            *logger = reinterpret_cast<void*>(&osdiLogMessage);
        }

        // Standard OSDI exports
        auto* num_desc_ptr = reinterpret_cast<uint32_t*>(getSymbol(hModule_, "OSDI_NUM_DESCRIPTORS"));
        auto* descriptor_size_ptr = reinterpret_cast<uint32_t*>(getSymbol(hModule_, "OSDI_DESCRIPTOR_SIZE"));
        auto* descriptors_base = reinterpret_cast<const char*>(getSymbol(hModule_, "OSDI_DESCRIPTORS"));

        if (num_desc_ptr && descriptor_size_ptr && descriptors_base) {
            uint32_t num = *num_desc_ptr;
            uint32_t descriptor_size = *descriptor_size_ptr;
            constexpr std::size_t supported_descriptor_size =
                offsetof(OsdiDescriptor, model_name);
            if (descriptor_size < supported_descriptor_size) {
                error_ = "OSDI descriptor size is unsupported in " + libraryPath +
                    " (module size " + std::to_string(descriptor_size) +
                    ", required prefix " + std::to_string(supported_descriptor_size) + ")";
                std::cerr << error_ << std::endl;
                return;
            }
            for (uint32_t i = 0; i < num; ++i) {
                OsdiDescriptor desc{};
                std::memcpy(
                    &desc,
                    descriptors_base + static_cast<size_t>(i) * descriptor_size,
                    supported_descriptor_size);
                available_models_.push_back(desc);
                metadata_.emplace_back(
                    available_models_.back(),
                    natures, nature_count_,
                    disciplines, discipline_count_,
                    attributes, attribute_count_);
            }
            loaded_ = true;
        } else {
            error_ = "OSDI library does not export OSDI_NUM_DESCRIPTORS/OSDI_DESCRIPTOR_SIZE/OSDI_DESCRIPTORS symbols: " + libraryPath;
            std::cerr << error_ << std::endl;
        }
    }

    ~OSDILoader() {
        if (hModule_) closeLibrary(hModule_);
    }

    const std::vector<OsdiDescriptor>& getAvailableModels() const {
        return available_models_;
    }

    const std::vector<OsdiDescriptorMetadata>& getAvailableMetadata() const {
        return metadata_;
    }

    bool isLoaded() const { return loaded_; }
    const std::string& getError() const { return error_; }
    const std::string& getPath() const { return library_path_; }
    uint32_t getVersionMajor() const { return version_major_; }
    uint32_t getVersionMinor() const { return version_minor_; }
    uint32_t getNatureCount() const { return nature_count_; }
    uint32_t getDisciplineCount() const { return discipline_count_; }
    uint32_t getAttributeCount() const { return attribute_count_; }

private:
#ifdef _WIN32
    using LibraryHandle = HMODULE;
    static LibraryHandle openLibrary(const std::string& path) {
        return LoadLibraryA(path.c_str());
    }
    static void* getSymbol(LibraryHandle handle, const char* name) {
        return reinterpret_cast<void*>(GetProcAddress(handle, name));
    }
    static void closeLibrary(LibraryHandle handle) {
        FreeLibrary(handle);
    }
    static std::string lastLoaderError() {
        return "Windows error " + std::to_string(GetLastError());
    }
#else
    using LibraryHandle = void*;
    static LibraryHandle openLibrary(const std::string& path) {
        return dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    }
    static void* getSymbol(LibraryHandle handle, const char* name) {
        return dlsym(handle, name);
    }
    static void closeLibrary(LibraryHandle handle) {
        dlclose(handle);
    }
    static std::string lastLoaderError() {
        const char* err = dlerror();
        return err ? std::string(err) : std::string("dynamic loader error");
    }
#endif

    LibraryHandle hModule_ = nullptr;
    bool loaded_ = false;
    std::string library_path_;
    std::string error_;
    uint32_t version_major_ = 0;
    uint32_t version_minor_ = 0;
    uint32_t nature_count_ = 0;
    uint32_t discipline_count_ = 0;
    uint32_t attribute_count_ = 0;
    std::vector<OsdiDescriptor> available_models_;
    std::vector<OsdiDescriptorMetadata> metadata_;
};

} // namespace gspice

#endif // GSPICE_OSDI_LOADER_HPP
