#ifndef GSPICE_PROVENANCE_HPP
#define GSPICE_PROVENANCE_HPP

#include <string>
#include <sstream>
#include <chrono>
#include <ctime>

namespace gspice {

struct ProvenanceInfo {
    std::string simulatorVersion = "0.1.0-beta";
    std::string gitCommitHash = "2026.07.23-academic-beta";
    std::string buildPlatform = "Windows / MSVC / OpenMP";
    std::string timestamp;
    std::string numericalPolicy = "reltol=1e-3 vntol=1e-6 abstol=1e-12";
    std::string modelHashes;

    std::string toJson() const {
        std::stringstream ss;
        ss << "{\n"
           << "  \"version\": \"" << simulatorVersion << "\",\n"
           << "  \"commit\": \"" << gitCommitHash << "\",\n"
           << "  \"platform\": \"" << buildPlatform << "\",\n"
           << "  \"numerical_policy\": \"" << numericalPolicy << "\"\n"
           << "}\n";
        return ss.str();
    }
};

} // namespace gspice

#endif // GSPICE_PROVENANCE_HPP
