#ifndef GSPICE_COMPACT_MODEL_HPP
#define GSPICE_COMPACT_MODEL_HPP

#include <algorithm>
#include <string>
#include <unordered_set>

namespace gspice {

class CompactModelRegistry {
public:
    static CompactModelRegistry& instance() {
        static CompactModelRegistry registry;
        return registry;
    }

    bool isNativeModelType(std::string type) const {
        normalize(type);
        return native_types_.find(type) != native_types_.end();
    }

    bool looksLikeCompactModel(std::string type) const {
        normalize(type);
        return type.find("PSP") != std::string::npos ||
               type.find("BSIM") != std::string::npos ||
               type.find("HICUM") != std::string::npos ||
               type.find("EKV") != std::string::npos;
    }

private:
    CompactModelRegistry() {
        native_types_.insert("NMOS");
        native_types_.insert("PMOS");
        native_types_.insert("N");
        native_types_.insert("P");
    }

    static void normalize(std::string& value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
    }

    std::unordered_set<std::string> native_types_;
};

} // namespace gspice

#endif // GSPICE_COMPACT_MODEL_HPP
