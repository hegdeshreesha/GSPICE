#ifndef GSPICE_COMPACT_MODEL_HPP
#define GSPICE_COMPACT_MODEL_HPP

#include <algorithm>
#include <cctype>
#include <string>
#include <unordered_set>

namespace gspice {

enum class CompactModelKind {
    PrimitiveMos,
    NativeCompact,
    ExternalCompact,
    Unknown
};

struct CompactModelInfo {
    CompactModelKind kind = CompactModelKind::Unknown;

    bool isPrimitiveMos() const {
        return kind == CompactModelKind::PrimitiveMos;
    }

    bool isNativeCompact() const {
        return kind == CompactModelKind::NativeCompact;
    }

    bool isExternalCompact() const {
        return kind == CompactModelKind::ExternalCompact;
    }
};

class CompactModelRegistry {
public:
    static CompactModelRegistry& instance() {
        static CompactModelRegistry registry;
        return registry;
    }

    CompactModelInfo classify(const std::string& type) const {
        if (isPrimitiveMosModelType(type)) return {CompactModelKind::PrimitiveMos};
        if (isNativeModelType(type)) return {CompactModelKind::NativeCompact};
        if (looksLikeCompactModel(type)) return {CompactModelKind::ExternalCompact};
        return {};
    }

    bool isNativeModelType(std::string type) const {
        normalize(type);
        return native_types_.find(type) != native_types_.end();
    }

    bool isPrimitiveMosModelType(std::string type) const {
        normalize(type);
        return primitive_mos_types_.find(type) != primitive_mos_types_.end();
    }

    bool looksLikeCompactModel(std::string type) const {
        normalize(type);
        return type.find("PSP") != std::string::npos ||
               type.find("BSIM") != std::string::npos ||
               type.find("HICUM") != std::string::npos ||
               type.find("EKV") != std::string::npos;
    }

    bool isPsp103ModelType(std::string type) const {
        normalize(type);
        return type == "PSP" || type == "PSP103" || type == "PSP103VA" ||
               type == "PSP103_VA" || type == "PSP103NQS" ||
               type == "PSPNQS103" || type == "PSPNQS103VA" ||
               type.rfind("PSP103", 0) == 0 ||
               type.rfind("PSPNQS103", 0) == 0;
    }

private:
    CompactModelRegistry() {
        primitive_mos_types_.insert("NMOS");
        primitive_mos_types_.insert("PMOS");
        primitive_mos_types_.insert("N");
        primitive_mos_types_.insert("P");
    }

    static void normalize(std::string& value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
    }

    std::unordered_set<std::string> native_types_;
    std::unordered_set<std::string> primitive_mos_types_;
};

} // namespace gspice

#endif // GSPICE_COMPACT_MODEL_HPP
