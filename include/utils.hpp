#ifndef GSPICE_UTILS_HPP
#define GSPICE_UTILS_HPP

#include <string>
#include <unordered_map>
#include <algorithm>
#include <cctype>

namespace gspice {

class Utils {
public:
    /**
     * Converts a SPICE-style value string (e.g., "1k", "10u", "5.0") to a double.
     */
    static double parseValue(std::string str) {
        if (str.empty()) return 0.0;

        // Remove trailing non-numeric characters except suffixes
        // (SPICE often allows "1kOhm" -> "1k")
        std::string suffix = "";
        size_t last_digit = str.find_last_of("0123456789.");
        if (last_digit != std::string::npos && last_digit < str.size() - 1) {
            suffix = str.substr(last_digit + 1);
            str = str.substr(0, last_digit + 1);
        }
        // Non-numeric tokens (e.g. "STEP=50p", "UIC") must not reach std::stod:
        // a leading letter means the caller handed us a keyword, not a value.
        if (str.empty() || (!std::isdigit(static_cast<unsigned char>(str[0])) &&
                            str[0] != '.' && str[0] != '+' && str[0] != '-')) {
            return 0.0;
        }

        double value = 0.0;
        try {
            value = std::stod(str);
        } catch (const std::exception&) {
            return 0.0;
        }
        if (suffix.empty()) return value;

        char s = std::tolower(suffix[0]);
        static const std::unordered_map<char, double> multipliers = {
            {'t', 1e12},  // Tera
            {'g', 1e9},   // Giga
            {'x', 1e6},   // Meg (SPICE uses X or MEG)
            {'k', 1e3},   // Kilo
            {'m', 1e-3},  // Milli
            {'u', 1e-6},  // Micro
            {'n', 1e-9},  // Nano
            {'p', 1e-12}, // Pico
            {'f', 1e-15}  // Femto
        };

        // Special case for 'meg'
        std::string low_suffix = suffix;
        std::transform(low_suffix.begin(), low_suffix.end(), low_suffix.begin(), ::tolower);
        if (low_suffix.substr(0, 3) == "meg") return value * 1e6;

        if (multipliers.count(s)) {
            return value * multipliers.at(s);
        }

        return value;
    }
};

} // namespace gspice

#endif // GSPICE_UTILS_HPP
