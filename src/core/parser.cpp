#include "parser.hpp"
#include "utils.hpp"
#include "devices/resistor.hpp"
#include "devices/capacitor.hpp"
#include "devices/diode.hpp"
#include "devices/voltage_source.hpp"
#include "devices/inductor.hpp"
#include "devices/mutual_inductor.hpp"
#include "devices/jfet.hpp"
#include "devices/mosvar_capacitor.hpp"
#include "devices/port.hpp"
#include "devices/mosfet.hpp"
#include "devices/gdi_device.hpp"
#include "devices/gmc_mosfet.hpp"
#include "devices/bsim3_model.hpp"
#include "devices/bsim3_parameters.hpp"
#include "devices/bsim4_model.hpp"
#include "devices/bsim4_parameters.hpp"
#include "devices/bsim_translator.hpp"
#include "devices/psp103_parameters.hpp"
#include "devices/psp103_model.hpp"
#include "devices/psp103_gsdi.hpp"
#include "devices/gmc_juncap_express.hpp"
#include "devices/probe.hpp"
#include "devices/current_source.hpp"
#include "devices/multi_port.hpp"
#include "devices/controlled_source.hpp"
#include "devices/bjt.hpp"
#include "devices/behavioral_source.hpp"
#include "expression.hpp"
#include "gsdi_device_adapter.hpp"
#include "gmc.hpp"
#ifdef GSPICE_HAVE_GMC_GENERATED
#include "gmc_hidden.hpp"
#include "gmc_va_probe.hpp"
#endif
#include <iostream>
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <set>
#include <vector>

namespace {

std::string toUpperCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

std::string toLowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string sanitizeIdentifier(std::string s) {
    for (char& c : s) {
        const unsigned char uc = static_cast<unsigned char>(c);
        if (!std::isalnum(uc) && c != '_' && c != '$') c = '_';
    }
    return s.empty() ? "node" : s;
}

std::string joinTokens(const std::vector<std::string>& tokens, size_t startIdx) {
    if (startIdx >= tokens.size()) return "";
    std::string out = tokens[startIdx];
    for (size_t i = startIdx + 1; i < tokens.size(); ++i) {
        out += " ";
        out += tokens[i];
    }
    return out;
}

std::string trimCopy(const std::string& text) {
    const auto start = text.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    const auto end = text.find_last_not_of(" \t\r\n");
    return text.substr(start, end - start + 1);
}

std::string stripQuotes(std::string text) {
    text = trimCopy(text);
    if (text.size() >= 2 && ((text.front() == '"' && text.back() == '"') || (text.front() == '\'' && text.back() == '\''))) {
        return text.substr(1, text.size() - 2);
    }
    return text;
}

std::string readEnvVar(const char* envName) {
#ifdef _MSC_VER
    char* buffer = nullptr;
    size_t length = 0;
    if (_dupenv_s(&buffer, &length, envName) != 0 || !buffer) {
        return "";
    }
    std::string value(buffer);
    std::free(buffer);
    return value;
#else
    const char* value = std::getenv(envName);
    return value ? std::string(value) : "";
#endif
}

bool envFlagEnabled(const char* envName) {
    std::string value = readEnvVar(envName);
    if (value.empty()) return false;
    value = toUpperCopy(value);
    return value == "1" || value == "YES" || value == "TRUE" || value == "ON";
}

bool verboseCompatWarnings() {
    return envFlagEnabled("GSPICE_VERBOSE_COMPAT_WARNINGS");
}

bool tryParseSpiceValue(const std::string& token, double& out);
bool isGroundName(const std::string& name);

bool parseVoltageProbeToken(const std::string& token, std::string& pos, std::string& neg) {
    std::string upper = toUpperCopy(token);
    if (upper.rfind("V(", 0) != 0 || token.back() != ')') return false;
    std::string payload = token.substr(2, token.size() - 3);
    size_t comma = payload.find(',');
    if (comma == std::string::npos) {
        pos = trimCopy(payload);
        neg = "0";
    } else {
        pos = trimCopy(payload.substr(0, comma));
        neg = trimCopy(payload.substr(comma + 1));
    }
    return !pos.empty() && !neg.empty();
}

bool parseInitialConditionToken(const std::string& token, std::string& node, double& value) {
    const auto eq = token.find('=');
    if (eq == std::string::npos || eq == 0 || eq + 1 >= token.size()) return false;
    std::string lhs = trimCopy(token.substr(0, eq));
    std::string rhs = trimCopy(token.substr(eq + 1));
    std::string neg;
    if (toUpperCopy(lhs).rfind("V(", 0) == 0 && lhs.back() == ')') {
        if (!parseVoltageProbeToken(lhs, node, neg)) return false;
        if (!isGroundName(neg)) return false;
    } else {
        node = lhs;
    }
    return !node.empty() && tryParseSpiceValue(stripQuotes(rhs), value);
}

bool parseSaveToken(const std::string& token, gspice::SaveSpec& save) {
    std::string pos;
    std::string neg;
    if (parseVoltageProbeToken(token, pos, neg)) {
        save.kind = "V";
        save.node_pos = pos;
        save.node_neg = neg;
        return true;
    }
    std::string cleaned = stripQuotes(trimCopy(token));
    if (cleaned.empty()) return false;
    if (toUpperCopy(cleaned).rfind("I(", 0) == 0) {
        save.kind = "I";
        save.node_pos = cleaned;
        return true;
    }
    save.kind = "V";
    save.node_pos = cleaned;
    save.node_neg = "0";
    return true;
}

std::string stripInlineComment(const std::string& text) {
    bool inQuote = false;
    char quote = '\0';
    for (size_t i = 0; i < text.size(); ++i) {
        char c = text[i];
        if ((c == '"' || c == '\'') && (i == 0 || text[i - 1] != '\\')) {
            if (!inQuote) {
                inQuote = true;
                quote = c;
            } else if (quote == c) {
                inQuote = false;
            }
        }
        if (!inQuote && (c == ';' || c == '$')) {
            return text.substr(0, i);
        }
    }
    return text;
}

std::string normalizeLine(const std::string& text) {
    return trimCopy(stripInlineComment(text));
}

std::filesystem::path resolveRelativePath(const std::filesystem::path& baseFile, const std::string& pathText) {
    std::filesystem::path p(stripQuotes(pathText));
    if (p.is_relative()) {
        p = baseFile.parent_path() / p;
    }
    return p.lexically_normal();
}

std::vector<std::string> tokenizeSimple(const std::string& line) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    while (ss >> token) tokens.push_back(token);
    return tokens;
}

std::string joinSimple(const std::vector<std::string>& tokens) {
    if (tokens.empty()) return "";
    std::string out = tokens[0];
    for (size_t i = 1; i < tokens.size(); ++i) {
        out += " ";
        out += tokens[i];
    }
    return out;
}

struct PreLine {
    std::string text;
    std::string source;
    int lineNo = 0;
};

struct SubcktDef {
    std::string name;
    std::vector<std::string> pins;
    std::unordered_map<std::string, std::string> params;
    std::vector<PreLine> body;
};

bool tryEvaluateParamExpression(
    const std::string& text,
    const std::unordered_map<std::string, std::string>& params,
    double& out,
    int depth = 0);

bool tryParseSpiceValue(const std::string& token, double& out);

std::pair<std::string, std::string> splitParameterToken(const std::string& token);
bool isStringValuedAssignmentKey(const std::string& key);

std::string remapPrimitiveLine(
    const PreLine& line,
    const std::string& instancePrefix,
    const std::unordered_map<std::string, std::string>& pinMap);

std::vector<PreLine> filterConditionalLines(
    const std::vector<PreLine>& lines,
    std::unordered_map<std::string, std::string>& params,
    std::vector<std::string>& errors,
    const std::string& context);

void replaceAll(std::string& text, const std::string& needle, const std::string& value) {
    if (needle.empty()) return;
    size_t pos = 0;
    while ((pos = text.find(needle, pos)) != std::string::npos) {
        text.replace(pos, needle.size(), value);
        pos += value.size();
    }
}

std::string formatNumericValue(double value) {
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

void addParamAlias(std::unordered_map<std::string, std::string>& params, const std::string& key, const std::string& value) {
    if (key.empty()) return;
    params[key] = value;
    params[toUpperCopy(key)] = value;
    params[toLowerCopy(key)] = value;
}

std::unordered_map<std::string, std::string> parseParameterAssignments(
    const std::vector<std::string>& tokens,
    size_t startIdx) {
    std::unordered_map<std::string, std::string> params;
    for (size_t i = startIdx; i < tokens.size(); ++i) {
        if (i + 2 < tokens.size() && tokens[i + 1] == "=") {
            addParamAlias(params, tokens[i], stripQuotes(tokens[i + 2]));
            i += 2;
            continue;
        }
        if (!tokens[i].empty() && tokens[i].back() == '=' && i + 1 < tokens.size()) {
            addParamAlias(params, tokens[i].substr(0, tokens[i].size() - 1), stripQuotes(tokens[i + 1]));
            ++i;
            continue;
        }
        size_t eq = tokens[i].find('=');
        if (eq == std::string::npos || eq == 0) continue;
        addParamAlias(params, tokens[i].substr(0, eq), stripQuotes(tokens[i].substr(eq + 1)));
    }
    return params;
}

bool lookupParamValue(
    const std::unordered_map<std::string, std::string>& params,
    const std::string& key,
    std::string& value) {
    auto it = params.find(key);
    if (it == params.end()) it = params.find(toUpperCopy(key));
    if (it == params.end()) it = params.find(toLowerCopy(key));
    if (it == params.end()) return false;
    value = it->second;
    return true;
}

std::string stripExpressionDelimiters(std::string text) {
    text = stripQuotes(trimCopy(text));
    if (text.size() >= 2 && text.front() == '{' && text.back() == '}') {
        text = trimCopy(text.substr(1, text.size() - 2));
    }
    text = stripQuotes(trimCopy(text));
    return text;
}

bool isExpressionFunctionName(const std::string& idUpper) {
    static const std::set<std::string> names = {
        "SIN", "COS", "TAN", "EXP", "LOG", "LN", "LOG10", "SQRT", "ABS",
        "POW", "MIN", "MAX", "LIMIT", "CLAMP", "IF", "SGN", "SIGN",
        "U", "STEP", "URAMP", "FLOOR", "CEIL", "CEILING", "ROUND",
        "GAUSS", "AGAUSS"
    };
    return names.count(idUpper) != 0;
}

bool isNumericSuffixStart(const std::string& text, size_t pos) {
    if (pos == 0) return false;
    const char prev = text[pos - 1];
    return std::isdigit(static_cast<unsigned char>(prev)) || prev == '.';
}

bool containsCircuitDependentExpression(const std::string& text) {
    const std::string upper = toUpperCopy(text);
    return upper.find("V(") != std::string::npos ||
           upper.find("I(") != std::string::npos ||
           upper.find("TIME") != std::string::npos;
}

std::string substituteParamIdentifiers(
    const std::string& expression,
    const std::unordered_map<std::string, std::string>& params,
    int depth) {
    std::string out;
    for (size_t i = 0; i < expression.size();) {
        const char c = expression[i];
        const bool identStart = std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '$';
        if (!identStart || isNumericSuffixStart(expression, i)) {
            out += c;
            ++i;
            continue;
        }

        const size_t start = i;
        while (i < expression.size()) {
            const char idc = expression[i];
            if (!(std::isalnum(static_cast<unsigned char>(idc)) || idc == '_' || idc == '$' || idc == ':' || idc == '.')) break;
            ++i;
        }
        const std::string id = expression.substr(start, i - start);
        size_t next = i;
        while (next < expression.size() && std::isspace(static_cast<unsigned char>(expression[next]))) ++next;
        const bool looksLikeFunction = next < expression.size() && expression[next] == '(' && isExpressionFunctionName(toUpperCopy(id));
        if (looksLikeFunction || toUpperCopy(id) == "PI" || toUpperCopy(id) == "E" || toUpperCopy(id) == "TIME" || toUpperCopy(id) == "T") {
            out += id;
            continue;
        }

        std::string rawValue;
        if (!lookupParamValue(params, id, rawValue)) {
            out += id;
            continue;
        }
        double evaluated = 0.0;
        if (tryEvaluateParamExpression(rawValue, params, evaluated, depth + 1)) {
            out += formatNumericValue(evaluated);
        } else {
            out += "(" + stripExpressionDelimiters(rawValue) + ")";
        }
    }
    return out;
}

bool tryEvaluateParamExpression(
    const std::string& text,
    const std::unordered_map<std::string, std::string>& params,
    double& out,
    int depth) {
    if (depth > 24) return false;
    const std::string stripped = stripExpressionDelimiters(text);
    if (stripped.empty()) return false;
    if (containsCircuitDependentExpression(stripped)) return false;
    if (tryParseSpiceValue(stripped, out)) return true;
    const std::string substituted = substituteParamIdentifiers(stripped, params, depth);
    if (tryParseSpiceValue(substituted, out)) return true;
    try {
        gspice::BehavioralExpression expr(substituted, [](const std::string&) { return -1; });
        gspice::VectorReal empty(0);
        out = expr.evaluate(empty, 0.0);
        return std::isfinite(out);
    } catch (const std::exception&) {
        return false;
    }
}

std::string resolvedParamString(
    const std::string& text,
    const std::unordered_map<std::string, std::string>& params) {
    double value = 0.0;
    if (tryEvaluateParamExpression(text, params, value)) return formatNumericValue(value);
    return stripExpressionDelimiters(text);
}

std::string resolveBracedNumericExpressions(
    const std::string& line,
    const std::unordered_map<std::string, std::string>& params) {
    std::string out;
    for (size_t i = 0; i < line.size();) {
        if (line[i] != '{') {
            out += line[i++];
            continue;
        }
        int depth = 0;
        size_t j = i;
        for (; j < line.size(); ++j) {
            if (line[j] == '{') ++depth;
            if (line[j] == '}') {
                --depth;
                if (depth == 0) break;
            }
        }
        if (j >= line.size()) {
            out += line.substr(i);
            break;
        }
        const std::string payload = line.substr(i + 1, j - i - 1);
        double value = 0.0;
        if (tryEvaluateParamExpression(payload, params, value)) {
            out += formatNumericValue(value);
        } else {
            out += line.substr(i, j - i + 1);
        }
        i = j + 1;
    }
    return out;
}

void resolveParameterMapExpressions(std::unordered_map<std::string, std::string>& params) {
    for (int pass = 0; pass < 8; ++pass) {
        bool changed = false;
        for (auto& item : params) {
            double value = 0.0;
            if (!tryEvaluateParamExpression(item.second, params, value)) continue;
            const std::string formatted = formatNumericValue(value);
            if (item.second != formatted) {
                item.second = formatted;
                changed = true;
            }
        }
        if (!changed) break;
    }
}

std::unordered_map<std::string, std::string> collectGlobalParams(
    const std::vector<PreLine>& input,
    std::vector<PreLine>& nonParamLines) {
    std::unordered_map<std::string, std::string> params;
    for (const auto& line : input) {
        auto tokens = tokenizeSimple(line.text);
        if (tokens.empty()) continue;
        std::string cmd = toUpperCopy(tokens[0]);
        if (cmd == ".PARAM" || cmd == ".PARAMS") {
            auto parsed = parseParameterAssignments(tokens, 1);
            for (const auto& [key, value] : parsed) {
                params[key] = value;
            }
        } else {
            nonParamLines.push_back(line);
        }
    }
    resolveParameterMapExpressions(params);
    return params;
}

std::string applyGlobalParams(
    std::string line,
    const std::unordered_map<std::string, std::string>& params) {
    for (const auto& [key, value] : params) {
        const std::string resolved = resolvedParamString(value, params);
        replaceAll(line, "{" + key + "}", resolved);
        replaceAll(line, "'" + key + "'", resolved);
    }
    line = resolveBracedNumericExpressions(line, params);
    auto tokens = tokenizeSimple(line);
    const bool isOptionsLine =
        !tokens.empty() && (toUpperCopy(tokens[0]) == ".OPTIONS" || toUpperCopy(tokens[0]) == ".OPTION" ||
                            toUpperCopy(tokens[0]) == ".OPT");
    if (!isOptionsLine) {
        for (size_t i = 0; i < tokens.size(); ++i) {
            auto& token = tokens[i];
            if (i + 2 < tokens.size() && tokens[i + 1] == "=") {
                double evaluated = 0.0;
                if (tryEvaluateParamExpression(tokens[i + 2], params, evaluated)) {
                    tokens[i + 2] = formatNumericValue(evaluated);
                }
                continue;
            }
            auto [key, value] = splitParameterToken(token);
            if (!key.empty()) {
                if (isStringValuedAssignmentKey(key)) continue;
                double evaluated = 0.0;
                if (tryEvaluateParamExpression(value, params, evaluated)) {
                    token = key + "=" + formatNumericValue(evaluated);
                }
                continue;
            }
            auto it = params.find(token);
            if (it != params.end()) token = resolvedParamString(it->second, params);
        }
    }
    return joinSimple(tokens);
}

bool tryParseSpiceValue(const std::string& token, double& out) {
    const std::string text = stripExpressionDelimiters(token);
    if (text.empty()) return false;
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (c == '*' || c == '/' || c == '^' || c == '(' || c == ')' ||
            c == '{' || c == '}' || c == '<' || c == '>' || c == '=' ||
            c == '?' || c == ':' || c == ',') {
            return false;
        }
        if ((c == '+' || c == '-') && i != 0) {
            const char prev = text[i - 1];
            if (prev != 'e' && prev != 'E') return false;
        }
    }
    try {
        out = gspice::Utils::parseValue(text);
        return true;
    } catch (...) {
        return false;
    }
}

bool isGroundName(const std::string& name) {
    std::string up = toUpperCopy(name);
    return up == "0" || up == "GND";
}

bool isParameterToken(const std::string& token) {
    return token.find('=') != std::string::npos;
}

bool isStringValuedAssignmentKey(const std::string& key) {
    const std::string upper = toUpperCopy(key);
    return upper == "MODE" || upper == "NOISEMODE" || upper == "NOISE_MODE" ||
           upper == "SOLVER" || upper == "METHOD" || upper == "FORMAT" ||
           upper == "PHASENOISE" || upper == "PHASE_NOISE" || upper == "JITTER";
}

bool parseSpiceBool(const std::string& value, bool fallback = true) {
    const std::string upper = toUpperCopy(stripQuotes(value));
    if (upper.empty()) return fallback;
    if (upper == "1" || upper == "YES" || upper == "TRUE" || upper == "ON") return true;
    if (upper == "0" || upper == "NO" || upper == "FALSE" || upper == "OFF") return false;
    return fallback;
}

bool isRfSweepKeyword(const std::string& token) {
    const std::string upper = toUpperCopy(token);
    return upper == "DEC" || upper == "OCT" || upper == "LIN" || upper == "VALUES";
}

void applyRfAnalysisOption(gspice::SimulationSettings& settings, const std::string& key, const std::string& value) {
    const std::string upper = toUpperCopy(key);
    if (value.empty()) return;
    if (upper == "FUND" || upper == "FUNDAMENTAL" || upper == "F0") {
        settings.f_fund.clear();
        settings.f_fund.push_back(gspice::Utils::parseValue(value));
    } else if (upper == "SIDEBANDS" || upper == "SIDEBAND" ||
               upper == "NHARMS" || upper == "N_HARMS" ||
               upper == "HARMS" || upper == "HARMONICS") {
        settings.n_harms = std::max(1, std::stoi(value));
    } else if (upper == "PHASENOISE" || upper == "PHASE_NOISE") {
        settings.pnoise_phase_noise = parseSpiceBool(value);
    } else if (upper == "JITTER") {
        settings.pnoise_jitter = parseSpiceBool(value);
    } else if (upper == "CARRIER" || upper == "CARRIERFREQ") {
        settings.pnoise_carrier = gspice::Utils::parseValue(value);
    } else if (upper == "NATIVE_HB_REQUIRED" || upper == "HB_NATIVE_REQUIRED" ||
               upper == "STRICT_HB" || upper == "SIGNOFF") {
        settings.hb_native_required = parseSpiceBool(value);
    }
}

void parseRfAssignmentOptions(gspice::SimulationSettings& settings,
                              const std::vector<std::string>& tokens,
                              size_t start) {
    for (size_t i = start; i < tokens.size(); ++i) {
        auto [key, value] = splitParameterToken(tokens[i]);
        if (!key.empty()) applyRfAnalysisOption(settings, key, value);
    }
}

std::pair<std::string, std::string> splitParameterToken(const std::string& token) {
    size_t eq = token.find('=');
    if (eq == std::string::npos || eq == 0) return {"", ""};
    return {token.substr(0, eq), stripQuotes(token.substr(eq + 1))};
}

void addParam(std::unordered_map<std::string, std::string>& params, const std::string& key, const std::string& value) {
    if (key.empty()) return;
    params[key] = value;
    params[toUpperCopy(key)] = value;
    params[toLowerCopy(key)] = value;
}

std::unordered_map<std::string, std::string> parseParameterTokens(
    const std::vector<std::string>& tokens,
    size_t startIdx) {
    return parseParameterAssignments(tokens, startIdx);
}

std::unordered_map<std::string, std::string> parseModelParamsFromText(std::string text) {
    for (char& c : text) {
        if (c == '(' || c == ')' || c == ',') c = ' ';
    }
    auto paramTokens = tokenizeSimple(text);
    return parseParameterTokens(paramTokens, 0);
}

double paramValue(
    const std::unordered_map<std::string, std::string>& params,
    const std::vector<std::string>& keys,
    double defaultValue) {
    for (const auto& key : keys) {
        auto it = params.find(key);
        if (it == params.end()) it = params.find(toUpperCopy(key));
        if (it == params.end()) it = params.find(toLowerCopy(key));
        if (it == params.end()) continue;
        double parsed = 0.0;
        if (tryEvaluateParamExpression(it->second, params, parsed)) return parsed;
    }
    return defaultValue;
}

bool parsePrimitiveValue(
    const std::vector<std::string>& tokens,
    size_t valueIdx,
    const std::vector<std::string>& paramKeys,
    double& out) {
const auto params = parseParameterTokens(tokens, valueIdx);
    const double keyed = paramValue(params, paramKeys, std::numeric_limits<double>::quiet_NaN());
    if (std::isfinite(keyed)) {
        out = keyed;
        return true;
    }
    auto [key, value] = splitParameterToken(tokens[valueIdx]);
    return tryEvaluateParamExpression(key.empty() ? tokens[valueIdx] : value, params, out);
}

void ensurePrimitiveMosGmcModels() {
    static const bool registered = [] {
        auto factory = [](const gspice::GmcModelDefinition& definition) -> std::unique_ptr<gspice::Device> {
            const auto value = [&](const std::vector<std::string>& keys, double fallback) {
                for (const auto& key : keys) {
                    auto it = definition.model_params.find(key);
                    if (it != definition.model_params.end()) {
                        double parsed = 0.0;
                        if (tryParseSpiceValue(it->second, parsed)) return parsed;
                    }
                }
                return fallback;
            };
            const auto instanceValue = [&](const std::vector<std::string>& keys, double fallback) {
                for (const auto& key : keys) {
                    auto it = definition.instance_params.find(key);
                    if (it != definition.instance_params.end()) {
                        double parsed = 0.0;
                        if (tryParseSpiceValue(it->second, parsed)) return parsed;
                    }
                }
                return fallback;
            };
            if (definition.nodes.size() != 4) return nullptr;
            const int type = (toUpperCopy(definition.type) == "PMOS" ||
                              toUpperCopy(definition.type) == "P") ? -1 : 1;
            auto legacy = std::make_unique<gspice::Mosfet>(
                definition.name,
                definition.nodes[0], definition.nodes[1],
                definition.nodes[2], definition.nodes[3],
                type,
                instanceValue({"W", "w"}, 1e-6),
                instanceValue({"L", "l"}, 1e-6),
                value({"VTO", "VT0", "VTH", "VTH0"}, 0.5),
                value({"KP", "BETA", "K"}, 100e-6),
                value({"LAMBDA", "LAMDA"}, 0.05),
                value({"GAMMA"}, 0.4),
                value({"PHI"}, 0.7));
            auto model = std::make_unique<gspice::GmcMosfetInstance>(
                type,
                instanceValue({"W", "w"}, 1e-6),
                instanceValue({"L", "l"}, 1e-6),
                value({"VTO", "VT0", "VTH", "VTH0"}, 0.5),
                value({"KP", "BETA", "K"}, 100e-6),
                value({"LAMBDA", "LAMDA"}, 0.05),
                value({"GAMMA"}, 0.4),
                value({"PHI"}, 0.7));
            return std::make_unique<gspice::GdiDevice>(
                definition.name, std::move(model), definition.nodes, std::move(legacy));
        };
        auto& registry = gspice::GmcRegistry::instance();
        for (const char* type : {"NMOS", "PMOS", "N", "P"}) {
            registry.registerModel(type, factory);
        }
        return true;
    }();
    (void)registered;
}

void ensureBsim3GmcModels() {
    static const bool registered = [] {
        auto factory = [](const gspice::GmcModelDefinition& definition) -> std::unique_ptr<gspice::Device> {
            if (definition.nodes.size() != 4) return nullptr;
            std::unordered_map<std::string, double> raw;
            for (const auto& [name, text] : definition.model_params) {
                double value = 0.0;
                if (!tryParseSpiceValue(text, value)) return nullptr;
                raw[name] = value;
            }
            const auto translated = gspice::BsimModelTranslator::translate(raw);
            if (!translated || translated.family != gspice::BsimFamily::Bsim3) return nullptr;
            const auto instanceValue = [&](const char* key, double fallback) {
                const auto it = definition.instance_params.find(key);
                if (it == definition.instance_params.end()) return fallback;
                double value = 0.0;
                return tryParseSpiceValue(it->second, value) ? value : fallback;
            };
            const double width = instanceValue("W", 1.0e-6);
            const double length = instanceValue("L", 1.0e-6);
            const auto prepared = gspice::Bsim3ParameterSet::from(translated.parameters).prepare(
                width, length, definition.temperature_c);
            if (!prepared.validation) return nullptr;
            const std::string typeName = toUpperCopy(definition.type);
            const int type = (typeName == "BSIM3_PMOS" || typeName == "PMOS") ? -1 : 1;
            auto device = std::make_unique<gspice::Bsim3Mosfet>(
                definition.name, 0, 1, 2, 3, type,
                prepared.width, prepared.length, prepared.vth0, prepared.kp,
                prepared.u0, prepared.vsat, prepared.k1, prepared.nfactor,
                prepared.toxe, prepared.kf, prepared.af, prepared.ef);
            return std::make_unique<gspice::GdiDevice>(
                definition.name,
                std::make_unique<gspice::GsdiDaeDeviceAdapter>(std::move(device), 4),
                definition.nodes);
        };
        auto& registry = gspice::GmcRegistry::instance();
        registry.registerModel("BSIM3_NMOS", factory);
        registry.registerModel("BSIM3_PMOS", factory);
        return true;
    }();
    (void)registered;
}

void ensureBsim4GmcModels() {
    static const bool registered = [] {
        auto factory = [](const gspice::GmcModelDefinition& definition) -> std::unique_ptr<gspice::Device> {
            if (definition.nodes.size() != 4) return nullptr;
            std::unordered_map<std::string, double> raw;
            for (const auto& [name, text] : definition.model_params) {
                double value = 0.0;
                if (!tryParseSpiceValue(text, value)) return nullptr;
                raw[name] = value;
            }
            const auto instanceValue = [&](const char* key, double fallback) {
                const auto it = definition.instance_params.find(key);
                if (it == definition.instance_params.end()) return fallback;
                double value = 0.0;
                return tryParseSpiceValue(it->second, value) ? value : fallback;
            };
            const auto prepared = gspice::Bsim4ParameterSet::from(raw).prepare(
                instanceValue("W", 1.0e-6), instanceValue("L", 1.0e-6),
                definition.temperature_c);
            if (!prepared.validation) return nullptr;
            const std::string typeName = toUpperCopy(definition.type);
            const int type = (typeName == "BSIM4_PMOS" || typeName == "PMOS") ? -1 : 1;
            auto device = std::make_unique<gspice::Bsim4Mosfet>(
                definition.name, 0, 1, 2, 3, prepared, type);
            return std::make_unique<gspice::GdiDevice>(
                definition.name,
                std::make_unique<gspice::GsdiDaeDeviceAdapter>(std::move(device), 4),
                definition.nodes);
        };
        auto& registry = gspice::GmcRegistry::instance();
        registry.registerModel("BSIM4_NMOS", factory);
        registry.registerModel("BSIM4_PMOS", factory);
        return true;
    }();
    (void)registered;
}

void ensureGmcGeneratedModels() {
#ifdef GSPICE_HAVE_GMC_GENERATED
    static const bool registered = [] {
        auto& registry = gspice::GmcRegistry::instance();
        registry.registerModel(
            "VA_PROBE",
            [](const gspice::GmcModelDefinition& definition)
                -> std::unique_ptr<gspice::Device> {
                gspice::GmcVaProbeModel model;
                const auto& descriptor = model.descriptor();
                if (definition.nodes.size() != static_cast<std::size_t>(descriptor.terminal_count)) {
                    return nullptr;
                }
                gspice::GsdiModelCard card;
                card.name = definition.name;
                card.type = definition.type;
                for (const auto& [name, value] : definition.model_params) {
                    double parsed = 0.0;
                    if (tryParseSpiceValue(value, parsed)) card.parameters[name] = parsed;
                }
                for (const auto& [name, value] : definition.instance_params) {
                    double parsed = 0.0;
                    if (tryParseSpiceValue(value, parsed)) card.parameters[name] = parsed;
                }
                return std::make_unique<gspice::GdiDevice>(
                    definition.name, model.createInstance(card, definition.nodes),
                    definition.nodes);
            });

        // Only nodes declared with role Internal are genuine hidden unknowns
        // needing their own MNA column; Collapsible nodes merge into a terminal.
        std::size_t internal_node_count = 0;
        for (const auto& node : gspice::GmcPspLikeModel{}.descriptor().nodes) {
            if (node.role == gspice::GsdiNodeRole::Internal) ++internal_node_count;
        }
        registry.registerModel(
            "PSP_LIKE",
            [internal_node_count](const gspice::GmcModelDefinition& definition)
                -> std::unique_ptr<gspice::Device> {
                gspice::GmcPspLikeModel model;
                const auto& descriptor = model.descriptor();
                if (definition.nodes.size() != static_cast<std::size_t>(descriptor.terminal_count) ||
                    definition.internal_nodes.size() != internal_node_count) {
                    return nullptr;
                }
                gspice::GsdiModelCard card;
                card.name = definition.name;
                card.type = definition.type;
                for (const auto& [name, value] : definition.model_params) {
                    double parsed = 0.0;
                    if (tryParseSpiceValue(value, parsed)) card.parameters[name] = parsed;
                }
                for (const auto& [name, value] : definition.instance_params) {
                    double parsed = 0.0;
                    if (tryParseSpiceValue(value, parsed)) card.parameters[name] = parsed;
                }
                std::vector<int> columns = definition.nodes;
                columns.insert(columns.end(), definition.internal_nodes.begin(),
                               definition.internal_nodes.end());
                return std::make_unique<gspice::GdiDevice>(
                    definition.name, model.createInstance(card, definition.nodes), columns,
                    gspice::GsdiCollapseMap::standard(descriptor));
            },
            internal_node_count);
        return true;
    }();
    (void)registered;
#endif
}

void ensurePsp103GmcModels() {
    static const bool registered = [] {
        auto factory = [](const gspice::GmcModelDefinition& definition)
            -> std::unique_ptr<gspice::Device> {
            if (definition.nodes.size() != 4) return nullptr;
            const auto model = gspice::Psp103ParameterSet::from(
                [&definition] {
                    gspice::GsdiParamMap values;
                    for (const auto& [name, value] : definition.model_params) {
                        double parsed = 0.0;
                        if (tryParseSpiceValue(value, parsed)) values[name] = parsed;
                    }
                    return values;
                }());
            gspice::GsdiParamMap instance;
            for (const auto& [name, value] : definition.instance_params) {
                double parsed = 0.0;
                if (tryParseSpiceValue(value, parsed)) instance[name] = parsed;
            }
            auto prepared = model.prepare(instance, definition.temperature_c);
            if (!prepared.valid()) return nullptr;
            // PSP103 is a pure 4-terminal model (D, G, S, B) with no hidden
            // nodes, so the device runs on local indices 0-3 and the route
            // mirrors Phase C/D (native evaluator -> GsdiDaeDeviceAdapter ->
            // collapse-aware GdiDevice). The descriptor carries the terminal
            // set and dense 4x4 pattern; standard() composes the identity
            // collapse plan (nothing collapses) so the DAE stamp layer is
            // numerically identical to the direct native path.
            auto device = std::make_unique<gspice::Psp103Mosfet>(
                definition.name, 0, 1, 2, 3, model, instance,
                definition.temperature_c);
            return std::make_unique<gspice::GdiDevice>(
                definition.name,
                std::make_unique<gspice::GsdiDaeDeviceAdapter>(
                    std::move(device), 4),
                definition.nodes,
                gspice::GsdiCollapseMap::standard(gspice::psp103GsdiDescriptor()));
        };
        auto& registry = gspice::GmcRegistry::instance();
        for (const char* type : {
                 "PSP", "PSP103", "PSP103VA", "PSP103_VA",
                 "PSP103NQS", "PSPNQS103", "PSPNQS103VA"}) {
            registry.registerModel(type, factory);
        }
        return true;
    }();
    (void)registered;
}

std::string psp103GmcType(const std::string& type) {
    return gspice::GmcRegistry::instance().hasModel(type) ? type : "PSP103VA";
}

bool isMosLevel103PspModel(const gspice::ModelCard& model) {
    const std::string type = toUpperCopy(model.type);
    if (type != "NMOS" && type != "PMOS" && type != "N" && type != "P") return false;
    return std::abs(paramValue(model.params, {"LEVEL"}, 1.0) - 103.0) < 1.0e-9;
}

bool isPsp103ModelCard(const gspice::ModelCard& model) {
    return gspice::CompactModelRegistry::instance().isPsp103ModelType(model.type) ||
           isMosLevel103PspModel(model);
}

std::string psp103GsdiType(const gspice::ModelCard& model) {
    return gspice::CompactModelRegistry::instance().isPsp103ModelType(model.type)
        ? model.type
        : "PSP103VA";
}

std::string psp103GmcType(const gspice::ModelCard& model) {
    return gspice::CompactModelRegistry::instance().isPsp103ModelType(model.type)
        ? psp103GmcType(model.type)
        : "PSP103VA";
}

std::unordered_map<std::string, std::string> psp103ModelParams(const gspice::ModelCard& model) {
    auto params = model.params;
    if (std::isfinite(paramValue(params, {"TYPE"}, std::numeric_limits<double>::quiet_NaN()))) {
        return params;
    }
    const std::string type = toUpperCopy(model.type);
    if (type == "PMOS" || type == "P") params["TYPE"] = "-1";
    if (type == "NMOS" || type == "N") params["TYPE"] = "1";
    return params;
}

std::string psp103IgnoredParameterSummary(const gspice::ModelCard& model) {
    std::set<std::string> ignored;
    for (const auto& [name, value] : model.params) {
        (void)value;
        if (!gspice::Psp103ParameterSet::handlesModelParameter(name)) {
            ignored.insert(toUpperCopy(name));
        }
    }
    if (ignored.empty()) return "";
    const std::size_t shown = std::min<std::size_t>(ignored.size(), 12);
    std::string summary = std::to_string(ignored.size()) + " unsupported PSP parameter(s) ignored by native evaluator: ";
    std::size_t i = 0;
    for (const auto& name : ignored) {
        if (i >= shown) break;
        if (i) summary += ", ";
        summary += name;
        ++i;
    }
    if (shown < ignored.size()) summary += ", ...";
    return summary;
}

double psp103RfGateResistance(
    const gspice::ModelCard& model,
    const std::unordered_map<std::string, std::string>& instanceParams) {
    const bool rfNamedModel = toUpperCopy(model.name).find("_RF") != std::string::npos;
    double rgo = paramValue(model.params, {"RGO", "RG", "RVPOLY"}, rfNamedModel ? 40.0 : 0.0);
    double rshg = paramValue(model.params, {"RSHG", "RSH", "RSHD"}, rfNamedModel ? 3.0 : 0.0);
    if (rfNamedModel && rgo <= 0.0) rgo = 40.0;
    if (rfNamedModel && rshg <= 0.0) rshg = 3.0;
    const double nf = std::max(paramValue(instanceParams, {"NF", "NG", "NFIN"}, 1.0), 1.0);
    const double m = std::max(paramValue(instanceParams, {"M", "MULT"}, 1.0), 1.0e-30);
    const double gateResistance = (std::max(rgo, 0.0) + std::max(rshg, 0.0) / nf) / m;
    return std::isfinite(gateResistance) && gateResistance > 0.0 ? gateResistance : 0.0;
}

std::string bsim4IgnoredParameterSummary(const gspice::ModelCard& model) {
    std::set<std::string> ignored;
    for (const auto& [name, value] : model.params) {
        (void)value;
        if (!gspice::Bsim4ImplementedParameters::isSupported(name)) {
            ignored.insert(toUpperCopy(name));
        }
    }
    if (ignored.empty()) return "";
    const std::size_t shown = std::min<std::size_t>(ignored.size(), 12);
    std::string summary = std::to_string(ignored.size()) + " unsupported BSIM4 model parameter(s) ignored by native evaluator: ";
    std::size_t i = 0;
    for (const auto& name : ignored) {
        if (i >= shown) break;
        if (i) summary += ", ";
        summary += name;
        ++i;
    }
    if (shown < ignored.size()) summary += ", ...";
    return summary;
}

std::string bsim3IgnoredParameterSummary(const gspice::ModelCard& model) {
    std::set<std::string> ignored;
    for (const auto& [name, value] : model.params) {
        (void)value;
        if (!gspice::Bsim3ImplementedParameters::isSupported(name)) {
            ignored.insert(toUpperCopy(name));
        }
    }
    if (ignored.empty()) return "";
    const std::size_t shown = std::min<std::size_t>(ignored.size(), 12);
    std::string summary = std::to_string(ignored.size()) + " unsupported BSIM3 model parameter(s) ignored by native evaluator: ";
    std::size_t i = 0;
    for (const auto& name : ignored) {
        if (i >= shown) break;
        if (i) summary += ", ";
        summary += name;
        ++i;
    }
    if (shown < ignored.size()) summary += ", ...";
    return summary;
}

void ensureJuncapExpressGmcModels() {
    static const bool registered = [] {
        auto factory = [](const gspice::GmcModelDefinition& definition)
            -> std::unique_ptr<gspice::Device> {
            if (definition.nodes.size() != 2) return nullptr;
            std::unordered_map<std::string, double> modelParams;
            for (const auto& [name, value] : definition.model_params) {
                double parsed = 0.0;
                if (tryParseSpiceValue(value, parsed)) modelParams[name] = parsed;
            }
            std::unordered_map<std::string, double> instanceParams;
            for (const auto& [name, value] : definition.instance_params) {
                double parsed = 0.0;
                if (tryParseSpiceValue(value, parsed)) instanceParams[name] = parsed;
            }
            auto instance = gspice::GmcJuncapExpressModule{}.createInstance(
                modelParams, instanceParams);
            if (!instance) return nullptr;
            return std::make_unique<gspice::GdiDevice>(
                definition.name, std::move(instance), definition.nodes);
        };
        gspice::GmcRegistry::instance().registerModel("JUNCAPEXP", factory);
        gspice::GmcRegistry::instance().registerModel("JUNCAP2", [](const gspice::GmcModelDefinition& definition)
            -> std::unique_ptr<gspice::Device> {
            if (definition.nodes.size() != 2) return nullptr;
            std::unordered_map<std::string, double> modelParams;
            for (const auto& [name, value] : definition.model_params) {
                double parsed = 0.0;
                if (tryParseSpiceValue(value, parsed)) modelParams[name] = parsed;
            }
            auto instance = gspice::GmcJuncap2Module{}.createInstance(modelParams, {});
            if (!instance) return nullptr;
            return std::make_unique<gspice::GdiDevice>(definition.name, std::move(instance), definition.nodes);
        });
        return true;
    }();
    (void)registered;
}

bool modelTypeMatches(const gspice::ModelCard* model, const std::vector<std::string>& types) {
    if (!model) return false;
    const std::string actual = toUpperCopy(model->type);
    for (const auto& type : types) {
        if (actual == toUpperCopy(type)) return true;
    }
    return false;
}

bool primitiveModelFallbackEnabled() {
    return envFlagEnabled("GSPICE_ALLOW_PRIMITIVE_MODEL_FALLBACK");
}

bool hasExtension(std::string path, const std::string& extension) {
    const auto first = path.find_first_not_of(" \t\r\n");
    const auto last = path.find_last_not_of(" \t\r\n");
    path = first == std::string::npos ? "" : path.substr(first, last - first + 1);
    if (path.size() >= 2 && ((path.front() == '"' && path.back() == '"') ||
                             (path.front() == '\'' && path.back() == '\''))) {
        path = path.substr(1, path.size() - 2);
    }
    std::filesystem::path fsPath(path);
    return toUpperCopy(fsPath.extension().string()) == toUpperCopy(extension);
}

std::string readWholeFile(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return "";
    std::ostringstream out;
    out << file.rdbuf();
    return out.str();
}

std::string extractJsonString(const std::string& text, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    const auto keyPos = text.find(needle);
    if (keyPos == std::string::npos) return "";
    const auto colon = text.find(':', keyPos + needle.size());
    if (colon == std::string::npos) return "";
    const auto quote = text.find('"', colon + 1);
    if (quote == std::string::npos) return "";
    const auto end = text.find('"', quote + 1);
    if (end == std::string::npos) return "";
    return text.substr(quote + 1, end - quote - 1);
}

int extractJsonInt(const std::string& text, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    const auto keyPos = text.find(needle);
    if (keyPos == std::string::npos) return 0;
    const auto colon = text.find(':', keyPos + needle.size());
    if (colon == std::string::npos) return 0;
    const auto start = text.find_first_of("0123456789", colon + 1);
    if (start == std::string::npos) return 0;
    const auto end = text.find_first_not_of("0123456789", start);
    return std::stoi(text.substr(start, end - start));
}

std::vector<std::string> extractGsdiParameterNames(const std::string& text) {
    std::vector<std::string> names;
    const std::string needle = "\"parameters\"";
    const auto keyPos = text.find(needle);
    if (keyPos == std::string::npos) return names;
    const auto arrayStart = text.find('[', keyPos + needle.size());
    const auto arrayEnd = text.find(']', arrayStart);
    if (arrayStart == std::string::npos || arrayEnd == std::string::npos) return names;
    std::string params = text.substr(arrayStart, arrayEnd - arrayStart + 1);
    std::size_t pos = 0;
    while ((pos = params.find("\"name\"", pos)) != std::string::npos) {
        const auto colon = params.find(':', pos + 6);
        const auto quote = params.find('"', colon + 1);
        const auto end = params.find('"', quote + 1);
        if (colon == std::string::npos || quote == std::string::npos || end == std::string::npos) break;
        names.push_back(params.substr(quote + 1, end - quote - 1));
        pos = end + 1;
    }
    return names;
}

bool gsdiAllowsParameter(const gspice::GsdiArtifactInfo& artifact, const std::string& name) {
    if (artifact.parameter_names.empty()) return true;
    const std::string wanted = toUpperCopy(name);
    for (const auto& declared : artifact.parameter_names) {
        if (toUpperCopy(declared) == wanted) return true;
    }
    return false;
}

std::string firstUnknownGsdiParameter(
    const gspice::GsdiArtifactInfo& artifact,
    const std::unordered_map<std::string, std::string>& params) {
    for (const auto& [name, value] : params) {
        (void)value;
        if (!gsdiAllowsParameter(artifact, name)) return name;
    }
    return "";
}

bool isLikelyCompactMosModelName(const std::string& modelName) {
    const std::string name = toUpperCopy(modelName);
    return name.find("SG13") != std::string::npos ||
           name.find("PSP") != std::string::npos ||
           name.find("BSIM") != std::string::npos ||
           name.find("HICUM") != std::string::npos ||
           name.find("EKV") != std::string::npos ||
           name.find("NFET") != std::string::npos ||
           name.find("PFET") != std::string::npos;
}

bool isLikelyCompactMosModelType(const std::string& modelType) {
    return gspice::CompactModelRegistry::instance().looksLikeCompactModel(modelType);
}

bool isUnsupportedCompactMosLevel(const gspice::ModelCard* model) {
    if (!model) return false;
    const double level = paramValue(model->params, {"LEVEL"}, 1.0);
    return level >= 49.0 && level <= 55.0;
}

std::string applyLocalParams(
    std::string line,
    const std::unordered_map<std::string, std::string>& params) {
    for (const auto& [key, value] : params) {
        const std::string resolved = resolvedParamString(value, params);
        replaceAll(line, "{" + key + "}", resolved);
        replaceAll(line, "'" + key + "'", resolved);
    }
    line = resolveBracedNumericExpressions(line, params);
    auto tokens = tokenizeSimple(line);
    for (size_t i = 0; i < tokens.size(); ++i) {
        auto& token = tokens[i];
        if (i + 2 < tokens.size() && tokens[i + 1] == "=") {
            double evaluated = 0.0;
            if (tryEvaluateParamExpression(tokens[i + 2], params, evaluated)) {
                tokens[i + 2] = formatNumericValue(evaluated);
            }
            continue;
        }
        auto [key, value] = splitParameterToken(token);
        if (!key.empty()) {
            if (isStringValuedAssignmentKey(key)) continue;
            double evaluated = 0.0;
            if (tryEvaluateParamExpression(value, params, evaluated)) {
                token = key + "=" + formatNumericValue(evaluated);
            }
        } else {
            std::string stripped = stripQuotes(token);
            auto it = params.find(stripped);
            if (it == params.end()) it = params.find(toUpperCopy(stripped));
            if (it != params.end()) token = resolvedParamString(it->second, params);
        }
    }
    return joinSimple(tokens);
}

bool isIhpLvMosWrapper(const std::string& subcktNameUpper, int& typeOut) {
    if (subcktNameUpper == "SG13_LV_NMOS") {
        typeOut = 1;
        return true;
    }
    if (subcktNameUpper == "SG13_LV_PMOS") {
        typeOut = -1;
        return true;
    }
    return false;
}

bool ihpR3CmcWrapperRshKey(const std::string& subcktNameUpper, std::string& key, double& nominal) {
    if (subcktNameUpper == "RSIL") {
        key = "RSH_RSIL";
        nominal = 7.0;
        return true;
    }
    if (subcktNameUpper == "RHIGH") {
        key = "RSH_RHIGH";
        nominal = 1360.0;
        return true;
    }
    if (subcktNameUpper == "RPPD") {
        key = "RSH_RPPD";
        nominal = 260.0;
        return true;
    }
    return false;
}

struct IhpR3WrapperDefaults {
    double xw = 0.0;
    double rzspec = 0.0;
    double tc1 = 3100e-6;
    double tc2 = 0.3e-6;
};

IhpR3WrapperDefaults ihpR3WrapperDefaults(const std::string& subcktNameUpper) {
    if (subcktNameUpper == "RSIL") return {0.01e-6, 4.5e-6, 3100e-6, 0.3e-6};
    if (subcktNameUpper == "RHIGH") return {-0.04e-6, 80e-6, -2300e-6, 2.1e-6};
    if (subcktNameUpper == "RPPD") return {0.006e-6, 35e-6, 170e-6, 0.4e-6};
    return {};
}

bool ihpTapWrapper(const std::string& subcktNameUpper) {
    return subcktNameUpper == "PTAP1" || subcktNameUpper == "NTAP1";
}

void addBuiltinPspGsdiArtifact(gspice::Netlist& netlist, const std::string& modelType) {
    if (netlist.findGsdiArtifact(modelType)) return;
    gspice::GsdiArtifactInfo artifact;
    artifact.model_type = modelType;
    artifact.path = "builtin:psp103";
    artifact.terminal_count = 4;
    netlist.addGsdiArtifact(artifact);
    netlist.addModelStatus("GSDI_BUILTIN: " + modelType + " PSP103 native evaluator");
}

#ifdef GSPICE_HAVE_GMC_GENERATED
// Embedded GMC-generated models (psp_like, va_probe) are baked into the
// gspice binary, so decks can reference them without an explicit .GSDI
// include. Register a builtin artifact derived from the model descriptor so
// the elaborator's terminal/parameter checks work the same as for loaded
// artifacts.
void addBuiltinGmcGeneratedGsdiArtifact(gspice::Netlist& netlist, const std::string& modelType) {
    if (netlist.findGsdiArtifact(modelType)) return;
    std::unique_ptr<gspice::GsdiModel> instance;
    const std::string upper = toUpperCopy(modelType);
    if (upper == "PSP_LIKE") {
        instance = std::make_unique<gspice::GmcPspLikeModel>();
    } else if (upper == "VA_PROBE") {
        instance = std::make_unique<gspice::GmcVaProbeModel>();
    }
    if (!instance) return;
    const auto& desc = instance->descriptor();
    gspice::GsdiArtifactInfo artifact;
    artifact.model_type = desc.model_type;
    artifact.path = "builtin:gmc-generated";
    artifact.terminal_count = desc.terminal_count;
    for (const auto& param : desc.parameters) {
        artifact.parameter_names.push_back(param.name);
    }
    netlist.addGsdiArtifact(artifact);
    netlist.addModelStatus("GSDI_BUILTIN: " + desc.model_type + " GMC-generated evaluator");
}
#endif

void appendUniquePath(std::vector<std::filesystem::path>& roots, const std::filesystem::path& path) {
    if (path.empty()) return;
    auto normalized = path.lexically_normal();
    if (std::find(roots.begin(), roots.end(), normalized) == roots.end()) {
        roots.push_back(normalized);
    }
}

std::string getEnvVar(const char* envName) {
    return readEnvVar(envName);
}

void appendEnvSearchRoots(std::vector<std::filesystem::path>& roots, const char* envName) {
    std::string value = getEnvVar(envName);
    if (value.empty()) return;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ';')) {
        item = trimCopy(item);
        if (!item.empty()) appendUniquePath(roots, item);
    }
}

bool primitiveIhpFallbackEnabled() {
    std::string value = getEnvVar("GSPICE_ALLOW_PRIMITIVE_IHP_FALLBACK");
    if (value.empty()) return false;
    value = toUpperCopy(value);
    return value == "1" || value == "YES" || value == "TRUE" || value == "ON";
}

std::string mapNodeToken(
    const std::string& token,
    const std::string& instancePrefix,
    const std::unordered_map<std::string, std::string>& pinMap) {
    auto it = pinMap.find(token);
    if (it != pinMap.end()) return it->second;
    auto itUpper = pinMap.find(toUpperCopy(token));
    if (itUpper != pinMap.end()) return itUpper->second;
    if (isGroundName(token) || isParameterToken(token) || token.empty()) return token;
    return instancePrefix + ":" + token;
}

std::vector<PreLine> readLogicalLines(
    const std::filesystem::path& filePath,
    bool skipTitle,
    std::vector<std::string>& errors) {
    std::vector<PreLine> lines;
    std::ifstream file(filePath.string());
    if (!file.is_open()) {
        errors.push_back("Could not open file: " + filePath.string());
        return lines;
    }

    std::string raw;
    int lineNo = 0;
    std::string pending;
    int pendingLine = 0;
    while (std::getline(file, raw)) {
        ++lineNo;
        if (skipTitle && lineNo == 1) continue;
        std::string line = normalizeLine(raw);
        if (line.empty() || line[0] == '*' || line[0] == '$') continue;
        if (!line.empty() && line[0] == '+') {
            pending += " ";
            pending += trimCopy(line.substr(1));
            continue;
        }
        if (!pending.empty()) {
            lines.push_back({pending, filePath.string(), pendingLine});
        }
        pending = line;
        pendingLine = lineNo;
    }
    if (!pending.empty()) {
        lines.push_back({pending, filePath.string(), pendingLine});
    }
    return lines;
}

std::vector<PreLine> loadNetlistFile(
    const std::filesystem::path& filePath,
    bool skipTitle,
    std::vector<std::string>& errors,
    std::set<std::string>& includeStack);

std::vector<PreLine> loadLibrarySection(
    const std::filesystem::path& filePath,
    const std::string& section,
    std::vector<std::string>& errors,
    std::set<std::string>& includeStack) {
    std::vector<PreLine> selected;
    auto lines = readLogicalLines(filePath, false, errors);
    const std::string wanted = toUpperCopy(section);
    bool inSection = false;
    for (const auto& line : lines) {
        auto tokens = tokenizeSimple(line.text);
        if (tokens.empty()) continue;
        std::string cmd = toUpperCopy(tokens[0]);
        if (cmd == ".LIB" && tokens.size() >= 2) {
            std::string name = toUpperCopy(stripQuotes(tokens.back()));
            if (name == wanted) inSection = true;
            continue;
        }
        if (cmd == ".ENDL" || cmd == ".ENDLIB") {
            if (inSection) break;
            continue;
        }
        if (inSection) {
            if ((cmd == ".INCLUDE" || cmd == ".INC") && tokens.size() >= 2) {
                auto includePath = resolveRelativePath(filePath, tokens[1]);
                auto included = loadNetlistFile(includePath, false, errors, includeStack);
                selected.insert(selected.end(), included.begin(), included.end());
                continue;
            }
            if (cmd == ".LIB" && tokens.size() >= 3) {
                auto libPath = resolveRelativePath(filePath, tokens[1]);
                auto nested = loadLibrarySection(libPath, tokens[2], errors, includeStack);
                selected.insert(selected.end(), nested.begin(), nested.end());
                continue;
            }
            selected.push_back(line);
        }
    }
    if (selected.empty()) {
        errors.push_back("Library section '" + section + "' not found in " + filePath.string());
    }
    return selected;
}

std::vector<PreLine> loadNetlistFile(
    const std::filesystem::path& filePath,
    bool skipTitle,
    std::vector<std::string>& errors,
    std::set<std::string>& includeStack) {
    std::vector<PreLine> output;
    std::filesystem::path normalized = filePath.lexically_normal();
    const std::string key = normalized.string();
    if (includeStack.count(key)) {
        errors.push_back("Recursive include detected: " + key);
        return output;
    }
    includeStack.insert(key);

    auto lines = readLogicalLines(normalized, skipTitle, errors);
    for (const auto& line : lines) {
        auto tokens = tokenizeSimple(line.text);
        if (tokens.empty()) continue;
        std::string cmd = toUpperCopy(tokens[0]);
        if ((cmd == ".INCLUDE" || cmd == ".INC") && tokens.size() >= 2) {
            auto includePath = resolveRelativePath(normalized, tokens[1]);
            auto included = loadNetlistFile(includePath, false, errors, includeStack);
            output.insert(output.end(), included.begin(), included.end());
            continue;
        }
        if (cmd == ".LIB" && tokens.size() >= 3) {
            auto libPath = resolveRelativePath(normalized, tokens[1]);
            auto selected = loadLibrarySection(libPath, tokens[2], errors, includeStack);
            output.insert(output.end(), selected.begin(), selected.end());
            continue;
        }
        output.push_back(line);
    }
    includeStack.erase(key);
    return output;
}

std::unordered_map<std::string, SubcktDef> collectSubckts(
    const std::vector<PreLine>& input,
    std::vector<PreLine>& topLevel,
    std::vector<std::string>& errors) {
    std::unordered_map<std::string, SubcktDef> subckts;
    bool inSubckt = false;
    SubcktDef current;
    for (const auto& line : input) {
        auto tokens = tokenizeSimple(line.text);
        if (tokens.empty()) continue;
        std::string cmd = toUpperCopy(tokens[0]);
        if (cmd == ".SUBCKT") {
            if (tokens.size() < 2) {
                errors.push_back("Invalid .SUBCKT at " + line.source + ":" + std::to_string(line.lineNo));
                continue;
            }
            inSubckt = true;
            current = SubcktDef{};
            current.name = tokens[1];
            for (size_t i = 2; i < tokens.size(); ++i) {
                if (isParameterToken(tokens[i])) {
                    auto [key, value] = splitParameterToken(tokens[i]);
                    addParam(current.params, key, value);
                } else {
                    current.pins.push_back(tokens[i]);
                }
            }
            continue;
        }
        if (cmd == ".ENDS" || cmd == ".ENDSUBCKT") {
            if (inSubckt) {
                subckts[toUpperCopy(current.name)] = current;
                inSubckt = false;
            }
            continue;
        }
        if (inSubckt && cmd == ".MODEL") {
            topLevel.push_back(line);
            current.body.push_back(line);
        } else if (inSubckt) {
            current.body.push_back(line);
        } else {
            topLevel.push_back(line);
        }
    }
    if (inSubckt) {
        errors.push_back("Unterminated .SUBCKT " + current.name);
    }
    return subckts;
}

size_t findSubcktNameIndex(const std::vector<std::string>& tokens) {
    for (size_t i = 1; i < tokens.size(); ++i) {
        if (isParameterToken(tokens[i])) return i > 1 ? i - 1 : i;
    }
    return tokens.empty() ? 0 : tokens.size() - 1;
}

std::string remapPrimitiveLine(
    const PreLine& line,
    const std::string& instancePrefix,
    const std::unordered_map<std::string, std::string>& pinMap) {
    auto tokens = tokenizeSimple(line.text);
    if (tokens.empty()) return line.text;
    tokens[0] = tokens[0] + "_" + instancePrefix;
    char first = static_cast<char>(std::toupper(tokens[0][0]));
    size_t firstNode = 1;
    size_t lastNodeExclusive = 1;
    if (first == 'R' || first == 'C' || first == 'L' || first == 'V' || first == 'I' || first == 'D' || first == 'W' || first == 'P') {
        lastNodeExclusive = std::min<size_t>(3, tokens.size());
    } else if (first == 'M') {
        lastNodeExclusive = std::min<size_t>(5, tokens.size());
    } else if (first == 'N') {
        lastNodeExclusive = std::min<size_t>(5, tokens.size());
    } else if (first == 'S') {
        lastNodeExclusive = tokens.size();
        for (size_t i = 1; i < tokens.size(); ++i) {
            if (isParameterToken(tokens[i])) {
                lastNodeExclusive = i;
                break;
            }
        }
    }
    for (size_t i = firstNode; i < lastNodeExclusive; ++i) {
        tokens[i] = mapNodeToken(tokens[i], instancePrefix, pinMap);
    }
    return joinSimple(tokens);
}

std::vector<PreLine> expandSubcktInstance(
    const PreLine& line,
    const std::unordered_map<std::string, SubcktDef>& subckts,
    std::vector<std::string>& errors,
    int depth);

std::vector<PreLine> expandLines(
    const std::vector<PreLine>& lines,
    const std::unordered_map<std::string, SubcktDef>& subckts,
    std::vector<std::string>& errors,
    int depth) {
    std::vector<PreLine> expanded;
    for (const auto& line : lines) {
        auto tokens = tokenizeSimple(line.text);
        if (!tokens.empty() && std::toupper(tokens[0][0]) == 'X') {
            auto sub = expandSubcktInstance(line, subckts, errors, depth + 1);
            expanded.insert(expanded.end(), sub.begin(), sub.end());
        } else {
            expanded.push_back(line);
        }
    }
    return expanded;
}

struct ConditionalFrame {
    bool parentActive = true;
    bool active = true;
    bool branchTaken = false;
    bool elseSeen = false;
};

bool conditionalStackActive(const std::vector<ConditionalFrame>& stack) {
    return stack.empty() ? true : stack.back().active;
}

std::vector<PreLine> filterConditionalLines(
    const std::vector<PreLine>& lines,
    std::unordered_map<std::string, std::string>& params,
    std::vector<std::string>& errors,
    const std::string& context) {
    std::vector<PreLine> filtered;
    std::vector<ConditionalFrame> stack;

    for (const auto& line : lines) {
        auto tokens = tokenizeSimple(line.text);
        if (tokens.empty()) continue;
        const std::string cmd = toUpperCopy(tokens[0]);
        if (cmd == ".PARAM" || cmd == ".PARAMS") {
            if (conditionalStackActive(stack)) {
                auto parsed = parseParameterAssignments(tokens, 1);
                for (const auto& [key, value] : parsed) params[key] = value;
                resolveParameterMapExpressions(params);
            }
            continue;
        }
        if (cmd == ".IF") {
            const bool parentActive = conditionalStackActive(stack);
            double condition = 0.0;
            bool conditionOk = false;
            if (!parentActive) {
                conditionOk = true;
            } else {
                conditionOk = tryEvaluateParamExpression(joinTokens(tokens, 1), params, condition);
            }
            if (!conditionOk) {
                errors.push_back("Could not evaluate .IF condition in " + context + " at " +
                                 line.source + ":" + std::to_string(line.lineNo) + ": " + line.text);
            }
            const bool branchActive = conditionOk && condition != 0.0;
            stack.push_back({parentActive, parentActive && branchActive, branchActive, false});
            continue;
        }
        if (cmd == ".ELIF" || cmd == ".ELSEIF") {
            if (stack.empty()) {
                errors.push_back("Unexpected " + cmd + " without .IF in " + context + " at " +
                                 line.source + ":" + std::to_string(line.lineNo));
                continue;
            }
            auto& frame = stack.back();
            if (frame.elseSeen) {
                errors.push_back(cmd + " after .ELSE in " + context + " at " +
                                 line.source + ":" + std::to_string(line.lineNo));
                frame.active = false;
                continue;
            }
            double condition = 0.0;
            bool branchActive = false;
            if (frame.parentActive && !frame.branchTaken) {
                if (!tryEvaluateParamExpression(joinTokens(tokens, 1), params, condition)) {
                    errors.push_back("Could not evaluate " + cmd + " condition in " + context + " at " +
                                     line.source + ":" + std::to_string(line.lineNo) + ": " + line.text);
                } else {
                    branchActive = condition != 0.0;
                }
            }
            frame.active = frame.parentActive && !frame.branchTaken && branchActive;
            frame.branchTaken = frame.branchTaken || branchActive;
            continue;
        }
        if (cmd == ".ELSE") {
            if (stack.empty()) {
                errors.push_back("Unexpected .ELSE without .IF in " + context + " at " +
                                 line.source + ":" + std::to_string(line.lineNo));
                continue;
            }
            auto& frame = stack.back();
            frame.active = frame.parentActive && !frame.branchTaken;
            frame.branchTaken = true;
            frame.elseSeen = true;
            continue;
        }
        if (cmd == ".ENDIF") {
            if (stack.empty()) {
                errors.push_back("Unexpected .ENDIF without .IF in " + context + " at " +
                                 line.source + ":" + std::to_string(line.lineNo));
                continue;
            }
            stack.pop_back();
            continue;
        }
        if (conditionalStackActive(stack)) {
            filtered.push_back(line);
        }
    }

    if (!stack.empty()) {
        errors.push_back("Unterminated .IF block in " + context);
    }
    return filtered;
}

std::vector<PreLine> expandSubcktInstance(
    const PreLine& line,
    const std::unordered_map<std::string, SubcktDef>& subckts,
    std::vector<std::string>& errors,
    int depth) {
    std::vector<PreLine> expanded;
    if (depth > 32) {
        errors.push_back("Subcircuit expansion depth exceeded at " + line.source + ":" + std::to_string(line.lineNo));
        return expanded;
    }
    auto tokens = tokenizeSimple(line.text);
    if (tokens.size() < 2) return expanded;
    size_t subcktIdx = findSubcktNameIndex(tokens);
    if (subcktIdx <= 1 || subcktIdx >= tokens.size()) {
        errors.push_back("Invalid subcircuit call: " + line.text);
        return expanded;
    }
    std::string subcktName = tokens[subcktIdx];
    auto it = subckts.find(toUpperCopy(subcktName));
    if (it == subckts.end()) {
        errors.push_back("Unknown subcircuit '" + subcktName + "' in line: " + line.text);
        return expanded;
    }
    const auto& def = it->second;
    if (subcktIdx - 1 < def.pins.size()) {
        errors.push_back("Subcircuit '" + subcktName + "' expects " + std::to_string(def.pins.size()) +
                         " pins but instance has " + std::to_string(subcktIdx - 1) + ": " + line.text);
        return expanded;
    }

    std::unordered_map<std::string, std::string> pinMap;
    for (size_t i = 0; i < def.pins.size(); ++i) {
        pinMap[def.pins[i]] = tokens[i + 1];
        pinMap[toUpperCopy(def.pins[i])] = tokens[i + 1];
    }
    std::string prefix = tokens[0];

    auto localParams = def.params;
    auto instanceParams = parseParameterTokens(tokens, subcktIdx + 1);
    for (const auto& [key, value] : instanceParams) {
        localParams[key] = value;
    }

    int ihpMosType = 0;
    if (isIhpLvMosWrapper(toUpperCopy(subcktName), ihpMosType) && def.pins.size() >= 4) {
        const std::string model = ihpMosType < 0 ? "sg13g2_lv_pmos_psp" : "sg13g2_lv_nmos_psp";
        std::string w = "1u";
        std::string l = "1u";
        std::string nf = "1";
        std::string mult = "1";
        std::string dta = "0";
        auto wit = localParams.find("w");
        if (wit == localParams.end()) wit = localParams.find("W");
        if (wit != localParams.end()) w = resolvedParamString(wit->second, localParams);
        auto lit = localParams.find("l");
        if (lit == localParams.end()) lit = localParams.find("L");
        if (lit != localParams.end()) l = resolvedParamString(lit->second, localParams);
        auto nfit = localParams.find("ng");
        if (nfit == localParams.end()) nfit = localParams.find("NG");
        if (nfit == localParams.end()) nfit = localParams.find("nf");
        if (nfit == localParams.end()) nfit = localParams.find("NF");
        if (nfit != localParams.end()) nf = resolvedParamString(nfit->second, localParams);
        auto mit = localParams.find("m");
        if (mit == localParams.end()) mit = localParams.find("M");
        if (mit != localParams.end()) mult = resolvedParamString(mit->second, localParams);
        auto dtait = localParams.find("trise");
        if (dtait == localParams.end()) dtait = localParams.find("TRISE");
        if (dtait == localParams.end()) dtait = localParams.find("dta");
        if (dtait == localParams.end()) dtait = localParams.find("DTA");
        if (dtait != localParams.end()) dta = resolvedParamString(dtait->second, localParams);
        std::vector<std::string> mosLine = {
            "NPSP_" + prefix,
            tokens[1],
            tokens[2],
            tokens[3],
            tokens[4],
            model,
            "W=" + stripQuotes(w),
            "L=" + stripQuotes(l),
            "NF=" + stripQuotes(nf),
            "M=" + stripQuotes(mult),
            "DTA=" + stripQuotes(dta)
        };
        expanded.push_back({joinSimple(mosLine), line.source, line.lineNo});
        return expanded;
    }

    auto activeBody = filterConditionalLines(def.body, localParams, errors, "subckt " + subcktName);

    if (ihpTapWrapper(toUpperCopy(subcktName)) && def.pins.size() >= 2) {
        const double resistance = std::max(paramValue(localParams, {"R"}, 262.8), 1e-12);
        expanded.push_back({
            ".GSPICEWARN IHP tap wrapper " + subcktName + " approximated as fixed resistor.",
            line.source,
            line.lineNo
        });
        expanded.push_back({
            "R_" + prefix + " " + tokens[1] + " " + tokens[2] + " " +
                formatNumericValue(resistance),
            line.source,
            line.lineNo
        });
        return expanded;
    }

    std::string rshKey;
    double nominalRsh = 0.0;
    if (ihpR3CmcWrapperRshKey(toUpperCopy(subcktName), rshKey, nominalRsh) && def.pins.size() >= 2) {
        const auto defaults = ihpR3WrapperDefaults(toUpperCopy(subcktName));
        const double rsh = paramValue(localParams, {rshKey}, nominalRsh);
        const double drawnW = paramValue(localParams, {"W"}, 0.5e-6);
        const double b = std::max(paramValue(localParams, {"B"}, 0.0), 0.0);
        const double kappa = std::max(paramValue(localParams, {"KAPPA"}, 1.85), 1e-12);
        const double ps = paramValue(localParams, {"PS"}, 0.18e-6);
        const double w = paramValue(localParams, {"WEFF"}, drawnW + defaults.xw);
        const double l = paramValue(localParams, {"LEFF"}, (b + 1.0) * paramValue(localParams, {"L"}, 0.5e-6) + (2.0 / kappa * w + ps) * b);
        const double m = std::max(paramValue(localParams, {"M"}, 1.0), 1e-30);
        const double resScale = paramValue(localParams, {"RES_RPARA"}, 1.0);
        const double postsim = paramValue(localParams, {"POSTSIM"}, 0.0);
        const double rqrc = paramValue(localParams, {"RQRC"}, 4.5e-6);
        const double rc = paramValue(localParams, {"RZ", "RC"}, defaults.rzspec / std::max(drawnW, 1e-30) - (postsim > 0.0 ? rqrc / std::max(drawnW, 1e-30) : 0.0));
        const double trise = paramValue(localParams, {"TRISE"}, 0.0);
        const double tc1 = paramValue(localParams, {"TC1"}, defaults.tc1);
        const double tc2 = paramValue(localParams, {"TC2"}, defaults.tc2);
        if (std::isfinite(rsh) && w > 0.0) {
            const double tempScale = std::max(1.0 + tc1 * trise + tc2 * trise * trise, 1e-12);
            const double resistance = std::max((rsh * resScale * l / w * tempScale + 2.0 * rc) / m, 1e-12);
            expanded.push_back({
                ".GSPICEWARN IHP r3_cmc wrapper " + subcktName +
                    " routed as native RSH*LEFF/WEFF with foundry wrapper width/contact/corner/temperature terms.",
                line.source,
                line.lineNo
            });
            expanded.push_back({
                "R_" + prefix + " " + tokens[1] + " " + tokens[2] +
                    " " + formatNumericValue(resistance),
                line.source,
                line.lineNo
            });
            return expanded;
        }
    }

    if (toUpperCopy(subcktName) == "SG13_HV_SVARICAP" && def.pins.size() >= 4) {
        const double l = paramValue(localParams, {"L"}, 600e-9);
        const double w = paramValue(localParams, {"W"}, 3e-6);
        const double nx = std::max(paramValue(localParams, {"NX"}, 1.0), 1.0);
        const double ny = std::max(paramValue(localParams, {"NY"}, 1.0), 1.0);
        const double toxo = std::max(paramValue(localParams, {"TOXO"}, 6.945e-9), 1e-12);
        const double fingers = nx * ny;
        const double epsOx = 3.453133e-11;
        const double cAccum = std::max(epsOx * l * w * fingers / toxo, 1e-18);
        const double cMin = std::max(0.15 * cAccum, 1e-18);
        const double vfb = paramValue(localParams, {"VFBO"}, -0.04009);
        const double slope = paramValue(localParams, {"MOSVAR_SLOPE", "SLOPE"}, 0.25);
        const double ck0 = paramValue(localParams, {"CK0"}, 4.267e-16);
        const double ckw = paramValue(localParams, {"CKW"}, 1.252e-10);
        const double ckwnx = paramValue(localParams, {"CKWNX"}, -7.948e-11);
        const double cCouple = std::max((ck0 + ckw * w + ckwnx * w / fingers) * fingers, 1e-18);
        const double rk0 = paramValue(localParams, {"RK0"}, 346.7);
        const double rkwnx = paramValue(localParams, {"RKWNX"}, 0.000631);
        const double rCouple = std::max((rk0 + rkwnx / std::max(w, 1e-30)) / fingers, 1e-6);
        const double rwell0 = paramValue(localParams, {"RWELL0"}, 32.85);
        const double rwellw = paramValue(localParams, {"RWELLW"}, -2.499e6);
        const double rwellnx = paramValue(localParams, {"RWELLNX"}, 4.759e6);
        const double rWell = std::max(((rwell0 + rwellw * w) + (rwellnx * w / fingers)) * fingers, 1e-6);
        const double rsubw0 = paramValue(localParams, {"RSUBW0"}, 0.2596);
        const double rsubwf = paramValue(localParams, {"RSUBWF"}, 0.0009212);
        const double rsubwexp = paramValue(localParams, {"RSUBWEXP"}, 0.6952);
        const double rSub = std::max((rsubw0 / std::sqrt(std::pow(nx, rsubwexp))) / (w + rsubwf), 1e-6);
        const double cic0 = std::max(paramValue(localParams, {"CIC0"}, 1e-18), 1e-24);
        const std::string w0 = prefix + "_W0";
        const std::string g2a = prefix + "_G2A";
        const std::string w1 = prefix + "_W1";
        expanded.push_back({
            ".GSPICEWARN IHP mosvar wrapper " + subcktName +
                " routed through native voltage-dependent MOSVAR capacitance plus foundry parasitic RC pieces.",
            line.source,
            line.lineNo
        });
        expanded.push_back({"C_" + prefix + "_IC1 " + tokens[1] + " " + w0 + " " + formatNumericValue(cic0), line.source, line.lineNo});
        expanded.push_back({"C_" + prefix + "_IC2 " + tokens[3] + " " + w0 + " " + formatNumericValue(cic0), line.source, line.lineNo});
        expanded.push_back({"R_" + prefix + "_WW0 " + tokens[2] + " " + w0 + " " + formatNumericValue(rWell), line.source, line.lineNo});
        expanded.push_back({"Y_" + prefix + "_MV1 " + tokens[1] + " " + w0 + " MOSVAR CACC=" + formatNumericValue(cAccum) + " CMIN=" + formatNumericValue(cMin) + " VFB=" + formatNumericValue(vfb) + " SLOPE=" + formatNumericValue(slope), line.source, line.lineNo});
        expanded.push_back({"Y_" + prefix + "_MV2 " + tokens[3] + " " + w0 + " MOSVAR CACC=" + formatNumericValue(cAccum) + " CMIN=" + formatNumericValue(cMin) + " VFB=" + formatNumericValue(vfb) + " SLOPE=" + formatNumericValue(slope), line.source, line.lineNo});
        expanded.push_back({".MODEL DAREA D IS=2.45e-17 N=4 CJO=1.444e-15", line.source, line.lineNo});
        expanded.push_back({"R_" + prefix + "_SUBW " + w1 + " " + tokens[4] + " " + formatNumericValue(rSub), line.source, line.lineNo});
        expanded.push_back({"D_" + prefix + "_SUBW " + w1 + " " + tokens[2] + " DAREA AREA=" + formatNumericValue(std::max(((nx * 0.38e-6) + (nx * l) + 1.11e-6) * (ny * (w + 0.97e-6)), 1e-30)), line.source, line.lineNo});
        expanded.push_back({"C_" + prefix + "_G1G2 " + tokens[1] + " " + g2a + " " + formatNumericValue(cCouple), line.source, line.lineNo});
        expanded.push_back({"R_" + prefix + "_G2A " + g2a + " " + tokens[3] + " " + formatNumericValue(rCouple), line.source, line.lineNo});
        return expanded;
    }

    std::vector<PreLine> bodyLines;
    for (const auto& body : activeBody) {
        auto bodyTokens = tokenizeSimple(body.text);
        if (bodyTokens.empty()) continue;
        PreLine paramBody = body;
        paramBody.text = applyLocalParams(paramBody.text, localParams);
        bodyTokens = tokenizeSimple(paramBody.text);
        const std::string bodyCmd = toUpperCopy(bodyTokens[0]);
        if (bodyCmd == ".MODEL") {
            continue;
        }
        if (std::toupper(bodyTokens[0][0]) == 'X') {
            size_t nestedIdx = findSubcktNameIndex(bodyTokens);
            for (size_t i = 1; i < nestedIdx; ++i) {
                bodyTokens[i] = mapNodeToken(bodyTokens[i], prefix, pinMap);
            }
            bodyTokens[0] = bodyTokens[0] + "_" + prefix;
            bodyLines.push_back({joinSimple(bodyTokens), paramBody.source, paramBody.lineNo});
        } else {
            bodyLines.push_back({remapPrimitiveLine(paramBody, prefix, pinMap), paramBody.source, paramBody.lineNo});
        }
    }
    return expandLines(bodyLines, subckts, errors, depth);
}

std::vector<double> parseParenNumberList(std::string text) {
    for (char& c : text) {
        if (c == ',' || c == '\t' || c == '(' || c == ')') c = ' ';
    }
    std::vector<double> values;
    std::stringstream ss(text);
    std::string tok;
    while (ss >> tok) {
        double v = 0.0;
        if (tryParseSpiceValue(tok, v)) {
            values.push_back(v);
        }
    }
    return values;
}

bool extractParenPayload(const std::string& sourceSpec, const std::string& keywordUpper, std::string& payloadOut) {
    std::string upper = toUpperCopy(sourceSpec);
    size_t keyPos = upper.find(keywordUpper);
    if (keyPos == std::string::npos) return false;
    size_t lpar = upper.find('(', keyPos);
    if (lpar == std::string::npos) return false;
    size_t rpar = upper.find(')', lpar + 1);
    if (rpar == std::string::npos) rpar = sourceSpec.size();
    payloadOut = sourceSpec.substr(lpar + 1, rpar - lpar - 1);
    return true;
}

void applyOptionToken(gspice::SimulationSettings& settings, const std::string& token) {
    std::string key;
    std::string value;
    size_t eq = token.find('=');
    if (eq == std::string::npos) {
        key = toUpperCopy(token);
        value = "";
    } else {
        key = toUpperCopy(token.substr(0, eq));
        value = token.substr(eq + 1);
    }

    double numeric = 0.0;
    const bool hasNumeric = !value.empty() && tryParseSpiceValue(value, numeric);
    auto truthy = [](std::string text) {
        std::transform(text.begin(), text.end(), text.begin(), ::toupper);
        return !(text == "0" || text == "NO" || text == "FALSE" || text == "OFF");
    };
    auto applyPreset = [&](std::string preset) {
        preset = toUpperCopy(preset);
        preset.erase(std::remove(preset.begin(), preset.end(), '_'), preset.end());
        preset.erase(std::remove(preset.begin(), preset.end(), '-'), preset.end());
        if (preset == "LIBERAL") preset = "LOW";
        else if (preset == "MODERATE") preset = "MEDIUM";
        else if (preset == "CONSERVATIVE") preset = "VERYHIGH";
        settings.tran_adaptive = true;
        if (preset == "LOW") {
            settings.reltol = 5e-3;
            settings.vntol = 10e-6;
            settings.abstol = 1e-12;
            settings.tran_lte_reltol = 2e-2;
            settings.tran_lte_abstol = 10e-6;
            settings.tran_max_iter = 40;
        } else if (preset == "MEDIUM") {
            settings.reltol = 1e-3;
            settings.vntol = 1e-6;
            settings.abstol = 1e-12;
            settings.tran_lte_reltol = 5e-3;
            settings.tran_lte_abstol = 1e-6;
            settings.tran_max_iter = 60;
        } else if (preset == "HIGH") {
            settings.reltol = 3e-4;
            settings.vntol = 300e-9;
            settings.abstol = 100e-15;
            settings.tran_lte_reltol = 1e-3;
            settings.tran_lte_abstol = 300e-9;
            settings.tran_trtol = 3.5;
            settings.tran_max_iter = 80;
        } else if (preset == "VERYHIGH") {
            settings.reltol = 1e-4;
            settings.vntol = 100e-9;
            settings.abstol = 10e-15;
            settings.tran_lte_reltol = 3e-4;
            settings.tran_lte_abstol = 100e-9;
            settings.chgtol = 1e-15;
            settings.tran_trtol = 1.0;
            settings.tran_method = "TRAPGEAR";
            settings.tran_max_iter = 120;
        }
    };
    auto applyNumericalPolicy = [&](std::string policy) {
        policy = toUpperCopy(policy);
        policy.erase(std::remove(policy.begin(), policy.end(), '_'), policy.end());
        policy.erase(std::remove(policy.begin(), policy.end(), '-'), policy.end());
        if (policy == "ROBUST" || policy == "CONSERVATIVE") {
            settings.source_stepping = true;
            settings.gmin_stepping = true;
            settings.line_search = true;
            settings.solver_singletons = true;
            settings.tran_adaptive = true;
            settings.tran_method = "GEAR2";
            settings.op_max_iter = std::max(settings.op_max_iter, 200);
            settings.tran_max_iter = std::max(settings.tran_max_iter, 160);
            settings.gmin = std::max(settings.gmin, 1e-12);
            settings.tran_lte_reltol = std::min(settings.tran_lte_reltol, 1e-3);
        } else if (policy == "FAST") {
            settings.source_stepping = false;
            settings.gmin_stepping = false;
            settings.line_search = false;
            settings.tran_adaptive = true;
            settings.tran_lte_reltol = std::max(settings.tran_lte_reltol, 5e-3);
        } else if (policy == "BALANCED" || policy == "DEFAULT") {
            settings.source_stepping = true;
            settings.gmin_stepping = true;
            settings.line_search = true;
            settings.solver_singletons = true;
            settings.tran_adaptive = true;
        }
    };
    if (key == "RELTOL" && hasNumeric && numeric > 0.0) {
        settings.reltol = numeric;
    }
    else if (key == "VNTOL" && hasNumeric && numeric > 0.0) settings.vntol = numeric;
    else if (key == "ABSTOL" && hasNumeric && numeric > 0.0) settings.abstol = numeric;
    else if (key == "GMIN" && hasNumeric && numeric >= 0.0) settings.gmin = numeric;
    else if ((key == "CSHUNT" || key == "CMIN" || key == "TRAN_CSHUNT") && hasNumeric && numeric >= 0.0) settings.cshunt = numeric;
    else if ((key == "THREADS" || key == "NUM_THREADS" || key == "NTHREADS" || key == "PARALLEL" || key == "CPUS") && hasNumeric && numeric > 0.0) {
        settings.num_threads = static_cast<int>(numeric);
    }
    else if ((key == "ACCURACY" || key == "ERRPRESET" || key == "ERR_PRESET") && !value.empty()) applyPreset(value);
    else if ((key == "NUMERICAL" || key == "CONVERGENCE" || key == "POLICY") && !value.empty()) applyNumericalPolicy(value);
    else if (key == "SIGNOFF" || key == "RF_SIGNOFF") {
        if (truthy(value)) {
            settings.hb_native_required = true;
            settings.max_pss_iter = std::max(settings.max_pss_iter, 20);
            settings.op_max_iter = std::max(settings.op_max_iter, 200);
            settings.tran_max_iter = std::max(settings.tran_max_iter, 160);
            settings.source_stepping = true;
            settings.gmin_stepping = true;
            settings.line_search = true;
            settings.solver_singletons = true;
            settings.tran_adaptive = true;
        }
    }
    else if (key == "TRTOL" && hasNumeric && numeric > 0.0) settings.tran_trtol = numeric;
    else if ((key == "TRAN_RELTOL" || key == "LTE_RELTOL") && hasNumeric && numeric > 0.0) settings.tran_lte_reltol = numeric;
    else if ((key == "TRABSTOL" || key == "TRAN_ABSTOL" || key == "LTE_ABSTOL") && hasNumeric && numeric > 0.0) settings.tran_lte_abstol = numeric;
    else if (key == "CHGTOL" && hasNumeric && numeric > 0.0) settings.chgtol = numeric;
    else if ((key == "MAXORD" || key == "MAXORDER") && hasNumeric) settings.tran_max_order = std::clamp(static_cast<int>(numeric), 1, 5);
    else if (key == "MAXSTEP" || key == "TRAN_MAXSTEP") {
        if (toUpperCopy(value) == "AUTO") {
            settings.tran_max_step_auto = true;
            settings.t_max_step = 0.0;
        } else if (hasNumeric && numeric > 0.0) {
            settings.tran_max_step_auto = false;
            settings.t_max_step = numeric;
        }
    }
    else if (key == "TMAX" && hasNumeric && numeric > 0.0) settings.t_max_step = numeric;
else if (key == "SAVE") {
        const std::string mode = toUpperCopy(value);
        if (mode == "NONE" || mode == "0" || mode == "FALSE" || mode == "OFF" || mode == "NO") {
            settings.save_none = true;
            settings.save_all = false;
            settings.saves.clear();
        } else if (mode == "ALL" || mode == "1" || mode == "TRUE" || mode == "ON" || mode == "YES") {
            settings.save_none = false;
            settings.save_all = true;
            settings.saves.clear();
        } else if (mode == "SELECTED") {
            settings.save_none = false;
            settings.save_all = false;
        }
    }
    else if ((key == "MINSTEP" || key == "TMIN" || key == "TRAN_MINSTEP") && hasNumeric && numeric > 0.0) settings.t_min_step = numeric;
    else if (key == "ADAPTIVE" || key == "TRAN_ADAPTIVE") settings.tran_adaptive = value.empty() ? true : truthy(value);
    else if (key == "SAVE_ADAPTIVE" || key == "SAVE_ADAPTIVE_STEPS" || key == "SAVEADAPTIVE") {
        settings.save_adaptive_steps = value.empty() ? true : truthy(value);
    }
    else if (key == "TRAN_PREDICTOR" || key == "PREDICTOR") settings.tran_predictor = value.empty() ? true : truthy(value);
    else if (key == "TRAN_ORDER_ADAPTIVE" || key == "ORDER_ADAPTIVE" || key == "ADAPTIVE_ORDER") {
        settings.tran_order_adaptive = value.empty() ? true : truthy(value);
    }
    else if (key == "TRAP_RINGING" || key == "TRAP_RINGING_CONTROL" || key == "TRAN_TRAP_RINGING") {
        settings.tran_trap_ringing = value.empty() ? true : truthy(value);
    }
    else if (key == "LTE_MODE" || key == "TRAN_LTE_MODE") {
        const std::string mode = toUpperCopy(value);
        if (mode == "PREDICTOR" || mode == "PC" || mode == "PREDICTORCORRECTOR") {
            settings.tran_lte_mode = "PREDICTOR";
        } else if (mode == "STEPDOUBLING" || mode == "DOUBLING" || mode == "ORACLE") {
            settings.tran_lte_mode = "STEPDOUBLING";
        }
    }
    else if ((key == "LTE_AUDIT_INTERVAL" || key == "TRAN_LTE_AUDIT_INTERVAL") && hasNumeric) {
        settings.tran_lte_audit_interval = std::clamp(static_cast<int>(numeric), 0, 1000000);
    }
    else if ((key == "SOLVER" || key == "LINEAR_SOLVER" || key == "SPARSE_SOLVER") && !value.empty()) {
        settings.solver_backend = toUpperCopy(value);
    }
    else if ((key == "ORDERING" || key == "MATRIX_ORDERING" || key == "SPARSE_ORDERING") && !value.empty()) {
        settings.solver_ordering = toUpperCopy(value);
    }
    else if (key == "SINGLETONS" || key == "SINGLETON_FILTER" || key == "SINGLETON_FILTERING") {
        settings.solver_singletons = value.empty() ? true : truthy(value);
    }
    else if (key == "SCALING" || key == "ROWSCALING" || key == "MATRIX_SCALING") {
        settings.solver_row_scaling = value.empty() ? true : truthy(value);
    }
    else if (key == "REFINEMENT" || key == "ITERATIVE_REFINEMENT") {
        if (hasNumeric) settings.solver_refinement_steps = std::clamp(static_cast<int>(numeric), 0, 3);
        else settings.solver_refinement_steps = value.empty() || truthy(value) ? 1 : 0;
    }
    else if (key == "SOURCESTEPPING" || key == "SRCSTEP" || key == "SOURCE_STEPPING") {
        settings.source_stepping = value.empty() ? true : truthy(value);
    }
    else if (key == "GMINSTEPPING" || key == "GMINSTEP" || key == "GMIN_STEPPING") {
        settings.gmin_stepping = value.empty() ? true : truthy(value);
    }
    else if (key == "LINESEARCH" || key == "DAMPING" || key == "NEWTON_DAMPING") {
        settings.line_search = value.empty() ? true : truthy(value);
    }
    else if (key == "NR_RESIDUALCHECK" || key == "RESIDUALCHECK" || key == "RESIDUAL_CHECK") {
        settings.nr_residual_check = value.empty() ? true : truthy(value);
    }
    else if (key == "NR_BYPASS" || key == "BYPASS" || key == "DEVICE_BYPASS") {
        settings.nr_bypass = value.empty() ? true : truthy(value);
    }
    else if ((key == "NR_BYPASS_TOL" || key == "BYPASS_TOL") && hasNumeric && numeric > 0.0) {
        settings.nr_bypass_tolerance = std::clamp(numeric, 1e-6, 1.0);
    }
    else if ((key == "NODESET_ITERS" || key == "NODESET_ITERATIONS") && hasNumeric && numeric >= 0.0) {
        settings.nodeset_iterations = std::clamp(static_cast<int>(numeric), 0, 20);
    }
    else if ((key == "NODESET_G" || key == "NODESET_CONDUCTANCE") && hasNumeric && numeric > 0.0) {
        settings.nodeset_conductance = numeric;
    }
    else if (key == "DAE_AUDIT" || key == "CHARGE_AUDIT" || key == "QJAC_AUDIT") {
        settings.dae_audit = value.empty() ? true : truthy(value);
    }
    else if ((key == "DAE_AUDIT_TOL" || key == "CHARGE_AUDIT_TOL") && hasNumeric && numeric > 0.0) {
        settings.dae_audit_tolerance = numeric;
    }
    else if (key == "FASTSPICE" || key == "FAST_SPICE" || key == "EVENT_DRIVEN") {
        settings.fastspice = value.empty() ? true : truthy(value);
    }
    else if (key == "MULTIRATE" || key == "MULTI_RATE" || key == "MULTI_TIMESTEP") {
        settings.multirate = value.empty() ? true : truthy(value);
    }
    else if (key == "PARALLEL_SOLVE" || key == "PARALLEL_BTF" || key == "PARALLEL_SOLVER") {
        settings.parallel_solve = value.empty() ? true : truthy(value);
    }
    else if (key == "TICER" || key == "MOR" || key == "INTERCONNECT_REDUCTION") {
        settings.ticer = value.empty() ? true : truthy(value);
    }
    else if ((key == "TICER_FMAX" || key == "TICERFMAX") && hasNumeric && numeric > 0.0) {
        settings.ticer_fmax = numeric;
    }
    else if ((key == "METHOD" || key == "TRAN_METHOD") && !value.empty()) {
        std::string method = toUpperCopy(value);
        method.erase(std::remove(method.begin(), method.end(), '_'), method.end());
        method.erase(std::remove(method.begin(), method.end(), '-'), method.end());
        if (method == "AUTO" || method == "BE" || method == "BACKWARDEULER" ||
            method == "TRAP" || method == "TRAPEZOIDAL" || method == "GEAR2" ||
            method == "GEAR" || method == "BDF" || method == "ADAMS" ||
            method == "ADAMSMOULTON" || method == "AM") {
            settings.tran_method = method;
        }
    }
    else if ((key == "ITL1" || key == "OP_MAX_ITER") && hasNumeric && numeric > 0.0) settings.op_max_iter = static_cast<int>(numeric);
    else if ((key == "ITL4" || key == "TRAN_MAX_ITER") && hasNumeric && numeric > 0.0) settings.tran_max_iter = static_cast<int>(numeric);
    else if ((key == "TEMP" || key == "TEMPERATURE") && hasNumeric) settings.temperature_c = numeric;
}

void parseSourceSpec(
    const std::string& sourceSpec,
    double& dcValue,
    double& acMagnitude,
    double& acPhaseDeg,
    bool& dcSeen,
    gspice::VoltageSource::WaveformType& waveformType,
    gspice::VoltageSource::PulseParams& pulse,
    gspice::VoltageSource::SinParams& sin,
    std::vector<double>& pwlTimes,
    std::vector<double>& pwlValues) {

    dcValue = 0.0;
    acMagnitude = 1.0;
    acPhaseDeg = 0.0;
    dcSeen = false;
    waveformType = gspice::VoltageSource::WaveformType::DC;
    pulse = gspice::VoltageSource::PulseParams{};
    sin = gspice::VoltageSource::SinParams{};
    pwlTimes.clear();
    pwlValues.clear();

    std::string scan = sourceSpec;
    for (char& c : scan) if (c == ',') c = ' ';
    std::stringstream ss(scan);
    std::vector<std::string> parts;
    std::string part;
    while (ss >> part) parts.push_back(part);

    for (size_t i = 0; i < parts.size(); ++i) {
        std::string pUpper = toUpperCopy(parts[i]);
        if (pUpper == "DC" && i + 1 < parts.size()) {
            double tmp = 0.0;
            if (tryParseSpiceValue(parts[i + 1], tmp)) {
                dcValue = tmp;
                dcSeen = true;
            }
            ++i;
            continue;
        }
        if (pUpper == "AC" && i + 1 < parts.size()) {
            double tmp = 0.0;
            if (tryParseSpiceValue(parts[i + 1], tmp)) acMagnitude = tmp;
            if (i + 2 < parts.size() && tryParseSpiceValue(parts[i + 2], tmp)) {
                acPhaseDeg = tmp;
                i += 2;
            } else {
                ++i;
            }
            continue;
        }
        if (!dcSeen && i == 0 && pUpper.find('(') == std::string::npos) {
            double tmp = 0.0;
            if (tryParseSpiceValue(parts[i], tmp)) {
                dcValue = tmp;
                dcSeen = true;
            }
        }
    }

    std::string payload;
    if (extractParenPayload(sourceSpec, "PULSE", payload)) {
        auto nums = parseParenNumberList(payload);
        if (nums.size() >= 2) {
            waveformType = gspice::VoltageSource::WaveformType::PULSE;
            pulse.v1 = nums[0];
            pulse.v2 = nums[1];
            if (nums.size() > 2) pulse.td = nums[2];
            if (nums.size() > 3) pulse.tr = nums[3];
            if (nums.size() > 4) pulse.tf = nums[4];
            if (nums.size() > 5) pulse.pw = nums[5];
            if (nums.size() > 6) pulse.per = nums[6];
        }
        return;
    }

    if (extractParenPayload(sourceSpec, "SIN", payload)) {
        auto nums = parseParenNumberList(payload);
        if (nums.size() >= 2) {
            waveformType = gspice::VoltageSource::WaveformType::SIN;
            sin.vo = nums[0];
            sin.va = nums[1];
            if (nums.size() > 2) sin.freq = nums[2];
            if (nums.size() > 3) sin.td = nums[3];
            if (nums.size() > 4) sin.theta = nums[4];
            if (nums.size() > 5) sin.phase_deg = nums[5];
        }
        return;
    }

    if (extractParenPayload(sourceSpec, "PWL", payload)) {
        auto nums = parseParenNumberList(payload);
        if (nums.size() >= 4) {
            waveformType = gspice::VoltageSource::WaveformType::PWL;
            for (size_t i = 0; i + 1 < nums.size(); i += 2) {
                pwlTimes.push_back(nums[i]);
                pwlValues.push_back(nums[i + 1]);
            }
        }
        return;
    }
}

bool applyMeasureParam(gspice::MeasureSpec& measure, const std::string& keyIn, const std::string& value) {
    std::string key = toUpperCopy(keyIn);
    if (key == "AT") {
        measure.at = gspice::Utils::parseValue(value);
        measure.has_at = true;
        return true;
    }
    if (key == "FROM") {
        measure.from = gspice::Utils::parseValue(value);
        measure.has_from = true;
        return true;
    }
    if (key == "TO") {
        measure.to = gspice::Utils::parseValue(value);
        measure.has_to = true;
        return true;
    }
    if (key == "VAL" || key == "VALUE") {
        measure.when_value = gspice::Utils::parseValue(value);
        measure.has_when_value = true;
        return true;
    }
    if (key == "RISE" || key == "FALL" || key == "CROSS") {
        measure.crossing = key;
        if (!value.empty() && toUpperCopy(value) != "LAST") {
            measure.crossing_count = std::max(1, static_cast<int>(gspice::Utils::parseValue(value)));
        }
        return true;
    }
    return false;
}

void parseMeasureParams(
    gspice::MeasureSpec& measure,
    const std::vector<std::string>& tokens,
    size_t startIdx,
    size_t endIdx) {
    endIdx = std::min(endIdx, tokens.size());
    for (size_t i = startIdx; i < endIdx; ++i) {
        auto [key, value] = splitParameterToken(tokens[i]);
        if (!key.empty()) {
            applyMeasureParam(measure, key, value);
            continue;
        }
        std::string upper = toUpperCopy(tokens[i]);
        if ((upper == "AT" || upper == "FROM" || upper == "TO" || upper == "VAL" ||
             upper == "VALUE" || upper == "RISE" || upper == "FALL" || upper == "CROSS") &&
            i + 1 < endIdx) {
            if (tokens[i + 1] == "=" && i + 2 < endIdx) {
                applyMeasureParam(measure, upper, tokens[i + 2]);
                i += 2;
            } else {
                applyMeasureParam(measure, upper, tokens[i + 1]);
                ++i;
            }
        }
    }
}

void parseMeasureParams(gspice::MeasureSpec& measure, const std::vector<std::string>& tokens, size_t startIdx) {
    parseMeasureParams(measure, tokens, startIdx, tokens.size());
}

std::string stripBehavioralExpressionDelimiters(std::string text) {
    text = trimCopy(stripQuotes(text));
    if (text.size() >= 2 && text.front() == '{' && text.back() == '}') {
        text = trimCopy(text.substr(1, text.size() - 2));
    }
    return text;
}

bool parseBehavioralSourceSpec(
    const std::vector<std::string>& tokens,
    size_t startIdx,
    gspice::BehavioralSource::Mode& mode,
    std::string& expression) {
    if (startIdx >= tokens.size()) return false;
    std::string spec = joinTokens(tokens, startIdx);
    spec = trimCopy(spec);
    if (spec.empty()) return false;
    auto upperSpec = toUpperCopy(spec);
    if (upperSpec.rfind("I=", 0) == 0) {
        mode = gspice::BehavioralSource::Mode::Current;
        expression = stripBehavioralExpressionDelimiters(spec.substr(2));
        return !expression.empty();
    }
    if (upperSpec.rfind("V=", 0) == 0) {
        mode = gspice::BehavioralSource::Mode::Voltage;
        expression = stripBehavioralExpressionDelimiters(spec.substr(2));
        return !expression.empty();
    }
    if (tokens.size() > startIdx + 2) {
        const std::string key = toUpperCopy(tokens[startIdx]);
        if ((key == "I" || key == "V") && tokens[startIdx + 1] == "=") {
            mode = key == "I" ? gspice::BehavioralSource::Mode::Current : gspice::BehavioralSource::Mode::Voltage;
            expression = stripBehavioralExpressionDelimiters(joinTokens(tokens, startIdx + 2));
            return !expression.empty();
        }
    }
    return false;
}

} // namespace

namespace gspice {

std::vector<std::string> Parser::tokenize(const std::string& line) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    while (ss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

Netlist Parser::parse(const std::string& filePath) {
    Netlist netlist;
    std::vector<std::string> preprocessErrors;
    std::set<std::string> includeStack;
    auto rawLines = loadNetlistFile(std::filesystem::path(filePath), true, preprocessErrors, includeStack);
    std::vector<PreLine> topLevelLines;
    auto subckts = collectSubckts(rawLines, topLevelLines, preprocessErrors);
    auto expandedLines = expandLines(topLevelLines, subckts, preprocessErrors, 0);
    std::vector<PreLine> paramFreeLines;
    auto globalParams = collectGlobalParams(expandedLines, paramFreeLines);
    paramFreeLines = filterConditionalLines(paramFreeLines, globalParams, preprocessErrors, "top-level deck");
    for (auto& line : paramFreeLines) {
        line.text = applyGlobalParams(line.text, globalParams);
    }
    for (const auto& error : preprocessErrors) {
        netlist.addError(error);
    }

    for (const auto& preLine : paramFreeLines) {
        std::string line = preLine.text;
        int lineNo = preLine.lineNo;

        auto tokens = tokenize(line);
        if (tokens.empty()) continue;

        std::string firstTokenUpper = toUpperCopy(tokens[0]);
        if (firstTokenUpper == "LOAD" || firstTokenUpper == "PRE_LOAD") {
            netlist.addError("Line " + std::to_string(lineNo) +
                             ": OSDI/OpenVAF-style load directives are disabled; use .GSDI with GMC-generated native models.");
            continue;
        }

        char firstChar = std::toupper(tokens[0][0]);
        
        if (firstChar == '.') {
            // Commands
            std::string cmd = tokens[0];
            std::transform(cmd.begin(), cmd.end(), cmd.begin(), ::toupper);
            
            if (cmd == ".END") {
                break;
            } else if (cmd == ".GSPICEWARN") {
                if (verboseCompatWarnings()) netlist.addWarning(joinTokens(tokens, 1));
            } else if (cmd == ".OP") {
                if (netlist.getSettings().type == "OP") {
                    SimulationSettings settings = netlist.getSettings();
                    settings.type = "OP";
                    netlist.setSettings(settings);
                } else {
                    netlist.addWarning(
                        "Line " + std::to_string(lineNo) +
                        ": .OP kept as operating-point initialization; active analysis remains ." +
                        netlist.getSettings().type);
                }
            } else if (cmd == ".OPTIONS" || cmd == ".OPTION" || cmd == ".OPT") {
                SimulationSettings settings = netlist.getSettings();
                for (size_t i = 1; i < tokens.size(); ++i) {
                    applyOptionToken(settings, tokens[i]);
                }
                netlist.setSettings(settings);
            } else if (cmd == ".TEMP" || cmd == ".TEMPERATURE") {
                if (tokens.size() < 2) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .TEMP line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.temperature_c = Utils::parseValue(tokens[1]);
                netlist.setSettings(settings);
            } else if (cmd == ".IC") {
                SimulationSettings settings = netlist.getSettings();
                for (size_t i = 1; i < tokens.size(); ++i) {
                    std::string node;
                    double value = 0.0;
                    if (!parseInitialConditionToken(tokens[i], node, value)) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .IC token ignored: " + tokens[i]);
                        continue;
                    }
                    settings.initial_conditions.push_back({netlist.getOrCreateNode(node), value});
                }
                netlist.setSettings(settings);
            } else if (cmd == ".NODESET") {
                SimulationSettings settings = netlist.getSettings();
                for (size_t i = 1; i < tokens.size(); ++i) {
                    std::string node;
                    double value = 0.0;
                    if (!parseInitialConditionToken(tokens[i], node, value)) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .NODESET token ignored: " + tokens[i]);
                        continue;
                    }
                    settings.nodesets.push_back({netlist.getOrCreateNode(node), value});
                }
                netlist.setSettings(settings);
            } else if (cmd == ".GLOBAL") {
                netlist.addWarning(
                    "Line " + std::to_string(lineNo) +
                    ": .GLOBAL accepted for compatibility; named nodes are already global in flat GSPICE decks.");
            } else if (cmd == ".SAVE" || cmd == ".PROBE" || cmd == ".PRINT" || cmd == ".PLOT") {
                SimulationSettings settings = netlist.getSettings();
                bool sawSave = false;
                for (size_t i = 1; i < tokens.size(); ++i) {
                    const std::string tokenUpper = toUpperCopy(tokens[i]);
                    if (tokenUpper == "ALL" || tokenUpper == "V(*)" || tokenUpper == "V(ALL)") {
                        settings.save_all = true;
                        settings.save_none = false;
                        settings.saves.clear();
                        sawSave = true;
                        continue;
                    }
                    SaveSpec save;
                    if (!parseSaveToken(tokens[i], save)) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid " + cmd +
                                           " token ignored: " + tokens[i]);
                        continue;
                    }
                    if (settings.save_all && toUpperCopy(save.kind) != "I") {
                        settings.save_all = false;
                        settings.save_none = false;
                        settings.saves.clear();
                    }
                    settings.saves.push_back(save);
                    sawSave = true;
                }
                if (!sawSave) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": empty " + cmd + " ignored.");
                }
                netlist.setSettings(settings);
            } else if (cmd == ".LIB" || cmd == ".INCLUDE" || cmd == ".INC") {
                netlist.addError("Line " + std::to_string(lineNo) + ": unresolved include/library directive: " + line);
            } else if (cmd == ".DC") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .DC line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "DC";
                settings.dc_sweeps.clear();
                for (size_t idx = 1; idx + 3 < tokens.size(); idx += 4) {
                    SweepSpec sweep;
                    sweep.source = tokens[idx];
                    sweep.start = Utils::parseValue(tokens[idx + 1]);
                    sweep.stop = Utils::parseValue(tokens[idx + 2]);
                    sweep.step = Utils::parseValue(tokens[idx + 3]);
                    if (sweep.step == 0.0) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": .DC sweep step cannot be zero: " + line);
                        continue;
                    }
                    settings.dc_sweeps.push_back(sweep);
                }
                if (settings.dc_sweeps.empty()) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .DC sweep groups ignored: " + line);
                    continue;
                }
                settings.dc_sweep_source = tokens[1];
                settings.dc_start = Utils::parseValue(tokens[2]);
                settings.dc_stop = Utils::parseValue(tokens[3]);
                settings.dc_step = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".STEP") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .STEP line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "STEP";
                settings.step_sweeps.clear();
                for (size_t idx = 1; idx + 3 < tokens.size(); idx += 4) {
                    SweepSpec sweep;
                    sweep.source = tokens[idx];
                    sweep.start = Utils::parseValue(tokens[idx + 1]);
                    sweep.stop = Utils::parseValue(tokens[idx + 2]);
                    sweep.step = Utils::parseValue(tokens[idx + 3]);
                    if (sweep.step == 0.0) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": .STEP step cannot be zero: " + line);
                        continue;
                    }
                    settings.step_sweeps.push_back(sweep);
                }
                netlist.setSettings(settings);
            } else if (cmd == ".MC" || cmd == ".MONTE" || cmd == ".MONTECARLO") {
                if (tokens.size() < 4) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid Monte Carlo line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "MC";
                settings.mc_runs = std::max(1, std::stoi(tokens[1]));
                settings.mc_source = tokens[2];
                std::string spec = joinTokens(tokens, 3);
                std::string payload;
                if (extractParenPayload(spec, "GAUSS", payload) || extractParenPayload(spec, "NORMAL", payload)) {
                    auto nums = parseParenNumberList(payload);
                    if (nums.size() >= 2) {
                        settings.mc_distribution = "GAUSSIAN";
                        settings.mc_mean = nums[0];
                        settings.mc_sigma = nums[1];
                    } else {
                        netlist.addError("Line " + std::to_string(lineNo) + ": Gaussian .MC requires mean and sigma: " + line);
                    }
                } else if (extractParenPayload(spec, "UNIFORM", payload) ||
                           extractParenPayload(spec, "UNIF", payload)) {
                    auto nums = parseParenNumberList(payload);
                    if (nums.size() >= 2) {
                        settings.mc_distribution = "UNIFORM";
                        settings.mc_lower = nums[0];
                        settings.mc_upper = nums[1];
                        if (settings.mc_upper < settings.mc_lower) {
                            std::swap(settings.mc_lower, settings.mc_upper);
                        }
                    } else {
                        netlist.addError("Line " + std::to_string(lineNo) + ": uniform .MC requires lower and upper bounds: " + line);
                    }
                } else {
                    settings.mc_distribution = "GAUSSIAN";
                    settings.mc_mean = Utils::parseValue(tokens[3]);
                    settings.mc_sigma = tokens.size() > 4 ? Utils::parseValue(tokens[4]) : 0.0;
                }
                for (size_t i = 4; i < tokens.size(); ++i) {
                    auto [key, value] = splitParameterToken(tokens[i]);
                    const std::string option = toUpperCopy(key);
                    if (option == "SEED") {
                        settings.mc_seed = static_cast<unsigned int>(std::stoul(value));
                    } else if (option == "LHS" || option == "LATINHYPERCUBE") {
                        const std::string upperValue = toUpperCopy(value);
                        settings.mc_latin_hypercube = value.empty() ||
                            !(upperValue == "0" || upperValue == "NO" || upperValue == "FALSE" || upperValue == "OFF");
                    }
                }
                netlist.setSettings(settings);
            } else if (cmd == ".CORNER" || cmd == ".CORNERS") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .CORNER line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                if (settings.type == "OP") settings.type = "CORNER";
                CornerSpec corner;
                corner.name = tokens[1];
                for (size_t i = 2; i < tokens.size(); ++i) {
                    auto [key, value] = splitParameterToken(tokens[i]);
                    if (key.empty()) {
                        netlist.addWarning(
                            "Line " + std::to_string(lineNo) +
                            ": ignoring malformed corner assignment '" + tokens[i] + "'");
                        continue;
                    }
                    corner.source_values.push_back({key, Utils::parseValue(value)});
                }
                if (corner.source_values.empty()) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": .CORNER has no source assignments: " + line);
                    continue;
                }
                settings.corners.push_back(corner);
                netlist.setSettings(settings);
            } else if (cmd == ".SPEC" || cmd == ".YIELD") {
                if (tokens.size() < 4) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .SPEC line ignored: " + line);
                    continue;
                }
                std::string outPos;
                std::string outNeg;
                if (!parseVoltageProbeToken(tokens[2], outPos, outNeg)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": only voltage .SPEC outputs are supported currently: " + line);
                    continue;
                }
                OutputSpec spec;
                spec.name = tokens[1];
                spec.node_pos = netlist.getOrCreateNode(outPos);
                spec.node_neg = netlist.getOrCreateNode(outNeg);
                for (size_t i = 3; i < tokens.size(); ++i) {
                    auto [key, value] = splitParameterToken(tokens[i]);
                    key = toUpperCopy(key);
                    if (key == "MIN" || key == "LOW") {
                        spec.min_value = Utils::parseValue(value);
                        spec.has_min = true;
                    } else if (key == "MAX" || key == "HIGH") {
                        spec.max_value = Utils::parseValue(value);
                        spec.has_max = true;
                    }
                }
                if (!spec.has_min && !spec.has_max) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": .SPEC has neither MIN nor MAX: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.output_specs.push_back(spec);
                netlist.setSettings(settings);
            } else if (cmd == ".TRAN") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .TRAN line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "TRAN";
                settings.t_step = Utils::parseValue(tokens[1]);
                settings.t_stop = Utils::parseValue(tokens[2]);
                if (settings.t_stop <= 0.0) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": .TRAN stop time must be positive (token '" +
                                       tokens[2] + "' did not parse as a time): " + line);
                    continue;
                }
                if (tokens.size() > 3) {
                    std::string opt = tokens[3];
                    std::transform(opt.begin(), opt.end(), opt.begin(), ::toupper);
                    if (opt == "UIC") {
                        settings.use_uic = true;
                    } else {
                        // SPICE-compatible .TRAN tstep tstop [tstart [tmax]].
                        settings.t_start = Utils::parseValue(tokens[3]);
                    }
                }
                if (tokens.size() > 4) {
                    std::string opt = tokens[4];
                    std::transform(opt.begin(), opt.end(), opt.begin(), ::toupper);
                    if (opt == "UIC") {
                        settings.use_uic = true;
                    } else if (!settings.tran_max_step_auto) {
                        settings.t_max_step = Utils::parseValue(tokens[4]);
                    }
                }
                if (tokens.size() > 5) {
                    std::string opt = tokens[5];
                    std::transform(opt.begin(), opt.end(), opt.begin(), ::toupper);
                    if (opt == "UIC") settings.use_uic = true;
                }
                netlist.setSettings(settings);
            } else if (cmd == ".TRANNOISE" || cmd == ".TRNOISE") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .TRANNOISE line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "TRAN";
                settings.transient_noise = true;
                settings.t_step = Utils::parseValue(tokens[1]);
                settings.t_stop = Utils::parseValue(tokens[2]);
                for (size_t i = 3; i < tokens.size(); ++i) {
                    auto [key, value] = splitParameterToken(tokens[i]);
                    key = toUpperCopy(key);
                    const size_t rawEq = tokens[i].find('=');
                    if (key == "SEED" && !value.empty()) {
                        settings.transient_noise_seed = static_cast<unsigned int>(std::stoul(value));
                    } else if (key == "SCALE" && !value.empty()) {
                        settings.transient_noise_scale = Utils::parseValue(value);
                    } else if (key == "FMAX" && !value.empty()) {
                        settings.transient_noise_fmax = Utils::parseValue(value);
                    } else if ((key == "NOISEMODE" || key == "MODE") && !value.empty()) {
                        const std::string rawValue = rawEq == std::string::npos ? value : stripQuotes(tokens[i].substr(rawEq + 1));
                        settings.transient_noise_mode = toUpperCopy(rawValue);
                    }
                }
                netlist.setSettings(settings);
            } else if (cmd == ".IC" || cmd == ".NODESET") {
                SimulationSettings settings = netlist.getSettings();
                for (size_t i = 1; i < tokens.size(); ++i) {
                    auto [key, val_str] = splitParameterToken(tokens[i]);
                    if (key.empty() && i + 2 < tokens.size() && tokens[i + 1] == "=") {
                        key = tokens[i];
                        val_str = tokens[i + 2];
                        i += 2;
                    }
                    if (!key.empty()) {
                        std::string node_name = key;
                        if (node_name.rfind("V(", 0) == 0 && node_name.back() == ')') {
                            node_name = node_name.substr(2, node_name.size() - 3);
                        }
                        int node_idx = netlist.getOrCreateNode(node_name);
                        double val = Utils::parseValue(val_str);
                        InitialConditionSpec ic;
                        ic.node = node_idx;
                        ic.value = val;
                        settings.initial_conditions.push_back(ic);
                    }
                }
                netlist.setSettings(settings);
            } else if (cmd == ".LOAD" || cmd == ".OSDI" || cmd == ".PRE_OSDI") {
                netlist.addError("Line " + std::to_string(lineNo) +
                                 ": external compact-model runtime loading is disabled; use .GSDI with GMC-generated native models.");
            } else if (cmd == ".GSDI" || cmd == ".PRE_GSDI") {
                if (tokens.size() < 2) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": .GSDI requires a .gsdi artifact path.");
                } else if (!hasExtension(tokens[1], ".gsdi")) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": .GSDI only accepts .gsdi artifacts, not " + tokens[1]);
                } else {
                    const auto artifactPath = resolveRelativePath(std::filesystem::path(preLine.source), tokens[1]);
                    const std::string artifact = readWholeFile(artifactPath);
                    if (artifact.empty()) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": could not read .GSDI artifact: " + artifactPath.string());
                        continue;
                    }
                    const std::string schema = extractJsonString(artifact, "schema");
                    const std::string modelType = extractJsonString(artifact, "model_type");
                    if (schema != "gsdi-artifact-v0" || modelType.empty()) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": invalid .GSDI artifact envelope: " + artifactPath.string());
                        continue;
                    }
                    GsdiArtifactInfo info;
                    info.model_type = modelType;
                    info.path = artifactPath.string();
                    info.terminal_count = extractJsonInt(artifact, "terminal_count");
                    info.parameter_names = extractGsdiParameterNames(artifact);
                    netlist.addGsdiArtifact(info);
                    netlist.addModelStatus("GSDI_LOADED: " + modelType + " from " + artifactPath.filename().string());
                }
            } else if (cmd == ".MODEL") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .MODEL line ignored: " + line);
                    continue;
                }
                ModelCard model;
                model.name = tokens[1];
                std::string typeToken = tokens[2];
                size_t paren = typeToken.find('(');
                if (paren != std::string::npos) {
                    model.type = typeToken.substr(0, paren);
                    std::string paramText = typeToken.substr(paren + 1);
                    if (tokens.size() > 3) paramText += " " + joinTokens(tokens, 3);
                    model.params = parseModelParamsFromText(paramText);
                } else {
                    model.type = typeToken;
                    std::string paramText = tokens.size() > 3 ? joinTokens(tokens, 3) : "";
                    model.params = parseModelParamsFromText(paramText);
                }
                auto modelContext = globalParams;
                for (const auto& [key, value] : model.params) {
                    modelContext[key] = value;
                }
                resolveParameterMapExpressions(modelContext);
                for (auto& [key, value] : model.params) {
                    double evaluated = 0.0;
                    if (tryEvaluateParamExpression(value, modelContext, evaluated)) {
                        value = formatNumericValue(evaluated);
                    }
                }
                if (isPsp103ModelCard(model)) {
                    addBuiltinPspGsdiArtifact(netlist, psp103GsdiType(model));
                    const std::string ignored = verboseCompatWarnings()
                        ? psp103IgnoredParameterSummary(model)
                        : "";
                    if (!ignored.empty()) {
                        netlist.addWarning(
                            "Line " + std::to_string(lineNo) + ": PSP model '" +
                            model.name + "' has " + ignored);
                    }
                }
#ifdef GSPICE_HAVE_GMC_GENERATED
                else {
                    addBuiltinGmcGeneratedGsdiArtifact(netlist, model.type);
                }
#endif
                netlist.addModelCard(model);
            } else if (cmd == ".AC") {
                if (tokens.size() < 2) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .AC line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "AC";
                settings.f_sweep_type = toUpperCopy(tokens[1]);
                settings.f_values.clear();
                if (settings.f_sweep_type == "VALUES") {
                    if (tokens.size() < 3) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .AC VALUES line ignored: " + line);
                        continue;
                    }
                    for (size_t j = 2; j < tokens.size(); ++j) settings.f_values.push_back(Utils::parseValue(tokens[j]));
                    settings.points_per_dec = static_cast<int>(settings.f_values.size());
                    settings.f_start = settings.f_values.front();
                    settings.f_stop = settings.f_values.back();
                } else {
                    if (tokens.size() < 5) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .AC line ignored: " + line);
                        continue;
                    }
                    settings.points_per_dec = std::stoi(tokens[2]);
                    settings.f_start = Utils::parseValue(tokens[3]);
                    settings.f_stop = Utils::parseValue(tokens[4]);
                }
                netlist.setSettings(settings);
            } else if (cmd == ".ACXF" || cmd == ".DCXF" || cmd == ".DCINC") {
                SimulationSettings settings = netlist.getSettings();
                settings.type = cmd.substr(1);
                if (cmd == ".DCINC") {
                    settings.f_sweep_type = "VALUES";
                    settings.f_values = {0.0};
                    settings.points_per_dec = 1;
                    settings.f_start = 0.0;
                    settings.f_stop = 0.0;
                    netlist.setSettings(settings);
                    continue;
                }
                if (tokens.size() >= 2) {
                    std::string outPos;
                    std::string outNeg;
                    if (parseVoltageProbeToken(tokens[1], outPos, outNeg)) {
                        settings.xf_out_pos = netlist.getOrCreateNode(outPos);
                        settings.xf_out_neg = netlist.getOrCreateNode(outNeg);
                    }
                }
                size_t sweep = (settings.xf_out_pos >= 0) ? 2 : 1;
                if (cmd == ".DCXF") {
                    settings.f_sweep_type = "VALUES";
                    settings.f_values = {0.0};
                    settings.points_per_dec = 1;
                    settings.f_start = 0.0;
                    settings.f_stop = 0.0;
                } else if (tokens.size() > sweep) {
                    settings.f_sweep_type = toUpperCopy(tokens[sweep]);
                    if (settings.f_sweep_type == "VALUES") {
                        for (size_t j = sweep + 1; j < tokens.size(); ++j) {
                            settings.f_values.push_back(Utils::parseValue(tokens[j]));
                        }
                        settings.points_per_dec = static_cast<int>(settings.f_values.size());
                        if (!settings.f_values.empty()) {
                            settings.f_start = settings.f_values.front();
                            settings.f_stop = settings.f_values.back();
                        }
                    } else if (tokens.size() > sweep + 3) {
                        settings.points_per_dec = std::stoi(tokens[sweep + 1]);
                        settings.f_start = Utils::parseValue(tokens[sweep + 2]);
                        settings.f_stop = Utils::parseValue(tokens[sweep + 3]);
                    }
                }
                netlist.setSettings(settings);
            } else if (cmd == ".PSS" || cmd == ".HB") {
                SimulationSettings settings = netlist.getSettings();
                settings.type = cmd.substr(1);
                if (cmd == ".PSS") settings.pss_requested = true;
                size_t i = 1;
                while (i < tokens.size()) {
                    const std::string tokenUpper = toUpperCopy(tokens[i]);
                    const size_t eq = tokens[i].find('=');
                    if (tokenUpper == "DRIVEN" || tokenUpper == "OSCILLATOR" ||
                        tokenUpper == "OSCILLATOR=YES" || tokenUpper == "USE_INITIAL_CONDITIONS=YES" ||
                        tokenUpper == "PSS_ADAPTIVE=YES" || tokenUpper == "PSS_CONTINUATION=YES") {
                        ++i;
                        continue;
                    }
                    if (eq != std::string::npos) {
                        const std::string key = toUpperCopy(tokens[i].substr(0, eq));
                        const std::string value = tokens[i].substr(eq + 1);
                        if (key == "TSTAB") settings.pss_tstab = Utils::parseValue(value);
                        else if (key == "TSTAB_PERIODS") settings.pss_tstab_periods = std::max(0, std::stoi(value));
                        else if (key == "PSS_RESIDUAL_GOAL" || key == "RESIDUAL_GOAL" ||
                                 key == "TOL" || key == "RELTOL") {
                            settings.pss_residual_goal = Utils::parseValue(value);
                        }
                        else if (key == "MAX_PSS_ITER" || key == "PSS_MAX_ITER" ||
                                 key == "MAXITER" || key == "MAX_ITERS" ||
                                 key == "HB_MAXITER" || key == "HB_MAX_ITER") {
                            settings.max_pss_iter = std::max(1, std::stoi(value));
                        }
                        else if (key == "PSS_CONTINUATION_STEPS") settings.max_pss_iter = std::max(settings.max_pss_iter, std::stoi(value) * 2);
                        else applyRfAnalysisOption(settings, key, value);
                    } else if (!tokens[i].empty() && std::all_of(tokens[i].begin(), tokens[i].end(), ::isdigit)) {
                        settings.n_harms = std::stoi(tokens[i]);
                    } else if (settings.f_fund.size() < 4) {
                        try {
                            settings.f_fund.push_back(Utils::parseValue(tokens[i]));
                        } catch (const std::exception&) {
                            netlist.addWarning("Line " + std::to_string(lineNo) + ": ignored unsupported " + cmd + " option: " + tokens[i]);
                        }
                    } else {
                        std::cerr << "Warning: GSPICE supports max 4 tones. Ignoring extra: " << tokens[i] << std::endl;
                    }
                    i++;
                }
                if (settings.f_fund.empty()) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": " + cmd + " requires at least 1 fundamental frequency.");
                }
                netlist.setSettings(settings);
            }
 else if (cmd == ".SP") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .SP line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "SP";
                settings.f_sweep_type = toUpperCopy(tokens[1]);
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".NOISE") {
                if (tokens.size() < 6) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .NOISE line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "NOISE";
                // tokens[1] is V(node)
                std::string outNode = tokens[1].substr(2, tokens[1].size()-3);
                settings.out_node = netlist.getOrCreateNode(outNode);
                size_t sweepIdx = 3;
                std::string sweepType = toUpperCopy(tokens[sweepIdx]);
                if (sweepType == "DEC" || sweepType == "OCT" || sweepType == "LIN") {
                    if (tokens.size() < 7) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .NOISE sweep line ignored: " + line);
                        continue;
                    }
                    settings.points_per_dec = std::stoi(tokens[4]);
                    settings.f_start = Utils::parseValue(tokens[5]);
                    settings.f_stop = Utils::parseValue(tokens[6]);
                } else {
                    settings.points_per_dec = std::stoi(tokens[3]);
                    settings.f_start = Utils::parseValue(tokens[4]);
                    settings.f_stop = Utils::parseValue(tokens[5]);
                }
                netlist.setSettings(settings);
            } else if (cmd == ".TF") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .TF line ignored: " + line);
                    continue;
                }
                std::string outPos;
                std::string outNeg;
                if (!parseVoltageProbeToken(tokens[1], outPos, outNeg)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": only voltage-output .TF is supported currently: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "TF";
                settings.tf_out_pos = netlist.getOrCreateNode(outPos);
                settings.tf_out_neg = netlist.getOrCreateNode(outNeg);
                settings.tf_input_source = tokens[2];
                netlist.setSettings(settings);
            } else if (cmd == ".SENS") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .SENS line ignored: " + line);
                    continue;
                }
                std::string outPos;
                std::string outNeg;
                if (!parseVoltageProbeToken(tokens[1], outPos, outNeg)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": only voltage-output .SENS is supported currently: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "SENS";
                settings.sens_out_pos = netlist.getOrCreateNode(outPos);
                settings.sens_out_neg = netlist.getOrCreateNode(outNeg);
                settings.sens_source = tokens[2];
                netlist.setSettings(settings);
            } else if (cmd == ".PZ") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .PZ line ignored: " + line);
                    continue;
                }
                std::string outPos;
                std::string outNeg;
                if (!parseVoltageProbeToken(tokens[1], outPos, outNeg)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": only voltage-output .PZ is supported currently: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "PZ";
                settings.tf_out_pos = netlist.getOrCreateNode(outPos);
                settings.tf_out_neg = netlist.getOrCreateNode(outNeg);
                settings.tf_input_source = tokens[2];
                settings.points_per_dec = 20;
                settings.f_start = 1.0;
                settings.f_stop = 1e12;
                if (tokens.size() >= 7) {
                    settings.points_per_dec = std::stoi(tokens[4]);
                    settings.f_start = Utils::parseValue(tokens[5]);
                    settings.f_stop = Utils::parseValue(tokens[6]);
                }
                netlist.setSettings(settings);
            } else if (cmd == ".MEAS" || cmd == ".MEASURE") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .MEASURE line ignored: " + line);
                    continue;
                }
                std::string outPos;
                std::string outNeg;
                MeasureSpec measure;
                measure.analysis = toUpperCopy(tokens[1]);
                measure.name = tokens[2];
                measure.op = toUpperCopy(tokens[3]);
                auto parseProbeIntoMeasure = [&](MeasureSpec& spec, const std::string& token) -> bool {
                    std::string probe = token;
                    auto [probeKeyLocal, probeValueLocal] = splitParameterToken(probe);
                    if (!probeKeyLocal.empty()) {
                        probe = probeKeyLocal;
                        spec.when_value = Utils::parseValue(probeValueLocal);
                        spec.has_when_value = true;
                    }
                    if (parseVoltageProbeToken(probe, outPos, outNeg)) {
                        spec.kind = "V";
                        spec.node_pos = netlist.getOrCreateNode(outPos);
                        spec.node_neg = netlist.getOrCreateNode(outNeg);
                        return true;
                    }
                    const std::string expr = stripQuotes(trimCopy(probe));
                    if (toUpperCopy(expr).rfind("I(", 0) != 0 || expr.back() != ')') return false;
                    spec.kind = "I";
                    spec.device_name = trimCopy(expr.substr(2, expr.size() - 3));
                    return true;
                };
                if (measure.op == "TRIG" || measure.op == "TRIGGER") {
                    size_t targetIdx = tokens.size();
                    for (size_t i = 5; i < tokens.size(); ++i) {
                        const std::string upper = toUpperCopy(tokens[i]);
                        if (upper == "TARG" || upper == "TARGET") {
                            targetIdx = i;
                            break;
                        }
                    }
                    if (targetIdx == tokens.size() || targetIdx + 1 >= tokens.size()) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": .MEASURE TRIG requires TARG probe: " + line);
                        continue;
                    }
                    if (!parseProbeIntoMeasure(measure, tokens[4])) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": unsupported .MEASURE TRIG expression: " + line);
                        continue;
                    }
                    parseMeasureParams(measure, tokens, 5, targetIdx);
                    MeasureSpec target;
                    if (!parseProbeIntoMeasure(target, tokens[targetIdx + 1])) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": unsupported .MEASURE TARG expression: " + line);
                        continue;
                    }
                    parseMeasureParams(target, tokens, targetIdx + 2, tokens.size());
                    if (!measure.has_when_value || !target.has_when_value) {
                        netlist.addError("Line " + std::to_string(lineNo) + ": .MEASURE TRIG/TARG requires VAL for both probes: " + line);
                        continue;
                    }
                    measure.op = "DELAY";
                    measure.has_target = true;
                    measure.target_kind = target.kind;
                    measure.target_device_name = target.device_name;
                    measure.target_node_pos = target.node_pos;
                    measure.target_node_neg = target.node_neg;
                    measure.target_has_when_value = target.has_when_value;
                    measure.target_when_value = target.when_value;
                    measure.target_crossing = target.crossing;
                    measure.target_crossing_count = target.crossing_count;
                    SimulationSettings settings = netlist.getSettings();
                    settings.measures.push_back(measure);
                    netlist.setSettings(settings);
                    continue;
                }
                std::string probeToken = tokens[4];
                auto [probeKey, probeValue] = splitParameterToken(probeToken);
                if (!probeKey.empty() && measure.op == "WHEN") {
                    probeToken = probeKey;
                    measure.when_value = Utils::parseValue(probeValue);
                    measure.has_when_value = true;
                }
                if (parseVoltageProbeToken(probeToken, outPos, outNeg)) {
                    measure.kind = "V";
                    measure.node_pos = netlist.getOrCreateNode(outPos);
                    measure.node_neg = netlist.getOrCreateNode(outNeg);
                } else {
                    const std::string expr = stripQuotes(trimCopy(probeToken));
                    if (toUpperCopy(expr).rfind("I(", 0) != 0 || expr.back() != ')') {
                        netlist.addError("Line " + std::to_string(lineNo) + ": unsupported .MEASURE expression: " + line);
                        continue;
                    }
                    measure.kind = "I";
                    measure.device_name = trimCopy(expr.substr(2, expr.size() - 3));
                }
                if (measure.op == "WHEN" && !measure.has_when_value && tokens.size() > 6 && tokens[5] == "=") {
                    measure.when_value = Utils::parseValue(tokens[6]);
                    measure.has_when_value = true;
                    parseMeasureParams(measure, tokens, 7);
                } else {
                parseMeasureParams(measure, tokens, 5);
                }
                if (measure.op == "WHEN" && !measure.has_when_value) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": .MEASURE WHEN requires a target value: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.measures.push_back(measure);
                netlist.setSettings(settings);
            } else if (cmd == ".FOUR") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .FOUR line ignored: " + line);
                    continue;
                }
                std::string outPos;
                std::string outNeg;
                if (!parseVoltageProbeToken(tokens[2], outPos, outNeg)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": only voltage-output .FOUR is supported currently: " + line);
                    continue;
                }
                FourSpec four;
                four.frequency = Utils::parseValue(tokens[1]);
                four.node_pos = netlist.getOrCreateNode(outPos);
                four.node_neg = netlist.getOrCreateNode(outNeg);
                for (size_t i = 3; i < tokens.size(); ++i) {
                    auto [key, value] = splitParameterToken(tokens[i]);
                    if (toUpperCopy(key) == "NHARMS" || toUpperCopy(key) == "HARMONICS") {
                        four.harmonics = std::max(1, std::stoi(value));
                    }
                }
                SimulationSettings settings = netlist.getSettings();
                if (settings.type == "OP") settings.type = "TRAN";
                settings.fours.push_back(four);
                netlist.setSettings(settings);
            } else if (cmd == ".STB") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .STB line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "STB";
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".PAC" || cmd == ".PXF" || cmd == ".PSSXF") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid " + cmd + " line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "PAC";
                if (cmd == ".PSSXF") settings.pss_requested = true;
                size_t i = 1;
                std::string outPos;
                std::string outNeg;
                if (i < tokens.size() && parseVoltageProbeToken(tokens[i], outPos, outNeg)) {
                    settings.out_node = netlist.getOrCreateNode(outPos);
                    ++i;
                }
                if (i < tokens.size() && !isRfSweepKeyword(tokens[i]) && tokens[i].find('=') == std::string::npos) {
                    settings.f_fund.clear();
                    settings.f_fund.push_back(Utils::parseValue(tokens[i]));
                    settings.pss_requested = true;
                    ++i;
                }
                if (i < tokens.size() && isRfSweepKeyword(tokens[i])) {
                    settings.f_sweep_type = toUpperCopy(tokens[i++]);
                }
                if (settings.f_sweep_type == "VALUES") {
                    settings.f_values.clear();
                    while (i < tokens.size() && tokens[i].find('=') == std::string::npos) {
                        settings.f_values.push_back(Utils::parseValue(tokens[i++]));
                    }
                    settings.points_per_dec = static_cast<int>(settings.f_values.size());
                    if (!settings.f_values.empty()) {
                        settings.f_start = settings.f_values.front();
                        settings.f_stop = settings.f_values.back();
                    }
                } else if (i + 2 < tokens.size()) {
                    settings.points_per_dec = std::stoi(tokens[i]);
                    settings.f_start = Utils::parseValue(tokens[i + 1]);
                    settings.f_stop = Utils::parseValue(tokens[i + 2]);
                    i += 3;
                }
                parseRfAssignmentOptions(settings, tokens, i);
                netlist.setSettings(settings);
            } else if (cmd == ".PSSPAC") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .PSSPAC line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.pss_requested = true;
                settings.type = "PAC";
                settings.f_fund.clear();
                settings.f_fund.push_back(Utils::parseValue(tokens[1]));
                settings.f_sweep_type = toUpperCopy(tokens[2]);
                settings.points_per_dec = std::stoi(tokens[3]);
                settings.f_start = Utils::parseValue(tokens[4]);
                settings.f_stop = tokens.size() > 5 ? Utils::parseValue(tokens[5]) : settings.f_start;
                parseRfAssignmentOptions(settings, tokens, 6);
                netlist.setSettings(settings);
            } else if (cmd == ".PNOISE") {
                if (tokens.size() < 2) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .PNOISE line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "PNOISE";
                size_t i = 1;
                std::string outPos;
                std::string outNeg;
                if (parseVoltageProbeToken(tokens[i], outPos, outNeg)) {
                    settings.out_node = netlist.getOrCreateNode(outPos);
                    ++i;
                } else {
                    settings.out_node = netlist.getOrCreateNode(tokens[i++]);
                }
                if (i < tokens.size() && !isRfSweepKeyword(tokens[i]) && tokens[i].find('=') == std::string::npos) {
                    settings.pnoise_input_source = tokens[i++];
                }
                if (i >= tokens.size() || tokens[i].find('=') != std::string::npos) {
                    settings.f_sweep_type = "DEC";
                    settings.points_per_dec = 50;
                    settings.f_start = 1.0;
                    settings.f_stop = 100000.0;
                    parseRfAssignmentOptions(settings, tokens, i);
                    if (!settings.f_fund.empty()) settings.pss_requested = true;
                    netlist.setSettings(settings);
                    continue;
                }
                std::string sweepType = toUpperCopy(tokens[i++]);
                settings.f_sweep_type = sweepType;
                if (sweepType == "DEC" || sweepType == "OCT" || sweepType == "LIN") {
                    if (i + 2 >= tokens.size()) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .PNOISE sweep line ignored: " + line);
                        continue;
                    }
                    settings.points_per_dec = std::stoi(tokens[i]);
                    settings.f_start = Utils::parseValue(tokens[i + 1]);
                    settings.f_stop = Utils::parseValue(tokens[i + 2]);
                    i += 3;
                } else if (sweepType == "VALUES") {
                    settings.f_values.clear();
                    while (i < tokens.size() && tokens[i].find('=') == std::string::npos) {
                        settings.f_values.push_back(Utils::parseValue(tokens[i++]));
                    }
                    settings.points_per_dec = static_cast<int>(settings.f_values.size());
                    if (!settings.f_values.empty()) {
                        settings.f_start = settings.f_values.front();
                        settings.f_stop = settings.f_values.back();
                    }
                }
                parseRfAssignmentOptions(settings, tokens, i);
                if (!settings.f_fund.empty()) settings.pss_requested = true;
                netlist.setSettings(settings);
            } else if (cmd == ".HBAC" || cmd == ".HBXF" || cmd == ".HBNOISE" || cmd == ".HBSP" || cmd == ".HBSTB" ||
                       cmd == ".PSSSP" || cmd == ".PSSP" || cmd == ".PSP" || cmd == ".PSSSTB" || cmd == ".PSTB") {
                SimulationSettings settings = netlist.getSettings();
                if (cmd == ".HBXF") settings.type = "HBAC";
                else if (cmd == ".PSTB") settings.type = "PSSSTB";
                else if (cmd == ".PSSP" || cmd == ".PSP") settings.type = "PSSSP";
                else settings.type = cmd.substr(1);

                settings.pss_requested = true;
                size_t i = 1;
                std::string outPos;
                std::string outNeg;
                if (i < tokens.size() && parseVoltageProbeToken(tokens[i], outPos, outNeg)) {
                    settings.out_node = netlist.getOrCreateNode(outPos);
                    ++i;
                }
                if (i < tokens.size()) {
                    if (!isRfSweepKeyword(tokens[i]) && tokens[i].find('=') == std::string::npos) {
                        settings.f_fund.clear();
                        while (i < tokens.size() && !isRfSweepKeyword(tokens[i]) &&
                               tokens[i].find('=') == std::string::npos) {
                            if (settings.f_fund.size() < 4) {
                                settings.f_fund.push_back(Utils::parseValue(tokens[i]));
                            } else {
                                std::cerr << "Warning: GSPICE supports max 4 tones. Ignoring extra: "
                                          << tokens[i] << std::endl;
                            }
                            ++i;
                        }
                    }
                }
                if (i < tokens.size()) {
                    settings.f_sweep_type = toUpperCopy(tokens[i]);
                    ++i;
                }
                if (settings.f_sweep_type == "VALUES") {
                    settings.f_values.clear();
                    for (; i < tokens.size(); ++i) {
                        if (tokens[i].find('=') != std::string::npos) break;
                        settings.f_values.push_back(Utils::parseValue(tokens[i]));
                    }
                    settings.points_per_dec = static_cast<int>(settings.f_values.size());
                    if (!settings.f_values.empty()) {
                        settings.f_start = settings.f_values.front();
                        settings.f_stop = settings.f_values.back();
                    }
                } else if (i + 2 < tokens.size()) {
                    settings.points_per_dec = std::stoi(tokens[i]);
                    settings.f_start = Utils::parseValue(tokens[i + 1]);
                    settings.f_stop = Utils::parseValue(tokens[i + 2]);
                    i += 3;
                }
                parseRfAssignmentOptions(settings, tokens, i);
                if (settings.f_fund.empty()) {
                    settings.f_fund.push_back(std::max(settings.f_start, 1.0));
                }
                netlist.setSettings(settings);
            } else {
                netlist.addError("Line " + std::to_string(lineNo) + ": unsupported directive; refusing to ignore active simulator syntax: " + line);
            }

        } else if (firstChar == 'R') {
            // Resistor: Rname N1 N2 Value or R=<value>
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid resistor line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = 0.0;
            if (!parsePrimitiveValue(tokens, 3, {"R", "RES", "VALUE"}, val)) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid resistor value: " + line);
                continue;
            }
            netlist.addDevice(std::make_unique<Resistor>(tokens[0], n1, n2, val));
        } else if (firstChar == 'C') {
            // Capacitor: Cname N1 N2 Value or C=<value>
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid capacitor line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
double val = 0.0;
            const ModelCard* capModel = netlist.findModelCard(tokens[3]);
            if (capModel && modelTypeMatches(capModel, {"C", "CAP"})) {
                const auto instanceParams = parseParameterTokens(tokens, 4);
                const double cj = paramValue(capModel->params, {"CJ", "C", "CAP"}, std::numeric_limits<double>::quiet_NaN());
                const double cjsw = paramValue(capModel->params, {"CJSW"}, 0.0);
                const double l = paramValue(instanceParams, {"L"}, 1.0);
                const double w = paramValue(instanceParams, {"W"}, 1.0);
                const double scale = paramValue(instanceParams, {"SCALE"}, 1.0);
                const double tempC = netlist.getSettings().temperature_c;
                const double tnom = paramValue(capModel->params, {"TNOM"}, 27.0);
                const double tc1 = paramValue(instanceParams, {"TC1"}, paramValue(capModel->params, {"TC1"}, 0.0));
                const double tc2 = paramValue(instanceParams, {"TC2"}, paramValue(capModel->params, {"TC2"}, 0.0));
                if (std::isfinite(cj)) {
                    const double deltaT = tempC - tnom;
                    const double tempScale = std::max(1.0 + tc1 * deltaT + tc2 * deltaT * deltaT, 1e-12);
                    val = (cj * l * w + cjsw * 2.0 * (l + w)) * scale * tempScale;
                    if (verboseCompatWarnings()) {
                        netlist.addWarning(
                            "Line " + std::to_string(lineNo) +
                            ": capacitor model '" + tokens[3] +
                            "' routed as (CJ*L*W + CJSW*perimeter)*SCALE with TC1/TC2 temperature scaling.");
                    }
                }
            } else {
                bool parsedOk = parsePrimitiveValue(tokens, 3, {"C", "CAP", "VALUE"}, val);
                (void)parsedOk;
                if (val == 0.0 && !parsedOk) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": invalid capacitor value: " + line);
                    continue;
                }
            }
            netlist.addDevice(std::make_unique<Capacitor>(tokens[0], n1, n2, val));
        } else if (firstChar == 'Y') {
            if (tokens.size() < 5 || toUpperCopy(tokens[3]) != "MOSVAR") {
                netlist.addError("Line " + std::to_string(lineNo) + ": unsupported Y element: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            const auto params = parseParameterTokens(tokens, 4);
            const double cacc = paramValue(params, {"CACC", "CACCUM", "C"}, std::numeric_limits<double>::quiet_NaN());
            if (!std::isfinite(cacc) || cacc <= 0.0) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid MOSVAR capacitance: " + line);
                continue;
            }
            const double cmin = paramValue(params, {"CMIN"}, 0.15 * cacc);
            const double vfb = paramValue(params, {"VFB", "VFBO"}, 0.0);
            const double slope = paramValue(params, {"SLOPE"}, 0.25);
            netlist.addDevice(std::make_unique<MosvarCapacitor>(tokens[0], n1, n2, cacc, cmin, vfb, slope));
        } else if (firstChar == 'L') {
            // Inductor: Lname N1 N2 Value
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid inductor line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Inductor>(tokens[0], n1, n2, val, -1));
        } else if (firstChar == 'T') {
            // Transmission line compatibility: Tname A B C D Z0=... TD=...
            // Implemented as a one-section passive LC approximation.
            if (tokens.size() < 7) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid transmission line: " + line);
                continue;
            }
            const auto params = parseParameterTokens(tokens, 5);
            const double z0 = paramValue(params, {"Z0", "ZO", "R", "IMPEDANCE"}, std::numeric_limits<double>::quiet_NaN());
            const double td = paramValue(params, {"TD", "DELAY", "TDELAY"}, std::numeric_limits<double>::quiet_NaN());
            if (!std::isfinite(z0) || z0 <= 0.0 || !std::isfinite(td) || td < 0.0) {
                netlist.addError("Line " + std::to_string(lineNo) + ": transmission line requires positive Z0 and non-negative TD: " + line);
                continue;
            }
            const int a = netlist.getOrCreateNode(tokens[1]);
            const int b = netlist.getOrCreateNode(tokens[2]);
            const int c = netlist.getOrCreateNode(tokens[3]);
            const int d = netlist.getOrCreateNode(tokens[4]);
            const std::string base = sanitizeIdentifier(tokens[0]);
            const double ctotal = td > 0.0 ? td / z0 : 1e-18;
            const double ltotal = td > 0.0 ? z0 * td : 1e-18;
            netlist.addDevice(std::make_unique<Capacitor>("C_" + base + "_IN", a, b, 0.5 * ctotal));
            netlist.addDevice(std::make_unique<Capacitor>("C_" + base + "_OUT", c, d, 0.5 * ctotal));
            if (b == d) {
                netlist.addDevice(std::make_unique<Inductor>("L_" + base + "_SER", a, c, ltotal, -1));
            } else {
                netlist.addDevice(std::make_unique<Inductor>("L_" + base + "_TOP", a, c, 0.5 * ltotal, -1));
                netlist.addDevice(std::make_unique<Inductor>("L_" + base + "_RET", d, b, 0.5 * ltotal, -1));
            }
            if (verboseCompatWarnings()) {
                netlist.addWarning(
                    "Line " + std::to_string(lineNo) +
                    ": transmission line '" + tokens[0] +
                    "' routed as a one-section LC approximation from Z0/TD.");
            }
        } else if (firstChar == 'K') {
            // Mutual coupling: Kname Lprimary Lsecondary coefficient
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid mutual inductor line: " + line);
                continue;
            }
            double coupling = 0.0;
            if (!parsePrimitiveValue(tokens, 3, {"K", "COUPLING", "VALUE"}, coupling)) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid mutual inductor coupling: " + line);
                continue;
            }
            netlist.addDevice(std::make_unique<MutualInductor>(tokens[0], tokens[1], tokens[2], coupling));
        } else if (firstChar == 'J') {
            // JFET/MESFET primitive: Jname D G S Model [AREA=...]
            if (tokens.size() < 5) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid JFET line: " + line);
                continue;
            }
            const ModelCard* modelCard = netlist.findModelCard(tokens[4]);
            if (!modelCard) {
                netlist.addError("Line " + std::to_string(lineNo) + ": JFET model card not found: " + line);
                continue;
            }
            const std::string modelType = toUpperCopy(modelCard->type);
            const bool pChannel = modelType == "PJF" || modelType == "PJFET" || modelType == "PJP" || modelType == "PMES";
            if (!pChannel && !(modelType == "NJF" || modelType == "NJFET" || modelType == "NJP" ||
                               modelType == "MES" || modelType == "MESFET" || modelType == "NMES")) {
                netlist.addError("Line " + std::to_string(lineNo) + ": unsupported JFET model type '" + modelCard->type + "': " + line);
                continue;
            }
            const auto instanceParams = parseParameterTokens(tokens, 5);
            const double area = std::max(paramValue(instanceParams, {"AREA", "M"}, 1.0), 0.0);
            const int type = pChannel ? -1 : 1;
            const double beta = paramValue(modelCard->params, {"BETA", "BET", "KP"}, 1e-4) * area;
            double vto = paramValue(modelCard->params, {"VTO", "VT0", "VP"}, pChannel ? 2.0 : -2.0);
            if (pChannel && vto < 0.0) vto = -vto;
            if (!pChannel && vto > 0.0) vto = -vto;
            const double lambda = paramValue(modelCard->params, {"LAMBDA", "L"}, 0.0);
            const double is = paramValue(modelCard->params, {"IS", "ISS"}, 1e-14) * area;
            const double n = paramValue(modelCard->params, {"N"}, 1.0);
            netlist.addDevice(std::make_unique<Jfet>(
                tokens[0],
                netlist.getOrCreateNode(tokens[1]),
                netlist.getOrCreateNode(tokens[2]),
                netlist.getOrCreateNode(tokens[3]),
                type,
                beta,
                vto,
                lambda,
                is,
                n));
        } else if (firstChar == 'B') {
            // Behavioral source: Bname N+ N- I={expr} or V={expr}.
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid behavioral source line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            BehavioralSource::Mode mode = BehavioralSource::Mode::Current;
            std::string expressionText;
            if (!parseBehavioralSourceSpec(tokens, 3, mode, expressionText)) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": behavioral source requires I={expr} or V={expr}: " + line);
                continue;
            }
            try {
                BehavioralExpression expression(
                    expressionText,
                    [&netlist](const std::string& nodeName) {
                        return netlist.getOrCreateNode(nodeName);
                    });
                netlist.addDevice(std::make_unique<BehavioralSource>(tokens[0], n1, n2, mode, expression, -1));
            } catch (const std::exception& ex) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": failed to parse behavioral source expression: " + std::string(ex.what()));
            }
        } else if (firstChar == 'G') {
            // VCCS: Gname N+ N- NC+ NC- value
            if (tokens.size() < 6) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid VCCS line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            int cp = netlist.getOrCreateNode(tokens[3]);
            int cn = netlist.getOrCreateNode(tokens[4]);
            double gm = Utils::parseValue(tokens[5]);
            netlist.addDevice(std::make_unique<VoltageControlledCurrentSource>(tokens[0], n1, n2, cp, cn, gm));
        } else if (firstChar == 'E') {
            // VCVS: Ename N+ N- NC+ NC- gain
            if (tokens.size() < 6) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid VCVS line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            int cp = netlist.getOrCreateNode(tokens[3]);
            int cn = netlist.getOrCreateNode(tokens[4]);
            double gain = Utils::parseValue(tokens[5]);
            netlist.addDevice(std::make_unique<VoltageControlledVoltageSource>(tokens[0], n1, n2, cp, cn, gain, -1));
        } else if (firstChar == 'F') {
            // CCCS: Fname N+ N- Vcontrol gain
            if (tokens.size() < 5) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid CCCS line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double gain = Utils::parseValue(tokens[4]);
            netlist.addDevice(std::make_unique<CurrentControlledCurrentSource>(tokens[0], n1, n2, tokens[3], gain));
        } else if (firstChar == 'H') {
            // CCVS: Hname N+ N- Vcontrol transresistance
            if (tokens.size() < 5) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid CCVS line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double transresistance = Utils::parseValue(tokens[4]);
            netlist.addDevice(std::make_unique<CurrentControlledVoltageSource>(tokens[0], n1, n2, tokens[3], transresistance, -1));
        } else if (firstChar == 'W') {
            // Stability Probe: Wname node_in node_out
            if (tokens.size() < 3) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid stability probe line: " + line);
                continue;
            }
            int nIn = netlist.getOrCreateNode(tokens[1]);
            int nOut = netlist.getOrCreateNode(tokens[2]);
            netlist.addDevice(std::make_unique<StabilityProbe>(tokens[0], nIn, nOut, -1));
        } else if (firstChar == 'P') {
            // Port: Pname N1 N2 PortNum Z0
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid port line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            int pNum = std::stoi(tokens[3]);
            double z0 = 50.0;
            if (tokens.size() > 4) z0 = Utils::parseValue(tokens[4]);
            netlist.addDevice(std::make_unique<Port>(tokens[0], n1, n2, pNum, z0));
        } else if (firstChar == 'S') {
            // S-parameter Multi-port Device: Sname N1 N2 ... Nn file="data.sNp"
            std::vector<int> nodes;
            size_t token_idx = 1;
            while (token_idx < tokens.size()) {
                bool is_num = !tokens[token_idx].empty() && (std::isdigit(tokens[token_idx][0]) || (tokens[token_idx][0] == '-' && tokens[token_idx].length() > 1 && std::isdigit(tokens[token_idx][1])));
                if (is_num) {
                    nodes.push_back(netlist.getOrCreateNode(tokens[token_idx]));
                } else {
                    break;
                }
                token_idx++;
            }
            std::string filename = "";
            for (; token_idx < tokens.size(); ++token_idx) {
                std::string tok = tokens[token_idx];
                std::transform(tok.begin(), tok.end(), tok.begin(), ::toupper);
                if (tok.substr(0, 5) == "FILE=") {
                    filename = tokens[token_idx].substr(5);
                    if (!filename.empty() && filename.front() == '"') filename = filename.substr(1, filename.size()-2);
                }
            }
            if (filename.empty()) {
                std::cerr << "Error: S-parameter device " << tokens[0] << " requires file=\"path\"" << std::endl;
                netlist.addError("Line " + std::to_string(lineNo) + ": S-parameter device requires file=\"path\": " + line);
                continue;
            }
            const int portCount = static_cast<int>(nodes.size());
            netlist.addDevice(std::make_unique<MultiPort>(tokens[0], nodes, portCount, filename));
        } else if (firstChar == 'X') {
            netlist.addError(
                "Line " + std::to_string(lineNo) +
                ": subcircuit instance is not implemented yet; refusing to ignore active device: " + line);
        } else if (firstChar == 'N') {
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid native compact-model N-device line: " + line);
                continue;
            }
            std::size_t modelIndex = 0;
            const ModelCard* model = nullptr;
            for (std::size_t i = 2; i < tokens.size(); ++i) {
                model = netlist.findModelCard(tokens[i]);
                if (model) {
                    modelIndex = i;
                    break;
                }
            }
            if (!model) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": model card not found for native compact device: " + line);
                continue;
            }
            const bool pspModel = isPsp103ModelCard(*model);
            const auto* artifact = netlist.findGsdiArtifact(pspModel ? psp103GsdiType(*model) : model->type);
            if (!artifact) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": .GSDI artifact for model type '" + model->type +
                    "' was not loaded before native compact device '" + tokens[0] + "'.");
                continue;
            }
            const std::size_t terminalCount = modelIndex - 1;
            if (artifact->terminal_count > 0 &&
                terminalCount != static_cast<std::size_t>(artifact->terminal_count)) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": native compact device '" + tokens[0] + "' has " +
                    std::to_string(terminalCount) + " terminal(s), but .GSDI model type '" +
                    model->type + "' declares " + std::to_string(artifact->terminal_count) + ".");
                continue;
            }
            auto instanceParams = parseParameterTokens(tokens, modelIndex + 1);
            const std::string unknownModelParam = firstUnknownGsdiParameter(*artifact, model->params);
            if (!unknownModelParam.empty()) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": .MODEL parameter '" + unknownModelParam + "' is not declared by .GSDI model type '" +
                    model->type + "'.");
                continue;
            }
            const std::string unknownInstanceParam = firstUnknownGsdiParameter(*artifact, instanceParams);
            if (!unknownInstanceParam.empty()) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": native compact device parameter '" + unknownInstanceParam +
                    "' is not declared by .GSDI model type '" + model->type + "'.");
                continue;
            }
            ensureGmcGeneratedModels();
            if (pspModel) {
                ensurePsp103GmcModels();
            }
            std::vector<int> deviceNodes;
            deviceNodes.reserve(terminalCount);
            for (std::size_t i = 1; i < modelIndex; ++i) {
                deviceNodes.push_back(netlist.getOrCreateNode(tokens[i]));
            }
            if (pspModel && deviceNodes.size() >= 4) {
                const double gateResistance = psp103RfGateResistance(*model, instanceParams);
                if (gateResistance > 0.0) {
                    const int externalGate = deviceNodes[1];
                    const int internalGate = netlist.createInternalNode(tokens[0] + ".rg", 0);
                    netlist.addDevice(std::make_unique<Resistor>(
                        "R_" + tokens[0] + "_RG", externalGate, internalGate, gateResistance));
                    deviceNodes[1] = internalGate;
                    netlist.addModelStatus("GMC PSP RF gate resistance: " + tokens[0]);
                }
            }
            std::vector<int> internalNodes;
            const std::size_t hiddenCount = GmcRegistry::instance().internalNodeCount(model->type);
            internalNodes.reserve(hiddenCount);
            for (std::size_t k = 0; k < hiddenCount; ++k) {
                internalNodes.push_back(netlist.createInternalNode(tokens[0], k));
            }
            GmcModelDefinition definition{
                tokens[0],
                pspModel ? psp103GmcType(*model) : model->type,
                pspModel ? psp103ModelParams(*model) : model->params,
                instanceParams,
                deviceNodes,
                netlist.getSettings().temperature_c,
                internalNodes
            };
            auto device = GmcRegistry::instance().create(definition);
            if (!device) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": no native GSDI/GMC model is registered for type '" + model->type +
                    "'; compile the Verilog-A with GMC to a .gsdi artifact and build/load that native model before running this deck.");
                continue;
            }
            netlist.addDevice(std::move(device));
        } else if (firstChar == 'M') {
            // MOSFET: Mname D G S B Model [W=..] [L=..]
            if (tokens.size() < 6) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid MOSFET line: " + line);
                continue;
            }
            int nD = netlist.getOrCreateNode(tokens[1]);
            int nG = netlist.getOrCreateNode(tokens[2]);
            int nS = netlist.getOrCreateNode(tokens[3]);
            int nB = netlist.getOrCreateNode(tokens[4]);
            
            double w = 1e-6, l = 1e-6;
            const ModelCard* modelCard = netlist.findModelCard(tokens[5]);
            int type = 1; // Default NMOS
            std::string modelNameUpper = toUpperCopy(tokens[5]);
            if (modelCard) {
                const std::string modelType = toUpperCopy(modelCard->type);
                if (modelType.find("PMOS") != std::string::npos || modelType == "P") type = -1;
                if (modelType.find("NMOS") != std::string::npos || modelType == "N") type = 1;
            } else if (modelNameUpper.find("PMOS") != std::string::npos) {
                type = -1;
            }

            double vth = type > 0 ? 0.5 : 0.5;
            double kp = 100e-6;
            double lambda = 0.05;
            double gamma = 0.4;
            double phi = 0.7;
            const auto instanceParams = parseParameterTokens(tokens, 6);
            w = paramValue(instanceParams, {"W"}, w);
            l = paramValue(instanceParams, {"L"}, l);
            if (modelCard && isPsp103ModelCard(*modelCard)) {
                const auto* artifact = netlist.findGsdiArtifact(psp103GsdiType(*modelCard));
                if (!artifact) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": .GSDI artifact for model type '" + modelCard->type +
                        "' was not loaded before PSP103 MOS device '" + tokens[0] + "'.");
                    continue;
                }
                if (artifact->terminal_count > 0 && artifact->terminal_count != 4) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": PSP103 MOS device '" + tokens[0] +
                        "' requires 4 terminals, but .GSDI model type '" +
                        modelCard->type + "' declares " +
                        std::to_string(artifact->terminal_count) + ".");
                    continue;
                }
                ensurePsp103GmcModels();
                std::array<int, 4> pspNodes = {nD, nG, nS, nB};
                const double gateResistance = psp103RfGateResistance(*modelCard, instanceParams);
                if (gateResistance > 0.0) {
                    const int internalGate = netlist.createInternalNode(tokens[0] + ".rg", 0);
                    netlist.addDevice(std::make_unique<Resistor>(
                        "R_" + tokens[0] + "_RG", nG, internalGate, gateResistance));
                    pspNodes[1] = internalGate;
                    netlist.addModelStatus("GMC PSP RF gate resistance: " + tokens[0]);
                }
                GmcModelDefinition definition{
                    tokens[0], psp103GmcType(*modelCard), psp103ModelParams(*modelCard),
                    instanceParams, {pspNodes[0], pspNodes[1], pspNodes[2], pspNodes[3]},
                    netlist.getSettings().temperature_c};
                auto device = GmcRegistry::instance().create(definition);
                if (!device) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": PSP103 model '" + tokens[5] +
                        "' was routed through native GSDI/GMC, but its electrical evaluator "
                        "is not implemented; refusing primitive Level-1 fallback.");
                } else {
                    netlist.addDevice(std::move(device));
                }
                continue;
            } else if (modelCard && isUnsupportedCompactMosLevel(modelCard) &&
                       static_cast<int>(paramValue(modelCard->params, {"LEVEL"}, 1.0)) == 49) {
                ensureBsim3GmcModels();
                GmcModelDefinition definition{
                    tokens[0], type > 0 ? "BSIM3_NMOS" : "BSIM3_PMOS",
                    modelCard->params, instanceParams, {nD, nG, nS, nB},
                    netlist.getSettings().temperature_c};
                auto device = GmcRegistry::instance().create(definition);
                if (!device) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": native BSIM3 translation/evaluation rejected this model; "
                        "refusing primitive fallback.");
                } else {
                    const std::string ignoredSummary = bsim3IgnoredParameterSummary(*modelCard);
                    if (verboseCompatWarnings() && !ignoredSummary.empty()) {
                        netlist.addWarning(
                            "Line " + std::to_string(lineNo) +
                            ": " + ignoredSummary);
                    }
                    netlist.addModelStatus("GMC BSIM3: " + tokens[0]);
                    netlist.addDevice(std::move(device));
                }
                continue;
            } else if (modelCard && isUnsupportedCompactMosLevel(modelCard) &&
                       static_cast<int>(paramValue(modelCard->params, {"LEVEL"}, 1.0)) == 54) {
                ensureBsim4GmcModels();
                GmcModelDefinition definition{
                    tokens[0], type > 0 ? "BSIM4_NMOS" : "BSIM4_PMOS",
                    modelCard->params, instanceParams, {nD, nG, nS, nB},
                    netlist.getSettings().temperature_c};
                auto device = GmcRegistry::instance().create(definition);
                if (!device) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": native BSIM4 translation/evaluation rejected this model; "
                        "refusing primitive fallback.");
                } else {
                    const std::string ignoredSummary = bsim4IgnoredParameterSummary(*modelCard);
                    if (verboseCompatWarnings() && !ignoredSummary.empty()) {
                        netlist.addWarning(
                            "Line " + std::to_string(lineNo) +
                            ": " + ignoredSummary);
                    }
                    netlist.addModelStatus("GMC BSIM4: " + tokens[0]);
                    netlist.addDevice(std::move(device));
                }
                continue;
            } else if (modelCard && modelTypeMatches(modelCard, {"NMOS", "PMOS", "N", "P"})) {
                if (isUnsupportedCompactMosLevel(modelCard) && !primitiveModelFallbackEnabled()) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": BSIM compact MOS level " + std::to_string(static_cast<int>(paramValue(modelCard->params, {"LEVEL"}, 1.0))) +
                        " is not implemented; refusing primitive Level-1 fallback.");
                    continue;
                }
                vth = paramValue(modelCard->params, {"VTO", "VT0", "VTH", "VTH0"}, vth);
                kp = paramValue(modelCard->params, {"KP", "BETA", "K"}, kp);
                lambda = paramValue(modelCard->params, {"LAMBDA", "LAMDA"}, lambda);
                gamma = paramValue(modelCard->params, {"GAMMA"}, gamma);
                phi = paramValue(modelCard->params, {"PHI"}, phi);
            } else if (modelCard) {
                if (isLikelyCompactMosModelType(modelCard->type) && !primitiveModelFallbackEnabled()) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": MOS model '" + tokens[5] + "' has compact/PDK type '" +
                        modelCard->type + "'. Refusing primitive Level-1 fallback. "
                        "Add the matching native GSDI/GMC model or set GSPICE_ALLOW_PRIMITIVE_MODEL_FALLBACK=1 only for debug smoke decks.");
                    continue;
                }
                netlist.addWarning(
                    "Line " + std::to_string(lineNo) +
                    ": MOS model '" + tokens[5] + "' has unsupported primitive type '" +
                    modelCard->type + "'; using GSPICE Level-1 fallback parameters.");
            } else if (isLikelyCompactMosModelName(tokens[5]) && !primitiveModelFallbackEnabled()) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": MOS model '" + tokens[5] + "' looks like a PDK compact model but no supported model card was loaded. "
                    "Refusing primitive Level-1 fallback; add the real native GSDI/GMC model first.");
                continue;
            }

            // Simple param parsing: look for W= and L=
            for (size_t i = 6; i < tokens.size(); ++i) {
                std::string upperTok = toUpperCopy(tokens[i]);
                if (upperTok.rfind("W=", 0) == 0) w = Utils::parseValue(stripQuotes(tokens[i].substr(2)));
                if (upperTok.rfind("L=", 0) == 0) l = Utils::parseValue(stripQuotes(tokens[i].substr(2)));
            }

            ensurePrimitiveMosGmcModels();
            GmcModelDefinition definition{
                tokens[0], type > 0 ? "NMOS" : "PMOS",
                {{"VTO", std::to_string(vth)}, {"KP", std::to_string(kp)},
                 {"LAMBDA", std::to_string(lambda)}, {"GAMMA", std::to_string(gamma)},
                 {"PHI", std::to_string(phi)}},
                {{"W", std::to_string(w)}, {"L", std::to_string(l)}},
                {nD, nG, nS, nB},
                netlist.getSettings().temperature_c
            };
            auto device = GmcRegistry::instance().create(definition);
            if (!device) {
                netlist.addError("Line " + std::to_string(lineNo) +
                                 ": native GMC primitive MOS factory is unavailable");
                continue;
            }
            netlist.addModelStatus("GMC primitive MOS: " + tokens[0] +
                                   " type=" + definition.type);
            netlist.addDevice(std::move(device));
        } else if (firstChar == 'V') {
            // Voltage Source: Vname N1 N2 [DC <value>] [AC <mag>] [PULSE(...)/SIN(...)/PWL(...)] or scalar value
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid voltage source line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            const std::string sourceSpec = joinTokens(tokens, 3);
            double dcValue = 0.0;
            double acMagnitude = 1.0;
            double acPhaseDeg = 0.0;
            bool dcSeen = false;
            VoltageSource::WaveformType wf = VoltageSource::WaveformType::DC;
            VoltageSource::PulseParams pulse;
            VoltageSource::SinParams sin;
            std::vector<double> pwlT;
            std::vector<double> pwlV;
            parseSourceSpec(sourceSpec, dcValue, acMagnitude, acPhaseDeg, dcSeen, wf, pulse, sin, pwlT, pwlV);
            if (!dcSeen) {
                if (wf == VoltageSource::WaveformType::PULSE) dcValue = pulse.v1;
                if (wf == VoltageSource::WaveformType::SIN) dcValue = sin.vo;
                if (wf == VoltageSource::WaveformType::PWL && !pwlV.empty()) dcValue = pwlV.front();
            }

            auto vsrc = std::make_unique<VoltageSource>(tokens[0], n1, n2, dcValue, -1);
            vsrc->setAcMagnitude(acMagnitude);
            vsrc->setAcPhaseDeg(acPhaseDeg);
            if (wf == VoltageSource::WaveformType::PULSE) vsrc->setPulse(pulse);
            if (wf == VoltageSource::WaveformType::SIN) vsrc->setSin(sin);
            if (wf == VoltageSource::WaveformType::PWL) vsrc->setPwl(pwlT, pwlV);
            netlist.addDevice(std::move(vsrc));
        } else if (firstChar == 'I') {
            // Current Source: Iname N1 N2 [DC <value>] [waveform...]
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid current source line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            const std::string sourceSpec = joinTokens(tokens, 3);
            double dcValue = 0.0;
            double acMagnitude = 1.0;
            double acPhaseDeg = 0.0;
            bool dcSeen = false;
            VoltageSource::WaveformType wf = VoltageSource::WaveformType::DC;
            VoltageSource::PulseParams pulse;
            VoltageSource::SinParams sin;
            std::vector<double> pwlT;
            std::vector<double> pwlV;
            parseSourceSpec(sourceSpec, dcValue, acMagnitude, acPhaseDeg, dcSeen, wf, pulse, sin, pwlT, pwlV);
            if (!dcSeen) {
                if (wf == VoltageSource::WaveformType::PULSE) dcValue = pulse.v1;
                if (wf == VoltageSource::WaveformType::SIN) dcValue = sin.vo;
                if (wf == VoltageSource::WaveformType::PWL && !pwlV.empty()) dcValue = pwlV.front();
            }
            auto isrc = std::make_unique<CurrentSource>(tokens[0], n1, n2, dcValue);
            isrc->setAcMagnitude(acMagnitude);
            isrc->setAcPhaseDeg(acPhaseDeg);
            if (wf == VoltageSource::WaveformType::PULSE) isrc->setPulse(pulse);
            if (wf == VoltageSource::WaveformType::SIN) isrc->setSin(sin);
            if (wf == VoltageSource::WaveformType::PWL) isrc->setPwl(pwlT, pwlV);
            netlist.addDevice(std::move(isrc));
        } else if (firstChar == 'Q') {
            // BJT: Qname C B E [S] model [area] [params...]
            if (tokens.size() < 5) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid BJT line: " + line);
                continue;
            }
            int nC = netlist.getOrCreateNode(tokens[1]);
            int nB = netlist.getOrCreateNode(tokens[2]);
            int nE = netlist.getOrCreateNode(tokens[3]);
            size_t modelIdx = 4;
            while (modelIdx < tokens.size() && !netlist.findModelCard(tokens[modelIdx])) {
                ++modelIdx;
            }
            if (modelIdx >= tokens.size()) {
                modelIdx = 4;
            }
            const ModelCard* modelCard = netlist.findModelCard(tokens[modelIdx]);
            int type = 1;
            double is = 1e-16;
            double bf = 100.0;
            double br = 1.0;
            double nf = 1.0;
            double nr = 1.0;
            double cje = 0.0;
            double cjc = 0.0;
            double tf = 0.0;
            double ikf = 0.0;
            double ikr = 0.0;
            double vaf = 0.0;
            double var = 0.0;
            double ise = 0.0;
            double isc = 0.0;
            double ne = 1.5;
            double nc = 2.0;
            if (modelCard) {
                const std::string modelType = toUpperCopy(modelCard->type);
                if (modelType == "PNP") type = -1;
                if (modelType == "NPN") type = 1;
                if (modelType != "NPN" && modelType != "PNP") {
                    netlist.addWarning(
                        "Line " + std::to_string(lineNo) +
                        ": BJT model '" + tokens[modelIdx] + "' has unsupported type '" +
                        modelCard->type + "'; using NPN fallback.");
                }
                is = paramValue(modelCard->params, {"IS"}, is);
                bf = paramValue(modelCard->params, {"BF", "BETA", "BETA_F"}, bf);
                br = paramValue(modelCard->params, {"BR", "BETA_R"}, br);
                nf = paramValue(modelCard->params, {"NF"}, nf);
                nr = paramValue(modelCard->params, {"NR"}, nr);
                cje = paramValue(modelCard->params, {"CJE", "CBE"}, cje);
                cjc = paramValue(modelCard->params, {"CJC", "CBC"}, cjc);
                tf = paramValue(modelCard->params, {"TF"}, tf);
                ikf = paramValue(modelCard->params, {"IKF", "IK"}, ikf);
                ikr = paramValue(modelCard->params, {"IKR"}, ikr);
                vaf = paramValue(modelCard->params, {"VAF", "VA", "BF_EARLY"}, vaf);
                var = paramValue(modelCard->params, {"VAR", "VB", "BR_EARLY"}, var);
                ise = paramValue(modelCard->params, {"ISE"}, ise);
                isc = paramValue(modelCard->params, {"ISC"}, isc);
                ne = paramValue(modelCard->params, {"NE"}, ne);
                nc = paramValue(modelCard->params, {"NC"}, nc);
            } else {
                netlist.addWarning(
                    "Line " + std::to_string(lineNo) +
                    ": BJT model '" + tokens[modelIdx] + "' was not found; using default NPN parameters.");
            }
            double area = 1.0;
            if (modelIdx + 1 < tokens.size() && tokens[modelIdx + 1].find('=') == std::string::npos) {
                double parsedArea = 1.0;
                if (tryParseSpiceValue(tokens[modelIdx + 1], parsedArea)) area = parsedArea;
            }
            auto instanceParams = parseParameterTokens(tokens, modelIdx + 1);
            area = paramValue(instanceParams, {"AREA", "M"}, area);
            netlist.addDevice(std::make_unique<Bjt>(
                tokens[0], nC, nB, nE, type, is, bf, br, nf, nr, area, cje, cjc, tf,
                ikf, ikr, vaf, var, ise, isc, ne, nc));
        } else if (firstChar == 'D') {
            // Diode: Dname N1 N2 [model] [area]
            if (tokens.size() < 3) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid diode line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            const ModelCard* modelCard = tokens.size() >= 4 ? netlist.findModelCard(tokens[3]) : nullptr;
            if (modelCard && (toUpperCopy(modelCard->type) == "JUNCAPEXP" ||
                              toUpperCopy(modelCard->type) == "JUNCAP2")) {
                const auto params = parseParameterTokens(tokens, 4);
                const double swjunexp = paramValue(modelCard->params, {"SWJUNEXP"}, 1.0);
                const bool express = toUpperCopy(modelCard->type) == "JUNCAPEXP" && swjunexp == 1.0;
                if (toUpperCopy(modelCard->type) == "JUNCAPEXP" && swjunexp != 0.0 && swjunexp != 1.0) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": JUNCAPEXP SWJUNEXP must be 0 or 1.");
                    continue;
                }
                ensureJuncapExpressGmcModels();
                GmcModelDefinition definition{
                    tokens[0], express ? "JUNCAPEXP" : "JUNCAP2", modelCard->params, params,
                    {n1, n2}, netlist.getSettings().temperature_c};
                auto device = GmcRegistry::instance().create(definition);
                if (!device) {
                    netlist.addError(
                        "Line " + std::to_string(lineNo) +
                        ": native JUNCAP factory could not create diode '" +
                        tokens[0] + "'.");
                } else {
                    netlist.addModelStatus(std::string(express ? "GMC JUNCAP Express: " : "GMC JUNCAP2: ") + tokens[0]);
                    netlist.addDevice(std::move(device));
                }
                continue;
            }
            double is = 1e-14;
            double n = 1.0;
            double cjo = 0.0;
            double rs = 0.0;
            double bv = 0.0;
            double ibv = 1e-10;
            double nbv = 1.0;
            double area = 1.0;
            if (tokens.size() >= 4) {
                if (modelCard && modelTypeMatches(modelCard, {"D", "DIODE"})) {
                    is = paramValue(modelCard->params, {"IS", "JS"}, is);
                    n = paramValue(modelCard->params, {"N", "NF"}, n);
                    cjo = paramValue(modelCard->params, {"CJO", "CJ0", "CJ"}, cjo);
                    rs = paramValue(modelCard->params, {"RS"}, rs);
                    bv = paramValue(modelCard->params, {"BV", "VJBR"}, bv);
                    ibv = paramValue(modelCard->params, {"IBV", "IJBR"}, ibv);
                    nbv = paramValue(modelCard->params, {"NBV", "NBR"}, nbv);
                } else if (modelCard) {
                    netlist.addWarning(
                        "Line " + std::to_string(lineNo) +
                        ": diode model '" + tokens[3] + "' has unsupported type '" +
                        modelCard->type + "'; using default diode parameters.");
                }
                if (tokens.size() >= 5 && tokens[4].find('=') == std::string::npos) {
                    double parsedArea = 1.0;
                    if (tryParseSpiceValue(tokens[4], parsedArea)) area = parsedArea;
                }
                auto instanceParams = parseParameterTokens(tokens, 4);
                area = paramValue(instanceParams, {"AREA", "M"}, area);
            }
            area = std::max(area, 1e-30);
            netlist.addDevice(std::make_unique<Diode>(
                tokens[0], n1, n2, is * area, n, cjo * area, rs / area, bv, ibv * area, nbv));
        } else {
            netlist.addError("Line " + std::to_string(lineNo) + ": unsupported element; refusing to ignore active device: " + line);
        }
    }

    return netlist;
}

} // namespace gspice
