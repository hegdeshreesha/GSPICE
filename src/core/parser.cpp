#include "parser.hpp"
#include "utils.hpp"
#include "devices/resistor.hpp"
#include "devices/capacitor.hpp"
#include "devices/diode.hpp"
#include "devices/voltage_source.hpp"
#include "devices/inductor.hpp"
#include "devices/port.hpp"
#include "devices/mosfet.hpp"
#include "devices/probe.hpp"
#include "devices/current_source.hpp"
#include "devices/multi_port.hpp"
#include "devices/osdi_device.hpp"
#include "devices/controlled_source.hpp"
#include "devices/bjt.hpp"
#include "devices/behavioral_source.hpp"
#include "expression.hpp"
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <sstream>
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
        "U", "STEP", "URAMP", "FLOOR", "CEIL", "CEILING", "ROUND"
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
    for (auto& token : tokens) {
        auto [key, value] = splitParameterToken(token);
        if (!key.empty()) {
            double evaluated = 0.0;
            if (tryEvaluateParamExpression(value, params, evaluated)) {
                token = key + "=" + formatNumericValue(evaluated);
            }
            continue;
        }
        auto it = params.find(token);
        if (it != params.end()) token = resolvedParamString(it->second, params);
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

bool modelTypeMatches(const gspice::ModelCard* model, const std::vector<std::string>& types) {
    if (!model) return false;
    const std::string actual = toUpperCopy(model->type);
    for (const auto& type : types) {
        if (actual == toUpperCopy(type)) return true;
    }
    return false;
}

bool primitiveModelFallbackEnabled() {
    const char* raw = std::getenv("GSPICE_ALLOW_PRIMITIVE_MODEL_FALLBACK");
    std::string value = raw ? std::string(raw) : "";
    std::transform(value.begin(), value.end(), value.begin(), ::toupper);
    return value == "1" || value == "YES" || value == "TRUE" || value == "ON";
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
    const std::string type = toUpperCopy(modelType);
    return type.find("PSP") != std::string::npos ||
           type.find("BSIM") != std::string::npos ||
           type.find("HICUM") != std::string::npos ||
           type.find("EKV") != std::string::npos ||
           type.find("OSDI") != std::string::npos;
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
    for (auto& token : tokens) {
        auto [key, value] = splitParameterToken(token);
        if (!key.empty()) {
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

bool hasOsdiDeviceInBody(const SubcktDef& def) {
    for (const auto& line : def.body) {
        auto tokens = tokenizeSimple(line.text);
        if (!tokens.empty() && std::toupper(tokens[0][0]) == 'N') return true;
    }
    return false;
}

std::vector<PreLine> expandActiveOsdiWrapperBody(
    const SubcktDef& def,
    const std::string& subcktName,
    const std::string& prefix,
    const std::unordered_map<std::string, std::string>& pinMap,
    std::unordered_map<std::string, std::string>& localParams,
    std::vector<std::string>& errors) {
    std::vector<PreLine> expanded;
    auto activeBody = filterConditionalLines(def.body, localParams, errors, "subckt " + subcktName);
    for (const auto& body : activeBody) {
        auto bodyTokens = tokenizeSimple(body.text);
        if (bodyTokens.empty()) continue;
        if (std::toupper(bodyTokens[0][0]) != 'N') continue;
        PreLine paramBody = body;
        paramBody.text = applyLocalParams(paramBody.text, localParams);
        expanded.push_back({remapPrimitiveLine(paramBody, prefix, pinMap), paramBody.source, paramBody.lineNo});
    }
    if (expanded.empty()) {
        errors.push_back("OSDI wrapper subcircuit '" + subcktName +
                         "' did not select an active compact-model device for instance " + prefix);
    }
    return expanded;
}

void appendUniquePath(std::vector<std::filesystem::path>& roots, const std::filesystem::path& path) {
    if (path.empty()) return;
    auto normalized = path.lexically_normal();
    if (std::find(roots.begin(), roots.end(), normalized) == roots.end()) {
        roots.push_back(normalized);
    }
}

std::string getEnvVar(const char* envName) {
#ifdef _WIN32
    char* buffer = nullptr;
    size_t length = 0;
    if (_dupenv_s(&buffer, &length, envName) != 0 || !buffer) {
        return "";
    }
    std::string value(buffer);
    free(buffer);
    return value;
#else
    const char* value = std::getenv(envName);
    return value ? std::string(value) : "";
#endif
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
    return !(value == "0" || value == "NO" || value == "FALSE" || value == "OFF");
}

std::vector<std::filesystem::path> osdiSearchRoots(const std::string& sourcePath) {
    std::vector<std::filesystem::path> roots;
    appendEnvSearchRoots(roots, "GSPICE_OSDI_DIR");
    appendEnvSearchRoots(roots, "NGSPICE_OSDI_DIR");
    std::filesystem::path source(sourcePath);
    appendUniquePath(roots, source.parent_path());
    appendUniquePath(roots, source.parent_path() / "osdi");
    appendUniquePath(roots, std::filesystem::current_path());
    appendUniquePath(roots, std::filesystem::current_path() / "osdi");
    appendUniquePath(roots, std::filesystem::current_path() / "build" / "osdi");
    return roots;
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
        if (hasOsdiDeviceInBody(def)) {
            return expandActiveOsdiWrapperBody(def, subcktName, prefix, pinMap, localParams, errors);
        }
        if (!primitiveIhpFallbackEnabled()) {
            errors.push_back("IHP wrapper subcircuit '" + subcktName + "' for instance " + prefix +
                             " does not contain a PSP/OSDI compact-model device. Refusing primitive MOS fallback; " +
                             "set GSPICE_ALLOW_PRIMITIVE_IHP_FALLBACK=1 only for placeholder smoke decks.");
            return expanded;
        }
        const std::string model = ihpMosType < 0 ? "PMOS" : "NMOS";
        std::string w = "1u";
        std::string l = "1u";
        auto wit = localParams.find("w");
        if (wit == localParams.end()) wit = localParams.find("W");
        if (wit != localParams.end()) w = resolvedParamString(wit->second, localParams);
        auto lit = localParams.find("l");
        if (lit == localParams.end()) lit = localParams.find("L");
        if (lit != localParams.end()) l = resolvedParamString(lit->second, localParams);
        std::vector<std::string> mosLine = {
            "MPSP_" + prefix,
            tokens[1],
            tokens[2],
            tokens[3],
            tokens[4],
            model,
            "W=" + stripQuotes(w),
            "L=" + stripQuotes(l)
        };
        expanded.push_back({
            ".GSPICEWARN IHP wrapper " + subcktName +
                " did not contain a PSP/OSDI device; using GSPICE's simple MOS fallback for compatibility only.",
            line.source,
            line.lineNo
        });
        expanded.push_back({joinSimple(mosLine), line.source, line.lineNo});
        return expanded;
    }

    auto activeBody = filterConditionalLines(def.body, localParams, errors, "subckt " + subcktName);
    std::vector<PreLine> bodyLines;
    for (const auto& body : activeBody) {
        auto bodyTokens = tokenizeSimple(body.text);
        if (bodyTokens.empty()) continue;
        PreLine paramBody = body;
        paramBody.text = applyLocalParams(paramBody.text, localParams);
        bodyTokens = tokenizeSimple(paramBody.text);
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
            settings.tran_max_iter = 80;
        } else if (preset == "VERYHIGH") {
            settings.reltol = 1e-4;
            settings.vntol = 100e-9;
            settings.abstol = 10e-15;
            settings.tran_lte_reltol = 3e-4;
            settings.tran_lte_abstol = 100e-9;
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
            settings.osdi_limiting_rhs = true;
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
    if (key == "RELTOL" && hasNumeric && numeric > 0.0) settings.reltol = numeric;
    else if (key == "VNTOL" && hasNumeric && numeric > 0.0) settings.vntol = numeric;
    else if (key == "ABSTOL" && hasNumeric && numeric > 0.0) settings.abstol = numeric;
    else if (key == "GMIN" && hasNumeric && numeric >= 0.0) settings.gmin = numeric;
    else if (key == "ACCURACY" && !value.empty()) applyPreset(value);
    else if ((key == "NUMERICAL" || key == "CONVERGENCE" || key == "POLICY") && !value.empty()) applyNumericalPolicy(value);
    else if (key == "TRTOL" && hasNumeric && numeric > 0.0) settings.tran_trtol = numeric;
    else if ((key == "TRAN_RELTOL" || key == "LTE_RELTOL") && hasNumeric && numeric > 0.0) settings.tran_lte_reltol = numeric;
    else if ((key == "TRABSTOL" || key == "TRAN_ABSTOL" || key == "LTE_ABSTOL") && hasNumeric && numeric > 0.0) settings.tran_lte_abstol = numeric;
    else if (key == "CHGTOL" && hasNumeric && numeric > 0.0) settings.chgtol = numeric;
    else if ((key == "MAXORD" || key == "MAXORDER") && hasNumeric) settings.tran_max_order = std::clamp(static_cast<int>(numeric), 1, 5);
    else if ((key == "MAXSTEP" || key == "TMAX" || key == "TRAN_MAXSTEP") && hasNumeric && numeric > 0.0) settings.t_max_step = numeric;
    else if ((key == "MINSTEP" || key == "TMIN" || key == "TRAN_MINSTEP") && hasNumeric && numeric > 0.0) settings.t_min_step = numeric;
    else if (key == "ADAPTIVE" || key == "TRAN_ADAPTIVE") settings.tran_adaptive = value.empty() ? true : truthy(value);
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
    else if (key == "OSDI_LIMITING_RHS" || key == "OSDILIMITINGRHS" || key == "OSDI_LIM_RHS" || key == "OSDI_LIMITING") {
        settings.osdi_limiting_rhs = value.empty() ? true : truthy(value);
    }
    else if (key == "OSDI_TRAN_JACOBIAN" || key == "OSDITRANJACOBIAN" || key == "OSDI_TRAN_JAC" || key == "OSDI_STANDARD_TRAN_JACOBIAN") {
        settings.osdi_tran_jacobian = value.empty() ? true : truthy(value);
    }
    else if (key == "OSDI_BIND_FULL_MODEL_PARAMS" || key == "OSDIFULLMODELPARAMS" || key == "OSDI_FULL_MODEL_PARAMS") {
        settings.osdi_bind_full_model_params = value.empty() ? true : truthy(value);
    }
    else if (key == "OSDI_INTERNAL_NODES" || key == "OSDIINTERNALNODES" || key == "OSDI_EXPAND_INTERNAL_NODES") {
        settings.osdi_internal_nodes = value.empty() ? true : truthy(value);
    }
    else if (key == "OSDI_SPICE_RHS" || key == "OSDISPICERHS" || key == "OSDI_NATIVE_RHS") {
        settings.osdi_spice_rhs = value.empty() ? true : truthy(value);
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
}

void parseSourceSpec(
    const std::string& sourceSpec,
    double& dcValue,
    double& acMagnitude,
    bool& dcSeen,
    gspice::VoltageSource::WaveformType& waveformType,
    gspice::VoltageSource::PulseParams& pulse,
    gspice::VoltageSource::SinParams& sin,
    std::vector<double>& pwlTimes,
    std::vector<double>& pwlValues) {

    dcValue = 0.0;
    acMagnitude = 1.0;
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
            ++i;
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
    return false;
}

void parseMeasureParams(gspice::MeasureSpec& measure, const std::vector<std::string>& tokens, size_t startIdx) {
    for (size_t i = startIdx; i < tokens.size(); ++i) {
        auto [key, value] = splitParameterToken(tokens[i]);
        if (!key.empty()) {
            applyMeasureParam(measure, key, value);
            continue;
        }
        std::string upper = toUpperCopy(tokens[i]);
        if ((upper == "AT" || upper == "FROM" || upper == "TO") && i + 1 < tokens.size()) {
            if (tokens[i + 1] == "=" && i + 2 < tokens.size()) {
                applyMeasureParam(measure, upper, tokens[i + 2]);
                i += 2;
            } else {
                applyMeasureParam(measure, upper, tokens[i + 1]);
                ++i;
            }
        }
    }
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
        if (firstTokenUpper == "OSDI" || firstTokenUpper == "PRE_OSDI" || firstTokenUpper == "LOAD") {
            if (tokens.size() < 2) {
                netlist.addError("Line " + std::to_string(lineNo) + ": " + tokens[0] + " requires a library path.");
                continue;
            }
            std::string osdiPath = stripQuotes(tokens[1]);
            const std::string osdiPathUpper = toUpperCopy(osdiPath);
            if (firstTokenUpper == "LOAD" &&
                osdiPathUpper.size() >= 3 &&
                osdiPathUpper.substr(osdiPathUpper.size() - 3) == ".VA") {
                netlist.addError("Line " + std::to_string(lineNo) +
                    ": LOAD of Verilog-A source requires OpenVAF compilation first. Compile to .osdi, then LOAD the .osdi file.");
                continue;
            }
            if (osdiPath.rfind("builtin:", 0) != 0) {
                osdiPath = resolveRelativePath(std::filesystem::path(preLine.source), osdiPath).string();
            }
            std::string loadError;
            if (!netlist.loadOsdiLibrary(osdiPath, loadError)) {
                netlist.addError("Line " + std::to_string(lineNo) + ": " + loadError);
            }
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
                netlist.addWarning(joinTokens(tokens, 1));
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
                    if (toUpperCopy(save.kind) == "I") {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": " + cmd +
                                           " current save '" + tokens[i] +
                                           "' is not written yet; voltage saves are supported now.");
                        continue;
                    }
                    if (settings.save_all) {
                        settings.save_all = false;
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
                    } else {
                        settings.t_max_step = Utils::parseValue(tokens[4]);
                    }
                }
                if (tokens.size() > 5) {
                    std::string opt = tokens[5];
                    std::transform(opt.begin(), opt.end(), opt.begin(), ::toupper);
                    if (opt == "UIC") settings.use_uic = true;
                }
                netlist.setSettings(settings);
            } else if (cmd == ".OSDI" || cmd == ".PRE_OSDI") {
                if (tokens.size() < 2) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": " + cmd + " requires a library path.");
                    continue;
                }
                std::string osdiPath = stripQuotes(tokens[1]);
                if (osdiPath.rfind("builtin:", 0) != 0) {
                    osdiPath = resolveRelativePath(std::filesystem::path(preLine.source), osdiPath).string();
                }
                std::string loadError;
                if (!netlist.loadOsdiLibrary(osdiPath, loadError)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": " + loadError);
                }
            } else if (cmd == ".LOAD") {
                if (tokens.size() < 2) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": .LOAD requires a library path.");
                    continue;
                }
                std::string osdiPath = stripQuotes(tokens[1]);
                const std::string osdiPathUpper = toUpperCopy(osdiPath);
                if (osdiPathUpper.size() >= 3 && osdiPathUpper.substr(osdiPathUpper.size() - 3) == ".VA") {
                    netlist.addError("Line " + std::to_string(lineNo) +
                        ": .LOAD of Verilog-A source requires OpenVAF compilation first. Compile to .osdi, then .LOAD the .osdi file.");
                    continue;
                }
                if (osdiPath.rfind("builtin:", 0) != 0) {
                    osdiPath = resolveRelativePath(std::filesystem::path(preLine.source), osdiPath).string();
                }
                std::string loadError;
                if (!netlist.loadOsdiLibrary(osdiPath, loadError)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": " + loadError);
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
                netlist.addModelCard(model);
                if (toUpperCopy(tokens[2]) == "OSDI" && tokens.size() > 3) {
                    std::string pathPart = tokens[3];
                    size_t start = pathPart.find('"');
                    size_t end = pathPart.find_last_of('"');
                    if (start != std::string::npos && end != std::string::npos && end > start) {
                        std::string path = pathPart.substr(start + 1, end - start - 1);
                        netlist.addOsdiModel(tokens[1], path);
                    }
                }
            } else if (cmd == ".AC") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .AC line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "AC";
                // tokens[1] is DEC/OCT/LIN
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".PSS" || cmd == ".HB") {
                SimulationSettings settings = netlist.getSettings();
                settings.type = cmd.substr(1);
                size_t i = 1;
                while (i < tokens.size()) {
                    bool is_int = !tokens[i].empty() && std::all_of(tokens[i].begin(), tokens[i].end(), ::isdigit);
                    if (is_int) {
                        settings.n_harms = std::stoi(tokens[i]);
                        break;
                    } else if (settings.f_fund.size() < 4) {
                        settings.f_fund.push_back(Utils::parseValue(tokens[i]));
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
                if (!parseVoltageProbeToken(tokens[4], outPos, outNeg)) {
                    netlist.addError("Line " + std::to_string(lineNo) + ": only voltage .MEASURE expressions are supported currently: " + line);
                    continue;
                }
                MeasureSpec measure;
                measure.analysis = toUpperCopy(tokens[1]);
                measure.name = tokens[2];
                measure.op = toUpperCopy(tokens[3]);
                measure.node_pos = netlist.getOrCreateNode(outPos);
                measure.node_neg = netlist.getOrCreateNode(outNeg);
                parseMeasureParams(measure, tokens, 5);
                SimulationSettings settings = netlist.getSettings();
                settings.measures.push_back(measure);
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
            } else if (cmd == ".PAC") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .PAC line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "PAC";
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".PNOISE") {
                if (tokens.size() < 6) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .PNOISE line ignored: " + line);
                    continue;
                }
                SimulationSettings settings = netlist.getSettings();
                settings.type = "PNOISE";
                std::string outNode = tokens[1].substr(2, tokens[1].size()-3);
                settings.out_node = netlist.getOrCreateNode(outNode);
                size_t sweepIdx = 3;
                std::string sweepType = toUpperCopy(tokens[sweepIdx]);
                if (sweepType == "DEC" || sweepType == "OCT" || sweepType == "LIN") {
                    if (tokens.size() < 7) {
                        netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .PNOISE sweep line ignored: " + line);
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
            } else if (cmd == ".HBAC" || cmd == ".HBNOISE" || cmd == ".HBSP" || cmd == ".HBSTB" || 
                       cmd == ".PSSSP" || cmd == ".PSSSTB") {
                SimulationSettings settings = netlist.getSettings();
                settings.type = cmd.substr(1);
                
                size_t i = 1;
                while (i < tokens.size()) {
                    bool is_int = !tokens[i].empty() && std::all_of(tokens[i].begin(), tokens[i].end(), ::isdigit);
                    if (is_int) {
                        settings.points_per_dec = std::stoi(tokens[i]);
                        if (i + 1 < tokens.size()) settings.f_start = Utils::parseValue(tokens[i+1]);
                        if (i + 2 < tokens.size()) settings.f_stop = Utils::parseValue(tokens[i+2]);
                        break;
                    } else {
                        settings.f_fund.push_back(Utils::parseValue(tokens[i]));
                    }
                    i++;
                }
                netlist.setSettings(settings);
            } else {
                netlist.addError("Line " + std::to_string(lineNo) + ": unsupported directive; refusing to ignore active simulator syntax: " + line);
            }

        } else if (firstChar == 'R') {
            // Resistor: Rname N1 N2 Value
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid resistor line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Resistor>(tokens[0], n1, n2, val));
        } else if (firstChar == 'C') {
            // Capacitor: Cname N1 N2 Value
            if (tokens.size() < 4) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid capacitor line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Capacitor>(tokens[0], n1, n2, val));
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
            netlist.addDevice(std::make_unique<MultiPort>(tokens[0], nodes, nodes.size(), filename));
        } else if (firstChar == 'X') {
            netlist.addError(
                "Line " + std::to_string(lineNo) +
                ": subcircuit instance is not implemented yet; refusing to ignore active device: " + line);
        } else if (firstChar == 'N') {
            // OSDI/OpenVAF-style compact model instance:
            // Nname D G S B modelName [instance params...]
            if (tokens.size() < 6) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid OSDI N-device line: " + line);
                continue;
            }
            std::vector<int> nodes = {
                netlist.getOrCreateNode(tokens[1]),
                netlist.getOrCreateNode(tokens[2]),
                netlist.getOrCreateNode(tokens[3]),
                netlist.getOrCreateNode(tokens[4])
            };
            const ModelCard* model = netlist.findModelCard(tokens[5]);
            if (!model) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": OSDI model card not found for '" + tokens[5] + "': " + line);
                continue;
            }
            const OsdiDescriptor* desc = netlist.findOsdiDescriptor(model->type);
            std::string autoLoadError;
            if (!desc) {
                auto roots = osdiSearchRoots(preLine.source);
                if (netlist.tryAutoLoadOsdiForModelType(model->type, roots, autoLoadError)) {
                    desc = netlist.findOsdiDescriptor(model->type);
                }
            }
            if (!desc) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": OSDI model type '" + model->type +
                    "' is not loaded. Add an OSDI/PRE_OSDI line for the compiled model before using: " +
                    line + (autoLoadError.empty() ? "" : " (" + autoLoadError + ")"));
                continue;
            }
            try {
                const auto& settings = netlist.getSettings();
                auto instanceParams = parseParameterTokens(tokens, 6);
                netlist.addDevice(std::make_unique<OSDIDevice>(
                    tokens[0], *desc, nodes, model->params, instanceParams,
                    settings.temperature_c,
                    settings.osdi_limiting_rhs,
                    settings.osdi_tran_jacobian,
                    settings.osdi_bind_full_model_params,
                    settings.osdi_spice_rhs));
                const char* descName = desc->name ? desc->name : desc->model_name;
                netlist.addModelStatus(
                    "OSDI_DEVICE instance=" + tokens[0] +
                    " model=" + tokens[5] +
                    " type=" + model->type +
                    " descriptor=" + std::string(descName ? descName : "<unnamed>") +
                    " osdi_nodes=" + std::to_string(desc->num_nodes) +
                    " internal_mode=" + std::string(settings.osdi_internal_nodes ? "bound-hidden" : "collapsed-only"));
            } catch (const std::exception& ex) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": failed to create OSDI device '" + tokens[0] + "': " + ex.what());
            }
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
            if (modelCard && modelTypeMatches(modelCard, {"NMOS", "PMOS", "N", "P"})) {
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
                        "Load the matching OSDI/OpenVAF model or set GSPICE_ALLOW_PRIMITIVE_MODEL_FALLBACK=1 only for debug smoke decks.");
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
                    "Refusing primitive Level-1 fallback; load the real OSDI/OpenVAF model first.");
                continue;
            }

            // Simple param parsing: look for W= and L=
            for (size_t i = 6; i < tokens.size(); ++i) {
                std::string upperTok = toUpperCopy(tokens[i]);
                if (upperTok.rfind("W=", 0) == 0) w = Utils::parseValue(stripQuotes(tokens[i].substr(2)));
                if (upperTok.rfind("L=", 0) == 0) l = Utils::parseValue(stripQuotes(tokens[i].substr(2)));
            }

            netlist.addDevice(std::make_unique<Mosfet>(tokens[0], nD, nG, nS, nB, type, w, l, vth, kp, lambda, gamma, phi));
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
            bool dcSeen = false;
            VoltageSource::WaveformType wf = VoltageSource::WaveformType::DC;
            VoltageSource::PulseParams pulse;
            VoltageSource::SinParams sin;
            std::vector<double> pwlT;
            std::vector<double> pwlV;
            parseSourceSpec(sourceSpec, dcValue, acMagnitude, dcSeen, wf, pulse, sin, pwlT, pwlV);
            if (!dcSeen) {
                if (wf == VoltageSource::WaveformType::PULSE) dcValue = pulse.v1;
                if (wf == VoltageSource::WaveformType::SIN) dcValue = sin.vo;
                if (wf == VoltageSource::WaveformType::PWL && !pwlV.empty()) dcValue = pwlV.front();
            }

            auto vsrc = std::make_unique<VoltageSource>(tokens[0], n1, n2, dcValue, -1);
            vsrc->setAcMagnitude(acMagnitude);
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
            bool dcSeen = false;
            VoltageSource::WaveformType wf = VoltageSource::WaveformType::DC;
            VoltageSource::PulseParams pulse;
            VoltageSource::SinParams sin;
            std::vector<double> pwlT;
            std::vector<double> pwlV;
            parseSourceSpec(sourceSpec, dcValue, acMagnitude, dcSeen, wf, pulse, sin, pwlT, pwlV);
            if (!dcSeen) {
                if (wf == VoltageSource::WaveformType::PULSE) dcValue = pulse.v1;
                if (wf == VoltageSource::WaveformType::SIN) dcValue = sin.vo;
                if (wf == VoltageSource::WaveformType::PWL && !pwlV.empty()) dcValue = pwlV.front();
            }
            auto isrc = std::make_unique<CurrentSource>(tokens[0], n1, n2, dcValue);
            isrc->setAcMagnitude(acMagnitude);
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
            if (tokens.size() >= 6 && !netlist.findModelCard(tokens[4]) && netlist.findModelCard(tokens[5])) {
                modelIdx = 5; // substrate node is accepted but ignored by this first-pass BJT.
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
                tokens[0], nC, nB, nE, type, is, bf, br, nf, nr, area, cje, cjc, tf));
        } else if (firstChar == 'D') {
            // Diode: Dname N1 N2 [model] [area]
            if (tokens.size() < 3) {
                netlist.addError("Line " + std::to_string(lineNo) + ": invalid diode line: " + line);
                continue;
            }
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double is = 1e-14;
            double n = 1.0;
            double cjo = 0.0;
            double area = 1.0;
            if (tokens.size() >= 4) {
                const ModelCard* modelCard = netlist.findModelCard(tokens[3]);
                if (modelCard && modelTypeMatches(modelCard, {"D", "DIODE"})) {
                    is = paramValue(modelCard->params, {"IS", "JS"}, is);
                    n = paramValue(modelCard->params, {"N", "NF"}, n);
                    cjo = paramValue(modelCard->params, {"CJO", "CJ0", "CJ"}, cjo);
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
            netlist.addDevice(std::make_unique<Diode>(tokens[0], n1, n2, is * area, n, cjo * area));
        } else {
            netlist.addError("Line " + std::to_string(lineNo) + ": unsupported element; refusing to ignore active device: " + line);
        }
    }

    return netlist;
}

} // namespace gspice
