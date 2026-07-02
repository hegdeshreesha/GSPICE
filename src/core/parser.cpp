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
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
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

void replaceAll(std::string& text, const std::string& needle, const std::string& value) {
    if (needle.empty()) return;
    size_t pos = 0;
    while ((pos = text.find(needle, pos)) != std::string::npos) {
        text.replace(pos, needle.size(), value);
        pos += value.size();
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
            for (size_t i = 1; i < tokens.size(); ++i) {
                size_t eq = tokens[i].find('=');
                if (eq == std::string::npos || eq == 0) continue;
                std::string key = tokens[i].substr(0, eq);
                std::string value = tokens[i].substr(eq + 1);
                params[key] = value;
                params[toUpperCopy(key)] = value;
            }
        } else {
            nonParamLines.push_back(line);
        }
    }
    return params;
}

std::string applyGlobalParams(
    std::string line,
    const std::unordered_map<std::string, std::string>& params) {
    for (const auto& [key, value] : params) {
        replaceAll(line, "{" + key + "}", value);
    }
    auto tokens = tokenizeSimple(line);
    for (auto& token : tokens) {
        auto it = params.find(token);
        if (it != params.end()) token = it->second;
    }
    return joinSimple(tokens);
}

bool tryParseSpiceValue(const std::string& token, double& out) {
    try {
        out = gspice::Utils::parseValue(token);
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
    std::unordered_map<std::string, std::string> params;
    for (size_t i = startIdx; i < tokens.size(); ++i) {
        if (i + 2 < tokens.size() && tokens[i + 1] == "=") {
            addParam(params, tokens[i], stripQuotes(tokens[i + 2]));
            i += 2;
            continue;
        }
        if (!tokens[i].empty() && tokens[i].back() == '=' && i + 1 < tokens.size()) {
            addParam(params, tokens[i].substr(0, tokens[i].size() - 1), stripQuotes(tokens[i + 1]));
            ++i;
            continue;
        }
        auto [key, value] = splitParameterToken(tokens[i]);
        addParam(params, key, value);
    }
    return params;
}

std::string applyLocalParams(
    std::string line,
    const std::unordered_map<std::string, std::string>& params) {
    for (const auto& [key, value] : params) {
        replaceAll(line, "{" + key + "}", value);
        replaceAll(line, "'" + key + "'", value);
    }
    auto tokens = tokenizeSimple(line);
    for (auto& token : tokens) {
        auto [key, value] = splitParameterToken(token);
        if (!key.empty()) {
            auto it = params.find(value);
            if (it == params.end()) it = params.find(toUpperCopy(value));
            if (it != params.end()) token = key + "=" + it->second;
        } else {
            std::string stripped = stripQuotes(token);
            auto it = params.find(stripped);
            if (it == params.end()) it = params.find(toUpperCopy(stripped));
            if (it != params.end()) token = it->second;
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
        const std::string model = ihpMosType < 0 ? "PMOS" : "NMOS";
        const std::string pspModel = ihpMosType < 0 ? "sg13g2_lv_pmos_psp" : "sg13g2_lv_nmos_psp";
        std::string w = "1u";
        std::string l = "1u";
        std::string ng = "1";
        std::string m = "1";
        auto wit = localParams.find("w");
        if (wit == localParams.end()) wit = localParams.find("W");
        if (wit != localParams.end()) w = wit->second;
        auto lit = localParams.find("l");
        if (lit == localParams.end()) lit = localParams.find("L");
        if (lit != localParams.end()) l = lit->second;
        auto ngit = localParams.find("ng");
        if (ngit == localParams.end()) ngit = localParams.find("NG");
        if (ngit != localParams.end()) ng = ngit->second;
        auto mit = localParams.find("m");
        if (mit == localParams.end()) mit = localParams.find("M");
        if (mit != localParams.end()) m = mit->second;

        if (hasOsdiDeviceInBody(def)) {
            std::vector<std::string> osdiLine = {
                "NPSP_" + prefix,
                tokens[1],
                tokens[2],
                tokens[3],
                tokens[4],
                pspModel,
                "w=" + stripQuotes(w),
                "l=" + stripQuotes(l),
                "nf=" + stripQuotes(ng),
                "mult=" + stripQuotes(m)
            };
            expanded.push_back({joinSimple(osdiLine), line.source, line.lineNo});
        } else {
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
                ".GSPICEWARN IHP PSP wrapper " + subcktName +
                    " is approximated with GSPICE's simple MOS model; results are not PSP/PDK accurate.",
                line.source,
                line.lineNo
            });
            expanded.push_back({joinSimple(mosLine), line.source, line.lineNo});
        }
        return expanded;
    }

    std::vector<PreLine> bodyLines;
    for (const auto& body : def.body) {
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
    if (key == "RELTOL" && hasNumeric && numeric > 0.0) settings.reltol = numeric;
    else if (key == "VNTOL" && hasNumeric && numeric > 0.0) settings.vntol = numeric;
    else if (key == "ABSTOL" && hasNumeric && numeric > 0.0) settings.abstol = numeric;
    else if (key == "GMIN" && hasNumeric && numeric >= 0.0) settings.gmin = numeric;
    else if (key == "ACCURACY" && !value.empty()) applyPreset(value);
    else if ((key == "TRTOL" || key == "TRAN_RELTOL" || key == "LTE_RELTOL") && hasNumeric && numeric > 0.0) settings.tran_lte_reltol = numeric;
    else if ((key == "TRABSTOL" || key == "TRAN_ABSTOL" || key == "LTE_ABSTOL") && hasNumeric && numeric > 0.0) settings.tran_lte_abstol = numeric;
    else if ((key == "MINSTEP" || key == "TMIN" || key == "TRAN_MINSTEP") && hasNumeric && numeric > 0.0) settings.t_min_step = numeric;
    else if (key == "ADAPTIVE" || key == "TRAN_ADAPTIVE") settings.tran_adaptive = value.empty() ? true : truthy(value);
    else if ((key == "METHOD" || key == "TRAN_METHOD") && !value.empty()) {
        std::string method = toUpperCopy(value);
        method.erase(std::remove(method.begin(), method.end(), '_'), method.end());
        method.erase(std::remove(method.begin(), method.end(), '-'), method.end());
        if (method == "AUTO" || method == "BE" || method == "BACKWARDEULER" ||
            method == "TRAP" || method == "TRAPEZOIDAL" || method == "GEAR2") {
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
        if (firstTokenUpper == "OSDI" || firstTokenUpper == "PRE_OSDI") {
            if (tokens.size() < 2) {
                netlist.addError("Line " + std::to_string(lineNo) + ": " + tokens[0] + " requires a library path.");
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
            } else if (cmd == ".LIB" || cmd == ".INCLUDE" || cmd == ".INC") {
                netlist.addWarning("Line " + std::to_string(lineNo) + ": model/include directive is not implemented yet: " + line);
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
            } else if (cmd == ".MODEL") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .MODEL line ignored: " + line);
                    continue;
                }
                ModelCard model;
                model.name = tokens[1];
                model.type = tokens[2];
                model.params = parseParameterTokens(tokens, 3);
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
                settings.points_per_dec = std::stoi(tokens[3]);
                settings.f_start = Utils::parseValue(tokens[4]);
                settings.f_stop = Utils::parseValue(tokens[5]);
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
                settings.points_per_dec = std::stoi(tokens[3]);
                settings.f_start = Utils::parseValue(tokens[4]);
                settings.f_stop = Utils::parseValue(tokens[5]);
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
                netlist.addWarning("Line " + std::to_string(lineNo) + ": unsupported directive ignored: " + line);
            }

        } else if (firstChar == 'R') {
            // Resistor: Rname N1 N2 Value
            if (tokens.size() < 4) continue;
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Resistor>(tokens[0], n1, n2, val));
        } else if (firstChar == 'C') {
            // Capacitor: Cname N1 N2 Value
            if (tokens.size() < 4) continue;
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Capacitor>(tokens[0], n1, n2, val));
        } else if (firstChar == 'L') {
            // Inductor: Lname N1 N2 Value
            if (tokens.size() < 4) continue;
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Inductor>(tokens[0], n1, n2, val, -1));
        } else if (firstChar == 'W') {
            // Stability Probe: Wname node_in node_out
            if (tokens.size() < 3) continue;
            int nIn = netlist.getOrCreateNode(tokens[1]);
            int nOut = netlist.getOrCreateNode(tokens[2]);
            netlist.addDevice(std::make_unique<StabilityProbe>(tokens[0], nIn, nOut, -1));
        } else if (firstChar == 'P') {
            // Port: Pname N1 N2 PortNum Z0
            if (tokens.size() < 4) continue;
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
                auto instanceParams = parseParameterTokens(tokens, 6);
                netlist.addDevice(std::make_unique<OSDIDevice>(tokens[0], *desc, nodes, model->params, instanceParams));
            } catch (const std::exception& ex) {
                netlist.addError(
                    "Line " + std::to_string(lineNo) +
                    ": failed to create OSDI device '" + tokens[0] + "': " + ex.what());
            }
        } else if (firstChar == 'M') {
            // MOSFET: Mname D G S B Model [W=..] [L=..]
            if (tokens.size() < 6) continue;
            int nD = netlist.getOrCreateNode(tokens[1]);
            int nG = netlist.getOrCreateNode(tokens[2]);
            int nS = netlist.getOrCreateNode(tokens[3]);
            int nB = netlist.getOrCreateNode(tokens[4]);
            
            double w = 1e-6, l = 1e-6;
            // Simple param parsing: look for W= and L=
            for (size_t i = 6; i < tokens.size(); ++i) {
                std::string upperTok = toUpperCopy(tokens[i]);
                if (upperTok.rfind("W=", 0) == 0) w = Utils::parseValue(stripQuotes(tokens[i].substr(2)));
                if (upperTok.rfind("L=", 0) == 0) l = Utils::parseValue(stripQuotes(tokens[i].substr(2)));
            }

            int type = 1; // Default NMOS
            std::string model = tokens[5];
            std::transform(model.begin(), model.end(), model.begin(), ::toupper);
            if (model.find("PMOS") != std::string::npos) type = -1;

            netlist.addDevice(std::make_unique<Mosfet>(tokens[0], nD, nG, nS, nB, type, w, l));
        } else if (firstChar == 'V') {
            // Voltage Source: Vname N1 N2 [DC <value>] [AC <mag>] [PULSE(...)/SIN(...)/PWL(...)] or scalar value
            if (tokens.size() < 4) continue;
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
            if (tokens.size() < 4) continue;
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
        } else if (firstChar == 'D') {
            // Diode: Dname N1 N2
            if (tokens.size() < 3) continue;
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            netlist.addDevice(std::make_unique<Diode>(tokens[0], n1, n2));
        } else {
            netlist.addWarning("Line " + std::to_string(lineNo) + ": unsupported element ignored: " + line);
        }
    }

    return netlist;
}

} // namespace gspice
