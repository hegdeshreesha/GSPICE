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
#include <iostream>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>

namespace {

std::string toUpperCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
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

bool tryParseSpiceValue(const std::string& token, double& out) {
    try {
        out = gspice::Utils::parseValue(token);
        return true;
    } catch (...) {
        return false;
    }
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

void parseSourceSpec(
    const std::string& sourceSpec,
    double& dcValue,
    double& acMagnitude,
    gspice::VoltageSource::WaveformType& waveformType,
    gspice::VoltageSource::PulseParams& pulse,
    gspice::VoltageSource::SinParams& sin,
    std::vector<double>& pwlTimes,
    std::vector<double>& pwlValues) {

    dcValue = 0.0;
    acMagnitude = 1.0;
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

    bool dcSeen = false;
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
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return netlist;
    }

    std::string line;
    // SPICE files often have a title as the first line
    std::getline(file, line); 
    int lineNo = 1;

    while (std::getline(file, line)) {
        ++lineNo;
        if (line.empty() || line[0] == '*' || line[0] == '$') continue;

        auto tokens = tokenize(line);
        if (tokens.empty()) continue;

        char firstChar = std::toupper(tokens[0][0]);
        
        if (firstChar == '.') {
            // Commands
            std::string cmd = tokens[0];
            std::transform(cmd.begin(), cmd.end(), cmd.begin(), ::toupper);
            
            if (cmd == ".END") {
                break;
            } else if (cmd == ".OP") {
                if (netlist.getSettings().type == "OP") {
                    SimulationSettings settings;
                    settings.type = "OP";
                    netlist.setSettings(settings);
                } else {
                    netlist.addWarning(
                        "Line " + std::to_string(lineNo) +
                        ": .OP kept as operating-point initialization; active analysis remains ." +
                        netlist.getSettings().type);
                }
            } else if (cmd == ".LIB" || cmd == ".INCLUDE" || cmd == ".INC") {
                netlist.addWarning("Line " + std::to_string(lineNo) + ": model/include directive is not implemented yet: " + line);
            } else if (cmd == ".TRAN") {
                if (tokens.size() < 3) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .TRAN line ignored: " + line);
                    continue;
                }
                SimulationSettings settings;
                settings.type = "TRAN";
                settings.t_step = Utils::parseValue(tokens[1]);
                settings.t_stop = Utils::parseValue(tokens[2]);
                if (tokens.size() > 3) {
                    std::string opt = tokens[3];
                    std::transform(opt.begin(), opt.end(), opt.begin(), ::toupper);
                    if (opt == "UIC") settings.use_uic = true;
                }
                netlist.setSettings(settings);
            } else if (cmd == ".MODEL") {
                // .MODEL modelName OSDI file="path"
                if (tokens.size() > 3 && tokens[2] == "OSDI") {
                     std::string pathPart = tokens[3];
                     size_t start = pathPart.find('\"');
                     size_t end = pathPart.find_last_of('\"');
                     if (start != std::string::npos && end != std::string::npos) {
                         std::string path = pathPart.substr(start+1, end-start-1);
                         netlist.addOsdiModel(tokens[1], path);
                     }
                }
            } else if (cmd == ".AC") {
                if (tokens.size() < 5) {
                    netlist.addWarning("Line " + std::to_string(lineNo) + ": invalid .AC line ignored: " + line);
                    continue;
                }
                SimulationSettings settings;
                settings.type = "AC";
                // tokens[1] is DEC/OCT/LIN
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".PSS" || cmd == ".HB") {
                SimulationSettings settings;
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
                SimulationSettings settings;
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
                SimulationSettings settings;
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
                SimulationSettings settings;
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
                SimulationSettings settings;
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
                SimulationSettings settings;
                settings.type = "PNOISE";
                std::string outNode = tokens[1].substr(2, tokens[1].size()-3);
                settings.out_node = netlist.getOrCreateNode(outNode);
                settings.points_per_dec = std::stoi(tokens[3]);
                settings.f_start = Utils::parseValue(tokens[4]);
                settings.f_stop = Utils::parseValue(tokens[5]);
                netlist.setSettings(settings);
            } else if (cmd == ".HBAC" || cmd == ".HBNOISE" || cmd == ".HBSP" || cmd == ".HBSTB" || 
                       cmd == ".PSSSP" || cmd == ".PSSSTB") {
                SimulationSettings settings;
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
                if (tokens[i].substr(0, 2) == "W=") w = Utils::parseValue(tokens[i].substr(2));
                if (tokens[i].substr(0, 2) == "L=") l = Utils::parseValue(tokens[i].substr(2));
            }

            int type = 1; // Default NMOS
            std::string model = tokens[5];
            std::transform(model.begin(), model.end(), model.begin(), ::toupper);
            if (model == "PMOS") type = -1;

            netlist.addDevice(std::make_unique<Mosfet>(tokens[0], nD, nG, nS, nB, type, w, l));
        } else if (firstChar == 'V') {
            // Voltage Source: Vname N1 N2 [DC <value>] [AC <mag>] [PULSE(...)/SIN(...)/PWL(...)] or scalar value
            if (tokens.size() < 4) continue;
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            const std::string sourceSpec = joinTokens(tokens, 3);
            double dcValue = 0.0;
            double acMagnitude = 1.0;
            VoltageSource::WaveformType wf = VoltageSource::WaveformType::DC;
            VoltageSource::PulseParams pulse;
            VoltageSource::SinParams sin;
            std::vector<double> pwlT;
            std::vector<double> pwlV;
            parseSourceSpec(sourceSpec, dcValue, acMagnitude, wf, pulse, sin, pwlT, pwlV);

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
            VoltageSource::WaveformType wf = VoltageSource::WaveformType::DC;
            VoltageSource::PulseParams pulse;
            VoltageSource::SinParams sin;
            std::vector<double> pwlT;
            std::vector<double> pwlV;
            parseSourceSpec(sourceSpec, dcValue, acMagnitude, wf, pulse, sin, pwlT, pwlV);
            if (wf != VoltageSource::WaveformType::DC) {
                netlist.addWarning(
                    "Line " + std::to_string(lineNo) + ": dynamic current-source waveform on " + tokens[0] +
                    " is parsed but transient waveform stamping is not implemented yet; using DC value.");
            }
            netlist.addDevice(std::make_unique<CurrentSource>(tokens[0], n1, n2, dcValue));
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
