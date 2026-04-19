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
#include <iostream>

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

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '*' || line[0] == '$') continue;

        auto tokens = tokenize(line);
        if (tokens.empty()) continue;

        char firstChar = std::toupper(tokens[0][0]);
        
        if (firstChar == '.') {
            // Commands
            std::string cmd = tokens[0];
            std::transform(cmd.begin(), cmd.end(), cmd.begin(), ::toupper);
            
            if (cmd == ".TRAN") {
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
                SimulationSettings settings;
                settings.type = "AC";
                // tokens[1] is DEC/OCT/LIN
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".PSS") {
                SimulationSettings settings;
                settings.type = "PSS";
                settings.f_fund = Utils::parseValue(tokens[1]);
                settings.n_harms = std::stoi(tokens[2]);
                netlist.setSettings(settings);
            } else if (cmd == ".SP") {
                SimulationSettings settings;
                settings.type = "SP";
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".NOISE") {
                SimulationSettings settings;
                settings.type = "NOISE";
                // tokens[1] is V(node)
                std::string outNode = tokens[1].substr(2, tokens[1].size()-3);
                settings.out_node = netlist.getOrCreateNode(outNode);
                settings.points_per_dec = std::stoi(tokens[3]);
                settings.f_start = Utils::parseValue(tokens[4]);
                settings.f_stop = Utils::parseValue(tokens[5]);
                netlist.setSettings(settings);
            } else if (cmd == ".HB") {
                SimulationSettings settings;
                settings.type = "HB";
                settings.f_fund = Utils::parseValue(tokens[1]);
                settings.n_harms = std::stoi(tokens[2]);
                netlist.setSettings(settings);
            } else if (cmd == ".STB") {
                SimulationSettings settings;
                settings.type = "STB";
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".PAC") {
                SimulationSettings settings;
                settings.type = "PAC";
                settings.points_per_dec = std::stoi(tokens[2]);
                settings.f_start = Utils::parseValue(tokens[3]);
                settings.f_stop = Utils::parseValue(tokens[4]);
                netlist.setSettings(settings);
            } else if (cmd == ".PNOISE") {
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
                settings.type = cmd.substr(1); // Remove dot
                settings.points_per_dec = std::stoi(tokens[1]);
                settings.f_start = Utils::parseValue(tokens[2]);
                settings.f_stop = Utils::parseValue(tokens[3]);
                netlist.setSettings(settings);
            }

        } else if (firstChar == 'R') {
            // Resistor: Rname N1 N2 Value
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Resistor>(tokens[0], n1, n2, val));
        } else if (firstChar == 'C') {
            // Capacitor: Cname N1 N2 Value
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Capacitor>(tokens[0], n1, n2, val));
        } else if (firstChar == 'L') {
            // Inductor: Lname N1 N2 Value
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            netlist.addDevice(std::make_unique<Inductor>(tokens[0], n1, n2, val, -1));
        } else if (firstChar == 'W') {
            // Stability Probe: Wname node_in node_out
            int nIn = netlist.getOrCreateNode(tokens[1]);
            int nOut = netlist.getOrCreateNode(tokens[2]);
            netlist.addDevice(std::make_unique<StabilityProbe>(tokens[0], nIn, nOut, -1));
        } else if (firstChar == 'P') {
            // Port: Pname N1 N2 PortNum Z0
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            int pNum = std::stoi(tokens[3]);
            double z0 = 50.0;
            if (tokens.size() > 4) z0 = Utils::parseValue(tokens[4]);
            netlist.addDevice(std::make_unique<Port>(tokens[0], n1, n2, pNum, z0));
        } else if (firstChar == 'M') {
            // MOSFET: Mname D G S B Model [W=..] [L=..]
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
            // Voltage Source: Vname N1 N2 Value
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            double val = Utils::parseValue(tokens[3]);
            // Branch current ID is assigned later during assembly, or we can handle it here
            // For now, let's keep it simple.
            netlist.addDevice(std::make_unique<VoltageSource>(tokens[0], n1, n2, val, -1)); 
        } else if (firstChar == 'D') {
            // Diode: Dname N1 N2
            int n1 = netlist.getOrCreateNode(tokens[1]);
            int n2 = netlist.getOrCreateNode(tokens[2]);
            netlist.addDevice(std::make_unique<Diode>(tokens[0], n1, n2));
        }
    }

    return netlist;
}

} // namespace gspice
