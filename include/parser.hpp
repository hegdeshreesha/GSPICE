#ifndef GSPICE_PARSER_HPP
#define GSPICE_PARSER_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "netlist.hpp"

namespace gspice {

class Parser {
public:
    static Netlist parse(const std::string& filePath);

private:
    static std::vector<std::string> tokenize(const std::string& line);
};

} // namespace gspice

#endif // GSPICE_PARSER_HPP
