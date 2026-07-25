#ifndef GSPICE_CIRCUIT_ELABORATOR_HPP
#define GSPICE_CIRCUIT_ELABORATOR_HPP

#include "simulator_context.hpp"
#include "parser.hpp"

namespace gspice {

class CircuitElaborator {
public:
    static bool elaborate(SimulatorContext& ctx, const std::string& netlistContent) {
        ctx.reset();
        // Write temporary file for parser
        std::string tmpPath = "tmp_elaborate.sp";
        {
            std::ofstream out(tmpPath);
            out << netlistContent;
        }
        ctx.netlist() = Parser::parse(tmpPath);
        std::remove(tmpPath.c_str());

        ctx.setNodeCount(ctx.netlist().getNumNodes());
        ctx.settings() = ctx.netlist().getSettings();
        return true;
    }
};

} // namespace gspice

#endif // GSPICE_CIRCUIT_ELABORATOR_HPP
