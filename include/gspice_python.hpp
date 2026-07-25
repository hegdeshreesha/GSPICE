#ifndef GSPICE_PYTHON_HPP
#define GSPICE_PYTHON_HPP

#include "simulator_context.hpp"
#include <cstddef>

namespace gspice {

struct SolutionBufferView {
    const double* dataPtr = nullptr;
    size_t sampleCount = 0;
    size_t variableCount = 0;
};

// Provides zero-copy buffer views into simulation results for Pybind11/NumPy
class PythonExportEngine {
public:
    static SolutionBufferView getSolutionBuffer(const SimulatorContext& ctx) {
        SolutionBufferView view;
        const auto& results = ctx.results().results();
        if (!results.empty()) {
            const auto& res = results.back();
            view.sampleCount = res.timeOrFreq.size();
            view.variableCount = res.variableNames.size();
            if (!res.values.empty() && !res.values[0].empty()) {
                view.dataPtr = res.values[0].data();
            }
        }
        return view;
    }
};

} // namespace gspice

#endif // GSPICE_PYTHON_HPP
