#ifndef GSPICE_SIMULATOR_CONTEXT_HPP
#define GSPICE_SIMULATOR_CONTEXT_HPP

#include "netlist.hpp"
#include "device.hpp"
#include "dae.hpp"
#include "device_arena.hpp"
#include "event_queue.hpp"
#include "result_database.hpp"

#include <memory>
#include <vector>
#include <string>

namespace gspice {

class SimulatorContext {
public:
    SimulatorContext() = default;
    ~SimulatorContext() = default;

    Netlist& netlist() { return netlist_; }
    const Netlist& netlist() const { return netlist_; }

    SimulationSettings& settings() { return settings_; }
    const SimulationSettings& settings() const { return settings_; }

    std::vector<std::unique_ptr<Device>>& devices() { return devices_; }
    const std::vector<std::unique_ptr<Device>>& devices() const { return devices_; }

    ResultDatabase& results() { return results_; }
    const ResultDatabase& results() const { return results_; }

    EventQueue& eventQueue() { return eventQueue_; }
    const EventQueue& eventQueue() const { return eventQueue_; }

    VectorReal& solution() { return solution_; }
    const VectorReal& solution() const { return solution_; }

    void setNodeCount(int count) {
        nodeCount_ = count;
        solution_ = VectorReal(count);
    }

    int nodeCount() const { return nodeCount_; }

    void reset() {
        devices_.clear();
        results_.clear();
        eventQueue_.clear();
        nodeCount_ = 0;
        solution_ = VectorReal();
        lastError_.clear();
    }

    void setLastError(const std::string& err) { lastError_ = err; }
    const std::string& getLastError() const { return lastError_; }

private:
    Netlist netlist_;
    SimulationSettings settings_;
    std::vector<std::unique_ptr<Device>> devices_;
    ResultDatabase results_;
    EventQueue eventQueue_;
    VectorReal solution_;
    int nodeCount_ = 0;
    std::string lastError_;
};

} // namespace gspice

#endif // GSPICE_SIMULATOR_CONTEXT_HPP
