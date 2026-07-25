#ifndef GSPICE_DEVICE_MANAGER_HPP
#define GSPICE_DEVICE_MANAGER_HPP

#include "device.hpp"
#include "device_arena.hpp"
#include "topology_engine.hpp"

#include <vector>
#include <memory>

namespace gspice {

class DeviceManager {
public:
    void addDevice(std::unique_ptr<Device> dev) {
        devices_.push_back(std::move(dev));
    }

    const std::vector<std::unique_ptr<Device>>& devices() const { return devices_; }
    std::vector<std::unique_ptr<Device>>& devices() { return devices_; }

    int count() const { return static_cast<int>(devices_.size()); }

    void prepareTopology(TopologyEngine& topology, int numNodes) {
        topology.reset(numNodes);
        // Identify zero-impedance branches (e.g. 0V DC voltage sources or ideal short switches)
        for (const auto& dev : devices_) {
            // Devices can register node merges here if applicable
        }
    }

    double collectBoundStep() const {
        double minBound = 1e30;
        for (const auto& dev : devices_) {
            // Query device step bounds if provided
        }
        return minBound < 1e20 ? minBound : 0.0;
    }

private:
    std::vector<std::unique_ptr<Device>> devices_;
    DeviceArena arena_;
};

} // namespace gspice

#endif // GSPICE_DEVICE_MANAGER_HPP
