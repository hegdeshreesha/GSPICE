#ifndef GSPICE_DEVICE_ARENA_HPP
#define GSPICE_DEVICE_ARENA_HPP

// ---------------------------------------------------------------------------
// DeviceStateArena — Contiguous device state memory pool (Pillar 2/5).
//
// Problem:
//   Each OSDIDevice currently stores prev_state_, next_state_, prev_react_,
//   prev2_react_, and prev_react_derivative_ as separate std::vector<double>
//   members. For a circuit with 10,000 PSP103.4 instances (each with ~50
//   state variables), this creates over 1,000,000 distinct heap allocations.
//   These allocations are scattered across the heap, destroying cache locality
//   when the Newton loop iterates over all devices in a tight loop.
//
// Solution:
//   Pre-allocate a single contiguous slab of memory at elaboration time.
//   Each device receives a pointer into a fixed sub-region of this slab.
//   Temporal access patterns are now cache-friendly: all device states for
//   iteration k are laid out consecutively in memory.
//
// Layout (for D devices, each with S state variables):
//
//   [device 0 state | device 1 state | ... | device D-1 state]
//    ↑ offset[0]      ↑ offset[1]            ↑ offset[D-1]
//
//   Multiple "pools" (e.g. prev, next, scratch) share the same offset[] map
//   but point to different contiguous slabs.
//
// Usage:
//   DeviceStateArena arena;
//   arena.allocate(device_idx, num_states);   // called for each device
//   arena.seal();                             // finalise layout
//   double* my_state = arena.ptr(pool_id, device_idx);
// ---------------------------------------------------------------------------

#include <cassert>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace gspice {

class DeviceStateArena {
public:
    // Number of parallel pools to allocate (e.g. prev, next, scratch).
    explicit DeviceStateArena(int num_pools = 3)
        : num_pools_(num_pools), sealed_(false), total_scalars_(0) {}

    // -----------------------------------------------------------------------
    // Elaboration phase — call once per device before seal().
    // -----------------------------------------------------------------------

    /// Register device `idx` as needing `count` scalar state variables.
    /// Must be called in device-index order (0, 1, 2, ...).
    void allocate(std::size_t device_idx, std::size_t count) {
        if (sealed_) throw std::logic_error("DeviceStateArena: allocate() after seal()");
        if (device_idx != offsets_.size()) {
            throw std::logic_error("DeviceStateArena: allocate() called out of order");
        }
        offsets_.push_back(total_scalars_);
        sizes_.push_back(count);
        total_scalars_ += count;
    }

    /// Finalise the layout and allocate all pools in a single calloc.
    void seal() {
        if (sealed_) return;
        sealed_ = true;
        const std::size_t pool_bytes = total_scalars_ * sizeof(double);
        // Allocate all pools in one contiguous block for maximum locality.
        data_.assign(static_cast<std::size_t>(num_pools_) * total_scalars_, 0.0);
    }

    // -----------------------------------------------------------------------
    // Runtime access — O(1) pointer arithmetic, no hash lookups.
    // -----------------------------------------------------------------------

    /// Raw pointer to pool `pool_id`'s block for device `device_idx`.
    double* ptr(int pool_id, std::size_t device_idx) {
        assert(sealed_);
        assert(pool_id >= 0 && pool_id < num_pools_);
        assert(device_idx < offsets_.size());
        return data_.data() +
               static_cast<std::size_t>(pool_id) * total_scalars_ +
               offsets_[device_idx];
    }

    const double* ptr(int pool_id, std::size_t device_idx) const {
        assert(sealed_);
        assert(pool_id >= 0 && pool_id < num_pools_);
        assert(device_idx < offsets_.size());
        return data_.data() +
               static_cast<std::size_t>(pool_id) * total_scalars_ +
               offsets_[device_idx];
    }

    /// Number of state scalars allocated for device `device_idx`.
    std::size_t size(std::size_t device_idx) const {
        assert(device_idx < sizes_.size());
        return sizes_[device_idx];
    }

    std::size_t numDevices()    const { return offsets_.size(); }
    std::size_t totalScalars()  const { return total_scalars_; }
    int         numPools()      const { return num_pools_; }
    bool        isSealed()      const { return sealed_; }

    /// Copy pool `src` into pool `dst` (e.g. commit accepted → previous).
    void copyPool(int src, int dst) {
        assert(sealed_);
        assert(src >= 0 && src < num_pools_);
        assert(dst >= 0 && dst < num_pools_);
        std::memcpy(data_.data() + static_cast<std::size_t>(dst) * total_scalars_,
                    data_.data() + static_cast<std::size_t>(src) * total_scalars_,
                    total_scalars_ * sizeof(double));
    }

    /// Zero a single pool.
    void zeroPool(int pool_id) {
        assert(sealed_);
        std::memset(data_.data() + static_cast<std::size_t>(pool_id) * total_scalars_,
                    0,
                    total_scalars_ * sizeof(double));
    }

private:
    int                   num_pools_;
    bool                  sealed_;
    std::size_t           total_scalars_;
    std::vector<std::size_t> offsets_; // per-device start in pool
    std::vector<std::size_t> sizes_;   // per-device scalar count
    std::vector<double>   data_;       // all pools, contiguous
};

} // namespace gspice

#endif // GSPICE_DEVICE_ARENA_HPP
