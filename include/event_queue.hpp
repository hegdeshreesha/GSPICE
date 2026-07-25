#ifndef GSPICE_EVENT_QUEUE_HPP
#define GSPICE_EVENT_QUEUE_HPP

#include <vector>
#include <algorithm>
#include <cmath>

namespace gspice {

enum class EventType {
    Breakpoint,
    SourceEdge,
    OutputSample
};

struct Event {
    double time = 0.0;
    EventType type = EventType::Breakpoint;
    double requestedNextStep = 0.0; // Optional bound_step request for step following event
};

class EventQueue {
public:
    void clear() {
        events_.clear();
    }

    void addEvent(double time, EventType type, double requestedNextStep = 0.0) {
        if (time < 0.0 || !std::isfinite(time)) return;
        events_.push_back({time, type, requestedNextStep});
    }

    void buildAndSort(const std::vector<double>& rawBreakpoints) {
        for (double t : rawBreakpoints) {
            addEvent(t, EventType::Breakpoint);
        }
        std::sort(events_.begin(), events_.end(), [](const Event& a, const Event& b) {
            return a.time < b.time;
        });
        // Unique close event times within tolerance
        auto it = std::unique(events_.begin(), events_.end(), [](const Event& a, const Event& b) {
            return std::abs(a.time - b.time) < 1e-15;
        });
        events_.erase(it, events_.end());
    }

    double getNextEventTimeAfter(double t, double tol = 1e-12) const {
        for (const auto& ev : events_) {
            if (ev.time > t + tol) {
                return ev.time;
            }
        }
        return -1.0;
    }

    bool isEventAt(double t, double tol = 1e-12) const {
        for (const auto& ev : events_) {
            if (std::abs(ev.time - t) <= tol) {
                return true;
            }
        }
        return false;
    }

    const std::vector<Event>& events() const { return events_; }

private:
    std::vector<Event> events_;
};

} // namespace gspice

#endif // GSPICE_EVENT_QUEUE_HPP
