#ifndef GSPICE_TRANSIENT_CONTROLLER_HPP
#define GSPICE_TRANSIENT_CONTROLLER_HPP

#include "event_queue.hpp"
#include "transient_control.hpp"
#include "dae.hpp"

#include <cmath>
#include <algorithm>

namespace gspice {

struct TransientControlOptions {
    double t_start = 0.0;
    double t_stop = 1e-3;
    double output_step = 1e-6;
    double max_step = 1e-5;
    double min_step = 1e-12;
    bool adaptive = true;
    double reltol = 1e-3;
    double vntol = 1e-6;
};

class TransientController {
public:
    explicit TransientController(TransientControlOptions opts)
        : opts_(opts), candidateStep_(opts.output_step) {}

    double getNextStep(double currentTime, const EventQueue& events) {
        double remaining = opts_.t_stop - currentTime;
        if (remaining <= 0.0) return 0.0;

        double cur = candidateStep_;
        cur = std::min(cur, remaining);

        double nextEventTime = events.getNextEventTimeAfter(currentTime);
        if (nextEventTime > 0.0) {
            double dist = nextEventTime - currentTime;
            if (dist > opts_.min_step && dist < cur) {
                cur = dist;
            }
        }
        return std::max(opts_.min_step, cur);
    }

    void noteStepAccepted(double currentTime, double usedStep, double lteError, bool isBreakpoint) {
        if (isBreakpoint) {
            // Restore step size immediately post-breakpoint landing
            candidateStep_ = opts_.adaptive ? std::min(opts_.max_step, opts_.output_step * 2.0) : opts_.output_step;
        } else if (opts_.adaptive) {
            double growthFactor = 1.2;
            if (lteError < 0.1) growthFactor = 1.5;
            else if (lteError > 0.8) growthFactor = 0.9;
            candidateStep_ = std::min(opts_.max_step, std::max(opts_.min_step, candidateStep_ * growthFactor));
        } else {
            candidateStep_ = opts_.output_step;
        }
    }

    void noteStepRejected() {
        candidateStep_ = std::max(opts_.min_step, candidateStep_ * 0.5);
    }

private:
    TransientControlOptions opts_;
    double candidateStep_;
};

} // namespace gspice

#endif // GSPICE_TRANSIENT_CONTROLLER_HPP
