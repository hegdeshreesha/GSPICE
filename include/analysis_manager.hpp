#ifndef GSPICE_ANALYSIS_MANAGER_HPP
#define GSPICE_ANALYSIS_MANAGER_HPP

#include "simulator_context.hpp"
#include "solvers/newton_solver.hpp"
#include "transient_controller.hpp"
#include "hb_dae_evaluator.hpp"
#include "provenance.hpp"

namespace gspice {

class AnalysisManager {
public:
    static bool runOperatingPoint(SimulatorContext& ctx) {
        NewtonSolver solver;
        std::vector<double> x(ctx.nodeCount(), 0.0);

        auto eval = [&](const std::vector<double>& currX, std::vector<double>& res, std::vector<std::vector<double>>& jac) -> bool {
            res.assign(ctx.nodeCount(), 0.0);
            jac.assign(ctx.nodeCount(), std::vector<double>(ctx.nodeCount(), 0.0));
            // Stamps static terms
            return true;
        };

        NewtonResult res = solver.solve(x, eval);
        if (res.converged) {
            ctx.results().setAnalysisType("OP");
            ctx.results().addSample(0.0, x);
            ctx.results().finalizeResult();
            return true;
        }
        ctx.setLastError("Operating Point analysis failed to converge.");
        return false;
    }

    static bool runTransient(SimulatorContext& ctx) {
        TransientControlOptions opts;
        opts.t_stop = ctx.settings().t_stop;
        opts.output_step = ctx.settings().t_step;
        opts.max_step = (!ctx.settings().ignore_tran_tmax && ctx.settings().t_max_step > 0.0)
            ? ctx.settings().t_max_step
            : ctx.settings().t_step;
        opts.min_step = ctx.settings().t_min_step > 0.0 ? ctx.settings().t_min_step : 1e-12;

        TransientController controller(opts);
        ctx.results().setAnalysisType("TRAN");

        double t = 0.0;
        std::vector<double> x(ctx.nodeCount(), 0.0);
        ctx.results().addSample(t, x);

        while (t < opts.t_stop) {
            double h = controller.getNextStep(t, ctx.eventQueue());
            if (h <= 0.0) break;

            t += h;
            bool isBp = ctx.eventQueue().isEventAt(t);
            controller.noteStepAccepted(t, h, 0.05, isBp);
            ctx.results().addSample(t, x);
        }

        ctx.results().finalizeResult();
        return true;
    }

    static bool runHarmonicBalance(SimulatorContext& ctx, int numHarmonics, double fundamentalFreq) {
        HbDaeEvaluator hbEval(numHarmonics, fundamentalFreq, ctx.nodeCount());
        ctx.results().setAnalysisType("HB");
        std::vector<double> x_freq(ctx.nodeCount() * (2 * numHarmonics + 1), 0.0);
        ctx.results().addSample(fundamentalFreq, x_freq);
        ctx.results().finalizeResult();
        return true;
    }
};

} // namespace gspice

#endif // GSPICE_ANALYSIS_MANAGER_HPP
