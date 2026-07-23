# GSPICE, VACASK, and Xyce comparison

Audit date: 2026-07-21 (America/Los_Angeles)

## Revisions audited

| Project | Upstream revision | Upstream location |
|---|---:|---|
| GSPICE | `046d69eb604afe91943652920832cd04b6c73d5c` (`main`, v1.2.0-beta baseline) | <https://github.com/hegdeshreesha/GSPICE> |
| VACASK | `ae9533635e26187bb57435189c09499014de59dd` (`main`) | <https://codeberg.org/arpadbuermen/VACASK> |
| Xyce | `d72b5846a0397ddf852a49305cb6f395457685ca` (`master`) | <https://github.com/Xyce/Xyce> |

The GSPICE column below describes the 1.3 academic-beta work in this tree, not
only the public baseline commit.

## Capability comparison

| Area | GSPICE 1.3 academic beta | VACASK main | Xyce master |
|---|---|---|---|
| Primary goal | Compact teaching/research core and Lumen backend | Extensible Verilog-A/OSDI simulator | Production-scale parallel SPICE-compatible simulator |
| Netlist language | Useful SPICE subset; strict rejection for unsupported active syntax | Native Spectre-like hierarchical language; SPICE converter under development | Broad SPICE compatibility plus XDM translation |
| Core analyses | OP, DC, TRAN, AC are beta on supported devices | OP, DC incremental/transfer, AC, stability, S-parameter, noise, TRAN, HB/HBAC, MC | Mature DC, TRAN, AC, noise plus HB, sensitivity, uncertainty and more |
| Periodic analyses | HB experimental; PSS/PAC/PNoise rejected | HB implemented; time-domain noise supported | HB and mature advanced-analysis infrastructure |
| Compact models | Experimental OSDI plus simplified built-ins | OSDI 0.4/OpenVAF-reloaded; PSP, BSIM3/4/BULK, VBIC and more | Large native/ADMS model library; some installer-only proprietary models |
| Nonlinear solve | Newton, residual check, damping, PN-junction/OSDI limiting, guarded device bypass, nodesets, adaptive source/gmin continuation | Residual convergence, bypass, several homotopies, nodesets | Mature solver-management and continuation infrastructure |
| Transient | Simulator-owned Q/Qdot history; transactional BE/trapezoid/Gear2 and variable-step BDF/Adams through order 5; polynomial predictor-corrector LTE, adaptive order, periodic step-doubling oracle, charge-aware LTE, breakpoints, trap-ringing damping | Predictor-corrector LTE, BE/trapezoid/Gear | Mature variable-step/order integration and restart infrastructure |
| Linear algebra | Optional KLU; internal sparse fallback | KLU | Trilinos ecosystem, KLU and distributed solvers |
| Parallelism | OpenMP stamping; one process | Primarily optimized single-process execution | MPI/distributed large-scale execution |
| Statistical analysis | One-source Gaussian/uniform MC, optional LHS; simple corners/spec yield | Parameter/model/instance sweeps and current MC/LHS framework | Sampling and uncertainty-quantification workflows |
| Output/API | ASCII RAW, CSV, CLI; no stable library API | ASCII/binary RAW, Python tools, linkable simulator library | Multiple output formats, C/Python/Matlab and coupling interfaces |
| Demonstrated scale | Small regression decks; no credible large-scale claim | Published 10k-transistor and multi-million-device multiplier scaling data | Designed and validated for circuits with millions of devices and MPI |
| Test/release maturity | 65 configured CTests on the audit host; native/OSDI DAE derivative-conservation audits, minimum-step starvation coverage, and two cross-simulator transient oracles | Broad system-test/demo tree and multi-OS CI | Large regression process; master is published after regression passes |
| License | Apache-2.0 | AGPL-3.0 | GPL-3.0 |

## Where GSPICE is currently better suited

- Small, understandable C++17 codebase for coursework and simulator research.
- Strict failure policy for unsupported active syntax and unimplemented
  analyses, reducing the chance of silently wrong teaching results.
- Lightweight build compared with Xyce's Trilinos-based stack.
- Direct Lumen integration, readable CSV/RAW export, deterministic sampling,
  and machine-readable capability maturity.
- Permissive Apache-2.0 licensing for academic and commercial embedding.

These are scope and integration advantages, not claims of numerical or model
parity.

## Where GSPICE is behind

VACASK is substantially ahead in OSDI fidelity, hierarchical parameterization,
arbitrary sweeps, compact-model breadth, transient/noise/HB algorithms,
runtime circuit alteration, binary RAW output, and demonstrated performance.
Its newest upstream series also adds Monte Carlo generators, Latin-hypercube
sampling, SPICE-to-VACASK conversion improvements, and cross-platform release
automation.

Xyce is substantially ahead in compact models, parser compatibility, solver
robustness, device limiting, sensitivity methods, uncertainty quantification,
parallel/distributed execution, restart/coupling APIs, test depth,
documentation, and validated large-circuit use.

## Reproducible transient baselines

Run `python tools/validate_transient.py --gspice <gspice> --xyce <Xyce>
--vacask <vacask> --source .` from the repository root. On the Windows audit
host, GSPICE 1.3 academic-beta was compared with Xyce 7.10 and the official
VACASK 0.3.4-rc1 binary:

| Deck | Reference | Maximum absolute difference | RMS difference |
|---|---:|---:|---:|
| 1 ms linear RC step, 501 GSPICE samples | analytic closed form | 2.77 uV | 1.67 uV |
| 1 ms linear RC step, 501 GSPICE samples | Xyce 7.10 | 44.2 uV | 8.52 uV |
| 1 ms linear RC step, 501 GSPICE samples | VACASK 0.3.4-rc1 | 43.8 uV | 26.1 uV |
| PSP103.4 CMOS inverter, 501 GSPICE samples | VACASK 0.3.4-rc1 | 7.77 mV | 1.38 mV |

These results establish agreement on two small public topologies only. The
VACASK executable predates the audited `main` revision, and the comparison does
not qualify temperature corners, large circuits, noise, AC compact-model
charge, discontinuous devices, or foundry model cards.

## How close GSPICE is to VACASK

Numerically, GSPICE is now close on the two transient baselines above. Its core
transient and nonlinear algorithm set now contains corresponding classes of
predictor-corrector control, order adaptation, bypass, nodesets, and
continuation, but those new paths do not yet have VACASK's validation depth.
As a simulator platform it is not close to feature or production parity:
VACASK has a substantially richer hierarchical language, OSDI 0.4 integration,
compact model library, broader homotopies, arbitrary
sweeps and runtime alteration, HB/HBAC, time-domain noise, Monte Carlo/LHS,
binary output, a library API, and demonstrated multi-million-device scale.
GSPICE should therefore be released as a small-circuit academic beta, not as a
VACASK replacement.

## Development direction

Do not attempt a one-week feature-count race. The credible path is:

1. Expand numeric golden comparisons against ngspice and Xyce.
2. Validate OSDI charge, Jacobian, limiting, AC, transient, and noise behavior.
3. Separate parser, analysis manager, nonlinear solver, device manager, and
   output layers so they can be tested independently.
4. Generalize sweeps/Monte Carlo to model and instance parameters.
5. Implement periodic engines only with convergence criteria and reference
   oracles; keep unsupported directives as hard failures until then.
6. Publish scale and accuracy data only after reproducible benchmark harnesses
   exist.
