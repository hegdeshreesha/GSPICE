# GSPICE Audit And Roadmap

This document is the current engineering baseline for making GSPICE a serious simulator backend for Lumen Circuit Studio.

## Immediate Findings

- Basic MNA for primitive DC circuits now works on sanity decks. A 1 V / 1 k / 1 k divider solves to `V(in)=1 V` and `V(out)=0.5 V`.
- GSPICE now preserves node names in output, so Lumen can map simulator vectors to schematic nets more safely.
- Singular matrices now fail with an explicit error instead of returning a fake zero solution.
- Unsupported subcircuit instances now fail explicitly. This is important because silently ignoring `X...` PDK devices makes simulation results meaningless.
- The current executable build issue was caused by duplicate `Path`/`PATH` process environment variables. Building succeeds when the process environment is cleaned before invoking CMake/MSBuild.

## Current Critical Gaps

- No real `.SUBCKT` expansion engine yet.
- `.LIB`, `.INCLUDE`, and foundry model library parsing are not implemented.
- PDK-grade MOS models are not implemented. The current primitive MOS is a simple Level-1 style model, not BSIM/PSP/HICUM/EKV class modeling.
- The solver path is named KLU, but currently falls back to dense Gaussian elimination. Real sparse KLU/SuiteSparse or an equivalent production sparse solver is still required.
- Transient integration is basic backward-Euler style stamping for capacitors/inductors. It needs proper timestep control, LTE estimation, breakpoints, Gear/trapezoidal options, and robust source breakpoint handling.
- Nonlinear solve lacks production convergence infrastructure: source stepping, gmin stepping, damped Newton/line search, device limiting, homotopy, and detailed convergence diagnostics.
- AC analysis only works where devices provide correct small-signal stamps. MOS AC stamping is currently omitted.
- RF/PSS/HB/PAC/STB paths are mostly prototypes or simplified demonstrations, not production analysis engines.
- Noise, distortion, sensitivity, transfer function, Monte Carlo, corners, measurements, and reliable raw-file output need significant work.
- Current sources parse dynamic waveforms but do not transient-stamp those waveforms yet.
- There is no regression suite with golden vectors against ngspice/Xyce/Spectre-style expected results.

## Reference Targets

ngspice provides a mature SPICE3-derived baseline with OP, DC sweep, AC, TRAN, noise, pole-zero, sensitivity, transfer function, and experimental PSS, plus broad device coverage and XSPICE mixed-signal capability.

Xyce provides a modern C++ simulator architecture with DC, transient, AC, noise, harmonic balance, sensitivity, uncertainty propagation, Verilog-A support, PDK translation support, and large-scale parallel simulation through Trilinos/KLU-class solver infrastructure.

## Required Roadmap

1. Make the primitive simulator honest and regression-tested.
   - Add golden tests for R, C, L, V, I, diode, and simple MOS circuits.
   - Add parser tests for units, comments, continuation lines, `.PARAM`, `.MODEL`, `.OP`, `.DC`, `.TRAN`, `.AC`, `.SAVE`, and `.PRINT`.
   - Emit stable ASCII RAW/CSV directly from GSPICE.

2. Build a real SPICE netlist front end.
   - Implement `.SUBCKT` definition/instantiation expansion.
   - Implement `.LIB`/`.INCLUDE` resolution with section/corner selection.
   - Implement parameter expressions and hierarchical parameter passing.
   - Add clear diagnostics for unsupported syntax.

3. Replace prototype solver behavior.
   - Integrate an actual sparse direct solver.
   - Add singularity detection, matrix conditioning diagnostics, and node/device contribution reports.
   - Add nonlinear convergence controls and robust Newton infrastructure.

4. Implement production transient.
   - Add trapezoidal and Gear methods.
   - Add adaptive timestepping, local truncation error control, and source breakpoints.
   - Add correct branch-current output for voltage sources, inductors, and selected terminals.

5. Add compact-model infrastructure.
   - Implement or integrate BSIM-class MOS support.
   - Make Verilog-A/OSDI support real and tested.
   - Add diode/BJT/JFET/MESFET/HBT/HICUM/PSP model support as needed by Lumen PDKs.

6. Add advanced analyses only after the core is reliable.
   - DC sweep, AC, noise, sensitivity, transfer function, pole-zero.
   - Harmonic balance and PSS with validated convergence behavior.
   - S-parameter and RF analyses with tested port normalization and Touchstone output.

7. Continuous validation.
   - Compare every supported analysis against ngspice and Xyce on small public decks.
   - Keep a Lumen smoke suite that verifies netlist generation, GSPICE execution, waveform file creation, and SigView loading.
   - Mark unsupported features as unsupported in the UI rather than silently falling back.

## Non-Negotiable Rule

GSPICE must never silently ignore an active device or source. If it cannot simulate a line correctly, it must fail loudly with a useful diagnostic.
