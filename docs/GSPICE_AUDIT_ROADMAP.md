# GSPICE Audit And Roadmap

This document is the current engineering baseline for making GSPICE a serious simulator backend for Lumen Circuit Studio.

## Immediate Findings

- Basic MNA for primitive DC circuits now works on sanity decks. A 1 V / 1 k / 1 k divider solves to `V(in)=1 V` and `V(out)=0.5 V`.
- GSPICE now preserves node names in output, so Lumen can map simulator vectors to schematic nets more safely.
- Singular matrices now fail with an explicit error instead of returning a fake zero solution.
- Unsupported subcircuit instances now fail explicitly. This is important because silently ignoring `X...` PDK devices makes simulation results meaningless.
- Transient independent voltage and current sources now evaluate PULSE/SIN/PWL waveforms at the actual simulation time, not from history length.
- Current sources now transient-stamp dynamic PULSE/SIN/PWL waveforms instead of degrading them to DC.
- Waveform sources without explicit `DC` now use the waveform initial value for operating-point initialization.
- `.OPTIONS` now supports core solver controls: `RELTOL`, `VNTOL`, `ABSTOL`, `GMIN`, `ITL1`, and `ITL4`.
- OP and transient Newton loops now use SPICE-style relative/absolute convergence checks instead of hard-coded absolute voltage deltas.
- The parser now handles continuation lines, `.INCLUDE`/`.INC`, `.LIB file section` extraction, basic `.SUBCKT` expansion, nested subcircuits, and simple global `.PARAM` substitution.
- IHP SG13G2 low-voltage MOS wrappers (`sg13_lv_nmos`, `sg13_lv_pmos`) can now be reduced to GSPICE's simple primitive MOS so Lumen inverter smoke decks run, but GSPICE emits a warning because this is not PSP/PDK-accurate simulation.
- The primitive MOS now includes first-order nonideal behavior: channel-length modulation, body-effect threshold shift, weak subthreshold current, simple gate/junction parasitic capacitance, and source/drain-to-bulk junction diode clamps.
- GSPICE now has first-pass OSDI netlist plumbing: `.OSDI`/`.PRE_OSDI` loader directives, `.MODEL` card storage, and `N...` compact-model instance parsing. A built-in OSDI smoke model verifies the loader/model/device path, but real OpenVAF `.osdi` ABI compatibility and parameter binding are still incomplete.
- GSPICE now auto-searches for missing OSDI model libraries through `GSPICE_OSDI_DIR`, `NGSPICE_OSDI_DIR`, the deck directory, and local `osdi` folders. The IHP PSP helper script `tools/build_ihp_psp_osdi.ps1` compiles `psp103.va` and `psp103_nqs.va` with OpenVAF when OpenVAF is installed.
- The current executable build issue was caused by duplicate `Path`/`PATH` process environment variables. Building succeeds when the process environment is cleaned before invoking CMake/MSBuild.

## Current Critical Gaps

- `.SUBCKT`, `.INCLUDE`, `.LIB`, and `.PARAM` support is only a first-pass compatibility layer. Hierarchical parameter passing, expression evaluation, conditional library content, and many model-library constructs are still incomplete.
- PDK-grade MOS models are not implemented. The current primitive MOS is an improved but still simplified Level-1 style model, not BSIM/PSP/HICUM/EKV class modeling. The IHP SG13G2 wrapper fallback is only a temporary smoke-test bridge.
- OSDI support is not yet ngspice-compatible. GSPICE can exercise its own descriptor ABI, but it still needs OpenVAF-generated `.osdi` loading, model/instance parameter binding, charge/capacitance Jacobian handling, and validated PSP examples.
- The solver path is named KLU, but currently falls back to dense Gaussian elimination. Real sparse KLU/SuiteSparse or an equivalent production sparse solver is still required.
- Transient integration is basic backward-Euler style stamping for capacitors/inductors. It needs proper timestep control, LTE estimation, breakpoints, Gear/trapezoidal options, and robust source breakpoint handling.
- Nonlinear solve has basic configurable tolerances, but still lacks production convergence infrastructure: source stepping, gmin stepping, damped Newton/line search, device limiting, homotopy, and detailed convergence diagnostics.
- AC analysis only works where devices provide correct small-signal stamps. MOS AC stamping is currently omitted.
- RF/PSS/HB/PAC/STB paths are mostly prototypes or simplified demonstrations, not production analysis engines.
- Noise, distortion, sensitivity, transfer function, Monte Carlo, corners, measurements, and reliable raw-file output need significant work.
- There is no regression suite with golden vectors against ngspice/Xyce/Spectre-style expected results.

## Reference Targets

ngspice provides a mature SPICE3-derived baseline with OP, DC sweep, AC, TRAN, noise, pole-zero, sensitivity, transfer function, and experimental PSS, plus broad device coverage, XSPICE mixed-signal capability, model libraries, subcircuits, parameters, behavioral sources, and Verilog-A/OSDI-oriented model flows. Reference: https://ngspice.sourceforge.io/docs/ngspice-manual.pdf

Xyce provides a modern C++ simulator architecture with SPICE compatibility, high-performance parallel execution, large-circuit capability, formal documentation, regression testing guidance, Verilog-A/ADMS workflows, and advanced analyses. Reference: https://xyce.sandia.gov/documentation-tutorials/

## Required Roadmap

1. Make the primitive simulator honest and regression-tested.
   - Add golden tests for R, C, L, V, I, diode, and simple MOS circuits.
   - Add parser tests for units, comments, continuation lines, `.PARAM`, `.MODEL`, `.OP`, `.DC`, `.TRAN`, `.AC`, `.SAVE`, and `.PRINT`.
   - Emit stable ASCII RAW/CSV directly from GSPICE.

2. Build a real SPICE netlist front end.
   - Expand the first-pass `.SUBCKT` implementation into full hierarchical parameter passing, defaults, and instance overrides.
   - Expand `.LIB`/`.INCLUDE` handling to cover real foundry model-library syntax, corner includes, nested includes, and conditional content.
   - Implement parameter expressions with arithmetic, units, functions, and dependency ordering.
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
   - Implement or integrate BSIM/PSP-class MOS support.
   - Make Verilog-A/OSDI support compatible with OpenVAF-generated `.osdi` libraries.
   - Bind `.MODEL` and instance parameters into OSDI model instances.
   - Stamp OSDI charges/capacitances for transient and small-signal analyses.
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
