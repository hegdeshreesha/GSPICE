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
- The parser now handles continuation lines, `.INCLUDE`/`.INC`, `.LIB file section` extraction, basic `.SUBCKT` expansion, nested subcircuits, `.IF`/`.ELSEIF`/`.ELSE`/`.ENDIF` conditionals, and numeric `.PARAM` expressions with units and common math functions.
- Numeric scalar parsing is now strict in the parser: expressions such as `1k*2` are evaluated as expressions instead of being truncated to `1k`.
- Unsupported active elements and unresolved simulator directives now fail loudly instead of being silently ignored. Output-selection directives such as `.SAVE`, `.PRINT`, `.PROBE`, and `.PLOT` are accepted as compatibility warnings while GSPICE writes all solved vectors.
- IHP SG13G2 low-voltage MOS wrappers (`sg13_lv_nmos`, `sg13_lv_pmos`) now prefer the real active PSP/OSDI wrapper branch and preserve calculated instance parameters such as geometry-dependent area/perimeter values. Primitive MOS fallback is disabled by default and requires `GSPICE_ALLOW_PRIMITIVE_IHP_FALLBACK=1` for placeholder smoke decks.
- The primitive MOS now includes first-order nonideal behavior: channel-length modulation, body-effect threshold shift, weak subthreshold current, simple gate/junction parasitic capacitance, and source/drain-to-bulk junction diode clamps.
- GSPICE now has first-pass OSDI netlist plumbing: `.OSDI`/`.PRE_OSDI` loader directives, `.MODEL` card storage, and `N...` compact-model instance parsing. A built-in OSDI smoke model verifies the loader/model/device path, but real OpenVAF `.osdi` ABI compatibility and parameter binding are still incomplete.
- GSPICE now auto-searches for missing OSDI model libraries through `GSPICE_OSDI_DIR`, `NGSPICE_OSDI_DIR`, the deck directory, and local `osdi` folders. The IHP PSP helper script `tools/build_ihp_psp_osdi.ps1` compiles `psp103.va` and `psp103_nqs.va` with OpenVAF when OpenVAF is installed.
- The current executable build issue was caused by duplicate `Path`/`PATH` process environment variables. Building succeeds when the process environment is cleaned before invoking CMake/MSBuild.
- SuiteSparse/KLU is now an optional real sparse backend through vcpkg, with deployment scripts and smoke tests that verify external KLU calls.
- `.DC` single-source sweeps are implemented and reuse the previous sweep solution as the next Newton initial guess.
- `.NOISE` now performs first-pass output-referred noise analysis by solving each equivalent device noise source through the AC matrix. Resistor thermal noise, diode shot noise, and primitive MOS channel noise are represented.
- Primitive MOS devices now stamp first-pass `gm/gds` small-signal AC behavior.
- Primitive diode and MOS instances now consume useful `.MODEL` parameters. Compact model syntax such as `.MODEL DNOM D(IS=... N=... CJO=...)` is parsed correctly.
- `G` VCCS and `E` VCVS controlled-source primitives are implemented.
- `F` CCCS and `H` CCVS current-controlled-source primitives are implemented with branch-current resolution.
- `.TF V(out) source` transfer-function analysis is implemented as a first-pass small-signal solve around the DC operating point.
- `.DC` now supports nested source sweeps such as `.DC V1 ... V2 ...`.
- DC operating-point solving now has automatic source stepping and gmin stepping recovery when direct Newton fails.
- DC operating-point solving now has first-pass line-search damping to reduce large Newton overshoot on difficult nonlinear decks.
- A primitive BJT (`Q`) device is implemented with `.MODEL NPN/PNP` parameter parsing, DC linearization, AC Jacobian stamping, and simple shot-noise hooks.
- `.STEP` is implemented for OP-style independent-source sweeps, including nested sweep groups.
- `.MC` is implemented as a first-pass Gaussian independent-source Monte Carlo analysis with deterministic seeding and summary statistics.
- `.SENS V(out) source` is implemented as a first-pass finite-difference DC source sensitivity analysis.
- `.PZ V(out) source` is implemented as a frequency sweep based pole/zero estimator around the DC operating point. It now reports multiple falling/rising gain-threshold estimates and a dominant pole estimate, but it is not yet a generalized eigenvalue pole-zero solver.
- `.MEASURE`/`.MEAS` supports first-pass transient voltage measurements for `FIND`, `MIN`, `MAX`, `AVG`, `RMS`, and `PP` with `AT`, `FROM`, and `TO` qualifiers.
- Behavioral `B` sources are implemented for `I={expression}` and `V={expression}` with nonlinear numerical Jacobian stamping. Expressions support `V(node)`, `V(node,ref)`, `I(branchDevice)`, `time`, arithmetic, comparisons, ternary `?:`, and common functions such as `sin`, `cos`, `exp`, `log`, `sqrt`, `abs`, `min`, `max`, `pow`, `limit`, and `if`.
- Corner sweeps are implemented for operating-point source assignments using `.CORNER name SRC=value ...`.
- Production specs are implemented with `.SPEC name V(node) MIN=value MAX=value`; OP/corner runs report pass/fail and Monte Carlo reports yield, mean, sigma, min, and max.
- OSDI compact-model fidelity improved: compact-model devices stamp transient, small-signal AC, and noise hooks when a compiled model exposes them. Simple model cards bind parseable parameters, while large foundry model cards are kept on the compiled model defaults except safe selector parameters such as `type` until full `.PARAM` expression evaluation is implemented.
- Model-fidelity reporting is now explicit: GSPICE prints `MODEL_STATUS: OSDI_LOADED` and `MODEL_STATUS: OSDI_DEVICE` lines, and refuses compact/PDK-looking MOS models from silently falling back to primitive Level-1 behavior unless an explicit debug environment override is set.
- Numerical-policy controls now include accuracy presets plus `NUMERICAL`/`CONVERGENCE`/`POLICY` presets such as `ROBUST`, which enable source stepping, GMIN stepping, line search, adaptive transient control, and higher Newton iteration limits.
- Production-deck parser compatibility now accepts `.TEMP`, `.GLOBAL`, `.NODESET`, and functional `.IC V(node)=value` seeds. `.IC` seeds the initial solution and is honored directly for transient runs with `UIC`.
- The CTest suite now includes syntax smoke tests plus an ASCII RAW numeric regression check for an RC step response.

## Current Critical Gaps

- `.SUBCKT`, `.INCLUDE`, `.LIB`, `.IF`, and `.PARAM` support is stronger but still not Spectre-class. Numeric parameter expressions, simple conditionals, nested includes, library sections, and subcircuit instance overrides work; advanced language constructs, complex quoting/tokenization, discipline-aware constructs, and many model-library corner cases remain incomplete.
- PDK-grade MOS/BJT models must come through compiled compact models. The primitive MOS remains an improved but still simplified Level-1 style model, not BSIM/PSP/HICUM/EKV class modeling. The IHP SG13G2 primitive fallback is explicit opt-in compatibility only, not a production simulation path.
- OSDI support now binds instance parameters, safe model selectors, and simple model-card parameters, and participates in DC, transient, AC, and noise when the compiled model exposes the needed hooks. It still needs full foundry `.PARAM` expression evaluation before all PSP/BSIM/HICUM model-card parameters can be safely bound, broader OpenVAF/ngspice ABI compatibility validation, charge-conservation stress decks, and validated golden examples.
- SuiteSparse/KLU is integrated for real and complex solves, but GSPICE still needs production-grade solver selection, scaling, condition diagnostics, iterative/preconditioned options, and larger stress tests.
- Transient integration has backward-Euler, trapezoidal/Gear2 options, adaptive step control, and source breakpoints, but still needs charge-conserving compact-model Jacobians, stronger trap-ringing controls, better event scheduling, and broad large-circuit validation.
- Nonlinear solve has configurable tolerances plus first-pass source stepping, gmin stepping, and line-search damping, but still lacks advanced homotopy/continuation, per-device convergence diagnostics, and robust device limiting.
- AC analysis only works where devices provide correct small-signal stamps. MOS has first-pass `gm/gds`; OSDI, controlled sources, passives, and basic sources are covered, but full compact-model small-signal behavior still needs validation.
- RF/PSS/HB/PAC/STB paths are mostly prototypes or simplified demonstrations, not production analysis engines.
- Noise has first-pass output-referred resistor/diode/primitive-MOS support, but flicker noise, correlation, compact-model noise hooks, input-referred noise, integrated noise, and PNoise remain incomplete.
- Transfer-function, sensitivity, pole/zero estimates, Monte Carlo, corners, production specs, behavioral sources, and transient measurements have first-pass implementations. Distortion, exact generalized-eigenvalue pole/zero extraction, analysis-wrapped `.STEP`, behavioral `.MEASURE` expressions, production process distributions, and reference-grade validation still need significant work.
- Numeric regression has started with RAW vector checks, but GSPICE still needs broad golden-vector comparisons against ngspice/Xyce/Spectre-style references.

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
