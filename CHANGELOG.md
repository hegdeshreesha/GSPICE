# Changelog

## 1.3.0-academic-beta - 2026-07-21

### Added

- CSV transient output selected by `--format csv` or a `.csv` output suffix.
- Gaussian and bounded-uniform source Monte Carlo sampling with deterministic
  seeds and optional one-dimensional Latin-hypercube sampling.
- Machine-readable `--capabilities` maturity report and explicit `--self-test`.
- Numeric/version/capability, CSV, Monte Carlo, false-success, DAE audit, and
  high-order integration regression tests, bringing the audit-host suite to 58
  tests.
- A device-neutral `F/Q` DAE contract shared by DC, transient, and AC for
  resistor, capacitor, inductor, diode, BJT, primitive MOS, and OSDI devices.
- Simulator-owned transactional `Q/Qdot` history and fixed-size opaque compact-
  model state frames with bounded rollback and no per-trial snapshot allocation.
- General variable-step BDF and implicit Adams-Moulton coefficients through
  order 5 with automatic startup order reduction.
- Runtime finite-difference derivative and charge-conservation auditing with
  `.OPTIONS DAE_AUDIT=1`.
- VACASK/Xyce comparison, authoritative limitations, release checklist,
  citation metadata, contribution/security policies, and Apache-2.0 license.
- Windows and Ubuntu CI builds that run the complete regression suite before
  uploading executables.

### Correctness and safety

- Unimplemented analyses are rejected instead of falling through to an
  operating-point result.
- PSS no longer reports convergence from a stub.
- STB without a stability probe and nonconverged experimental HB now fail.
- Unknown CLI options, invalid formats, missing option arguments, and multiple
  input files now produce usage errors.
- The binary version is generated from the CMake project version instead of a
  stale hard-coded string.
- Adaptive two-half-step acceptance now records the accepted midpoint in both
  solution and DAE histories, preventing time/state history misalignment.
- Multistep solution and Q/Qdot histories restart at source breakpoints, and
  `METHOD=AUTO` uses backward-Euler startup followed by trapezoidal integration
  for substantially lower error on smooth waveforms while retaining a damped
  Gear2/BDF2 branch for compact models.
- Diode, BJT, and primitive MOS transient storage is charge-conserving and is
  assembled from the same `Q` Jacobian used by small-signal AC.
- Corrected OSDI Jacobian row/column orientation, contribution-presence flag
  interpretation, and residual/Jacobian sign convention across DC, transient,
  and AC; a state-preserving PSP audit now verifies the ABI mapping.

### Compatibility notes

- Invoking GSPICE with no netlist now prints help; use `--self-test` to run the
  built-in OSDI smoke test.
- PSS/PAC/PNoise/SP and periodic derivative directives are hard errors until a
  validated execution engine exists.
