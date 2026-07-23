# GSPICE

GSPICE is Lumen's experimental C++17 SPICE-like simulator for teaching,
research, reproducible small-circuit experiments, and native integration with
Lumen Circuit Studio. Version 1.3 is an **academic beta**. It is not a
replacement for Xyce, ngspice, VACASK, or a commercial signoff simulator.

## What works

- Primitive R, C, L, V, I, diode, BJT, simplified MOS, controlled-source, and
  behavioral-source devices.
- Operating point, nested DC sweeps, transient, small-signal AC, first-pass
  noise, transfer function, source sensitivity, corners, specifications, and
  source Monte Carlo analyses.
- Backward Euler, trapezoidal, Gear2, and general variable-step BDF/implicit
  Adams-Moulton integration through order 5; adaptive steps, source
  breakpoints, nonuniform polynomial predictor-corrector LTE, adaptive order,
  automatic trap-ringing damping, transactional rejected-step rollback,
  charge-aware LTE control, and adaptive convergence recovery.
- Hierarchical subcircuit expansion, includes, library sections, conditionals,
  parameter expressions, model cards, initial conditions, and save selection
  for the supported syntax subset.
- Optional SuiteSparse/KLU, OpenMP-assisted stamping, and experimental OSDI
  compact-model loading.
- SPICE ASCII RAW and CSV transient output.
- Deterministic Gaussian or uniform source Monte Carlo sampling, including
  optional Latin-hypercube sampling.

Run `gspice --capabilities` for the binary's machine-readable maturity report.
GSPICE refuses unsupported active syntax and analyses instead of returning a
substitute operating-point result.

## Build

Prerequisites are CMake 3.10+, a C++17 compiler, and OpenMP. SuiteSparse/KLU is
optional for development builds and required by the Windows vcpkg release
preset.

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

For a release build with SuiteSparse/KLU, set `VCPKG_ROOT` to a bootstrapped
vcpkg checkout and use the checked-in manifest preset:

```powershell
$env:VCPKG_ROOT = "C:\path\to\vcpkg"
cmake --preset windows-vcpkg-release
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
```

The preset enables `GSPICE_REQUIRE_KLU`, so configuration fails instead of
silently producing an internal-solver release. Verify the resulting executable
with `gspice --capabilities`; `sparse_backend` must be `SuiteSparse-KLU`.

## Use

```powershell
gspice circuit.sp --threads 1
gspice circuit.sp --output waveform.raw
gspice circuit.sp --output waveform.csv --format csv
gspice --version
gspice --capabilities
```

Monte Carlo examples:

```spice
.MC 100 V1 GAUSS(1.0, 0.05) SEED=7
.MC 100 V1 UNIFORM(0.9, 1.1) SEED=7 LHS=1
```

Monte Carlo currently varies one independent source during an OP analysis; it
does not yet vary arbitrary model or instance parameters.

## Transient reference validation

The public harness compares a linear RC transient with an analytic solution,
Xyce, and VACASK, plus a PSP103.4 inverter with VACASK:

```powershell
python tools/validate_transient.py `
  --gspice build/Release/gspice.exe `
  --xyce "C:/Program Files/XyceNF_7.10/bin/Xyce.exe" `
  --vacask "C:/path/to/vacask.exe" `
  --source .
```

`TRTOL` is a multiplier on the LTE tolerance (normally 1 or larger), while
`LTE_RELTOL`, `TRABSTOL`, `CHGTOL`, and `MAXORD` control the underlying
transient error test. Sub-unity `TRTOL` values emit a warning, and repeated LTE
failure at the minimum timestep terminates with an actionable diagnostic rather
than crawling indefinitely. `MAXORD` supports orders 1 through 5 for explicit
`METHOD=BDF` and `METHOD=ADAMS`; `GEAR2` remains the fixed order-2 alias.
`METHOD=AUTO` uses backward-Euler startup and restarts at waveform
breakpoints, then selects trapezoidal integration for smooth native devices,
temporarily changes to Gear2/BDF2 when ringing is detected, or stays damped
when a compact model requests it. `LTE_MODE=PREDICTOR` is the default and uses
a periodic step-doubling audit (`LTE_AUDIT_INTERVAL=16` by default);
`LTE_MODE=STEPDOUBLING` selects the validation oracle for every eligible step.
`TRAN_ORDER_ADAPTIVE`, `TRAP_RINGING`, `NR_BYPASS`, `NODESET_ITERS`, and
`NODESET_G` expose the associated controls. OSDI hidden internal-node expansion
and limiting-RHS application are available as opt-in qualification paths via
`OSDI_INTERNAL_NODES=1` and `OSDI_LIMITING_RHS=1`.
`DAE_AUDIT=1` independently finite-differences migrated device `F` and `Q`
Jacobians after the operating point and checks terminal-charge conservation;
`DAE_AUDIT_TOL` controls its relative tolerance.

## Release status

- [Comparison with VACASK and Xyce](docs/COMPARISON_VACASK_XYCE.md)
- [DAE and charge architecture](docs/DAE_ARCHITECTURE.md)
- [Complete known limitations](docs/LIMITATIONS.md)
- [Academic beta release checklist](docs/ACADEMIC_RELEASE_CHECKLIST.md)
- [Beta testing guide](BETA_TESTING.md)
- [Engineering audit and roadmap](docs/GSPICE_AUDIT_ROADMAP.md)

Please cite GSPICE using `CITATION.cff`. Contributions are welcome under
`CONTRIBUTING.md`.

## License

Apache License 2.0. See `LICENSE`.
