# GSPICE 1.3 academic-beta limitations

This is the authoritative known-limitations list for the academic beta. A
successful run means only that the requested deck stayed within an implemented
path and its numerical solve completed. It is not a signoff guarantee.

## Scientific-use boundary

- Not qualified for tapeout, safety-critical design, foundry signoff, device
  characterization, or publication-quality claims without independent
  reference validation.
- Validation is concentrated on small decks. Reproducible RC analytic transient
  comparison now exists, but coverage
  across operating regions, temperature, tolerances, stiff systems,
  pathological topologies, and large matrices is sparse.
- No formal accuracy specification, backward-compatibility policy, or stable
  numerical reproducibility guarantee exists across compilers/platforms.
- Experimental analyses and models must be cross-checked against independent
  references or measured data.

## Netlist and language

- Implements a SPICE subset, not a complete production simulator language.
- Subcircuits, includes, library sections, conditionals, continuation lines,
  and parameter expressions cover common cases but not all quoting, escaping,
  scoping, recursion, function, vector, string, table, or corner-library rules.
- No general user-defined functions, discipline/nature system, digital or
  mixed-signal language, XSPICE code models, or arbitrary runtime alteration.
- `.SAVE` output supports node voltages and device currents for branch-current
  devices plus DAE/static semiconductor current probes. General expressions and
  all dynamic element current histories are incomplete.
- Multiple analysis blocks and analysis-wrapped nested sweeps are not a general
  control language.
- Error source locations can be imprecise after nested preprocessing.

## Devices and models

- Primitive MOS is a simplified Level-1-style teaching/debug model. It is not
  BSIM, PSP, EKV, HiSIM, or foundry-grade.
- Native BSIM4 level-54 DC has a clean-room BSIM4.8.3 parity gate. Native
  BSIM3 level-49 has an ngspice weak/moderate/strong DC parity gate; broader
  geometry, temperature, capacitance, and transient charge parity still need
  deeper BSIM3 equation coverage.
- Primitive BJT, diode, and noise implementations are first-pass models with
  limited temperature, charge, breakdown, high-injection, flicker, and
  correlation behavior.
- Broad transmission-line families, broad magnetic-core modeling,
  semiconductor switches, IBIS, TCAD devices, and many standard SPICE model
  levels are absent. Primitive `NJF`/`PJF` JFET devices, MESFET aliases routed
  through the primitive JFET approximation, one-section `Txxx` Z0/TD line
  approximations, and linear two-winding `Kxxx L1 L2 k` mutual inductance are
  supported for OP/AC/transient through the DAE path.
- native compact-model compatibility is experimental. Reactive residuals participate in
  transient integration and charge-aware LTE checks. Recognized model
  parameters are resolved and bound strictly; limiting flags participate in
  convergence; epoch-scoped evaluation bypass and hidden uncollapsed internal
  unknowns are implemented. Limiting-RHS use and internal-node expansion remain
  opt-in while broader native compact-model ABI, foundry-card, collapse-topology, charge
  conservation, noise, and operating-point-variable validation is completed.
- IHP/PSP flows depend on compatible external native compact models and environment
  configuration; primitive fallback is intentionally refused by default.
- No bundled, validated cross-platform compact-model package is released.

## Analyses

- OP/DC recovery uses adaptive source continuation, adaptive logarithmic gmin
  descent, line search, temporary `.NODESET` constraints, and worst-unknown
  diagnostics. Device limiting coverage and hard large-circuit homotopy cases
  still need broader qualification.
- TRAN supports BE, trapezoidal, Gear2, and general variable-step BDF and
  implicit Adams-Moulton formulas through order 5. Trial steps are
  transactional, rejected model/Q/Qdot state rolls back, and
  voltage/current/charge LTE uses `TRTOL` and `CHGTOL`. The normal adaptive
  path uses nonuniform polynomial predictor-corrector error estimation with
  `p-1`/`p`/`p+1` order selection; periodic step doubling supplies an
  independent oracle. `METHOD=AUTO` detects alternating weakly damped slopes
  and temporarily switches from trapezoidal to Gear2. These paths have unit,
  smoke and analytic checks, but still lack broad
  convergence-order qualification on stiff nonlinear compact-model circuits,
  dense event scheduling, checkpoint/restart, and large-circuit validation.
- AC is correct only for devices with complete, validated small-signal stamps.
- Noise is output-referred and first-pass. PNOISE/HBNOISE execution exists for
  the implemented periodic small-signal path, but input-referred/integrated
  noise, comprehensive flicker/correlation, and transient noise remain limited
  and need external oracle sweeps before signoff use.
- `.TF` supports a voltage output and independent source input subset.
- `.SENS` is finite-difference source sensitivity, not general direct/adjoint
  device/model sensitivity.
- `.PZ` is a frequency-sweep threshold estimator, not a generalized-eigenvalue
  pole-zero solver.
- `.MEASURE` covers transient voltage probes and device-current probes for
  `FIND`, `MIN`, `MAX`, `PP`, `AVG`, `RMS`, `INTEG`, `DERIV`, and simple
  `WHEN probe=value` crossings with `RISE`/`FALL`/`CROSS` counts. Common
  transient trigger/target delay measurements are supported. General
  expressions, AC/DC modes, and broad SPICE compatibility are absent.
- `.STEP`, `.CORNER`, and `.MC` operate mainly on independent-source OP runs.
  Arbitrary instance/model/option distributions and correlations are absent.
- Monte Carlo supports Gaussian or bounded-uniform variation of one source,
  deterministic seeds, and one-dimensional LHS. It has no correlation matrix,
  process/mismatch hierarchy, Sobol sequence, importance sampling, or confidence
  interval reporting.
- HB has a native collocation Newton refinement path with shooting
  initialization and signoff-required mode, but it is still prototype-grade
  without production continuation/preconditioning coverage.
- PSS, PAC/PXF, PNOISE, SP, HBAC/HBNOISE/HBSP/HBSTB/HBXF, PSSSP, and PSSSTB
  now execute for supported small decks and have smoke/signoff regressions.
  They remain prototype RF engines until each has independent oracle sweeps
  across nonlinear devices, sidebands, ports, noise, and convergence stress
  cases.
- STB is experimental and requires the specific implemented probe topology.
- `.FOUR` supports transient voltage-probe harmonic tables over the final
  period. Distortion, envelope, transient sensitivity, and general uncertainty
  quantification are absent.

## Numerics, scale, and performance

- Development builds may omit SuiteSparse/KLU and use the less externally
  validated internal sparse solver. The vcpkg release preset requires KLU and
  fails configuration if the dependency cannot be found or linked.
- Row scaling and one-pass iterative refinement are enabled by default. There
  is still no condition-number estimate, column equilibration, iterative/
  preconditioned solver family, distributed matrix, or MPI support.
- OpenMP parallelism is limited to selected stamping loops and can be disabled
  for small workloads; thread count does not imply speedup.
- Symbolic reuse and singleton reduction exist, but large, ill-conditioned,
  highly coupled, and memory-stress circuits lack benchmark coverage.
- No restart/checkpoint files, crash recovery, deterministic parallel reduction
  contract, or resource-limit controls.

## Output, APIs, and tooling

- ASCII RAW and CSV transient output are supported; binary RAW, HDF5,
  Touchstone analysis output, rich metadata, compression, and streaming APIs
  are incomplete or absent.
- CSV uses a simple invariant numeric representation and does not yet carry
  units, provenance, tolerances, model hashes, or rejected-step history.
- There is no stable public C/C++ library ABI, Python package, remote API, or
  simulator-coupling contract.
- CLI exit codes distinguish usage, parse, and simulation failures but are not
  yet a formally versioned API.
- macOS is not part of the initial release CI matrix.

## Release and project maturity

- Version 1.3 is an academic beta with a small contributor/test base.
- CI covers Windows and Ubuntu builds/tests; hardware/compiler diversity is
  still narrow.
- The configured suite contains smoke, parser, RAW, and analytic checks, but many remain smoke/regex tests
  rather than independent high-precision numeric oracles.
- Fuzzing, sanitizers, static analysis, supply-chain attestations, signed
  artifacts, SBOMs, reproducible builds, and a vulnerability-response history
  are not established.
- Apache-2.0 covers GSPICE code only. External native compact models, PDKs, SuiteSparse,
  OpenMP runtimes, and redistributed dependencies retain their own licenses.
