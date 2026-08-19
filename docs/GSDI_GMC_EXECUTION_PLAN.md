# GSPICE GSDI + GMC Compact-Model Execution Plan

Consolidated engineering plan for replacing the external OSDI/OpenVAF-style
compact-model stack with in-house components:

- **GSDI** — GSPICE Standard Device Interface (runtime model ABI) replaces OSDI.
- **GMC** — the build-time BSIM/PSP-required Verilog-A-subset compiler
  (`tools/gmc/*.py`) replaces OpenVAF.

Target models, in order: **BSIM4.8.3, BSIM3, PSP103** (plus the existing
JUNCAP/JUNCAP2 native modules). External simulators and ngspice remain benchmark
oracles only.

---

## 1. Context and Goal

The ecosystem is being mapped against external references (PSP103.4 + OSDI 0.4 +
OpenVAF-reloaded; C6288 benchmark: 10 112 transistors, 57.98 s / 138.5 MB vs
ngspice 71.81 s, Xyce 151.57 s) to define competitive targets for GSPICE.

### Long-term target runtime flow (simplified)

```
parser (.MODEL, M-elements)
  -> GmcRegistry factories
  -> GsdiModel / GsdiInstance          <- GMC-generated headers emit these directly
  -> GdiDevice (collapse-aware stamps)  <- GsdiDaeDeviceAdapter only bridges
                                           hand-written native Device models
  -> Device DAE  F(x) + dQ(x)/dt
  -> OP / DC / TRAN / AC / NOISE / PSS / PNOISE
```

`GsdiDaeDeviceAdapter` remains useful for bridging native `Device` models into
GSDI, but generated GMC models should eventually emit `GsdiModel` /
`GsdiInstance` directly and not stop at `GdiInstanceBase`.

Every model exposes a local node space (terminals + internal + collapsible
nodes) mapped onto the global MNA matrix via `GsdiCollapseMap` (pair collapse
`V(x,y) <+ 0` and ground collapse), so hidden nodes (e.g. series resistance)
either merge into a kept terminal column or become genuine extra unknowns.

---

## 2. Architectural Principles (non-negotiable)

1. **GSDI replaces OSDI, GMC replaces OpenVAF.** No OSDI or OpenVAF code,
   headers, struct layouts, comments, tests, or implementation organization
   are copied into this repository. Concepts such as model descriptors,
   setup/evaluate split, F/Q decomposition, and Jacobian patterns are standard
   simulator architecture and inform the design without being "concept
   property" — nothing is copied in any concrete form. The external reference simulator is a
   benchmark reference only.
2. **Build-time generated C++ headers, no `dlopen`.** No external compact-model
   compilers and no run-time loading of model binaries. Models are compiled to
   native C++ at build time (CMake custom command, e.g.
   `gmc_va_probe.va -> gmc_va_probe.hpp -> executable`).
3. **Fail-closed compact-model enablement.** Unsupported or unverified models
   are rejected with a clear error; no silent fallback to primitive devices.
4. **One shared `F(x) + dQ(x)/dt` DAE path.** All compact models, native or
   generated, flow through the same device-neutral DAE contract; no
   per-analysis companion-model hacks.
5. **External simulators / ngspice as benchmark oracles only.** ngspice (BSD) is the numeric
   oracle for reference gates; other simulators are benchmark references only.
6. **Target order: BSIM4.8.3, BSIM3, PSP103.** BSIM4.8.3 first — it is the
   broadest overlap with the ecosystem benchmarks; then BSIM3 (shared
   architecture, faster to back-port); PSP103 last via the same pipeline.
7. **License boundary on BSIM sources.** BSIM4.8.3
   (`C:\EDA\BSIM4_4.8.3_standard_05192025`) is Educational Community License
   v2.0 (ECL-2.0, educational / non-commercial). `b4par.c` etc. are used as a
   **numeric oracle only**; the GSPICE implementation is a clean-room
   reimplementation and never copies code from it into the Apache-2.0 tree.
8. **Test gate per step.** Every implementation slice ships with specific
   in-order focused tests; the full suite runs green before any merge/release
   (see §5).

---

## 3. Execution Flow

| Layer | Component | Responsibility |
|---|---|---|
| Parser | `src/core/parser.cpp` | `.MODEL` cards, M-element elaboration, factory registration |
| Factories | `GmcRegistry` (`include/gmc.hpp`) | creates models/devices from model definitions |
| Runtime ABI | `include/gsdi.hpp`, `include/gdi.hpp` | descriptors, collapse maps, eval contract, dual numbers |
| Generated bridge | GMC-emitted `GsdiModel` / `GsdiInstance` | direct target for generated models |
| Native bridge | `include/gsdi_device_adapter.hpp` | `GsdiDaeDeviceAdapter` for hand-written native `Device` models only |
| Global adapter | `include/devices/gdi_device.hpp` | local index -> global MNA stamping, collapse-aware |
| Compiler | `tools/gmc/*.py` | BSIM/PSP-required Verilog-A subset -> C++ headers (parser, AST, emitter) |

---

## 4. Phases

### Phase A — GSDI v1.0 (runtime model interface) — largely complete

- `GsdiModelDescriptor`: nodes (roles Terminal / Internal / Collapsible),
  parameters (model vs instance, min/max), opvars, static `jacobian_pattern`,
  analysis capability flags.
- `GsdiCollapseMap`: merge-chain resolution (`valueSource`), `standard()` plan
  from descriptor, fail-closed on cycles / out-of-range partners.
- Evaluation contract `GsdiEvalRequest` / `GsdiEvalResult` (static and dynamic
  residual + Jacobian, conservation groups, noise sources, opvars, limiting /
  bypass flags).
- Dual-number arithmetic `GmcDual<N>`; `GmcDual4` alias preserved (generated
  code checks `std::is_same_v<Scalar, GmcDual4>`), `GmcDual8` for BSIM-class
  models; `gmcTanh`, `gmcAbs`, `gmcMin`, `gmcMax`, `gmcLim`.
- `GdiDevice` collapse-aware stamping (static / dynamic / AC / noise);
  collapsed values read the kept column, ground-merged rows/cols map to the
  reference index (-1) and are emitted so charge-conservation groups close
  (stampers no-op negative indices; matches the Capacitor/Diode/BJT/BSIM
  emission contract).
- Gates: `gmc_dual`, `gsdi_metadata`, `gsdi_collapse` — passing.

### Phase B — GMC v1.0 (build-time compiler)

- BSIM/PSP-required Verilog-A subset: parser, AST, C++ emitter
  (`tools/gmc/`). Full Verilog-A is out of scope — too broad and unnecessary.
- Emit the model descriptor **and** evaluation bodies using `GmcDual<N>` sized
  to ports + internal nodes.
- **Emit GSDI-native descriptor metadata for generated `.va` models, and route
  the generated model through `GsdiModel` / `GsdiInstance` directly** (instead
  of stopping at `GdiInstanceBase`). This is the best next implementation step.
- Generated header consumed via CMake custom command (pattern:
  `gmc_compile.py --input <model>.va --output <build>/gmc_va_probe.hpp`).
- **Deferred: no optimizer.** Do not build code optimization yet. First make
  the BSIM/PSP equations correct; optimize after profiling (Phase F).
- Gates: `gmc_veriloga_parser`, `gmc_generated_va`, `gmc_hidden_nodes`,
  `gmc_hidden_internal_elaboration`; generated models must pass the collapse and
  dual finite-difference gates.

### Phase C — BSIM4.8.3 production enablement

- Clean-room implementation of BSIM4.8.3 under GMC/GSDI, guided by `b4par.c`
  as numeric oracle (ECL-2.0 boundary: oracle only, no code copying).
- Parser: enable production registration of BSIM4 compact models; internal-node
  MNA allocation for non-collapsible hidden nodes.
- Validation: reference gates vs ngspice across parameters / DC / transient /
  noise / charge conservation; retire the experimental factories.
- Gates: full `bsim4_*` suite, reference probes, smoke gates. **This phase is
  the current big blocker.**

### Phase D — BSIM3 production enablement

- Back-port the BSIM4 pipeline to BSIM3 (shared architecture) once BSIM4.8.3
  gates are green.
- Validation: ngspice reference gates across the `bsim3_*` suites.
- Gates: `bsim3_*` suites + reference probes + smoke gates.

### Phase E — PSP103 production

- Recast `psp103_core` under the GSDI descriptor + collapse plan, using the
  same pipeline as C/D: the native evaluator now routes through
  `GsdiDaeDeviceAdapter` + collapse-aware `GdiDevice` (local indices 0-3,
  `GsdiCollapseMap::standard` from `psp103GsdiDescriptor()` in
  `include/devices/psp103_gsdi.hpp`). PSP103 is a pure 4-terminal model
  (D/G/S/B); the descriptor declares terminals only, so `standard()` composes
  the identity plan (no hidden `Td`/`Tint` nodes in the native core — the
  collapsible/internal-node machinery remains exercised by the GMC-generated
  `psp_like`/`gmc_hidden` gates).
- Gate: `psp103_gsdi_descriptor` — descriptor metadata (4 terminals, dense
  4x4 jacobian pattern, op/tran/ac support), identity collapse plan, and
  numeric identity of the routed vs direct evaluator path (static + dynamic
  residuals/Jacobians).
- Validation: ngspice reference gates first; benchmark parity only after
  correctness gates pass.
- Reference gates (both green, `temp=21`, `psp103_*-2.mod` CMC ref data):
  - `psp103_idvg_e2e_parity` — NMOS IdVg vs `psp_ref_idvg.raw`, 31 pts,
    max rel err 0.097%.
  - `psp103_pmos_idvg_e2e_parity` — PMOS IdVg vs `psp_ref_pmosidvg.raw`
    (ngspice `nmp1` instance convention; gspice uses `M`-prefix device),
    31 pts, max rel err 0.097%.

### Phase F — Performance (after correctness)

- Benchmark against external simulator / ngspice / Xyce numbers **only after the BSIM/PSP
  correctness gates are green**. C6288 (10 112 transistors) parity lives here,
  not in earlier phases.
- Optimize after profiling: KLU backend and sparse stamp optimization, DAE
  AUDIT gates, fail-closed robustness checks (bypass, limiting, unsupported
  models). Any optimizer work for generated code belongs here, driven by
  profiler output.

---

## 5. Test-Gate Policy

- **Focused gates per slice:** each slice runs its targeted suites in-order
  (e.g. `gmc_dual`, `gsdi_metadata`, `gsdi_collapse`; then `bsim4_*`; ...).
- **Full suite before merge/release** — "full suite green every slice" is the
  ideal, but on this Windows sandbox MSBuild sometimes needs elevated
  FileTracker access; a full-suite run may be skipped mid-slice when the
  sandbox blocks it, and is mandatory only at merge/release boundaries.

Current test inventory (measured):

- `ctest -C Debug` in `build/`: **136/136 pass** (100%), including all GMC/GSDI
  gates (`gmc_dual`, `gmc_generated_va`, `gmc_hidden_nodes`,
  `gmc_hidden_internal_elaboration`, `gmc_deep_va`, `gmc_veriloga_parser`,
  `gsdi_metadata`, `gsdi_collapse`, `gmc_juncap_express`, ...), the BSIM dual
  Jacobian/parameter gates (`bsim3_dual_jacobian`, `bsim4_parameters`), and the
  PSP103 IdVg parity gates.
- `ctest -C Release` in `build/`: **136/136 pass**.

Windows build note: the build must run inside the VS x64 developer
environment (`VsDevCmd.bat -arch=x64`); a bare `cmake --build` fails with
MSVC C1083 (missing `INCLUDE` environment).

---

## 6. Current Status (factual, against the current tree)

**Done:**

- GSDI metadata/collapse infrastructure exists: `GsdiModelDescriptor`,
  `GsdiCollapseMap` (chains, `standard()`), `GmcDual<N>` (with `GmcDual4`
  alias preserved for generated code and `GmcDual8` for BSIM-class models),
  plus the `is_gmc_dual_v` dispatch trait used by generated code.
- Collapse-aware `GdiDevice` stamping in place; gates `gmc_dual`,
  `gsdi_metadata`, `gsdi_collapse` pass.
- GMC pipeline works end-to-end for a tiny model:
  `.va -> gmc_compile.py -> .hpp -> compiled executable`
  (`gmc_va_probe` pattern; `gmc_generated_va` gate passes).
- **GMC now emits GSDI-native code**: the generated header contains a
  `GsdiModelDescriptor` (model type, terminal nodes, parameters, dense
  `jacobian_pattern`, `supports_*` derived from `ddt` presence) plus
  `GsdiModel` / `GsdiInstance` classes evaluated with `GmcDual<N>`. The
  generated instance is exercised both directly and through the `GdiDevice`
  adapter in the `gmc_generated_va` gate; no generated code uses
  `GdiInstanceBase` / `GdiModuleBase` anymore.
- VA-subset hidden-node syntax implemented (`terminal <n>;` and
  `collapsible <node> with <terminal>;` in `tools/gmc/gmc_parser.py`), carried
  through the AST (`ModuleNode.terminal_count`, `collapsible_pairs`) and the IR
  dump, and emitted as `GsdiNodeRole::Terminal / Internal / Collapsible`
  descriptor entries with `collapse_partner`, `terminalCount()` /
  `internalNodeCount()`, and `GmcDual<N>` sized to the full local node count.
  Fail-closed validation rejects collapsible partners that are not external
  terminals and terminal counts that exceed the port list. Gate `gmc_hidden_nodes`
  (model `psp_like` in `tests/gmc_hidden.va` with internal node `Td` and
  collapsible node `Tint -> s`) verifies descriptor metadata, the standard
  collapse plan, local and globally-stamped residuals/Jacobians, and the
  chain-rule cancellation of the collapsed pair's dynamic stamps; `gmc_veriloga_parser`
  covers the new syntax plus the failure cases.
- **VA-subset expression depth for the deep-intrinsic slice is in place** (Phase
  B expansion): the parser supports the ternary `?:` operator and the `**`
  power operator, `$`-prefixed system identifiers (`$temperature`, `$vt`) and
  `$`-prefixed intrinsics (e.g. `$ln`); the emitter resolves `sin` / `cos` /
  `tan` / `sinh` / `cosh` / `log10` / `min` / `max` / `pow` (plus the earlier
  `exp` / `ln` / `log` / `sqrt` / `abs` / `tanh` / `limexp`) through a
  `MATH_INTRINSICS` table into dual-dispatch lambdas that switch between
  `GmcDual<N>` AD routines and plain `std::` doubles. Unknown functions fail
  closed with `ValueError`; `ddt` outside a top-level branch contribution is
  rejected. The emitted `GsdiInstance` exposes `temperature_` / `vt_` (with the
  Boltzmann factor) so `$temperature` / `$vt` resolve to instance state. Gate
  `gmc_deep_va` (model `deep_probe` in `tests/gmc_deep.va`) cross-checks the full
  intrinsic/ternary surface via direct evaluation and the `GdiDevice` stamp path
  against a plain-double mirror model using central finite differences.
- **User-defined `analog function`s are supported** (the single biggest missing
  construct for the BSIM3/BSIM4 sources, which define ~15+ helper functions
  each, e.g. `exp_lim`, `fetlim`, `pnjlim`, `limvgs`): the parser handles
  `analog function real <name>;` headers, `input`/`integer`/`real`
  declarations, `begin/end` or bare bodies, and assignment-to-function-name
  returns; the emitter emits each function as a generic capture-by-reference
  lambda (callee-first order via topological sort, `Scalar` deduced from the
  first argument) so calls resolve identically for plain doubles and
  `GmcDual<N>` AD. Functions are pure: only arguments, locals, intrinsics, and
  other user functions may be referenced — anything else fails closed
  (`recursive analog function`, `out-of-scope identifier`, non-`real` return).
  The `gmc_deep_va` model now calls a two-function chain (`gclamp` ->
  `gsoft`) and the gate re-evaluates at a clamp-bound corner.
- **`case` statements and `integer` variables are supported** (BSIM4 uses 13
  `case` selectors — `geo`, `rgeo`, `dioMod`, `rgateMod`, `tnoiMod`,
  `fnoiMod`, ... — and 273+ module-scope `integer` declarations): the parser
  handles `case (expr) ... endcase` with comma-separated value lists
  (`0, 10:`) and a trailing `default:` item, and `integer` variable
  declarations at module scope; the emitter lowers each `case` to a chained
  `if / else if / else` ladder comparing the selector (`Scalar`/`GmcDual`)
  against each item value, with `ddt` / user-function walkers covering nested
  case bodies. The `gmc_deep_va` model gains an `integer sel` selector whose
  `case` item folds `1e-9*v` into the mirror model's expected current, so the
  FD cross-check exercises the ladder through the `GdiDevice` stamp path.
- **A Verilog-A text preprocessor is implemented** (`tools/gmc/gmc_preprocessor.py`)
  and runs before the GMC lexer/parser. It resolves `` `define `` / `` `undef `` /
  `` `ifdef `` / `` `ifndef `` / `` `else `` / `` `endif `` / `` `include `` /
  `` `resetall `` with lazy rescanning macro expansion: object-like macros,
  function-like macros with parenthesis-depth + string-aware argument
  splitting (VBIC-style string args containing commas and empty strings),
  trailing-backslash line continuations, and multi-line macro invocations
  (BSIM statement style). `` `include `` resolves relative to the including
  file then against `--include` dirs, with cycle/nesting guards; the parser
  skips `nature`/`discipline`/`connectrules` sections from `disciplines.vams`.
  `gmc_compile.py` reports include dependencies.
- **The parser now consumes the remaining real-model constructs**: `&&`/`||`/
  `!` (plus `&`/`|`/`%` tokenized, parsed, and failed closed at emission),
  quoted strings as single tokens, `$strobe`/`$display`/`$discontinuity`/
  `$bound_step`/`$fatal`/... as no-op statements, `for`/`while` loops with
  `begin`/`end` bodies (emitted as C++ `for`/`while` with GmcDual-compatible
  `gmc_nz` conditions), `@(initial_step)` event blocks (collapsed into the
  per-evaluation pass with a comment), `white_noise`/`flicker_noise`
  contributions (parsed, fail closed at emission), return-type-less
  `analog function` declarations (default `real`), and non-literal parameter
  defaults (`parameter real dsub = drout;` — parsed into `default_expr`, fail
  closed at emission). Smoke check (manual, AGPL sources stay out of repo):
  preprocessed `bsim3v3.va` (711 parameters, 1218 variables, 17 functions,
  1349 analog statements), `bsim4v8.va` (897/2414/1) and `vbic_1p3.va` now
  **parse fully to IR**; preprocess still leaves zero residual backticks.
- **The parser and emitter now survive the full BSIM3/BSIM4 analog bodies**:
  bare `begin`/`end` grouping blocks (BSIM's `if (...) begin begin ... end end`
  style) are consumed without disturbing control-flow nesting, so the whole
  analog body reaches the IR instead of truncating at a stray `end`;
  `electrical` node declarations pull names that are not module ports into the
  port list as `Internal` nodes (bsim3: `di`,`si`; bsim4:
  `di`,`si`,`gi`,`gm`,`bi`,`sbulk`,`dbulk`), so `V(x,y)` probes and `I(x,y) <+`
  contributions stamp internal-node rows instead of `vn[-1]` UB. Emitter
  additions: `$limit(x, ...)` folds to its first argument (all other args are
  domain bounds); the charge-partition idiom `C * ddt(q)` with a
  time-invariant multiplier hoists to `ddt(C * q)`; `white_noise`/
  `flicker_noise` emit through the GSDI noise-source vector and mark generated
  descriptors as noise-capable. Repro smoke (manual, AGPL): `bsim3v3.va` and
  `bsim4v8.va` now **emit the full analog body to C++** and the generated
  `GmcBsim3Model`/`GmcBsim4Model` **compile clean** on MSVC 2026 and evaluate
  through `GsdiInstance::evaluate` with correct terminal/internal residual and
  Jacobian sizes (bsim3: 6×7 matrix over 4 terminals + 2 internal; bsim4:
  11 nodes). Registering these factories with the elaborator and the
  ngspice-clean cross-check are Phase C/D.
- **Elaborator MNA allocation for non-collapsible internal nodes is in place**:
  `GmcRegistry::registerModel(type, factory, internal_node_count)` +
  `internalNodeCount(type)`; `GmcModelDefinition.internal_nodes` (trailing field,
  keeps legacy aggregate initializers valid); `Netlist::createInternalNode`
  allocates `%internal.<device>.<k>` columns from the normal node-ID pool. The
  parser's M-element path allocates terminal nodes first, then hidden columns,
  and the `PSP_LIKE` factory (behind `GSPICE_HAVE_GMC_GENERATED`, embedding the
  generated headers into the `gspice` binary) builds a `GdiDevice` with
  `GsdiCollapseMap::standard`. Gate `gmc_hidden_internal_elaboration` runs the
  `tests/decks/gmc_hidden_internal.sp` deck (`NPSP1 1 2 0 0 psp1`) and checks the
  internal `Td` node solves to 4.8 V with VDD=5/VGG=1.
- **`GdiDevice` emits reference-index rows for ground-merged terminals** so
  charge-conservation groups close (previously dropped, which unbalanced the
  `smoke_dae_audit` group pool for any device with a grounded terminal); the
  audit now reports `charge=1.56e-14 charge_J=0` on the DAE audit deck.

**Open / blockers:**

Current update: BSIM3 level-49 and BSIM4 level-54 model cards now route
through native GSDI/GMC without experimental environment flags and have
no-warning smoke gates. Full independent reference parity across broad
DC/AC/tran/noise grids remains Phase C/D work before they can be called fully
validated. JUNCAP2 TAT and reverse-breakdown branches now have native smoothed
implementations with derivative and parser/GSDI smoke coverage.

- **BSIM reference gates are the big blocker** — BSIM3/BSIM4 production
  enablement is still disabled and the ngspice reference gates fail until the
  clean-room models land (Phase C/D).
- Parser/elaborator: internal-node MNA allocation exists only for models
  registered with an internal-node count (the `PSP_LIKE` path); BSIM4/BSIM3
  factories with 6/11 local nodes (4 terminals + hidden `di`/`si`/...) are not
  yet registered into the elaborator — Phase C/D work.
- Generated bsim4 evaluates with NaN Jacobian entries on the internal-node
  columns at the default bias probe (degenerate AD terms, not an
  emit/compile defect) — numeric cross-check vs ngspice is Phase D.
- PSP103 descriptor + collapse migration: native evaluator recast through
  `GsdiDaeDeviceAdapter` + `GdiDevice` with `psp103GsdiDescriptor()` (Phase E);
  NMOS/PMOS IdVg ngspice reference parity gates stay green (0.097% max rel
  err) and gate `psp103_gsdi_descriptor` passes. PSP103's core is a pure
  4-terminal model, so the descriptor declares terminals only.
- **Debug suite is fully green (136/136)**: the last four Debug-only WIP gates
  were closed — `gsdi_collapse` (test data misindexed `x_vals` against
  `nodes={5,7}`; Release "pass" was vacuous since `assert` compiles out under
  NDEBUG), `gmc_deep_va` (off-corner probe at v=-0.5 made `log10(v*10)` NaN in
  both mirror and generated model; probe moved to v=0.1 inside the `gsoft`
  lower-clamp region), `bsim4_parameters` (geometry-binning expectations
  ignored the `/Leff` division; `hot`/`hotJunction` cases omitted `AT`, so the
  default `AT=3.3e4` at 100 C clamped VSAT to 0 and failed validation; plus a
  real model fix — `junctionTemperatureScale` was applied to `IS/JS/JSW`
  before the card values were read and binned, so it was discarded, and it
  used `XTI/EG` before those were read; the scale now applies after the binned
  reads), and `bsim3_dual_jacobian` (FD/AD cross-check ran exactly on the
  `vsb = vs - vb = 0` kink where `bsim3Positive` is non-differentiable; bias
  moved off the kink into the smooth region).

---

## 7. Timeline (hidden by default)

<details>
<summary>Timeline (hidden)</summary>

Phase sequencing and schedule are intentionally withheld from this document.

Ordering is strict: **A -> B -> C (BSIM4.8.3) -> D (BSIM3) -> E (PSP103) -> F
(performance)**, with focused test gates between every slice and a full-suite
gate at merge/release boundaries. Correctness precedes performance; C6288
benchmark work is Phase F only.

</details>
