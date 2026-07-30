# GSPICE vs VACASK bottleneck audit

Audit date: 2026-07-30

Upstream VACASK project: <https://codeberg.org/arpadbuermen/VACASK>

Local VACASK checkout: `C:\EDA\VACASK`

Local VACASK revision: `764f69317d840cbe99d0cf93164f0e7b5ee38515`

## Executive summary

GSPICE is not primarily blocked by KLU on the small reference decks. The dominant PSP inverter bottleneck is transient step rejection and the repeated Newton work caused by LTE control around compact-model switching events. The second bottleneck is matrix/device assembly overhead: GSPICE still stamps into sparse triplets and collapses/builds CCS around solver calls, while VACASK binds devices to persistent matrix value slots and refactors a stable KLU matrix in place.

VACASK is ahead architecturally because it treats circuit elaboration, matrix sparsity, sweep alteration, and nonlinear solve as persistent objects. GSPICE has recently added KLU symbolic/numeric reuse, but it still rebuilds too much surrounding state per iteration and per sweep/job.

## Reference runs on this host

Commands used:

- `C:\EDA\GSPICE\build-vcpkg\gspice.exe`
- `C:\EDA\Tools\vacask_0.3.4.rc1\vacask_0.3.4.rc1_windows-x86_64\bin\vacask.exe`

### RC transient baseline

GSPICE quiet deck, RAW output suppressed with `--save none`:

- Wall time: `0.072535 s`
- Accepted/rejected: `61 / 21`
- Newton solves/iterations: `138 / 211`
- Solver calls: `237`
- KLU time: `0.001012 s`

VACASK `--no-output`:

- Wall time: `0.077677 s`
- VACASK-reported elapsed simulation time: `0.0010756 s`

Conclusion: the tiny RC deck is dominated by executable startup, parsing, and general overhead. It is not a meaningful KLU benchmark.

### PSP103.4 inverter baseline

GSPICE reference deck, no `--adaptive-maxstep`:

- Wall time: `0.348148 s`
- Accepted/rejected: `388 / 371`
- LTE rejected: `369`
- Newton solves/iterations: `1015 / 2163`
- Stamp/solve in Newton summary: `0.06 s / 0.06 s`
- KLU calls/time: `2815 / 0.021270 s`
- Solver total: `0.043221 s`

VACASK `--no-output`:

- Wall time: `0.136670 s`
- VACASK-reported elapsed simulation time: `0.0109083 s`

GSPICE PSP with `LTE_AUDIT_INTERVAL=0`:

- Wall time: `0.278785 s`
- Accepted/rejected: `390 / 384`
- LTE rejected: `382`
- Newton solves/iterations: `956 / 2020`
- Stamp/solve in Newton summary: `0.04 s / 0.05 s`
- KLU calls/time: `2956 / 0.020346 s`

Conclusion: disabling LTE audits helps, but rejection remains high. The bottleneck is the LTE controller and repeated nonlinear work, not symbolic KLU.

## Architectural differences that matter

### Matrix assembly

VACASK:

- Builds a `SparsityMap`, enumerates it, performs KLU analysis, and then device instances use bound matrix value pointers.
- `KluMatrixCore::refactor()` reuses the same KLU symbolic/numeric object with `klu_refactor()`.
- The matrix object exposes `valuePtr()`/`valueArray()` so load code can write directly to pre-enumerated storage.

GSPICE:

- Uses `SparseMatrix::add()` triplets plus thread triplet buffers.
- Every solve calls `matrix.getEntries()`, builds CCS, optionally column-scales, and then calls KLU.
- KLU symbolic/numeric reuse is present, but CCS construction and entry collapse still happen around it.

Expected bottleneck on larger decks: GSPICE will spend more time in stamping/collapse/fill than VACASK because device stamps are not direct writes into stable matrix arrays.

### Transient LTE/rejection policy

VACASK:

- Uses predictor-corrector LTE control and Gear/trapezoidal/BE methods.
- Benchmarks show very low rejection on several analog baselines and explicit discussion of residual-check vs speed tradeoffs.

GSPICE:

- PSP reference run shows `371` rejected steps for `388` accepted steps.
- Most rejections are LTE rejections, not convergence rejections.
- `--adaptive-maxstep` made PSP worse by allowing a larger ceiling: accepted/rejected changed to `294 / 562`.

Expected bottleneck: step-size proposal grows too aggressively or LTE estimation is too conservative/noisy around PSP transitions, causing many rejected steps and repeated device evaluation.

### Nonlinear/device bypass

VACASK:

- Has continuation bypass and inactive element bypass controls documented in the benchmark methodology.
- Its C6288 notes show continuation bypass improves time/iteration from about `20.4 ms` to `15.6 ms` with residual checks disabled.

GSPICE:

- Reports `bypassed_devices=0` on the PSP reference run.
- OSDI/device bypass exists in code paths, but it is not firing on this baseline.

Expected bottleneck on digital/mostly-inactive circuits: GSPICE evaluates too many unchanged devices every Newton iteration.

### Sweeps and alteration

VACASK:

- Can sweep parameters/options/instance/model values and alter without reloading/re-elaborating the whole circuit.

GSPICE:

- SimENV currently expands variable sweeps into separate GSPICE jobs/decks because GSPICE does not support `.STEP PARAM`.
- This is compatible and correct, but not fast for deep sweeps.

Expected bottleneck on ADE workflows: per-point process/deck startup and full parse/elaboration dominate small and medium sweep points.

## Recommended GSPICE work plan

1. Add a first-class profiling mode:
   - Report parse/elaboration, matrix clear/fill/collapse, CCS build, device evaluation, KLU refactor, KLU solve, output serialization, LTE estimation, and rejected-step replay time.

2. Replace triplet stamping hot paths:
   - Build a persistent sparsity map after elaboration.
   - Bind each device contribution to direct matrix/RHS slots.
   - Keep triplets only for topology-changing paths and debug.

3. Tune transient control on compact-model edges:
   - Add rejection histograms by reason and time window.
   - Clamp growth after recent LTE rejection.
   - Make LTE audit optional by analysis preset and avoid duplicate expensive solve paths unless requested.
   - Revisit absolute LTE scaling for small capacitive nodes and PSP charge variables.

4. Make OSDI continuation bypass visible and useful:
   - Report eligible/bypassed/device-evaluated counts.
   - Enable safe first-iteration bypass after accepted transient points for time-invariant devices.
   - Add regression decks where bypass should trigger.

5. Add native variable/model/instance sweep execution:
   - Parse `.STEP PARAM`.
   - Reuse parsed netlist and elaboration where topology is unchanged.
   - Reuse operating points between nearby sweep points, like VACASK continuation mode.

6. Create comparable benchmarks:
   - Port VACASK RC, graetz, ring, PSP inverter, and a smaller multiplier to valid GSPICE SPICE decks.
   - Run with output disabled, KLU required, single thread, same stop/maxstep/reltol/vntol/abstol policy.

## Current ranking of likely bottlenecks

1. PSP transient LTE rejection and repeated Newton replay.
2. Matrix assembly/collapse/CCS build around KLU.
3. OSDI/device evaluation without effective continuation bypass.
4. Full process/deck restart for variable sweeps.
5. Output/progress overhead on small jobs.
6. KLU numeric solve itself, only after the above are fixed or on larger matrices.
