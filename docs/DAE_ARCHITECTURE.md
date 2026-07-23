# GSPICE DAE and Charge Architecture

## Purpose

GSPICE uses a device-neutral modified-nodal differential-algebraic
equation (DAE) contract:

```text
F(x, t) + dQ(x, t)/dt = 0
```

`F` contains instantaneous flow contributions. `Q` contains conserved storage
quantities such as terminal charge and branch flux. Devices provide both
residuals and their Jacobians. Each analysis applies its own mathematical
operator to the same device data:

| Analysis | Residual or matrix |
|---|---|
| Operating point | `F`, `dF/dx` |
| Transient | `F + dQ/dt`, `dF/dx + alpha*dQ/dx` |
| AC | `dF/dx + j*omega*dQ/dx` |
| Harmonic balance | `F` and `Q` at time collocation points |

This prevents transient, AC, and future harmonic-balance implementations from
drifting into different physical charge models.

## Clean-room boundary

VACASK, Xyce, ngspice, published numerical-analysis literature, compact-model
standards, and public simulator behavior may be used to identify requirements
and compare numerical results. GSPICE implementation work follows these rules:

1. Do not copy source code, pseudocode, internal names, comments, tests, or file
   organization from another simulator.
2. Derive algorithms independently from mathematical definitions, standards,
   and published numerical methods.
3. Write GSPICE-specific interfaces and tests before migrating complex models.
4. Use other simulators as black-box numerical references. Keep reference decks
   and tolerance decisions attributable and reproducible.
5. Review contributions for license compatibility before merging them into the
   Apache-2.0 codebase.

Architectural inspiration is allowed; source-level imitation is not.

## Current implementation

`include/dae.hpp` defines the independent DAE request, result, and assembly
types. Resistor, capacitor, inductor, diode, BJT, primitive MOS, and OSDI
devices expose this contract. DC, transient, and AC prefer DAE evaluation and
retain a compatibility fallback for devices that have not migrated. All three
analyses therefore consume the same `F`, `Q`, `dF/dx`, and `dQ/dx` data for
migrated devices.

`include/transient_state_store.hpp` provides contiguous numeric and opaque
transactional state stores. The transient engine owns accepted `Q` and `Qdot`
history, per-device opaque model state, candidate frames, bounded speculative
commits, and cursor-based rollback. Trial steps do not allocate polymorphic
device snapshots.

`include/integration_formula.hpp` derives variable-step BDF and implicit
Adams-Moulton coefficients from time history. `METHOD=BDF` and `METHOD=ADAMS`
support orders 1 through 5, with automatic startup order reduction. The
simulator restarts multistep solution and `Q/Qdot` history at discontinuity
breakpoints. `include/transient_control.hpp` provides nonuniform polynomial
prediction, an error multiplier derived from the selected predictor/corrector
pair, and projected `p-1`/`p`/`p+1` order selection. Step doubling remains a
periodic independent LTE oracle and is also used during startup and around
breakpoints. `METHOD=AUTO` starts with backward Euler, uses trapezoidal
integration on smooth native-device trajectories, temporarily switches to
Gear2/BDF2 after alternating weakly damped slopes, and returns to trapezoidal
after a stable interval. Compact models that request damping remain on Gear2.

OSDI evaluation uses epoch-scoped bypass caching. The cache key includes the
analysis request, time, integration coefficient, evaluation epoch, and device
inputs; it is invalidated by accepted-state changes, restore, and internal-node
binding. Residual verification always performs a full high-precision
evaluation. OSDI limiting flags prevent false Newton convergence. Optional
hidden internal-unknown expansion binds only uncollapsed topology roots into
MNA without exposing them as user-visible voltage nodes.

`include/dae_audit.hpp` independently finite-differences `F` and `Q`, compares
them with the supplied Jacobians, and checks per-group terminal charge and
Jacobian conservation. It is available in tests and at runtime with
`.OPTIONS DAE_AUDIT=1` for audit-safe devices.

## Invariants

Every migrated charge model must demonstrate:

- terminal charge conservation;
- an analytic or finite-difference match between `Q` and `dQ/dx`;
- consistency between AC susceptance and the transient charge derivative;
- no accepted-state mutation during a rejected candidate step;
- convergence-order tests under changing timestep ratios;
- comparison against at least two independent reference implementations when a
  compatible compact model is available.

## Completed migration sequence

1. Migrate linear passives and establish analytical DAE tests.
2. Bind the transactional state store to devices without changing legacy state.
3. Migrate OSDI evaluation into the common DAE result while retaining a guarded
   compatibility path.
4. Replace device-owned charge histories with simulator-owned `Q` history.
5. Add general variable-step BDF and Adams integration formulas.
6. Migrate diode junction charge, then charge-conserving BJT and MOS models.
7. Make AC consume the same stored-quantity Jacobian.

Remaining work includes making harmonic balance consume the contract, broader
convergence-order qualification on stiff compact-model circuits, and a
restart/checkpoint representation for the complete transient state.

## Going beyond current reference simulators

The intended differentiators are correctness features, not a larger checkbox
list:

- direct terminal-charge LTE as a second guard alongside solution LTE;
- runtime charge-conservation and derivative diagnostics for audit-safe models;
- deterministic transactional rollback across model state and simulator state;
- derivative auditing by finite differences in a developer mode;
- cross-analysis consistency tests (`dQ/dx` versus AC and transient behavior);
- bypass invalidation that includes timestep, integration coefficients,
  limiting state, model parameters, and discontinuities;
- reproducible reference validation with stored error budgets.
