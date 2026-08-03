# Independent Compact Model Architecture

GSPICE must remain a self-contained simulator core. Do not depend on external
model compilers, generated plugin ABIs, or copied simulator implementation code
for compact-device support.

## Principles

- Implement compact models natively inside GSPICE.
- Use permissively licensed specifications, papers, and measurement decks only
  as references; do not copy implementation code.
- Fail closed when a deck requests an unsupported compact model.
- Keep model equations behind the existing DAE `F/Q` device contract so OP, DC,
  transient, AC, noise, PSS, and PNOISE share one path.
- Make performance a first-class design constraint: analytic derivatives,
  stable sparse patterns, cached model state, and structure-of-arrays storage
  for large instance groups.

## Replacement Layers

1. `CompactModelRegistry`
   - Maps `.MODEL` type names to native model factories.
   - Reports unsupported model types with clear diagnostics.

2. `CompactModel`
   - Owns parsed model-card parameters, defaults, temperature transforms, and
     validation.
   - Creates shared immutable model data for instances.

3. `CompactInstance`
   - Owns geometry and instance parameters.
   - Evaluates terminal currents, charges, Jacobians, operating-point variables,
     and noise sources.

4. `CompactDevice`
   - Adapts `CompactInstance` to the existing `Device` interface.
   - Stamps DC, transient DAE `F/Q`, AC small signal, and noise without a plugin
     boundary.

## Parser Policy

- `.MODEL nch PSP103VA ...` is accepted only when `PSP103VA` exists in the
  native registry.
- PDK-grade model names must never fall back to primitive Level-1 MOS unless an
  explicit debug-only environment variable is set.
- External compiled-model directives and N-device plugin instances are
  unsupported in the independent architecture.

## Performance Plan

- Pre-parse numeric and expression parameters once per model card.
- Group native compact-device instances by model type for cache-friendly
  evaluation.
- Reuse sparse matrix structure across Newton iterations and accepted transient
  steps.
- Keep debug finite-difference audits out of release hot paths.
- Add model-owned noise callbacks so normal noise and PNOISE reuse the same
  source definitions.

## Implementation Order

1. Remove external compiled-model headers, loader, emulator, device, scripts,
   docs, and tests.
2. Add the native compact-model registry and unsupported-model diagnostics.
3. Move the existing primitive MOS path behind `CompactDevice`.
4. Add native PSP/BSIM-class skeletons with parameter parsing but fail for
   unimplemented equation blocks.
5. Implement one model family at a time with analytic OP/DC/TRAN/AC/noise
   regressions.
