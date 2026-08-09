# GSPICE Production Readiness

This checklist is for the native GSPICE + Lumen IHP flow.

## Required Gate

Run:

```powershell
python C:\EDA\GSPICE\tools\production_signoff.py `
  --source C:\EDA\GSPICE `
  --build C:\EDA\GSPICE\build-klu `
  --gspice C:\EDA\GSPICE\build-klu\Release\gspice.exe
```

Optional Lumen deck gate:

```powershell
python C:\EDA\GSPICE\tools\production_signoff.py `
  --lumen-deck C:\Users\hegde\Downloads\scratch\simenv_Dummy1_single\20260808_113041\input.sp `
  --lumen-raw C:\Users\hegde\Downloads\scratch\simenv_Dummy1_single\20260808_113041\production_signoff.raw
```

The gate requires:

- `SuiteSparse-KLU` capability.
- `psp103_native`, `ihp_psp_rf_full`, `ihp_passive_wrapper_full`.
- `bsim3_gsdi`, `bsim4_gsdi`.
- `juncap2_tat`, `juncap2_avalanche`.
- Broad native GSDI/GMC regression pass.
- Optional Lumen deck completion with zero warnings and RAW output.

The optional Lumen deck gate also runs a verbose compatibility scan. Normal
production mode remains quiet, but verbose warnings are reported as `risks`
because they may indicate places where the native model is faster than an
OSDI/OpenVAF flow by evaluating a smaller model surface.

The same deck is also scanned by `tools/rank_psp103_ignored.py`. The generated
`build-klu/psp103_ignored_ranking.json` groups both ignored PSP103 parameters
and accepted-but-approximate PSP103 parameters into critical/high/medium/low
buckets such as noise, self-heating, junction leakage, RF parasitics,
geometry/binning, and charge/capacitance. Treat critical/high items as
production parity blockers for designs that exercise those effects.

For strict model-fidelity signoff, add:

```powershell
  --strict-fidelity
```

This fails the signoff if verbose compact-model compatibility warnings are
present, if the PSP103 ignored-parameter ranking contains critical/high
occurrences, or if the deck uses accepted PSP103 parameters that the native
evaluator still handles approximately.

## VACASK Comparison Bridge

VACASK is an optional external oracle. To check whether a Lumen/GSPICE deck can
be converted into VACASK syntax without touching the factory PDK, run:

```powershell
python C:\EDA\GSPICE\tools\vacask_lumen_bridge.py `
  --deck C:\Users\hegde\Downloads\scratch\simenv_Dummy1_single\20260808_113041\input.sp `
  --workdir C:\EDA\GSPICE\build-klu\vacask_bridge
```

Add `--run-vacask` after conversion succeeds. A rejected conversion exits with
the standard skip code `77` unless `--require-vacask` is set, and writes a JSON
report with the converter stdout/stderr tails. This keeps VACASK validation
actionable without making normal GSPICE runs depend on VACASK.

## Lumen Production Defaults

- GSPICE timeout is bounded by default at 900 seconds.
- Simulation Cockpit exposes a timeout control; `Auto` uses the default.
- KLU preference is persisted and requests `SOLVER=KLU`.
- `MAXSTEP=AUTO` is conservative and caps the internal transient step to the
  `.TRAN` print step. Use an explicit larger `MAXSTEP` only when throughput is
  more important than edge/phase fidelity.
- Normal runs keep compatibility diagnostics off; verbose diagnostics are opt-in.
- Run manifests record command, artifacts, effective timeout, save policy, and GSPICE capability summary.

## Remaining Signoff Gap

The remaining gaps are independent full-reference BSIM parity, full official
PSP parameter-surface parity, and arbitrary OSDI runtime model loading. The
current native PSP/BSIM GSDI paths are tested and normal-mode warning-clean, but
verbose diagnostics and the PSP scanner can still expose official model-card
parameters that the native evaluator either does not consume or only
approximates. Treat very large speedups over OSDI/OpenVAF as valid only after
the verbose scan, the accepted-approximation scan, and an external reference
sweep all pass for the target design.
