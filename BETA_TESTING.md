# GSPICE Beta Testing Guide

Date: 2026-07-20

GSPICE is Lumen's native experimental SPICE-like simulator. This beta focuses on correctness, diagnostics, reproducibility, and comparison against independent references. It is not yet a signoff replacement for mature or commercial simulators.

The 1.3 academic beta also adds CSV output, Gaussian/uniform source Monte Carlo
with optional Latin-hypercube sampling, and `--capabilities`. Unimplemented
periodic analyses now fail explicitly instead of returning substitute results.

## Beta Goal

The beta should establish that GSPICE can reliably run supported small decks, write readable RAW/log output, and fail loudly when a deck requires unsupported model syntax or compact-model behavior.

## Build

```powershell
cd C:\EDA\GSPICE
cmake --build build --config Release
```

If `cmake` is not on `PATH`, use:

```powershell
& "C:\Program Files\CMake\bin\cmake.exe" --build C:\EDA\GSPICE\build --config Release
```

If Windows reports duplicate `Path` / `PATH` variables during CMake/MSBuild, launch from a cleaned process environment:

```powershell
if (Test-Path Env:PATH) { $env:Path = $env:PATH; Remove-Item Env:PATH -ErrorAction SilentlyContinue }
cmake --build C:\EDA\GSPICE\build --config Release
```

## Run

```powershell
C:\EDA\GSPICE\build\Release\gspice.exe --threads 1 path\to\input.sp
```

Use `--threads` to select the requested worker count. For first beta validation, use `--threads 1` unless you are specifically testing parallel execution.

## Recommended Smoke Decks

Start with:

- `tests/decks/regression_rc_step.sp`
- `tests/decks/save_voltage_selection.sp`
- `tests/decks/behavioral_voltage.sp`
- `tests/decks/controlled_sources.sp`
- `tests/decks/measure_tran.sp`
- `tests/decks/error_ihp_primitive_fallback.sp` for compact-model fallback checks

## Validation

Focused beta check:

```powershell
ctest --test-dir C:\EDA\GSPICE\build -C Release -R "(error_unsupported_compact_model|error_ihp_primitive_fallback|regression_save_voltage_selection|smoke_tran_predictor_options|smoke_behavioral|smoke_controlled|smoke_measure|smoke_pz|smoke_tf|smoke_noise|smoke_step|smoke_mc)" --output-on-failure
```

Full local check:

```powershell
ctest --test-dir C:\EDA\GSPICE\build -C Release --output-on-failure
```

## Compact Model Policy

GSPICE must not silently downgrade PDK-grade active devices to primitive MOS behavior. External compiled-model plugins are disabled in the Apache build; unsupported compact models must fail with a diagnostic until native compact-model implementations are available.

## Known Limitations

- Native PSP/BSIM/HICUM-class compact models are not production-ready yet.
- Internal-node expansion for native compact models is not implemented yet.
- Full foundry `.PARAM` expression binding is incomplete.
- RF/PSS/HB/STB and exact pole-zero extraction are not production grade.
- PSS, PAC, PNoise, SP execution, and periodic derivative analyses are
  explicitly unsupported in 1.3; HB and STB remain experimental.
- Primitive MOS is only a debug/simple-device model, not a production compact model.

## Reference Simulator Direction

GSPICE aims to combine independently implemented native compact-model rigor with:

- strict unsupported-feature diagnostics,
- reference comparison against public analytic and measurement-backed decks,
- reproducible run manifests through Lumen,
- model-fidelity reporting per active instance,
- and a growing regression suite for every advertised capability.
