# GSPICE Beta Testing Guide

Date: 2026-07-20

GSPICE is Lumen's native experimental SPICE-like simulator. This beta focuses on correctness, diagnostics, reproducibility, and comparison against Ngspice/Xyce. It is not yet a signoff replacement for Spectre, Ngspice, Xyce, or VACASK.

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
- `tests/decks/osdi_psp_options.sp` when OSDI files are available

## Validation

Focused beta check:

```powershell
ctest --test-dir C:\EDA\GSPICE\build -C Release -R "(smoke_osdi_options_deck|smoke_osdi_psp_inverter|smoke_osdi_limiting_rhs_option|smoke_osdi_model_status|regression_save_voltage_selection|smoke_tran_predictor_options|smoke_behavioral|smoke_controlled|smoke_measure|smoke_pz|smoke_tf|smoke_noise|smoke_step|smoke_mc)" --output-on-failure
```

Full local check:

```powershell
ctest --test-dir C:\EDA\GSPICE\build -C Release --output-on-failure
```

## Compact Model Policy

GSPICE must not silently downgrade PDK-grade active devices to primitive MOS behavior. If a deck needs PSP/BSIM/HICUM-class modeling, it should use compiled compact models through OSDI/OpenVAF-style libraries, or fail with a diagnostic explaining what is missing.

Advanced OSDI controls are intentionally explicit:

- `OSDI_LIMITING_RHS`
- `OSDI_SPICE_RHS`
- `OSDI_TRAN_JACOBIAN`
- `OSDI_INTERNAL_NODES`
- `OSDI_BIND_FULL_MODEL_PARAMS`

Some controls are still experimental and should be validated against Ngspice/Xyce before trusting results.

## Known Limitations

- PSP/OSDI support is under active validation.
- Internal-node expansion is experimental.
- Full foundry `.PARAM` expression binding is incomplete.
- RF/PSS/HB/STB and exact pole-zero extraction are not production grade.
- PSS, PAC, PNoise, SP execution, and periodic derivative analyses are
  explicitly unsupported in 1.3; HB and STB remain experimental.
- Primitive MOS is only a debug/simple-device model, not a production compact model.

## Better Than VACASK Direction

GSPICE aims to combine VACASK-inspired Verilog-A/OSDI rigor with:

- strict unsupported-feature diagnostics,
- reference comparison against Ngspice/Xyce,
- reproducible run manifests through Lumen,
- model-fidelity reporting per active instance,
- and a growing regression suite for every advertised capability.
