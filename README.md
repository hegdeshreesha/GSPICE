# GSPICE

**Lumen's Experimental Native Circuit Simulator Core**

GSPICE is an experimental C++ SPICE-like simulator being developed as the native simulation engine for Lumen Circuit Studio. The current codebase is useful for primitive-device bring-up and Lumen integration work, but it is not yet a replacement for production simulators such as Cadence Spectre, ngspice, or Xyce.

## Current Features

- Primitive MNA stamping for R, C, L, V, I, diode, simple MOS, ports, probes, and simple S-parameter devices.
- `.OP`, `.TRAN`, and `.AC` paths for supported primitive circuits.
- OpenMP-assisted device stamping.
- Named-node stdout output for Lumen waveform import.
- Explicit diagnostics for unsupported subcircuit instances and unsupported model/include directives.

## Major Gaps

- No real `.SUBCKT` expansion yet.
- No `.LIB` / `.INCLUDE` foundry model parsing yet.
- No production BSIM/PSP/HICUM-class compact model support yet.
- The current KLU bridge still falls back to a dense solver.
- RF/PSS/HB/noise/stability analyses are prototype-level and need validation.

See `docs/GSPICE_AUDIT_ROADMAP.md` for the current audit and roadmap.

## Installation And Build

GSPICE uses CMake.

### Prerequisites

- CMake 3.10+
- C++17 compiler, such as MSVC, GCC, or Clang
- OpenMP

### Build

```powershell
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

If MSBuild fails on Windows with a duplicate `Path` / `PATH` environment error, run CMake from a clean process environment.

## Usage

```powershell
gspice my_circuit.sp --threads 1
gspice --help
gspice --version
```

## License

GSPICE is licensed under the Apache License, Version 2.0.
