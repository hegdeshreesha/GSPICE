# GSPICE
**The Next-Generation Open-Source Circuit Simulator Core**

GSPICE is a high-performance, object-oriented circuit simulator built from the ground up to be a "Super Form" successor to Ngspice and Xyce. It is designed for modern semiconductor processes and RF design.

## Key Features
*   🚀 **High Performance:** Built in modern C++ with a native Sparse Matrix architecture.
*   🔋 **Universal Models:** Native **OSDI (Open Simulation Device Interface)** host, allowing dynamic loading of professional Verilog-A models (BSIM, PSP, HICUM).
*   🛰️ **RF Ready:** Built-in support for advanced RF analyses:
    *   **PSS** (Periodic Steady State via Shooting Method)
    *   **HB** (Harmonic Balance Engine)
    *   **SP** (S-Parameter characterization)
*   🛡️ **Stable Solvers:** Robust non-linear Newton-Raphson with G-min stepping and damping.
*   ⚡ **Scalable:** Integrated bridge for **KLU (SuiteSparse)**, the industry gold-standard for circuit simulation.

## Supported Analyses
*   `.OP`: DC Operating Point
*   `.TRAN`: Transient Analysis (Backward Euler)
*   `.AC`: Small-signal Frequency Domain
*   `.SP`: Scattering Parameters (RF)
*   `.NOISE`: Adjoint-matrix Noise Density
*   `.PSS`: Periodic Steady State
*   `.HB`: Harmonic Balance

## Installation & Build
GSPICE uses CMake for a cross-platform build experience.

### Prerequisites
*   CMake 3.10+
*   C++17 Compiler (MSVC, GCC, or Clang)

### Build Instructions
```powershell
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

## Usage
```powershell
gspice my_circuit.sp [options]
```
Use `gspice --help` for full command line options.

## License
Licensed under the **Apache License, Version 2.0**.
