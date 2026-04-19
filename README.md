# GSPICE
**The Next-Generation Parallel Open-Source Circuit Simulator Core**

GSPICE is a high-performance, object-oriented circuit simulator built from the ground up to be a "Super Form" successor to Ngspice and Xyce. It combines modern C++ architecture with industrial-grade stability and advanced RF analysis capabilities.

## 🚀 Key Features
*   **Parallel Core:** Fully multi-threaded device evaluation and matrix stamping using OpenMP. Scalable up to 16 threads for massive designs.
*   **Universal Models (OSDI):** Native **Open Simulation Device Interface** host. Dynamically load professional Verilog-A models (BSIM4, PSP, HICUM) compiled with OpenVAF.
*   **Industrial Stability:** Verified against **Level-50 Industrial MOSFET models**. Robust convergence using G-min stepping and Newton-Raphson damping.
*   **Scalable Math:** High-efficiency **Sparse Matrix** architecture with a bridge for **KLU (SuiteSparse)**, the industry gold-standard for circuit simulation.
*   **Professional CLI:** Full-featured command-line interface with thread control, help systems, and standard SPICE output formatting.

## 🛰️ Advanced Analysis Suite
GSPICE supports 15+ industrial analysis types:
*   **Standard:** `.OP`, `.TRAN`, `.AC`, `.NOISE`
*   **RF Core:** `.PSS` (Shooting Method), `.HB` (Harmonic Balance), `.SP` (S-Parameters)
*   **RF Advanced:** `.PAC`, `.PNOISE`, `.HBAC`, `.HBNOISE`, `.HBSP`
*   **Stability:** `.STB`, `.HBSTB`, `.PSSSTB` (Using the Tian Loop-Break Algorithm)

## 🛠️ Installation & Build
GSPICE uses CMake for a professional cross-platform build experience.

### Prerequisites
*   CMake 3.10+
*   C++17 Compiler (MSVC 2022, GCC, or Clang)
*   OpenMP 2.0+

### Build Instructions
```powershell
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

## 💻 Usage
```powershell
# Run a simulation on 8 threads
gspice my_circuit.sp --threads 8

# Specify a custom output file
gspice lna_test.sp -o lna_results.raw

# View help and version
gspice --help
gspice --version
```

## 📄 License
GSPICE is licensed under the **Apache License, Version 2.0**.
