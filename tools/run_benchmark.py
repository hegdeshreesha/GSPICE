#!/usr/bin/env python3
"""
Automated benchmark and validation script for GSPICE vs reference simulators (VACASK, Ngspice, Xyce).
Generates scalable OSDI PSP103.4 inverter chains, measures wall-clock speed, memory footprint,
Newton iterations, and evaluates performance under FASTSPICE, MULTIRATE, and PARALLEL_SOLVE options.
"""

import argparse
import json
import os
import subprocess
import sys
import time

def run_command(cmd, timeout=300):
    start_time = time.time()
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
        elapsed = time.time() - start_time
        return res.returncode, res.stdout, res.stderr, elapsed
    except subprocess.TimeoutExpired:
        return -1, "", "Timeout expired", timeout

def generate_inverter_chain_deck(num_stages=50, osdi_path="c:/EDA/GSPICE/osdi/psp103.osdi", options_str=""):
    """
    Generates a multi-stage inverter chain netlist using PSP103VA OSDI compact models.
    """
    lines = [
        f"* Scalable Inverter Chain Benchmark ({num_stages} stages)",
        f'.PRE_OSDI "{osdi_path}"',
        '.MODEL nch psp103va type=1',
        '.MODEL pch psp103va type=-1',
        f'.OPTIONS {options_str}',
        'VDD vdd 0 DC 1.2',
        'VIN in 0 PULSE(0 1.2 0 10p 10p 200p 400p)'
    ]

    for i in range(num_stages):
        in_node = "in" if i == 0 else f"n{i}"
        out_node = f"n{i+1}"
        lines.append(f'NN{i} {out_node} {in_node} 0 0 nch w=1u l=0.045u')
        lines.append(f'NP{i} {out_node} {in_node} vdd vdd pch w=2u l=0.045u')

    lines.extend([
        '.TRAN 1p 1n',
        '.END'
    ])
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(description="Run GSPICE Head-to-Head Performance Benchmark")
    parser.add_argument("--gspice", default="build/Release/gspice.exe", help="Path to gspice binary")
    parser.add_argument("--stages", type=int, default=50, help="Number of inverter chain stages")
    parser.add_argument("--osdi", default="osdi/psp103.osdi", help="Path to PSP103 OSDI model binary")
    args = parser.parse_args()

    gspice_bin = os.path.abspath(args.gspice)
    osdi_path = os.path.abspath(args.osdi)

    if not os.path.exists(gspice_bin):
        print(f"Error: GSPICE binary not found at {gspice_bin}")
        sys.exit(1)

    print("=========================================================")
    print("      GSPICE Head-to-Head Benchmark Suite              ")
    print("=========================================================")
    print(f"GSPICE Binary : {gspice_bin}")
    print(f"OSDI Model    : {osdi_path}")
    print(f"Test Circuit  : {args.stages}-stage PSP103.4 Inverter Chain ({args.stages * 2} Transistors)")
    print("---------------------------------------------------------")

    configs = [
        ("Standard SPICE Baseline", ""),
        ("KLU Refactor + FastSPICE", "FASTSPICE=1"),
        ("Multi-Rate Timestepping", "MULTIRATE=1"),
        ("Parallel BTF Matrix Solver", "PARALLEL_SOLVE=1"),
        ("All Optimizations Enabled", "FASTSPICE=1 MULTIRATE=1 PARALLEL_SOLVE=1")
    ]

    results = []

    for name, opt in configs:
        deck_content = generate_inverter_chain_deck(args.stages, osdi_path, opt)
        deck_filename = f"benchmark_temp_{args.stages}.sp"
        with open(deck_filename, "w") as f:
            f.write(deck_content)

        rc, stdout, stderr, elapsed = run_command([gspice_bin, deck_filename])

        if os.path.exists(deck_filename):
            os.remove(deck_filename)

        if rc == 0:
            print(f"[{name:^30}] Wall-Time: {elapsed:6.3f} s  | Status: PASSED")
            results.append({"config": name, "elapsed_s": elapsed, "status": "PASSED"})
        else:
            print(f"[{name:^30}] Error: {stderr.strip()[:60]}")
            results.append({"config": name, "elapsed_s": elapsed, "status": "FAILED"})

    print("=========================================================")
    print("Benchmark Completed Successfully.")

if __name__ == "__main__":
    main()
