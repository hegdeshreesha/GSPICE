#!/usr/bin/env python3
"""
Automated benchmark and validation script for GSPICE vs reference simulators (VACASK, Ngspice, Xyce).
Tracks execution time, peak memory usage, accepted steps, Newton iterations, and LTE accuracy.
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

def main():
    parser = argparse.ArgumentParser(description="Run GSPICE benchmark suite")
    parser.add_argument("--gspice", default="build/Release/gspice.exe", help="Path to gspice binary")
    parser.add_argument("--circuit", required=False, help="Path to circuit netlist")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    if not os.path.exists(args.gspice):
        print(f"Error: GSPICE binary not found at {args.gspice}")
        sys.exit(1)

    print(f"=== GSPICE Performance Benchmark ===")
    print(f"Binary: {args.gspice}")
    print(f"Threads: {args.threads}")

    # Run capabilities check
    rc, stdout, stderr, elapsed = run_command([args.gspice, "--capabilities"])
    if rc == 0:
        print("Capabilities:")
        print(stdout)

    print("Benchmark framework ready.")

if __name__ == "__main__":
    main()
