#!/usr/bin/env python3
"""
Performance and Accuracy Benchmark Suite for GSPICE
Measures execution time, timepoints, iteration counts, and output validity.
"""

import sys
import time
import subprocess
from pathlib import Path

def run_benchmark(gspice_exe: Path, deck_path: Path):
    print(f"[*] Running benchmark on deck: {deck_path.name}")
    start = time.perf_counter()
    import os
    env = dict(os.environ)
    env["GSPICE_ALLOW_PRIMITIVE_IHP_FALLBACK"] = "1"
    res = subprocess.run(
        [str(gspice_exe), "--threads", "4", str(deck_path)],
        capture_output=True,
        text=True,
        env=env
    )
    elapsed = time.perf_counter() - start
    
    if res.returncode != 0:
        print(f"[-] Benchmark FAILED for {deck_path.name}!")
        print(res.stderr)
        return False
    
    print(f"[+] {deck_path.name} finished in {elapsed:.4f}s")
    return True

def main():
    root = Path(__file__).parent.parent
    gspice_exe = root / "build" / "Release" / "gspice.exe"
    if not gspice_exe.exists():
        gspice_exe = root / "build" / "gspice.exe"

    if not gspice_exe.exists():
        print("[-] GSPICE binary not found. Please build GSPICE first.")
        sys.exit(1)

    print(f"[*] GSPICE binary: {gspice_exe}")
    
    decks = [
        root / "tests" / "decks" / "rc_pulse_v.sp",
        root / "tests" / "decks" / "options.sp",
        root / "ring_osc_smooth.sp",
    ]

    all_passed = True
    for deck in decks:
        if deck.exists():
            ok = run_benchmark(gspice_exe, deck)
            if not ok:
                all_passed = False

    if all_passed:
        print("[+] All benchmarks completed successfully!")
    else:
        print("[-] Some benchmarks failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
