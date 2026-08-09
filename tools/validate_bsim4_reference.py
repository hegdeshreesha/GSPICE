"""Compare the native BSIM4 DC probe against the clean-room reference oracle.

The oracle (tools/bsim4_reference_oracle.py) is a faithful transcription of
the official BSIM4.8.3 CMC source (see its module docstring for the exact
source locations).  ngspice is intentionally NOT used as the reference: the
BSIM4.6.x bundled with ngspice ignores PDIBL1/PDIBL2/ETAB even when
version=4.5, so it cannot validate the DIBL pathway.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default="python")
    parser.add_argument("--gspice-probe", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--relative-tolerance", type=float, default=0.01)
    parser.add_argument("--absolute-floor", type=float, default=1.0e-15)
    args = parser.parse_args()
    oracle = pathlib.Path(args.source) / "tools" / "bsim4_reference_oracle.py"
    decks = ["a", "b", "c", "d"]
    worst = 0.0
    failures = 0
    total_points = 0
    for deck in decks:
        reference = subprocess.run(
            [args.python, str(oracle), "--deck", deck, "--json"], check=True,
            capture_output=True, text=True).stdout.split()
        reference = [float(value) for value in reference]
        output = subprocess.run([args.gspice_probe, deck], check=True,
                                capture_output=True, text=True).stdout
        native = [float(line) for line in output.splitlines() if line.strip()]
        total_points += len(native)
        if len(reference) != len(native):
            print(f"BSIM4 deck '{deck}' length mismatch: "
                  f"oracle={len(reference)} native={len(native)}")
            failures += 1
            continue
        for index, (expected, actual) in enumerate(zip(reference, native)):
            delta = abs(actual - expected)
            relative = delta / max(abs(expected), args.absolute_floor)
            worst = max(worst, relative)
            if delta > args.absolute_floor and relative > args.relative_tolerance:
                failures += 1
                print(f"deck={deck} point={index} oracle={expected:.6e} "
                      f"native={actual:.6e} relative_error={relative:.3e}")
    print(f"BSIM4 reference validation: decks={len(decks)} "
          f"points={total_points} failures={failures} worst_relative={worst:.3e}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())