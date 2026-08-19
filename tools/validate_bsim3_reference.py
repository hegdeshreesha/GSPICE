"""Compare native BSIM3 probe currents with an independent ngspice deck."""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys
import tempfile


def parse_ngspice(path: pathlib.Path) -> list[float]:
    values = []
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"^\d+\s+([-+0-9.eE]+)\s+[-+0-9.eE]+\s+([-+0-9.eE]+)", line.strip())
        if match:
            values.append(float(match.group(2)))
    return values


def parse_probe(output: str) -> list[float]:
    return [float(line.split()[1]) for line in output.splitlines() if len(line.split()) == 2]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ngspice", required=True)
    parser.add_argument("--gspice", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--relative-tolerance", type=float, default=0.05)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-12)
    args = parser.parse_args()

    source = pathlib.Path(args.source)
    deck = source / "tests" / "decks" / "bsim3_reference_dc.sp"
    with tempfile.TemporaryDirectory(prefix="gspice_bsim3_ngspice_") as tmp:
        log = pathlib.Path(tmp) / "bsim3_ngspice_validation.log"
        subprocess.run([args.ngspice, "-b", "-o", str(log), str(deck)], check=True)
        reference = parse_ngspice(log)
    probe = subprocess.run([args.gspice], check=True, capture_output=True, text=True)
    native = parse_probe(probe.stdout)
    if len(reference) != len(native):
        print(f"BSIM3 reference length mismatch: ngspice={len(reference)} native={len(native)}")
        return 1
    errors = []
    for index, (expected, actual) in enumerate(zip(reference, native)):
        # ngspice reports I(VDS), opposite to the positive drain-terminal
        # current convention used by the native evaluator.
        absolute = abs(abs(actual) - abs(expected))
        scale = max(abs(expected), args.absolute_tolerance)
        relative = absolute / scale
        if absolute > args.absolute_tolerance and relative > args.relative_tolerance:
            errors.append((index, expected, actual, relative))
    for index, expected, actual, relative in errors:
        print(f"point={index} ngspice={expected:.6e} native={actual:.6e} relative_error={relative:.3e}")
    print(f"BSIM3 reference validation: points={len(reference)} failures={len(errors)}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
