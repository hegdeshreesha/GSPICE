#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def run_case(exe, deck):
    with tempfile.TemporaryDirectory() as tmp:
        csv_path = Path(tmp) / "out.csv"
        proc = subprocess.run(
            [exe, "--threads", "1", "--format", "csv", "--output", str(csv_path), deck],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            print(proc.stdout)
            raise SystemExit(proc.returncode)
        with csv_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            raise AssertionError("CSV has no transient rows")
        return proc.stdout, rows[-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--deck", required=True)
    parser.add_argument("--signal", required=True)
    parser.add_argument("--measure", required=True)
    parser.add_argument("--min", type=float, required=True)
    args = parser.parse_args()

    output, last = run_case(args.exe, args.deck)
    if args.signal not in last:
        raise AssertionError(f"CSV missing {args.signal}; columns={','.join(last)}")
    csv_value = float(last[args.signal])
    match = re.search(rf"MEASURE {re.escape(args.measure)} = (\S+)", output)
    if not match:
        raise AssertionError(f"stdout missing MEASURE {args.measure}")
    measure_value = float(match.group(1))
    if csv_value <= args.min:
        raise AssertionError(f"{args.signal}={csv_value:.12g} <= {args.min:.12g}")
    if measure_value <= args.min:
        raise AssertionError(f"{args.measure}={measure_value:.12g} <= {args.min:.12g}")
    print(f"ok {args.signal}={csv_value:.12g} {args.measure}={measure_value:.12g}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"current probe validation failed: {exc}", file=sys.stderr)
        raise
