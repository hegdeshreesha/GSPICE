#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
from pathlib import Path


def parse_diode_dc(log_text):
    rows = []
    for line in log_text.splitlines():
        parts = line.split()
        if len(parts) < 3 or not parts[0].isdigit():
            continue
        try:
            sweep = float(parts[1])
            vout = float(parts[2])
        except ValueError:
            continue
        rows.append({"sweep": sweep, "V(out)": vout})
    if not rows:
        raise RuntimeError("ngspice diode DC output had no numeric rows")
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ngspice", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--deck", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()

    log = Path(args.workdir) / f"{args.case}.log"
    proc = subprocess.run(
        [args.ngspice, "-b", "-o", str(log), args.deck],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        if log.exists():
            print(log.read_text(encoding="utf-8", errors="replace"))
        raise SystemExit(proc.returncode)
    text = log.read_text(encoding="utf-8", errors="replace")
    if args.case != "ngspice_diode_dc":
        raise RuntimeError(f"unsupported ngspice oracle case: {args.case}")
    rows = parse_diode_dc(text)
    with Path(args.csv).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["sweep", "V(out)"])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
