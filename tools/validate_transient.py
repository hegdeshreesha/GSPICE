#!/usr/bin/env python3
"""Reproducible transient comparison against analytic, Xyce, and VACASK oracles."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def run(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    completed = subprocess.run(
        command, cwd=cwd, env=env, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=180, check=False)
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n{completed.stdout}")


def normalized(name: str) -> str:
    return name.strip().strip('"{}').lower().replace(" ", "")


def read_csv_signal(path: Path, signal: str) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        keys = [normalized(item) for item in header]
        time_index = next((i for i, key in enumerate(keys) if key in {"time", "index"}), None)
        wanted = normalized(signal)
        signal_index = next((i for i, key in enumerate(keys) if key == wanted), None)
        if time_index is None or signal_index is None:
            raise RuntimeError(f"missing time/{signal} in {path}: {header}")
        samples = []
        for row in reader:
            if len(row) <= max(time_index, signal_index):
                continue
            try:
                samples.append((float(row[time_index]), float(row[signal_index])))
            except ValueError:
                continue
    if len(samples) < 3:
        raise RuntimeError(f"too few samples in {path}")
    return samples


def read_ascii_raw(path: Path, signal: str) -> list[tuple[float, float]]:
    variables: dict[str, int] = {}
    values = False
    samples: list[tuple[float, float]] = []
    time_index = signal_index = variable_count = None
    pending: list[float] = []
    for original in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = original.strip()
        if not values and line.lower() == "values:":
            values = True
            time_index = variables.get("time")
            signal_index = variables.get(normalized(signal))
            variable_count = len(variables)
            continue
        if not values:
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                variables[normalized(parts[1])] = int(parts[0])
            continue
        parts = line.split()
        try:
            if parts and parts[0].isdigit():
                pending = [float(item) for item in parts[1:]]
            else:
                pending.extend(float(item) for item in parts)
        except ValueError:
            continue
        if (variable_count is not None and len(pending) >= variable_count and
                time_index is not None and signal_index is not None):
            samples.append((pending[time_index], pending[signal_index]))
            pending = []
    if len(samples) < 3:
        raise RuntimeError(f"too few samples or missing signal {signal} in {path}")
    return samples


def interpolate(samples: list[tuple[float, float]], target: float) -> float:
    if target <= samples[0][0]:
        return samples[0][1]
    if target >= samples[-1][0]:
        return samples[-1][1]
    lo, hi = 0, len(samples) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if samples[mid][0] <= target:
            lo = mid
        else:
            hi = mid
    t0, y0 = samples[lo]
    t1, y1 = samples[hi]
    alpha = (target - t0) / max(t1 - t0, 1e-300)
    return y0 + alpha * (y1 - y0)


def metrics(candidate: list[tuple[float, float]], reference: list[tuple[float, float]]) -> tuple[float, float, int]:
    errors = []
    for time, value in candidate:
        errors.append(value - interpolate(reference, time))
    return max(map(abs, errors)), math.sqrt(sum(e * e for e in errors) / len(errors)), len(errors)


def analytic_samples(candidate: list[tuple[float, float]]) -> list[tuple[float, float]]:
    delay = 1e-3
    rise = 1e-6
    tau = 1e-3
    ramp_end = (rise - tau + tau * math.exp(-rise / tau)) / rise
    samples = []
    for time, _ in candidate:
        elapsed = time - delay
        if elapsed <= 0.0:
            value = 0.0
        elif elapsed < rise:
            value = (elapsed - tau + tau * math.exp(-elapsed / tau)) / rise
        else:
            value = 1.0 + (ramp_end - 1.0) * math.exp(-(elapsed - rise) / tau)
        samples.append((time, value))
    return samples


def locate_output(directory: Path, preferred: str, suffix: str) -> Path:
    direct = directory / preferred
    if direct.exists():
        return direct
    matches = sorted(directory.glob(f"*{suffix}"), key=lambda path: path.stat().st_mtime, reverse=True)
    if not matches:
        raise RuntimeError(f"no {suffix} output found in {directory}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gspice", required=True, type=Path)
    parser.add_argument("--xyce", type=Path)
    parser.add_argument("--vacask", type=Path)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--max-abs", type=float, default=2e-3)
    parser.add_argument("--psp-max-abs", type=float, default=2e-2)
    args = parser.parse_args()
    args.gspice = args.gspice.resolve()
    if args.xyce:
        args.xyce = args.xyce.resolve()
    if args.vacask:
        args.vacask = args.vacask.resolve()
    reference_dir = args.source.resolve() / "tests" / "reference"

    with tempfile.TemporaryDirectory(prefix="gspice_transient_reference_") as name:
        work = Path(name)
        gspice_csv = work / "gspice.csv"
        run([str(args.gspice), "--threads", "1", "--format", "csv", "--output",
             str(gspice_csv), str(reference_dir / "rc_step_gspice.sp")], work)
        gspice = read_csv_signal(gspice_csv, "v(out)")

        max_abs, rms, count = metrics(gspice, analytic_samples(gspice))
        print(f"analytic RC: max_abs={max_abs:.6e} rms={rms:.6e} samples={count}")
        failures = int(max_abs > args.max_abs)

        if args.xyce:
            shutil.copy2(reference_dir / "rc_step_xyce.cir", work / "rc_step_xyce.cir")
            run([str(args.xyce), "rc_step_xyce.cir"], work)
            xyce = read_csv_signal(locate_output(work, "rc_step_xyce.csv", ".csv"), "v(out)")
            max_abs, rms, count = metrics(gspice, xyce)
            print(f"Xyce RC:     max_abs={max_abs:.6e} rms={rms:.6e} samples={count}")
            failures += int(max_abs > args.max_abs)
        else:
            print("Xyce RC:     SKIPPED (no executable supplied)")

        if args.vacask:
            shutil.copy2(reference_dir / "rc_step_vacask.sim", work / "rc_step_vacask.sim")
            run([str(args.vacask), "rc_step_vacask.sim"], work, os.environ.copy())
            vacask = read_ascii_raw(locate_output(work, "tran_ref.raw", ".raw"), "out")
            max_abs, rms, count = metrics(gspice, vacask)
            print(f"VACASK RC:   max_abs={max_abs:.6e} rms={rms:.6e} samples={count}")
            failures += int(max_abs > args.max_abs)

            shutil.copy2(reference_dir / "psp_inverter_vacask.sim", work / "psp_inverter_vacask.sim")
            run([str(args.vacask), "psp_inverter_vacask.sim"], work, os.environ.copy())
            vacask_psp = read_ascii_raw(locate_output(work, "psp_ref.raw", ".raw"), "out")
            gspice_psp_csv = work / "gspice_psp.csv"
            gspice_env = os.environ.copy()
            gspice_env.setdefault("GSPICE_OSDI_DIR", str(args.source.resolve() / "osdi"))
            run([str(args.gspice), "--threads", "1", "--format", "csv", "--output",
                 str(gspice_psp_csv), str(reference_dir / "psp_inverter_gspice.sp")], work, gspice_env)
            gspice_psp = read_csv_signal(gspice_psp_csv, "v(out)")
            max_abs, rms, count = metrics(gspice_psp, vacask_psp)
            print(f"VACASK PSP:  max_abs={max_abs:.6e} rms={rms:.6e} samples={count}")
            failures += int(max_abs > args.psp_max_abs)
        else:
            print("VACASK RC:   SKIPPED (no executable supplied)")
            print("VACASK PSP:  SKIPPED (no executable supplied)")

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
