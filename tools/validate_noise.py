#!/usr/bin/env python3
"""Physics sanity checks for GSPICE noise and transient-noise analyses.

These are not golden-waveform regressions. They check properties that should
remain true as the implementation evolves:

* RC thermal noise integrates to roughly kT/C.
* Transient noise is deterministic for a fixed seed and changes with seed.
* OSDI flicker noise has stronger low-frequency density than high-frequency
  density, and transient colored-noise accounting exercises flicker sources.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def run(command: list[str], cwd: Path, timeout: int = 180) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n{completed.stdout}"
        )
    return completed


def write(path: Path, text: str) -> Path:
    path.write_text(text.strip() + "\n", encoding="utf-8")
    return path


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        rows: list[dict[str, float]] = []
        for row in reader:
            rows.append({key.strip().strip('"'): float(value) for key, value in row.items()})
    if not rows:
        raise RuntimeError(f"no CSV rows in {path}")
    return rows


def read_transient_csv(path: Path, signal: str) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        keys = [item.strip().strip('"') for item in header]
        if "time" not in keys or signal not in keys:
            raise RuntimeError(f"missing time/{signal} in {path}: {header}")
        ti = keys.index("time")
        si = keys.index(signal)
        samples = [(float(row[ti]), float(row[si])) for row in reader if len(row) > max(ti, si)]
    if len(samples) < 3:
        raise RuntimeError(f"too few transient samples in {path}")
    return samples


def validate_rc_thermal(gspice: Path, work: Path) -> None:
    deck = write(
        work / "rc_thermal_noise.sp",
        """
        * RC thermal noise should integrate to kT/C.
        V1 in 0 DC 0 AC 0
        R1 in out 1k
        C1 out 0 1p
        .NOISE V(out) V1 DEC 50 1 1e11
        .END
        """,
    )
    out = work / "rc_thermal_noise.csv"
    run([str(gspice), "--threads", "1", "--format", "csv", "-o", str(out), str(deck)], work)
    rows = read_csv(out)
    measured = rows[-1]["integrated_psd"]
    expected = 1.380649e-23 * 298.15 / 1e-12
    ratio = measured / expected
    if not 0.85 <= ratio <= 1.15:
        raise RuntimeError(
            f"RC thermal noise failed: integrated={measured:.6e}, "
            f"kT/C={expected:.6e}, ratio={ratio:.3f}"
        )
    print(f"RC thermal noise: integrated={measured:.6e}, kT/C={expected:.6e}, ratio={ratio:.3f}")


def transient_seed_deck(seed: int) -> str:
    return f"""
    * Transient-noise repeatability deck.
    V1 in 0 PULSE(0 1 0 1n 1n 5n 10n)
    R1 in out 1k
    C1 out 0 1p
    .SAVE V(out)
    .OPTIONS TRAN_PROGRESS_INTERVAL=0
    .TRANNOISE 0.25n 8n FMAX=1e9 SEED={seed} SCALE=50 NOISEMODE=SDE OVERSAMPLE=2
    .END
    """


def validate_transient_seed_behavior(gspice: Path, work: Path) -> None:
    samples = []
    for name, seed in [("a", 1234), ("b", 1234), ("c", 4321)]:
        deck = write(work / f"trannoise_seed_{name}.sp", transient_seed_deck(seed))
        out = work / f"trannoise_seed_{name}.csv"
        run([str(gspice), "--threads", "1", "--format", "csv", "-o", str(out), str(deck)], work)
        samples.append(read_transient_csv(out, "V(out)"))
    same = max(abs(x[1] - y[1]) for x, y in zip(samples[0], samples[1]))
    diff = max(abs(x[1] - y[1]) for x, y in zip(samples[0], samples[2]))
    if same > 1e-15:
        raise RuntimeError(f"same-seed transient noise is not repeatable: max diff={same:.3e}")
    if diff < 1e-15:
        raise RuntimeError(f"different transient-noise seeds did not change waveform: max diff={diff:.3e}")
    print(f"Transient noise seeds: same_seed_max={same:.3e}, different_seed_max={diff:.3e}")


def validate_osdi_flicker(gspice: Path, source: Path, work: Path) -> None:
    osdi = source / "osdi" / "psp103.osdi"
    if not osdi.exists():
        print(f"SKIP: OSDI PSP model not found: {osdi}")
        sys.exit(77)
    deck = write(
        work / "osdi_flicker_noise.sp",
        f"""
        * OSDI PSP flicker trend validation.
        .PRE_OSDI "{osdi}"
        .MODEL pch psp103va type=-1
        VDD vdd 0 DC 2 AC 0
        VG gate 0 DC 0
        RLOAD out 0 100k
        N1 out gate vdd vdd pch w=1u l=0.13u nf=1 mult=1
        .NOISE V(out) VDD VALUES 1 10 100 1k
        .END
        """,
    )
    out = work / "osdi_flicker_noise.csv"
    run([str(gspice), "--threads", "1", "--format", "csv", "-o", str(out), str(deck)], work)
    rows = read_csv(out)
    low = rows[0]["onoise_psd"]
    high = rows[-1]["onoise_psd"]
    flicker_count = rows[0]["flicker"]
    if flicker_count < 1:
        raise RuntimeError("OSDI flicker validation did not report flicker noise sources")
    if not low > high:
        raise RuntimeError(f"OSDI flicker trend failed: onoise_psd(1Hz)={low:.6e}, high={high:.6e}")

    tran_deck = source / "tests" / "decks" / "osdi_psp_trannoise.sp"
    completed = run([str(gspice), "--threads", "1", str(tran_deck)], source)
    match = re.search(r"Transient noise summary:.*flicker=([1-9][0-9]*).*colored=([1-9][0-9]*)", completed.stdout)
    if not match:
        raise RuntimeError("OSDI transient colored-noise summary did not show flicker/colored sources")
    print(
        f"OSDI flicker trend: onoise_psd(1Hz)={low:.6e}, onoise_psd(1kHz)={high:.6e}, "
        f"flicker_sources={int(flicker_count)}"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gspice", required=True, type=Path)
    parser.add_argument("--source", default=Path.cwd(), type=Path)
    parser.add_argument(
        "--case",
        choices=["all", "rc", "seed", "osdi"],
        default="all",
        help="Validation case to run.",
    )
    args = parser.parse_args()
    source = args.source.resolve()
    gspice = args.gspice.resolve()
    with tempfile.TemporaryDirectory(prefix="gspice_noise_validate_") as tmp:
        work = Path(tmp)
        if args.case in {"all", "rc"}:
            validate_rc_thermal(gspice, work)
        if args.case in {"all", "seed"}:
            validate_transient_seed_behavior(gspice, work)
        if args.case in {"all", "osdi"}:
            validate_osdi_flicker(gspice, source, work)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
