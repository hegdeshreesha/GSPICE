"""Compare GSPICE results against a EXTERNAL_ORACLE-produced reference CSV.

EXTERNAL_ORACLE stays an oracle here, not a runtime dependency.  The script always runs
GSPICE for a named validation case and compares it with either:

* a precomputed EXTERNAL_ORACLE CSV (--external_oracle-csv), or
* a EXTERNAL_ORACLE command template (--external_oracle-command / EXTERNAL_ORACLE_COMMAND) that writes CSV.

The command template may use {deck}, {csv}, {case}, and {workdir}.  This keeps
the local EXTERNAL_ORACLE netlist-conversion details outside GSPICE and avoids making
normal builds depend on a second simulator.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import pathlib
import re
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass


REPO = pathlib.Path(__file__).resolve().parent.parent
SKIP = 77


@dataclass(frozen=True)
class Case:
    deck: str
    analysis: str
    columns: tuple[str, ...]
    compare: tuple[str, ...]


CASES: dict[str, Case] = {
    "ihp_lv_dc": Case(
        deck="tests/decks/ihp_lv_mos_dc_native.sp",
        analysis="dc",
        columns=("sweep", "V(vdd)", "V(gate)", "V(drain)"),
        compare=("V(drain)",),
    ),
    "ihp_lv_ac": Case(
        deck="tests/decks/ihp_lv_mos_ac_native.sp",
        analysis="ac",
        columns=("frequency", "V(vdd)", "V(gate)", "V(drain)"),
        compare=("V(drain).real", "V(drain).imag"),
    ),
    "ihp_lv_tran": Case(
        deck="tests/decks/ihp_lv_inverter_tran_native.sp",
        analysis="raw",
        columns=(),
        compare=("V(out)", "V(drain)", "V(net2)"),
    ),
    "ihp_lv_noise": Case(
        deck="tests/decks/ihp_lv_mos_noise_native.sp",
        analysis="noise",
        columns=("frequency", "onoise_sqrt", "onoise_psd", "noise_sources"),
        compare=("onoise_psd",),
    ),
    "ngspice_diode_dc": Case(
        deck="tests/decks/external_oracle_diode_dc.sp",
        analysis="dc",
        columns=("sweep", "V(in)", "V(out)"),
        compare=("V(out)",),
    ),
}


def run(args: list[str], *, cwd: pathlib.Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=cwd, check=False, capture_output=True, text=True)


def parse_float(text: str) -> float:
    return float(text.strip())


def parse_complex(text: str) -> tuple[float, float]:
    match = re.match(r"\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)", text.strip())
    if not match:
        raise ValueError(f"not a complex tuple: {text!r}")
    return float(match.group(1)), float(match.group(2))


def parse_dc_stdout(stdout: str, case: Case) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    pattern = re.compile(r"^\s*([-+0-9.eE]+)\s*\|\s*(.*)$")
    for line in stdout.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        values = [parse_float(match.group(1))]
        values.extend(parse_float(part) for part in match.group(2).split())
        if len(values) == len(case.columns):
            rows.append(dict(zip(case.columns, values)))
    return rows


def parse_ac_stdout(stdout: str, case: Case) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    pattern = re.compile(r"^\s*([-+0-9.eE]+)\s*\|\s*(.*)$")
    for line in stdout.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        complexes = re.findall(r"\([^)]*\)", match.group(2))
        if len(complexes) != len(case.columns) - 1:
            continue
        row: dict[str, float] = {"frequency": parse_float(match.group(1))}
        for name, text in zip(case.columns[1:], complexes):
            real, imag = parse_complex(text)
            row[f"{name}.real"] = real
            row[f"{name}.imag"] = imag
        rows.append(row)
    return rows


def parse_noise_stdout(stdout: str, case: Case) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    pattern = re.compile(r"^\s*([-+0-9.eE]+)\s*\|\s*(.*)$")
    for line in stdout.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        values = [parse_float(match.group(1))]
        values.extend(parse_float(part) for part in match.group(2).split())
        if len(values) == len(case.columns):
            rows.append(dict(zip(case.columns, values)))
    return rows


def parse_raw(path: pathlib.Path) -> list[dict[str, float]]:
    variables: list[str] = []
    rows: list[dict[str, float]] = []
    in_variables = False
    in_values = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line == "Variables:":
            in_variables = True
            continue
        if line == "Values:":
            in_variables = False
            in_values = True
            continue
        if in_variables:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                variables.append(parts[1])
            continue
        if in_values:
            parts = line.split()
            if len(parts) == len(variables) + 1 and parts[0].isdigit():
                parts = parts[1:]
            if len(parts) == len(variables):
                rows.append({name: parse_float(value) for name, value in zip(variables, parts)})
    return rows


def run_gspice(gspice: pathlib.Path, case: Case, source: pathlib.Path) -> list[dict[str, float]]:
    deck = source / case.deck
    with tempfile.TemporaryDirectory(prefix="gspice_external_oracle_") as tmp:
        raw = pathlib.Path(tmp) / "gspice.raw"
        result = run([str(gspice), "--threads", "1", "-o", str(raw), str(deck)], cwd=source)
        if result.returncode != 0:
            sys.stderr.write(result.stdout)
            sys.stderr.write(result.stderr)
            raise RuntimeError(f"GSPICE failed for {deck} with exit code {result.returncode}")
        if case.analysis == "dc":
            rows = parse_dc_stdout(result.stdout, case)
        elif case.analysis == "ac":
            rows = parse_ac_stdout(result.stdout, case)
        elif case.analysis == "noise":
            rows = parse_noise_stdout(result.stdout, case)
        elif raw.exists():
            rows = parse_raw(raw)
        else:
            rows = []
        if not rows:
            raise RuntimeError(f"no GSPICE rows parsed for {deck}")
        return rows


def read_csv(path: pathlib.Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise RuntimeError(f"CSV has no header: {path}")
        return [{key: parse_float(value) for key, value in row.items() if value != ""} for row in reader]


def write_csv(path: pathlib.Path, rows: list[dict[str, float]]) -> None:
    names: list[str] = []
    for row in rows:
        for key in row:
            if key not in names:
                names.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def resolve_external_oracle_csv(args: argparse.Namespace, case_name: str, deck: pathlib.Path) -> pathlib.Path | None:
    if args.external_oracle_csv:
        path = pathlib.Path(args.external_oracle_csv)
        return path / f"{case_name}.csv" if path.is_dir() else path

    template = args.external_oracle_command or os.environ.get("EXTERNAL_ORACLE_COMMAND")
    exe = args.external_oracle_exe or os.environ.get("EXTERNAL_ORACLE_EXE")
    if not template and not exe:
        return None

    workdir = pathlib.Path(tempfile.mkdtemp(prefix="external_oracle_oracle_"))
    csv_path = workdir / f"{case_name}.csv"
    if template:
        command = template.format(
            deck=str(deck),
            csv=str(csv_path),
            case=case_name,
            workdir=str(workdir),
        )
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
    else:
        result = run([str(exe), str(deck), "--csv", str(csv_path)], cwd=workdir)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"EXTERNAL_ORACLE command failed with exit code {result.returncode}")
    if not csv_path.exists():
        raise RuntimeError(f"EXTERNAL_ORACLE did not create expected CSV: {csv_path}")
    return csv_path


def choose_compare_columns(case: Case, gspice_rows: list[dict[str, float]], external_oracle_rows: list[dict[str, float]]) -> list[str]:
    available = set(gspice_rows[0]) & set(external_oracle_rows[0])
    chosen = [name for name in case.compare if name in available]
    if chosen:
        return chosen
    ignored = {"time", "frequency", "sweep"}
    chosen = sorted(name for name in available if name not in ignored)
    if not chosen:
        raise RuntimeError("no common numeric columns to compare")
    return chosen


def compare(
    case_name: str,
    case: Case,
    gspice_rows: list[dict[str, float]],
    external_oracle_rows: list[dict[str, float]],
    reltol: float,
    abstol: float,
) -> int:
    if len(gspice_rows) != len(external_oracle_rows):
        raise RuntimeError(f"row-count mismatch: gspice={len(gspice_rows)} external_oracle={len(external_oracle_rows)}")
    columns = choose_compare_columns(case, gspice_rows, external_oracle_rows)
    failures = 0
    worst = 0.0
    worst_text = ""
    for index, (g_row, v_row) in enumerate(zip(gspice_rows, external_oracle_rows)):
        for name in columns:
            actual = g_row[name]
            expected = v_row[name]
            delta = abs(actual - expected)
            relative = delta / max(abs(expected), abstol)
            if relative > worst:
                worst = relative
                worst_text = f"row={index} column={name} gspice={actual:.12e} external_oracle={expected:.12e}"
            if delta > abstol and relative > reltol:
                failures += 1
                print(f"FAIL {case_name}: row={index} column={name} "
                      f"gspice={actual:.12e} external_oracle={expected:.12e} rel={relative:.3e}")
    print(f"EXTERNAL_ORACLE parity {case_name}: rows={len(gspice_rows)} columns={','.join(columns)} "
          f"failures={failures} worst_relative={worst:.3e} {worst_text}")
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, default=REPO)
    parser.add_argument("--gspice", type=pathlib.Path)
    parser.add_argument("--case", choices=sorted(CASES), default="ihp_lv_dc")
    parser.add_argument("--list-cases", action="store_true")
    parser.add_argument("--external_oracle-csv")
    parser.add_argument("--external_oracle-command")
    parser.add_argument("--external_oracle-exe")
    parser.add_argument("--emit-gspice-csv", type=pathlib.Path)
    parser.add_argument("--require-external_oracle", action="store_true")
    parser.add_argument("--relative-tolerance", type=float, default=0.05)
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-12)
    args = parser.parse_args()

    if args.list_cases:
        for name, case in sorted(CASES.items()):
            print(f"{name}: {case.deck} ({case.analysis})")
        return 0
    if not args.gspice:
        parser.error("--gspice is required unless --list-cases is used")

    source = args.source.resolve()
    case = CASES[args.case]
    deck = source / case.deck
    gspice_rows = run_gspice(args.gspice.resolve(), case, source)
    if args.emit_gspice_csv:
        write_csv(args.emit_gspice_csv, gspice_rows)
        print(f"wrote GSPICE CSV: {args.emit_gspice_csv}")

    external_oracle_csv = resolve_external_oracle_csv(args, args.case, deck)
    if not external_oracle_csv:
        message = "EXTERNAL_ORACLE oracle not configured; set EXTERNAL_ORACLE_COMMAND, EXTERNAL_ORACLE_EXE, or pass --external_oracle-csv"
        if args.require_external_oracle:
            print(message, file=sys.stderr)
            return 2
        print(message + " (skipped)")
        return SKIP
    external_oracle_rows = read_csv(external_oracle_csv)
    return compare(args.case, case, gspice_rows, external_oracle_rows,
                   args.relative_tolerance, args.absolute_tolerance)


if __name__ == "__main__":
    sys.exit(main())
