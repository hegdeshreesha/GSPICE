#!/usr/bin/env python3
"""Prepare a VACASK comparison run for a Lumen/GSPICE deck.

This bridge keeps VACASK optional and external.  It copies the input deck to a
scratch directory, tries the bundled Ngspice-to-VACASK converter, optionally
runs VACASK on the converted deck, and emits a JSON status report.  It never
edits PDK files.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


DEFAULT_VACASK = (
    Path("C:/EDA/Tools/vacask_0.3.4.rc1/vacask_0.3.4.rc1_windows-x86_64/bin/vacask.exe")
)
DEFAULT_VACASK_PYTHON = (
    Path("C:/EDA/Tools/vacask_0.3.4.rc1/vacask_0.3.4.rc1_windows-x86_64/lib/python")
)
SKIP = 77
INCLUDE_TOKEN = __import__("re").compile(r"(?im)^\s*\.(?:include|lib)\s+\"?([^\"\s]+)\"?")


def run(cmd: list[str], cwd: Path, timeout: int, env: dict[str, str] | None = None) -> dict:
    started = time.time()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            timeout=timeout,
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        return {
            "command": cmd,
            "return_code": proc.returncode,
            "elapsed_s": round(time.time() - started, 3),
            "stdout": proc.stdout,
            "stderr": proc.stderr,
            "timed_out": False,
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "command": cmd,
            "return_code": 124,
            "elapsed_s": round(time.time() - started, 3),
            "stdout": exc.stdout or "",
            "stderr": exc.stderr or "",
            "timed_out": True,
        }


def converter_env(vacask_python: Path) -> dict[str, str]:
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(vacask_python) if not existing else str(vacask_python) + os.pathsep + existing
    return env


def inferred_source_paths(deck: Path) -> list[Path]:
    found: list[Path] = [deck.parent]
    text = deck.read_text(encoding="utf-8", errors="replace")
    for match in INCLUDE_TOKEN.finditer(text):
        raw = Path(match.group(1).strip())
        candidate = raw if raw.is_absolute() else deck.parent / raw
        parent = candidate.parent.resolve()
        if parent not in found:
            found.append(parent)
    return found


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, default=Path("C:/EDA/GSPICE/build-klu/vacask_bridge"))
    parser.add_argument("--vacask", type=Path, default=DEFAULT_VACASK)
    parser.add_argument("--vacask-python", type=Path, default=DEFAULT_VACASK_PYTHON)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--run-vacask", action="store_true")
    parser.add_argument("--require-vacask", action="store_true")
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()

    deck = args.deck.resolve()
    workdir = args.workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    copied = workdir / deck.name
    converted = workdir / (deck.stem + ".sim")
    shutil.copy2(deck, copied)

    report = {
        "status": "UNKNOWN",
        "deck": str(deck),
        "workdir": str(workdir),
        "copied_deck": str(copied),
        "converted_deck": str(converted),
        "steps": {},
        "notes": [],
    }

    if not args.vacask_python.exists():
        report["status"] = "SKIP"
        report["notes"].append(f"VACASK Python path not found: {args.vacask_python}")
    else:
        command = [args.python, "-m", "ng2vc"]
        for source_path in inferred_source_paths(deck):
            command.extend(["-sp", str(source_path)])
        command.extend([str(copied), str(converted)])
        convert = run(
            command,
            workdir,
            args.timeout,
            converter_env(args.vacask_python.resolve()),
        )
        report["steps"]["convert"] = {
            "return_code": convert["return_code"],
            "elapsed_s": convert["elapsed_s"],
            "timed_out": convert["timed_out"],
            "stdout_tail": convert["stdout"][-2000:],
            "stderr_tail": convert["stderr"][-2000:],
        }
        if convert["return_code"] != 0 or not converted.exists():
            report["status"] = "SKIP"
            report["notes"].append(
                "VACASK converter rejected this deck. This is expected for some "
                "Lumen/GSPICE syntax; use the report tails to extend conversion coverage."
            )
        else:
            report["status"] = "CONVERTED"

    if report["status"] == "CONVERTED" and args.run_vacask:
        if not args.vacask.exists():
            report["status"] = "SKIP"
            report["notes"].append(f"VACASK executable not found: {args.vacask}")
        else:
            vacask = run([str(args.vacask.resolve()), str(converted)], workdir, args.timeout)
            report["steps"]["vacask"] = {
                "return_code": vacask["return_code"],
                "elapsed_s": vacask["elapsed_s"],
                "timed_out": vacask["timed_out"],
                "stdout_tail": vacask["stdout"][-4000:],
                "stderr_tail": vacask["stderr"][-4000:],
            }
            report["status"] = "PASS" if vacask["return_code"] == 0 else "FAIL"

    report_path = args.report.resolve() if args.report else workdir / "vacask_lumen_bridge.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

    if report["status"] in {"PASS", "CONVERTED"}:
        return 0
    if report["status"] == "SKIP" and not args.require_vacask:
        return SKIP
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
