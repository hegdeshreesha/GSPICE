#!/usr/bin/env python3
"""Production readiness gate for the native GSPICE + Lumen IHP flow.

The script is intentionally conservative: it checks that the selected GSPICE
binary advertises the production-candidate native features, runs the broad
GSDI/GMC regression slice, and can optionally execute a real Lumen deck with a
zero-warning requirement.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path


REQUIRED_FEATURES = (
    "psp103_native",
    "ihp_psp_rf_full",
    "ihp_passive_wrapper_full",
    "bsim3_gsdi",
    "bsim4_gsdi",
    "juncap2_tat",
    "juncap2_avalanche",
)

BROAD_NATIVE_REGEX = (
    "gmc_veriloga_parser|gmc_generated_va|gmc_hidden|gmc_deep|gsdi|"
    "bsim|juncap|feature_registry|ihp|psp103|cli_capabilities"
)


def run_command(cmd: list[str], cwd: Path, timeout: int) -> dict:
    started = time.time()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
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


def load_capabilities(gspice: Path, cwd: Path) -> tuple[dict, dict]:
    result = run_command([str(gspice), "--capabilities"], cwd, 30)
    payload: dict = {}
    if result["return_code"] == 0:
        try:
            payload = json.loads(result["stdout"])
        except json.JSONDecodeError as exc:
            result["stderr"] += f"\nCould not parse capabilities JSON: {exc}"
            result["return_code"] = 1
    return payload, result


def check_capabilities(payload: dict) -> list[str]:
    failures: list[str] = []
    if payload.get("sparse_backend") != "SuiteSparse-KLU":
        failures.append("GSPICE must advertise sparse_backend=SuiteSparse-KLU")
    features = payload.get("features", {})
    for name in REQUIRED_FEATURES:
        info = features.get(name, {})
        if info.get("available") is not True:
            failures.append(f"required feature is unavailable: {name}")
        if info.get("maturity") not in {"tested", "validated"}:
            failures.append(f"required feature is not tested/validated: {name}")
    return failures


def run_lumen_deck(gspice: Path, deck: Path, output: Path, timeout: int) -> dict:
    if output.exists():
        try:
            output.unlink()
        except OSError:
            output = output.with_name(f"{output.stem}_{int(time.time())}{output.suffix}")
    result = run_command(
        [str(gspice), str(deck), "-o", str(output), "--format", "raw"],
        deck.parent,
        timeout,
    )
    combined = f"{result['stdout']}\n{result['stderr']}"
    result["warning_count"] = len(re.findall(r"(?i)\bwarning\b", combined))
    result["raw_exists"] = output.exists()
    result["raw_size"] = output.stat().st_size if output.exists() else 0
    result["raw_path"] = str(output)
    result["completed"] = "Simulation Completed Successfully" in combined
    return result


def run_lumen_fidelity_scan(gspice: Path, deck: Path, timeout: int) -> dict:
    env = os.environ.copy()
    env["GSPICE_VERBOSE_COMPAT_WARNINGS"] = "1"
    started = time.time()
    try:
        proc = subprocess.run(
            [str(gspice), str(deck), "--save", "none"],
            cwd=str(deck.parent),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
            env=env,
        )
        combined = f"{proc.stdout}\n{proc.stderr}"
        warnings = [
            line.strip()
            for line in combined.splitlines()
            if re.search(r"(?i)(unsupported .*parameter|ignored by native evaluator|approximated)", line)
        ]
        return {
            "return_code": proc.returncode,
            "elapsed_s": round(time.time() - started, 3),
            "compat_warning_count": len(warnings),
            "examples": warnings[:12],
            "timed_out": False,
        }
    except subprocess.TimeoutExpired as exc:
        combined = f"{exc.stdout or ''}\n{exc.stderr or ''}"
        warnings = [
            line.strip()
            for line in combined.splitlines()
            if re.search(r"(?i)(unsupported .*parameter|ignored by native evaluator|approximated)", line)
        ]
        return {
            "return_code": 124,
            "elapsed_s": round(time.time() - started, 3),
            "compat_warning_count": len(warnings),
            "examples": warnings[:12],
            "timed_out": True,
        }


def run_psp103_ignored_ranking(source: Path, gspice: Path, deck: Path, timeout: int) -> dict:
    report_path = source / "build-klu" / "psp103_ignored_ranking.json"
    result = run_command(
        [
            sys.executable,
            str(source / "tools" / "rank_psp103_ignored.py"),
            "--deck",
            str(deck),
            "--gspice",
            str(gspice),
            "--timeout",
            str(timeout),
            "--json",
            str(report_path),
        ],
        source,
        timeout,
    )
    payload: dict = {
        "return_code": result["return_code"],
        "elapsed_s": result["elapsed_s"],
        "timed_out": result["timed_out"],
        "report": str(report_path),
    }
    if report_path.exists():
        try:
            full = json.loads(report_path.read_text(encoding="utf-8"))
            payload["summary"] = full.get("summary", {})
        except json.JSONDecodeError as exc:
            payload["parse_error"] = str(exc)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GSPICE production readiness checks.")
    parser.add_argument("--source", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--build", default="")
    parser.add_argument("--gspice", default="")
    parser.add_argument("--lumen-deck", default="")
    parser.add_argument("--lumen-raw", default="")
    parser.add_argument("--skip-ctest", action="store_true")
    parser.add_argument("--strict-fidelity", action="store_true")
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--report", default="")
    args = parser.parse_args()

    source = Path(args.source).resolve()
    build = Path(args.build).resolve() if args.build else source / "build-klu"
    gspice = Path(args.gspice).resolve() if args.gspice else build / "Release" / "gspice.exe"
    report_path = Path(args.report).resolve() if args.report else source / "build-klu" / "production_signoff.json"

    report: dict = {
        "status": "UNKNOWN",
        "source": str(source),
        "build": str(build),
        "gspice": str(gspice),
        "required_features": list(REQUIRED_FEATURES),
        "checks": {},
        "failures": [],
        "risks": [],
    }

    capabilities, cap_result = load_capabilities(gspice, source)
    report["checks"]["capabilities"] = {
        "return_code": cap_result["return_code"],
        "elapsed_s": cap_result["elapsed_s"],
        "sparse_backend": capabilities.get("sparse_backend", ""),
    }
    if cap_result["return_code"] != 0:
        report["failures"].append("capabilities command failed")
    report["failures"].extend(check_capabilities(capabilities))

    if not args.skip_ctest:
        ctest = run_command(
            [
                "ctest",
                "--test-dir",
                str(build),
                "-C",
                "Release",
                "-R",
                BROAD_NATIVE_REGEX,
                "--output-on-failure",
            ],
            source,
            max(args.timeout, 240),
        )
        report["checks"]["ctest"] = {
            "return_code": ctest["return_code"],
            "elapsed_s": ctest["elapsed_s"],
            "timed_out": ctest["timed_out"],
        }
        if ctest["return_code"] != 0:
            report["failures"].append("broad native CTest suite failed")

    if args.lumen_deck:
        deck = Path(args.lumen_deck).resolve()
        raw = Path(args.lumen_raw).resolve() if args.lumen_raw else deck.parent / "production_signoff.raw"
        lumen = run_lumen_deck(gspice, deck, raw, args.timeout)
        report["checks"]["lumen_deck"] = {
            "return_code": lumen["return_code"],
            "elapsed_s": lumen["elapsed_s"],
            "warning_count": lumen["warning_count"],
            "raw_exists": lumen["raw_exists"],
            "raw_size": lumen["raw_size"],
            "raw_path": lumen["raw_path"],
            "completed": lumen["completed"],
        }
        if lumen["return_code"] != 0 or not lumen["completed"]:
            report["failures"].append("Lumen deck did not complete")
        if lumen["warning_count"] != 0:
            report["failures"].append("Lumen deck emitted warnings")
        if not lumen["raw_exists"] or lumen["raw_size"] <= 0:
            report["failures"].append("Lumen deck did not write RAW output")

        fidelity = run_lumen_fidelity_scan(gspice, deck, args.timeout)
        report["checks"]["lumen_fidelity_scan"] = fidelity
        if fidelity["compat_warning_count"]:
            message = (
                f"Lumen deck has {fidelity['compat_warning_count']} verbose compact-model "
                "fidelity warning(s); normal mode is quiet, but strict OSDI/OpenVAF parity is not proven"
            )
            report["risks"].append(message)
            if args.strict_fidelity:
                report["failures"].append(message)

        psp_ranking = run_psp103_ignored_ranking(source, gspice, deck, args.timeout)
        report["checks"]["psp103_ignored_ranking"] = psp_ranking
        summary = psp_ranking.get("summary", {})
        if summary.get("critical_or_high", 0):
            report["risks"].append(
                f"PSP103 ignored-parameter ranking has {summary['critical_or_high']} "
                "critical/high ignored or approximated occurrence(s); review psp103_ignored_ranking.json"
            )
            if args.strict_fidelity:
                report["failures"].append(
                    "strict fidelity requires zero critical/high PSP103 ignored or approximated occurrences"
                )
        if summary.get("accepted_approximate", 0):
            report["risks"].append(
                f"PSP103 scanner found {summary['accepted_approximate']} accepted-but-approximate "
                "parameter occurrence(s); normal mode is quiet, but full OpenVAF equation parity is not proven"
            )
            if args.strict_fidelity:
                report["failures"].append(
                    "strict fidelity requires zero accepted-but-approximate PSP103 parameter occurrences"
                )

    report["status"] = "PASS" if not report["failures"] else "FAIL"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
