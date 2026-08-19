#!/usr/bin/env python3
"""Rank ignored PSP103 model-card parameters by expected fidelity impact.

The tool is intentionally conservative.  It does not claim numerical error;
it tells us which unsupported official PSP parameters are most likely to
explain a speed/fidelity difference versus an external PSP103 reference run.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MODEL_BLOCK = re.compile(r"(?im)^\s*\.model\s+(\S+)\s+(\S+)(?:[^\n]*(?:\n\s*\+.*)*)")
PARAM_TOKEN = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*=")
INCLUDE_TOKEN = re.compile(r"(?im)^\s*\.(?:include|lib)\s+\"?([^\"\s]+)\"?")


BUCKETS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("critical", "noise", ("NOI", "FNOI", "FN", "ALPNOI", "FNT", "KF", "AF", "EF")),
    ("critical", "self-heating", ("RTH", "CTH", "STRTH", "SH", "STRT")),
    ("high", "junction-leakage-capacitance",
     ("CJ", "CJO", "CJOR", "IDSAT", "VBR", "VBI", "FBBT", "CTAT", "CSR", "PBR",
      "PHIG", "XJUN", "SWJUN", "TRJ", "JUN", "AVL")),
    ("high", "rf-series-parasitics",
     ("RGO", "RSHG", "RBULK", "RWELL", "RJUN", "RVPOLY", "RG", "RS", "RD")),
    ("high", "gate-leakage-gidl",
     ("IG", "AIG", "BIG", "CIG", "GIDL", "AGIDL", "BGIDL", "CGIDL")),
    ("medium", "geometry-binning",
     ("LP", "WP", "A1", "A2", "A3", "DLS", "DWS", "DL", "DW", "LVAR", "WVAR")),
    ("medium", "temperature-scaling",
     ("TC", "TCE", "TBG", "TR", "ST", "XT", "UTE", "KT")),
    ("medium", "charge-capacitance",
     ("Q", "CG", "CAP", "CFR", "COV", "NOV", "LOV", "WOV")),
    ("low", "selector-or-corner-control", ("SW", "VERSION", "LEVEL", "TYPE")),
)

SEVERITY_ORDER = {"critical": 0, "high": 1, "medium": 2, "low": 3, "unknown": 4}

APPROXIMATE_PARAMS: dict[str, tuple[str, str, str]] = {
    # Native PSP accepts these aliases to keep IHP decks quiet, but the full
    # OpenVAF/CMC sub-equation surface is not yet reproduced one-for-one.
    "SWIMPACT": ("high", "impact-ionization", "switch accepted; full PSP impact-ionization branch is not equation-complete"),
    "SWNQS": ("high", "non-quasi-static", "switch accepted; NQS state equations are not implemented"),
    "SWNUD": ("high", "non-uniform-doping", "switch accepted; NUD correction is approximated"),
    "SWEDGE": ("high", "edge-transistor", "switch accepted; edge-transistor current/charge/noise is approximated"),
    "SWDELVTAC": ("medium", "ac-threshold-shift", "switch accepted; AC delta-VT correction is approximated"),
    "IGINVLW": ("high", "gate-leakage", "length/width gate-leakage correction accepted but not used as an independent term"),
    "IGOVW": ("high", "gate-leakage", "overlap gate-leakage width correction accepted but not used as an independent term"),
    "IGOVDW": ("high", "gate-leakage", "drain-overlap gate-leakage width correction accepted but not used as an independent term"),
    "CJORGAT": ("high", "juncap2-sidewall", "gate-side junction capacitance accepted but folded into the conservative native junction model"),
    "CJORGATD": ("high", "juncap2-sidewall", "drain gate-side junction capacitance accepted but folded into the conservative native junction model"),
    "CJORSTI": ("high", "juncap2-sidewall", "STI-side junction capacitance accepted but folded into the conservative native junction model"),
    "CJORSTID": ("high", "juncap2-sidewall", "drain STI-side junction capacitance accepted but folded into the conservative native junction model"),
    "PBRGAT": ("high", "juncap2-sidewall", "gate-side breakdown potential accepted but approximated"),
    "PBRGATD": ("high", "juncap2-sidewall", "drain gate-side breakdown potential accepted but approximated"),
    "PBRSTI": ("high", "juncap2-sidewall", "STI-side breakdown potential accepted but approximated"),
    "PBRSTID": ("high", "juncap2-sidewall", "drain STI-side breakdown potential accepted but approximated"),
    "XJUNGAT": ("high", "juncap2-sidewall", "gate-side grading accepted but approximated"),
    "XJUNGATD": ("high", "juncap2-sidewall", "drain gate-side grading accepted but approximated"),
    "XJUNSTI": ("high", "juncap2-sidewall", "STI-side grading accepted but approximated"),
    "XJUNSTID": ("high", "juncap2-sidewall", "drain STI-side grading accepted but approximated"),
    "FBBTRBOT": ("high", "juncap2-bbt", "bottom BBT temperature branch accepted but approximated"),
    "FBBTRBOTD": ("high", "juncap2-bbt", "drain bottom BBT temperature branch accepted but approximated"),
    "FBBTRGAT": ("high", "juncap2-bbt", "gate-side BBT branch accepted but approximated"),
    "FBBTRGATD": ("high", "juncap2-bbt", "drain gate-side BBT branch accepted but approximated"),
    "FBBTRSTI": ("high", "juncap2-bbt", "STI-side BBT branch accepted but approximated"),
    "FBBTRSTID": ("high", "juncap2-bbt", "drain STI-side BBT branch accepted but approximated"),
    "CTATBOT": ("high", "juncap2-tat", "bottom TAT branch accepted but native equation is smoothed/approximate"),
    "CTATBOTD": ("high", "juncap2-tat", "drain bottom TAT branch accepted but native equation is smoothed/approximate"),
    "CTATGAT": ("high", "juncap2-tat", "gate-side TAT branch accepted but approximated"),
    "CTATGATD": ("high", "juncap2-tat", "drain gate-side TAT branch accepted but approximated"),
    "CTATSTI": ("high", "juncap2-tat", "STI-side TAT branch accepted but approximated"),
    "CTATSTID": ("high", "juncap2-tat", "drain STI-side TAT branch accepted but approximated"),
    "CSRHGAT": ("high", "juncap2-srh", "gate-side SRH branch accepted but approximated"),
    "CSRHGATD": ("high", "juncap2-srh", "drain gate-side SRH branch accepted but approximated"),
    "CSRHSTI": ("high", "juncap2-srh", "STI-side SRH branch accepted but approximated"),
    "CSRHSTID": ("high", "juncap2-srh", "drain STI-side SRH branch accepted but approximated"),
    "FNTEDGEO": ("critical", "edge-noise", "edge flicker-noise parameter accepted but exact PSP edge-noise network is not implemented"),
    "EFEDGEO": ("critical", "edge-noise", "edge noise exponent accepted but exact PSP edge-noise network is not implemented"),
    "BETEDGEW": ("high", "edge-transistor", "edge beta geometry term accepted but approximated"),
    "FBETEDGE": ("high", "edge-transistor", "edge beta factor accepted but approximated"),
    "LPEDGE": ("high", "edge-transistor", "edge length parameter accepted but approximated"),
    "WEDGE": ("high", "edge-transistor", "edge width parameter accepted but approximated"),
}


def native_psp_read_names() -> set[str]:
    text = (REPO / "include" / "devices" / "psp103_parameters.hpp").read_text(
        encoding="utf-8", errors="replace")
    block = re.search(r"handled\s*=\s*\{(.*?)\}", text, re.S)
    if not block:
        raise RuntimeError("could not find PSP103 handled parameter set")
    return {name.upper() for name in re.findall(r'"([^"]+)"', block.group(1))}


def official_psp_names() -> set[str]:
    sys.path.insert(0, str(REPO / "tools"))
    import psp103_reference_oracle  # pylint: disable=import-error,import-outside-toplevel

    return {name.upper() for name in psp103_reference_oracle.extract_parameter_table()}


def classify(name: str) -> tuple[str, str]:
    upper = name.upper()
    for severity, bucket, prefixes in BUCKETS:
        if any(upper.startswith(prefix) for prefix in prefixes):
            return severity, bucket
    return "unknown", "unclassified"


def resolve_sources(deck: Path, seen: set[Path] | None = None) -> list[Path]:
    if seen is None:
        seen = set()
    deck = deck.resolve()
    if deck in seen or not deck.exists():
        return []
    seen.add(deck)
    sources = [deck]
    text = deck.read_text(encoding="utf-8", errors="replace")
    for match in INCLUDE_TOKEN.finditer(text):
        raw = match.group(1).strip()
        include = Path(raw)
        if not include.is_absolute():
            include = deck.parent / include
        sources.extend(resolve_sources(include, seen))
    return sources


def scan_source(deck: Path, native: set[str], official: set[str]) -> list[dict]:
    text = deck.read_text(encoding="utf-8", errors="replace")
    findings: list[dict] = []
    for match in MODEL_BLOCK.finditer(text):
        model_name = match.group(1)
        model_type = match.group(2)
        is_psp = "psp" in model_type.lower()
        if not is_psp:
            continue
        tokens = {token.upper() for token in PARAM_TOKEN.findall(match.group(0))}
        ignored = sorted(name for name in tokens if name not in native and name != "LEVEL")
        for name in ignored:
            severity, bucket = classify(name)
            findings.append({
                "deck": str(deck),
                "model": model_name,
                "model_type": model_type,
                "parameter": name,
                "official": name in official,
                "severity": severity,
                "bucket": bucket,
            })
        approximated = sorted(name for name in tokens if name in native and name in APPROXIMATE_PARAMS)
        for name in approximated:
            severity, bucket, reason = APPROXIMATE_PARAMS[name]
            findings.append({
                "deck": str(deck),
                "model": model_name,
                "model_type": model_type,
                "parameter": name,
                "official": name in official,
                "severity": severity,
                "bucket": bucket,
                "status": "accepted_approximate",
                "reason": reason,
            })
    return findings


def run_gspice_verbose(gspice: Path, deck: Path, timeout: int) -> dict:
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
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
        return {
            "return_code": proc.returncode,
            "elapsed_s": round(time.time() - started, 3),
            "stdout": proc.stdout,
            "stderr": proc.stderr,
            "timed_out": False,
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "return_code": 124,
            "elapsed_s": round(time.time() - started, 3),
            "stdout": exc.stdout or "",
            "stderr": exc.stderr or "",
            "timed_out": True,
        }


def summarize(findings: list[dict]) -> dict:
    by_bucket: dict[str, dict] = {}
    by_model: dict[str, dict] = {}
    for item in findings:
        bucket = item["bucket"]
        severity = item["severity"]
        param = item["parameter"]
        entry = by_bucket.setdefault(bucket, {
            "severity": severity,
            "count": 0,
            "parameters": set(),
            "official_count": 0,
        })
        entry["count"] += 1
        entry["parameters"].add(param)
        if item["official"]:
            entry["official_count"] += 1

        model = item["model"]
        model_entry = by_model.setdefault(model, {
            "model_type": item["model_type"],
            "count": 0,
            "critical_or_high": 0,
        })
        model_entry["count"] += 1
        if severity in {"critical", "high"}:
            model_entry["critical_or_high"] += 1

    buckets = []
    for name, entry in by_bucket.items():
        buckets.append({
            "bucket": name,
            "severity": entry["severity"],
            "count": entry["count"],
            "official_count": entry["official_count"],
            "parameters": sorted(entry["parameters"]),
        })
    buckets.sort(key=lambda row: (SEVERITY_ORDER[row["severity"]], -row["count"], row["bucket"]))
    models = [{"model": name, **entry} for name, entry in by_model.items()]
    models.sort(key=lambda row: (-row["critical_or_high"], -row["count"], row["model"]))
    ignored_count = sum(1 for item in findings if item.get("status", "ignored") == "ignored")
    approximate_count = sum(1 for item in findings if item.get("status") == "accepted_approximate")
    return {
        "total_ignored": ignored_count,
        "accepted_approximate": approximate_count,
        "official_ignored": sum(1 for item in findings if item["official"] and item.get("status", "ignored") == "ignored"),
        "critical_or_high": sum(1 for item in findings if item["severity"] in {"critical", "high"}),
        "critical_or_high_approximate": sum(1 for item in findings if item.get("status") == "accepted_approximate" and item["severity"] in {"critical", "high"}),
        "buckets": buckets,
        "models": models,
    }


def printable(report: dict) -> str:
    lines = [
        f"PSP103 ignored-parameter ranking: total={report['summary']['total_ignored']} "
        f"accepted_approximate={report['summary']['accepted_approximate']} "
        f"official={report['summary']['official_ignored']} "
        f"critical_or_high={report['summary']['critical_or_high']}"
    ]
    for bucket in report["summary"]["buckets"]:
        params = ", ".join(bucket["parameters"][:16])
        if len(bucket["parameters"]) > 16:
            params += ", ..."
        lines.append(
            f"  [{bucket['severity']}] {bucket['bucket']}: "
            f"{bucket['count']} occurrence(s), {bucket['official_count']} official; {params}"
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck", type=Path, action="append", default=[],
                        help="SPICE deck or model file to scan. Can be repeated.")
    parser.add_argument("--gspice", type=Path,
                        help="Optional GSPICE executable; enables verbose runtime scan.")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--fail-on-critical", action="store_true")
    parser.add_argument("--fail-on-approximate", action="store_true",
                        help="Fail if accepted PSP103 parameters are only approximated by native PSP.")
    args = parser.parse_args()

    if not args.deck:
        parser.error("at least one --deck is required")

    native = native_psp_read_names()
    official = official_psp_names()
    findings: list[dict] = []
    runtime: list[dict] = []
    for deck in args.deck:
        resolved = deck.resolve()
        for source in resolve_sources(resolved):
            findings.extend(scan_source(source, native, official))
        if args.gspice:
            result = run_gspice_verbose(args.gspice.resolve(), resolved, args.timeout)
            runtime.append({
                "deck": str(resolved),
                "return_code": result["return_code"],
                "elapsed_s": result["elapsed_s"],
                "timed_out": result["timed_out"],
                "verbose_warning_count": len(re.findall(
                    r"unsupported PSP parameter\(s\) ignored by native evaluator",
                    f"{result['stdout']}\n{result['stderr']}")),
            })

    report = {
        "native_parameter_count": len(native),
        "official_parameter_count": len(official),
        "summary": summarize(findings),
        "runtime": runtime,
        "findings": findings,
    }
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(printable(report))
    if args.fail_on_critical and report["summary"]["critical_or_high"]:
        return 1
    if args.fail_on_approximate and report["summary"]["accepted_approximate"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
