"""Audit native compact-model parameter coverage.

bsim4 mode   : official BSIM4.8.3 parameter set (from bsim4def.h structs)
               vs names the native evaluator actually reads.
psp103 mode  : implemented PSP103 allow-list vs the full deck corpus.
deck-scan    : for every corpus deck that instantiates a compact model,
               report which model-card parameters are ignored by the native
               evaluator (set-but-ignored, official name) and which are typos
               (neither official nor implemented).

Exit code 1 if any deck has a set-but-ignored parameter; 0 otherwise.
"""

from __future__ import annotations

import argparse
import importlib.util
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
OFFICIAL_HEADER = REPO.parent / "BSIM4_4.8.3_standard_05192025" / "code" / "bsim4def.h"
DEVICES = REPO / "include" / "devices"
DECKS = REPO / "tests" / "decks"

GET_BRACES = re.compile(r"\bget\(\s*\{(.*?)\}\s*[),]", re.S)
QUOTED = re.compile(r'"(?:[^"\\]|\\.)*"')
BLOCK_QUOTED = re.compile(r'"(?:[^"\\]|\\.)*"')
OFFICIAL_FIELD = re.compile(r"\bdouble\s+BSIM4([a-z][A-Za-z0-9_]*)\s*;")
MODEL_BLOCK = re.compile(r"(?m)^\s*\.model\s+(\S+)\s+(\S+)[\s\S]*?(?=^\s*\.(?!model)|\Z)")
LEVEL_TOKEN = re.compile(r"\blevel\s*[=:]?\s*(\d+)", re.I)

FAMILIES = [
    (r"^(W|L)[a-z]", "geometry-multiplier (l*/w*)"),
    (r"^(DLC|DWC|DLL|DLWL|DWLC|DWB|DWG|DWR|DCC|DMC|DLCAD|DWCAD|DLOV|DWOV)",
     "geometry-deviation"),
    (r"^(AGIDL|BGIDL|CGIDL)", "GIDL"),
    (r"^CGD|^CGS", "overlap-capacitance"),
    (r"^CJ|^PB|^MJ|^XJ", "junction"),
    (r"^NOI|^NLEV|^TNOIMOD", "noise"),
    (r"^EF$|^AF$|^KF$|^EM$|^NOIMOD$", "flicker"),
    (r"^XPART$|^CAPMOD$", "charge-switch"),
    (r"^(QOV|QSUB0|QACC|QBD|QBS)", "charge-groups"),
    (r"^UA1|^UB1|^UC1|^UD1", "temp-mobility"),
    (r"^(KT1|KT2|UTE)(L)?", "vth-temp"),
    (r"^(PDIBL|DR|PCLM|PVAG|FPROUT|PDITS|PDISC)", "rout"),
    (r"^(A0|AGS|A1|A2|A3|B0|B1)$", "post-mobility"),
    (r"^(BK|BG)", "dual-k"),
    (r"^(K0|K1|K2|KETA|K3|K3B|W0|VBN)", "body"),
    (r"^(ETA|NFACTOR|VTH0|VFB|VNTH|LT0)", "subthreshold-core"),
    (r"^(DVT|DSUB|CDSC)", "short-channel"),
    (r"^(N0|N1|N2|N3|N4|N5|N6|N7|N8)", "mobility-basics"),
    (r"^(UO|UA|UB|UC|UD|UP|VSAT|MU0)", "mobility-basics"),
    (r"^(TOX|TOXE|TOXM|D_TOX|NTOX)", "oxide"),
    (r"^(XJ|XJ0|NDEP|NSD|PHIN|PHI)", "physics"),
    (r"^(VTH0|VTH|VT0)", "threshold"),
]
DEFAULT_FAMILY = "other"


def extract_cpp_read_set(patterns: tuple[str, ...] = ("bsim4*.hpp",)) -> set[str]:
    names: set[str] = set()
    for pattern in patterns:
        for path in DEVICES.glob(pattern):
            text = path.read_text(encoding="utf-8", errors="replace")
            for match in GET_BRACES.finditer(text):
                for quoted in QUOTED.findall(match.group(1)):
                    names.add(quoted.strip('"').upper())
    return names


def extract_official_bsim4() -> set[str]:
    if not OFFICIAL_HEADER.exists():
        return set()
    text = OFFICIAL_HEADER.read_text(encoding="utf-8", errors="replace")
    found: set[str] = set()
    for match in OFFICIAL_FIELD.finditer(text):
        found.add(match.group(1).upper())
    return found


def extract_psp103_implemented() -> set[str]:
    path = DEVICES / "psp103_parameters.hpp"
    text = path.read_text(encoding="utf-8", errors="replace")
    block = re.search(r"handled\s*=\s*\{(.*?)\}", text, re.S)
    names: set[str] = set()
    for quoted in BLOCK_QUOTED.findall(block.group(1)):
        names.add(quoted.strip('"').upper())
    return names


def family_of(name: str) -> str:
    for pattern, label in FAMILIES:
        if re.match(pattern, name, re.I):
            return label
    return DEFAULT_FAMILY


def grouped_missing(based: set[str], given: set[str]) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {}
    for name in sorted(given - based):
        grouped.setdefault(family_of(name), []).append(name)
    return grouped


INSTANCE_TOKENS = {
    "L", "W", "M", "NF", "NFIN", "AD", "AS", "PD", "PS", "NRD", "NRS", "RD",
    "RS", "SA", "SB", "SD", "MULT", "DELVTO", "LDT", "LD", "WD", "LINT", "WINT",
}


def load_psp_official() -> set[str] | None:
    try:
        spec = importlib.util.spec_from_file_location(
            "psp103_reference_oracle",
            REPO / "tools" / "psp103_reference_oracle.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return set(module.extract_parameter_table())
    except Exception:
        return None


def deck_scan(implemented: set[str], official: set[str]) -> list[tuple[str, str, str, str]]:
    findings: list[tuple[str, str, str, str]] = []
    sources = sorted(DECKS.glob("*.sp")) + sorted(DECKS.glob("*.mod"))
    for path in sources:
        text = path.read_text(encoding="utf-8", errors="replace")
        for block in MODEL_BLOCK.finditer(text):
            model_name = block.group(1).strip()
            model_type = block.group(2).strip()
            body = block.group(0)
            is_psp = "psp" in model_type.lower()
            if "psp" not in model_type.lower() and (
                    not model_type.upper().startswith(("NMOS", "PMOS"))
                    or not LEVEL_TOKEN.search(body)):
                continue
            tokens = set(re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*=", body))
            psp_table = load_psp_official() if is_psp else None
            for token in sorted(tokens):
                upper = token.upper()
                if upper in {"LEVEL"} or upper in implemented or upper in INSTANCE_TOKENS:
                    continue
                if psp_table is not None and upper in psp_table:
                    status = "IGNORED-OFFICIAL"
                else:
                    status = "IGNORED-OFFICIAL" if upper in official else "UNKNOWN"
                findings.append((path.name, model_name, upper, status))
    return findings


def print_grouped(title: str, grouped: dict[str, list[str]]):
    print(f"\n{title}")
    for family in sorted(grouped):
        if family == DEFAULT_FAMILY:
            continue
        print(f"  [{family}] {', '.join(grouped[family])}")
    if DEFAULT_FAMILY in grouped:
        print(f"  [other] {', '.join(grouped[DEFAULT_FAMILY])}")


def allowlist_block_exists() -> set[str] | None:
    path = DEVICES / "bsim4_parameters.hpp"
    text = path.read_text(encoding="utf-8", errors="replace")
    block = re.search(r"Bsim4ImplementedParameters.*?supported\s*=\s*\{(.*?)\};", text, re.S)
    if not block:
        return None
    out: set[str] = set()
    for quoted in BLOCK_QUOTED.findall(block.group(1)):
        out.add(quoted.strip('"').upper())
    return out


def extract_allowlist_declared() -> set[str] | None:
    path = DEVICES / "bsim4_parameters.hpp"
    text = path.read_text(encoding="utf-8", errors="replace")
    block = re.search(r"Bsim4ImplementedParameters.*?supported\s*=\s*\{(.*?)\};", text, re.S)
    if not block:
        return None
    out: set[str] = set()
    for quoted in BLOCK_QUOTED.findall(block.group(1)):
        out.add(quoted.strip('"').upper())
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["bsim4", "psp103", "deck-scan", "check-allowlist"],
                        default="bsim4")
    args = parser.parse_args()

    if args.mode == "check-allowlist":
        read_set = extract_cpp_read_set()
        declared = extract_allowlist_declared()
        if declared is None:
            print("check-allowlist: Bsim4ImplementedParameters allow-list block not found.")
            return 1
        missing = read_set - declared
        stale = declared - read_set
        for name in sorted(missing):
            print(f"read by get() but MISSING from allow-list: {name}")
        for name in sorted(stale):
            print(f"declared in allow-list but not read by get(): {name}")
        if missing or stale:
            print("check-allowlist: DRIFT between native read-set and parser allow-list.")
            return 1
        print(f"check-allowlist: allow-list matches read-set ({len(read_set)} names).")
        return 0

    if args.mode == "bsim4":
        official = extract_official_bsim4()
        implemented = extract_cpp_read_set()
        print(f"official BSIM4.8.3 parameters (bsim4def.h): {len(official)}")
        print(f"parameters read by native evaluator: {len(implemented)}")
        missing = grouped_missing(official, implemented)
        covered = len(official) - sum(len(v) for v in missing.values())
        print(f"implemented-and-official: {covered} ({covered / len(official):.1%})")
        print_grouped("official params NOT implemented by native evaluator:", missing)
        return 0

    if args.mode == "psp103":
        implemented = extract_psp103_implemented()
        print(f"PSP103 parameters handled by native evaluator: {len(implemented)}")
        print(" ".join(sorted(implemented)))
        return 0

    implemented = extract_cpp_read_set(("bsim4*.hpp", "bsim3*.hpp")) | extract_psp103_implemented()
    official = OFFICIAL_HEADER.exists() and extract_official_bsim4() or set()
    findings = deck_scan(implemented, official)
    if not findings:
        print("deck-scan: no set-but-ignored compact-model parameters found.")
        return 0
    for deck, model, param, status in findings:
        print(f"{deck}: .model {model} sets {param} -> {status}")
    print(f"\ndeck-scan: {len(findings)} set-but-ignored parameter(s).")
    return 1


if __name__ == "__main__":
    sys.exit(main())