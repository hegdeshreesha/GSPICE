"""PSP103 reference oracle (standalone transcription, no OSDI/OpenVAF).

Source of truth: the official PSP103 Verilog-A (103.6.0, NXP/CEA/TUD/ASU, ECL-2.0)
mirrored under tools/reference_sources/psp103_103.6/.

Modes
-----
--list-params            : dump the official parameter table extracted from
                           PSP103_module.include (941+ declarations).
--mod-card <file>        : audit a .mod card against the official VA table.
                           Reports unrecognized (typos) and official-but-not-yet
                           -in-C++ names alongside what the native evaluator
                           honors (Psp103ParameterSet::handlesModelParameter).
--points                 : evaluate a bias grid via the transcription.
                           Not yet implemented (transcription is staged:
                           geometry -> surface potential -> currents).

Nagenda (arithmetic transcription, staged):
   stage1  geometry/prepare: binning (LPLIC), DW/DL, temperature scaling
   stage2  surface-potential iteration (SPCalculation)
   stage3  DC currents (RDC/THESAT/ALP + gate/GIDL leakage)
   stage4  charges (source/drain/overlap/NQS)
   stage5  noise
Each stage mirrors the VA macros one-to-one (PSP103_module.include).
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
VA_DIR = REPO / "tools" / "reference_sources" / "psp103_103.6"
MODULE = VA_DIR / "PSP103_module.include"
COMMON_MACRODEFS = VA_DIR / "Common103_macrodefs.include"

PARAM_TABLE_RE = re.compile(
    r"`(?:MPRnb|MPRcz|MPRco|MPRcc|IPRnb|IPRco|IPRcc|MPIcc|MPRnc|MPRnco)\s*"
    r"\(\s*([A-Z][A-Za-z0-9_]*)\s*,\s*([-+0-9.eE]+)\s*,")


def extract_parameter_table() -> dict[str, float]:
    table: dict[str, float] = {}
    for path in (MODULE, VA_DIR / "JUNCAP200_parlist.include",
                 VA_DIR / "PSP103_binpars.include"):
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for match in PARAM_TABLE_RE.finditer(text):
            name, default = match.group(1), match.group(2)
            try:
                table[name] = float(default)
            except ValueError:
                continue
    return table


def parse_model_card(path: pathlib.Path) -> dict[str, float]:
    card: dict[str, float] = {}
    text = path.read_text(encoding="utf-8", errors="replace")
    for match in re.finditer(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([-+0-9.eE]+)", text):
        card[match.group(1).upper()] = float(match.group(2))
    return card


def native_read_names() -> set[str]:
    text = (REPO / "include" / "devices" / "psp103_parameters.hpp").read_text(
        encoding="utf-8", errors="replace")
    block = re.search(r"handled\s*=\s*\{(.*?)\}", text, re.S)
    names: set[str] = set()
    for quoted in re.findall(r'"([^"]*)"', block.group(1)):
        names.add(quoted.upper())
    return names


def run_list_params(table: dict[str, float]) -> int:
    print(f"official PSP103 VA parameter-table entries: {len(table)}")
    for name in sorted(table):
        text = f"{table[name]:.12g}".rstrip("0").rstrip(".")
        print(f"{name} = {text or '0'}")
    return 0


def run_check_card(path: pathlib.Path) -> int:
    table = extract_parameter_table()
    native = native_read_names()
    card = parse_model_card(path)
    not_in_va = sorted(set(card) - set(table))
    in_va_but_ignored = sorted(set(card) & set(table) - native)
    print(f"model card: {path.name}  ({len(card)} parameters)")
    print(f"official VA table: {len(table)} parameters; "
          f"native Psp103ParameterSet reads {len(native)}")
    if not_in_va:
        print("\nNOT in official PSP103 table (typos / custom extensions):")
        for name in not_in_va:
            print(f"  {name} = {card[name]:g}")
    if in_va_but_ignored:
        print(f"\nOFFICIAL but ignored by native evaluator ({len(in_va_but_ignored)}):")
        for name in in_va_but_ignored:
            print(f"  {name} = {card[name]:g}")
    used = sorted(set(card) & native)
    print(f"\nhonored by native evaluator ({len(used)}):")
    print("  " + " ".join(used))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-params", action="store_true")
    parser.add_argument("--check-card", type=pathlib.Path)
    parser.add_argument("--points", action="store_true")
    args = parser.parse_args()

    if args.points:
        print("error: --points transcription not implemented yet "
              "(stages: geometry -> SP -> DC; see module docstring)")
        return 1
    if args.check_card:
        return run_check_card(args.check_card)
    return run_list_params(extract_parameter_table())


if __name__ == "__main__":
    sys.exit(main())