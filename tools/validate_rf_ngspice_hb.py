#!/usr/bin/env python3
import argparse
import math
import pathlib
import re
import subprocess
import tempfile


def run(args):
    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if proc.returncode != 0:
        print(proc.stdout)
        raise SystemExit(proc.returncode)
    return proc.stdout


def parse_gspice_second_harmonic(output):
    match = re.search(r"^2\s+2\.000000000e\+03 \|[^\n]*\| (\S+)", output, re.MULTILINE)
    if not match:
        raise AssertionError("could not find GSPICE HB second-harmonic V(out) row")
    return float(match.group(1))


def parse_ngspice_fourier_second_harmonic(log_text):
    for line in log_text.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0] == "2":
            return 0.5 * float(fields[2])
    raise AssertionError("could not find ngspice .four second-harmonic row")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gspice", required=True)
    parser.add_argument("--ngspice", required=True)
    parser.add_argument("--deck", required=True)
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-2)
    args = parser.parse_args()

    gspice_h2 = parse_gspice_second_harmonic(run([args.gspice, "--threads", "1", args.deck]))
    with tempfile.TemporaryDirectory(prefix="gspice_rf_ngspice_") as tmp:
        deck = pathlib.Path(tmp) / "hb_diode_fourier.cir"
        log = pathlib.Path(tmp) / "ngspice.log"
        deck.write_text(
            """* generated independent transient/Fourier oracle for the HB diode deck
V1 in 0 SIN(0.55 0.08 1k)
R1 in out 100
D1 out 0 DHB
.MODEL DHB D(IS=1e-14 N=1 CJO=0.1p)
.TRAN 0.2u 20m 15m
.FOUR 1k V(out)
.END
""",
            encoding="utf-8",
        )
        run([args.ngspice, "-b", "-o", str(log), str(deck)])
        ngspice_h2 = parse_ngspice_fourier_second_harmonic(
            log.read_text(encoding="utf-8", errors="replace")
        )

    rel_error = abs(gspice_h2 - ngspice_h2) / max(abs(ngspice_h2), 1.0e-30)
    print(
        "RF ngspice HB oracle: "
        f"gspice_h2={gspice_h2:.12g} ngspice_h2={ngspice_h2:.12g} "
        f"rel_error={rel_error:.6g} limit={args.relative_tolerance:.6g}"
    )
    if not math.isfinite(rel_error) or rel_error > args.relative_tolerance:
        raise AssertionError("RF ngspice HB oracle mismatch")


if __name__ == "__main__":
    main()
