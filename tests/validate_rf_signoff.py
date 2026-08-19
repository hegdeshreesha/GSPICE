#!/usr/bin/env python3
import argparse
import json
import math
import re
import subprocess
import sys


def run_deck(exe, deck):
    proc = subprocess.run(
        [exe, "--threads", "1", deck],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        raise SystemExit(proc.returncode)
    return proc.stdout


def parse_hb_table(output):
    header = None
    rows = {}
    for line in output.splitlines():
        if line.startswith("harmonic freq_hz"):
            header = line
            continue
        if not header:
            continue
        fields = [part.strip() for part in line.split("|")]
        if len(fields) < 2:
            continue
        head = fields[0].split()
        if len(head) != 2 or not head[0].isdigit():
            continue
        harmonic = int(head[0])
        values = []
        for part in fields[1:]:
            cols = part.split()
            if len(cols) >= 2:
                values.append((float(cols[0]), float(cols[1])))
        rows[harmonic] = values
    return rows


def assert_between(name, value, lo, hi):
    if not (lo <= value <= hi):
        raise AssertionError(f"{name}={value:.12g} outside [{lo:.12g}, {hi:.12g}]")
    print(f"ok {name}={value:.12g}")


DEFAULT_BOUNDS = {
    "hb_linear_rc": {
        "hb_linear_rc_ratio": None,
    },
    "hb_diode": {
        "hb_diode_second_harmonic": [1.8e-3, 2.1e-3],
    },
    "hb_multitone": {
        "hb_multitone_1khz_vout": [0.499, 0.501],
        "hb_multitone_1p5khz_vout": [0.499, 0.501],
    },
    "hbac_diode": {
        "hbac_diode_vout": [0.84, 0.88],
    },
    "hbnoise_diode": {
        "hbnoise_diode_psd": [1.0e-18, 1.8e-18],
    },
    "pnoise_phase_jitter": {
        "pnoise_phase_1khz_dbc_per_hz": [-166.0, -163.0],
        "pnoise_jitter_1khz_s_per_rtHz": [8.0e-13, 1.1e-12],
    },
    "hbsp": {
        "hbsp_s11_real": [-1e-8, 1e-8],
        "hbsp_s11_imag": [-1e-8, 1e-8],
    },
    "psssp": {
        "psssp_s11_real": [-1e-8, 1e-8],
        "psssp_s11_imag": [-1e-8, 1e-8],
    },
    "hbstb": {
        "hbstb_mag": [0.999, 1.001],
        "hbstb_phase_deg": [179.0, 181.0],
    },
    "pssstb": {
        "pssstb_mag": [0.999, 1.001],
        "pssstb_phase_deg": [179.0, 181.0],
    },
    "psspac_linear": {
        "psspac_linear_vout_1khz": [0.9998, 1.0001],
    },
    "hb_psp103": {
        "hb_psp103_drain_fundamental": [4.5e-3, 4.9e-3],
        "hb_psp103_drain_second_harmonic": [2.8e-4, 3.3e-4],
    },
}


def load_bounds(path):
    if not path:
        return dict(DEFAULT_BOUNDS)
    with open(path, "r", encoding="utf-8") as stream:
        data = json.load(stream)
    bounds = dict(DEFAULT_BOUNDS)
    bounds.update(data.get("cases", data))
    return bounds


def metric_bounds(bounds, case, metric):
    values = bounds[case][metric]
    if values is None:
        expected = 1.0 / math.sqrt(1.0 + (2.0 * math.pi * 1e3 * 1e3 * 1e-9) ** 2)
        return [expected - 2e-5, expected + 2e-5]
    if isinstance(values, dict):
        reference = float(values["reference"])
        abs_tol = float(values.get("abs_tol", 0.0))
        rel_tol = float(values.get("rel_tol", 0.0))
        tolerance = abs_tol + rel_tol * abs(reference)
        return [reference - tolerance, reference + tolerance]
    return values


def check_metric(bounds, case, metric, value):
    lo, hi = metric_bounds(bounds, case, metric)
    assert_between(metric, value, lo, hi)


def validate_hb_linear_rc(output, bounds):
    rows = parse_hb_table(output)
    vin = rows[1][0][0]
    vout = rows[1][1][0]
    ratio = vout / vin
    check_metric(bounds, "hb_linear_rc", "hb_linear_rc_ratio", ratio)


def validate_hb_diode(output, bounds):
    rows = parse_hb_table(output)
    second = rows[2][1][0]
    check_metric(bounds, "hb_diode", "hb_diode_second_harmonic", second)


def validate_hb_multitone(output, bounds):
    rows = parse_hb_table(output)
    h1 = rows[2][1][0]
    h2 = rows[3][1][0]
    check_metric(bounds, "hb_multitone", "hb_multitone_1khz_vout", h1)
    check_metric(bounds, "hb_multitone", "hb_multitone_1p5khz_vout", h2)


def validate_hbac_diode(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| V\(in\)=\S+ V\(out\)=(\S+)", output)
    if not match:
        raise AssertionError("could not find HBAC diode V(out) row")
    gain = float(match.group(1))
    check_metric(bounds, "hbac_diode", "hbac_diode_vout", gain)


def validate_hbnoise_diode(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| \S+ (\S+) 128", output)
    if not match:
        raise AssertionError("could not find HBNOISE diode PSD row")
    psd = float(match.group(1))
    check_metric(bounds, "hbnoise_diode", "hbnoise_diode_psd", psd)


def validate_pnoise_phase_jitter(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| \S+ \S+ \d+ (\S+) (\S+)", output)
    if not match:
        raise AssertionError("could not find PNOISE phase/jitter row")
    phase = float(match.group(1))
    jitter = float(match.group(2))
    check_metric(bounds, "pnoise_phase_jitter", "pnoise_phase_1khz_dbc_per_hz", phase)
    check_metric(bounds, "pnoise_phase_jitter", "pnoise_jitter_1khz_s_per_rtHz", jitter)


def validate_hbsp(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| (\S+) (\S+)", output)
    if not match:
        raise AssertionError("could not find HBSP S11 row")
    real = float(match.group(1))
    imag = float(match.group(2))
    check_metric(bounds, "hbsp", "hbsp_s11_real", real)
    check_metric(bounds, "hbsp", "hbsp_s11_imag", imag)


def validate_psssp(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| (\S+) (\S+)", output)
    if not match:
        raise AssertionError("could not find PSSSP S11 row")
    real = float(match.group(1))
    imag = float(match.group(2))
    check_metric(bounds, "psssp", "psssp_s11_real", real)
    check_metric(bounds, "psssp", "psssp_s11_imag", imag)


def validate_hbstb(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| Mag: (\S+) Phase: (\S+)", output)
    if not match:
        raise AssertionError("could not find HBSTB row")
    mag = float(match.group(1))
    phase = float(match.group(2))
    check_metric(bounds, "hbstb", "hbstb_mag", mag)
    check_metric(bounds, "hbstb", "hbstb_phase_deg", phase)


def validate_pssstb(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| Mag: (\S+) Phase: (\S+)", output)
    if not match:
        raise AssertionError("could not find PSSSTB row")
    mag = float(match.group(1))
    phase = float(match.group(2))
    check_metric(bounds, "pssstb", "pssstb_mag", mag)
    check_metric(bounds, "pssstb", "pssstb_phase_deg", phase)


def validate_psspac_linear(output, bounds):
    match = re.search(r"1\.000000000e\+03 \| V\(in\)=\S+ V\(out\)=(\S+)", output)
    if not match:
        raise AssertionError("could not find PSSPAC linear V(out) row")
    gain = float(match.group(1))
    check_metric(bounds, "psspac_linear", "psspac_linear_vout_1khz", gain)


def validate_hb_psp103(output, bounds):
    rows = parse_hb_table(output)
    fundamental = rows[1][2][0]
    second = rows[2][2][0]
    check_metric(bounds, "hb_psp103", "hb_psp103_drain_fundamental", fundamental)
    check_metric(bounds, "hb_psp103", "hb_psp103_drain_second_harmonic", second)


VALIDATORS = {
    "hb_linear_rc": validate_hb_linear_rc,
    "hb_diode": validate_hb_diode,
    "hb_multitone": validate_hb_multitone,
    "hbac_diode": validate_hbac_diode,
    "hbnoise_diode": validate_hbnoise_diode,
    "pnoise_phase_jitter": validate_pnoise_phase_jitter,
    "hbsp": validate_hbsp,
    "psssp": validate_psssp,
    "hbstb": validate_hbstb,
    "pssstb": validate_pssstb,
    "psspac_linear": validate_psspac_linear,
    "hb_psp103": validate_hb_psp103,
}

REQUIRES_NATIVE_HB = {
    "hb_linear_rc",
    "hb_diode",
    "hb_multitone",
    "hbac_diode",
    "hbnoise_diode",
    "hbsp",
    "hbstb",
    "hb_psp103",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--deck", required=True)
    parser.add_argument("--case", required=True, choices=sorted(VALIDATORS))
    parser.add_argument("--golden")
    args = parser.parse_args()
    output = run_deck(args.exe, args.deck)
    if args.case in REQUIRES_NATIVE_HB and "Native HB Converged" not in output:
        print(output)
        raise AssertionError("native HB did not converge")
    VALIDATORS[args.case](output, load_bounds(args.golden))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"RF signoff validation failed: {exc}", file=sys.stderr)
        raise
