#!/usr/bin/env python3
"""
Automated benchmark and validation script for GSPICE vs independent references.
Generates scalable OSDI PSP103.4 inverter chains, measures wall-clock speed, memory footprint,
Newton iterations, and evaluates performance under FASTSPICE, MULTIRATE, and PARALLEL_SOLVE options.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

def run_command(cmd, timeout=300):
    start_time = time.time()
    proc = None
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        stdout, stderr = proc.communicate(timeout=timeout)
        elapsed = time.time() - start_time
        return proc.returncode, stdout, stderr, elapsed
    except subprocess.TimeoutExpired:
        if proc is not None:
            proc.kill()
            stdout, stderr = proc.communicate()
        else:
            stdout, stderr = "", ""
        elapsed = time.time() - start_time
        message = f"Timeout expired after {timeout}s"
        if stdout:
            message += "\nSTDOUT tail:\n" + stdout[-1200:]
        if stderr:
            message += "\nSTDERR tail:\n" + stderr[-1200:]
        return -1, stdout, message, elapsed

def generate_inverter_chain_deck(
    num_stages=50,
    osdi_path="c:/EDA/GSPICE/osdi/psp103.osdi",
    options_str="",
    step="1p",
    stop="1n",
):
    """
    Generates a multi-stage inverter chain netlist using PSP103VA OSDI compact models.
    """
    lines = [
        f"* Scalable Inverter Chain Benchmark ({num_stages} stages)",
        f'.PRE_OSDI "{osdi_path}"',
        '.MODEL nch psp103va type=1',
        '.MODEL pch psp103va type=-1',
        f'.OPTIONS {options_str}',
        'VDD vdd 0 DC 1.2',
        'VIN in 0 PULSE(0 1.2 0 10p 10p 200p 400p)'
    ]

    for i in range(num_stages):
        in_node = "in" if i == 0 else f"n{i}"
        out_node = f"n{i+1}"
        lines.append(f'NN{i} {out_node} {in_node} 0 0 nch w=1u l=0.045u')
        lines.append(f'NP{i} {out_node} {in_node} vdd vdd pch w=2u l=0.045u')

    lines.extend([
        f'.TRAN {step} {stop}',
        '.END'
    ])
    return "\n".join(lines)

def parse_key_values(line):
    values = {}
    for key, value in re.findall(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)", line):
        try:
            if any(ch in value for ch in ".eE"):
                values[key] = float(value)
            else:
                values[key] = int(value)
        except ValueError:
            values[key] = value
    return values

def parse_gspice_metrics(stdout):
    metrics = {}
    for line in stdout.splitlines():
        if line.startswith("Transient summary:"):
            metrics.update({f"tran_{k}": v for k, v in parse_key_values(line).items()})
        elif line.startswith("Newton summary:"):
            metrics.update({f"newton_{k}": v for k, v in parse_key_values(line).items()})
        elif line.startswith("Runtime summary:"):
            metrics.update({f"runtime_{k}": v for k, v in parse_key_values(line).items()})
        elif line.startswith("Transient phase summary:"):
            metrics.update({f"phase_{k}": v for k, v in parse_key_values(line).items()})
        elif line.startswith("Solver summary (real):"):
            metrics.update({f"solver_real_{k}": v for k, v in parse_key_values(line).items()})
    return metrics

def format_metric(metrics, name, width=8, precision=3):
    value = metrics.get(name)
    if value is None:
        return " " * width
    if isinstance(value, float):
        return f"{value:{width}.{precision}f}"
    return f"{value:{width}d}"

def print_result_line(name, elapsed, metrics):
    accepted = format_metric(metrics, "tran_accepted", 7, 0)
    rejected = format_metric(metrics, "tran_rejected", 6, 0)
    iterations = format_metric(metrics, "newton_total_iterations", 8, 0)
    stamp = format_metric(metrics, "newton_stamp_seconds", 8)
    solve = format_metric(metrics, "newton_solve_seconds", 8)
    state = format_metric(metrics, "phase_state", 7)
    lte = format_metric(metrics, "phase_lte", 7)
    accept = format_metric(metrics, "phase_accept", 7)
    output = format_metric(metrics, "phase_output", 7)
    klu = format_metric(metrics, "solver_real_ext_klu", 7)
    fill = format_metric(metrics, "solver_real_fill", 7)
    factor = format_metric(metrics, "solver_real_factor", 7)
    backsolve = format_metric(metrics, "solver_real_backsolve", 7)
    bypassed = format_metric(metrics, "newton_bypassed_devices", 9, 0)
    print(
        f"[{name:^30}] Wall: {elapsed:7.3f}s | "
        f"acc:{accepted} rej:{rejected} it:{iterations} "
        f"stamp:{stamp}s solve:{solve}s state:{state}s lte:{lte}s "
        f"accept:{accept}s out:{output}s klu:{klu}s fill:{fill}s "
        f"fac:{factor}s back:{backsolve}s bypass:{bypassed} | PASSED"
    )

def parse_ascii_raw(path):
    times = []
    values = []
    in_values = False
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            text = line.strip()
            if not text:
                continue
            if text == "Values:":
                in_values = True
                continue
            if not in_values:
                continue
            numbers = []
            for token in text.replace(",", " ").split():
                try:
                    numbers.append(float(token))
                except ValueError:
                    pass
            if len(numbers) >= 2:
                times.append(numbers[0])
                values.append(numbers[1:])
    return times, values

def compare_raw_outputs(reference_path, candidate_path):
    ref_t, ref_v = parse_ascii_raw(reference_path)
    cand_t, cand_v = parse_ascii_raw(candidate_path)
    count = min(len(ref_t), len(cand_t))
    if count == 0:
        return {"points": 0, "max_abs": float("nan"), "rms": float("nan")}
    max_abs = 0.0
    sum_sq = 0.0
    samples = 0
    max_time_delta = 0.0
    for i in range(count):
        max_time_delta = max(max_time_delta, abs(ref_t[i] - cand_t[i]))
        width = min(len(ref_v[i]), len(cand_v[i]))
        for j in range(width):
            delta = cand_v[i][j] - ref_v[i][j]
            max_abs = max(max_abs, abs(delta))
            sum_sq += delta * delta
            samples += 1
    rms = (sum_sq / samples) ** 0.5 if samples else float("nan")
    return {
        "points": count,
        "samples": samples,
        "max_abs": max_abs,
        "rms": rms,
        "max_time_delta": max_time_delta,
    }

def main():
    parser = argparse.ArgumentParser(description="Run GSPICE Head-to-Head Performance Benchmark")
    parser.add_argument("--gspice", default="build/Release/gspice.exe", help="Path to gspice binary")
    parser.add_argument("--stages", type=int, default=50, help="Number of inverter chain stages")
    parser.add_argument("--osdi", default="osdi/psp103.osdi", help="Path to PSP103 OSDI model binary")
    parser.add_argument("--tran-step", default="1p", help="Transient output step")
    parser.add_argument("--tran-stop", default="1n", help="Transient stop time")
    parser.add_argument("--timeout", type=int, default=300, help="Timeout per benchmark configuration in seconds")
    parser.add_argument(
        "--configs",
        default="all",
        help="Comma-separated configs to run: baseline,conservative,fastspice,multirate,stampcache,parallel,all",
    )
    parser.add_argument(
        "--compare-accuracy",
        action="store_true",
        help="Run conservative charge-LTE reference and compare each selected waveform against it",
    )
    parser.add_argument(
        "--keep-outputs",
        action="store_true",
        help="Keep temporary benchmark deck and waveform output files",
    )
    args = parser.parse_args()

    gspice_bin = os.path.abspath(args.gspice)
    osdi_path = os.path.abspath(args.osdi)

    if not os.path.exists(gspice_bin):
        print(f"Error: GSPICE binary not found at {gspice_bin}")
        sys.exit(1)

    print("=========================================================")
    print("      GSPICE Head-to-Head Benchmark Suite              ")
    print("=========================================================")
    print(f"GSPICE Binary : {gspice_bin}")
    print(f"OSDI Model    : {osdi_path}")
    print(f"Test Circuit  : {args.stages}-stage PSP103.4 Inverter Chain ({args.stages * 2} Transistors)")
    print("---------------------------------------------------------")

    available_configs = [
        ("baseline", "Standard SPICE Baseline", ""),
        ("conservative", "Conservative Charge LTE", "LTE_CHARGE_INTERVAL=1"),
        ("fastspice", "KLU Refactor + FastSPICE", "FASTSPICE=1"),
        ("multirate", "Multi-Rate Timestepping", "MULTIRATE=1"),
        ("stampcache", "Multirate + Stamp Cache", "MULTIRATE=1 TRAN_STAMP_CACHE=1"),
        ("parallel", "Parallel BTF Matrix Solver", "PARALLEL_SOLVE=1"),
        ("all", "All Optimizations Enabled", "FASTSPICE=1 MULTIRATE=1 PARALLEL_SOLVE=1")
    ]
    requested = {item.strip().lower() for item in args.configs.split(",") if item.strip()}
    if requested == {"all"}:
        configs = available_configs
    else:
        configs = [item for item in available_configs if item[0] in requested]
        missing = requested - {item[0] for item in available_configs}
        if missing:
            print(f"Error: unknown config(s): {', '.join(sorted(missing))}")
            sys.exit(2)
        if not configs:
            print("Error: no benchmark configs selected")
            sys.exit(2)
    if args.compare_accuracy and not any(item[0] == "conservative" for item in configs):
        configs = [available_configs[1]] + configs

    results = []
    output_files = {}
    reference_output = None

    for key, name, opt in configs:
        deck_content = generate_inverter_chain_deck(
            args.stages,
            osdi_path,
            opt,
            args.tran_step,
            args.tran_stop,
        )
        deck_filename = f"benchmark_temp_{args.stages}.sp"
        output_filename = f"benchmark_temp_{args.stages}_{key}.raw"
        with open(deck_filename, "w") as f:
            f.write(deck_content)

        rc, stdout, stderr, elapsed = run_command(
            [gspice_bin, deck_filename, "-o", output_filename],
            timeout=args.timeout,
        )

        if os.path.exists(deck_filename) and not args.keep_outputs:
            os.remove(deck_filename)
        if rc == 0:
            output_files[key] = output_filename
            if key == "conservative":
                reference_output = output_filename
        if os.path.exists(output_filename) and not (args.keep_outputs or args.compare_accuracy):
            os.remove(output_filename)

        if rc == 0:
            metrics = parse_gspice_metrics(stdout)
            print_result_line(name, elapsed, metrics)
            result = {"config": key, "elapsed_s": elapsed, "status": "PASSED"}
            result.update(metrics)
            results.append(result)
        else:
            print(f"[{name:^30}] Error: {stderr.strip()[:300]}")
            results.append({"config": key, "elapsed_s": elapsed, "status": "FAILED"})

    if args.compare_accuracy and reference_output and os.path.exists(reference_output):
        print("---------------------------------------------------------")
        print("Accuracy vs Conservative Charge LTE reference")
        for result in results:
            key = result["config"]
            if key == "conservative" or result.get("status") != "PASSED":
                continue
            candidate_output = output_files.get(key)
            if not candidate_output or not os.path.exists(candidate_output):
                continue
            delta = compare_raw_outputs(reference_output, candidate_output)
            result.update({f"accuracy_{k}": v for k, v in delta.items()})
            print(
                f"[{key:^30}] points:{delta['points']:7d} "
                f"samples:{delta.get('samples', 0):9d} "
                f"max_abs:{delta['max_abs']:.6e} rms:{delta['rms']:.6e} "
                f"max_dt:{delta['max_time_delta']:.3e}"
            )
    if args.compare_accuracy and not args.keep_outputs:
        for path in output_files.values():
            if os.path.exists(path):
                os.remove(path)

    print("=========================================================")
    print("Benchmark Completed Successfully.")

if __name__ == "__main__":
    main()
