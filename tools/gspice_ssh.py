#!/usr/bin/env python3
"""
GSPICE Remote Simulation using SSH + SCP

Usage:
  python gspice_ssh.py input.sp \\
      --host <hostname> --user <username> --remote-gspice <path> \\
      [--output result.raw] [--key ssh_key_file]
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile


def run_cmd(cmd, desc, print_on_error=True):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if r.returncode != 0 and print_on_error:
            if r.stdout: print(r.stdout, file=sys.stdout)
            if r.stderr: print(r.stderr, file=sys.stderr)
        return r
    except FileNotFoundError:
        print(f"ERROR: '{cmd[0]}' not found. Is SSH installed?", file=sys.stderr)
        sys.exit(5)
    except subprocess.TimeoutExpired:
        print("ERROR: command timed out", file=sys.stderr)
        sys.exit(5)


def create_remote_dir(ssh_prefix, host, remote_dir):
    win_cmd = f"powershell -Command \"if (!(Test-Path '{remote_dir}')) {{ New-Item -ItemType Directory -Force -Path '{remote_dir}' | Out-Null }}\""
    r = run_cmd(ssh_prefix + [host, win_cmd], "create remote dir (win)", print_on_error=False)
    if r.returncode == 0:
        return True
    r = run_cmd(ssh_prefix + [host, f"mkdir -p {remote_dir}"], "create remote dir (nix)", print_on_error=False)
    return r.returncode == 0


def check_remote_file_exists(ssh_prefix, host, remote_file):
    win_cmd = f"powershell -Command \"if (Test-Path '{remote_file}') {{ 'EXISTS' }} else {{ 'MISSING' }}\""
    r = run_cmd(ssh_prefix + [host, win_cmd], "check output (win)", print_on_error=False)
    if "EXISTS" in r.stdout:
        return True
    r = run_cmd(ssh_prefix + [host, f"test -f {remote_file} && echo EXISTS || echo MISSING"], "check output (nix)", print_on_error=False)
    return "EXISTS" in r.stdout


def cleanup_remote_dir(ssh_prefix, host, remote_dir):
    win_cmd = f"powershell -Command \"Remove-Item -Recurse -Force '{remote_dir}' -ErrorAction SilentlyContinue\""
    r = run_cmd(ssh_prefix + [host, win_cmd], "cleanup (win)", print_on_error=False)
    if r.returncode != 0:
        run_cmd(ssh_prefix + [host, f"rm -rf {remote_dir}"], "cleanup (nix)", print_on_error=False)


def process_and_deploy_netlist_dependencies(input_sp, host, scp_prefix, ssh_prefix, remote_dir):
    """
    Parse input.sp for .lib, .include, .inc directives referencing local files.
    Upload referenced model files/directories to Machine B and update netlist paths.
    Returns (path_to_upload, is_temporary_file).
    """
    try:
        with open(input_sp, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except Exception as exc:
        print(f"WARNING: Could not read netlist for dependency scanning: {exc}", file=sys.stderr)
        return input_sp, False

    pattern = re.compile(r'^\s*\.(?:lib|include|inc)\s+["\']?([^"\'\s\n]+)["\']?', re.MULTILINE | re.IGNORECASE)
    matches = pattern.findall(content)
    if not matches:
        return input_sp, False

    uploaded_dirs = set()
    path_replacements = {}

    remote_models_dir = f"{remote_dir}/models"
    create_remote_dir(ssh_prefix, host, remote_models_dir)

    for ref_path in matches:
        clean_path = ref_path.strip('"\'')
        if os.path.isfile(clean_path):
            parent_dir = os.path.abspath(os.path.dirname(clean_path))
            if parent_dir not in uploaded_dirs:
                uploaded_dirs.add(parent_dir)
                print(f"Deploying PDK model directory '{parent_dir}' to Machine B ...")
                files_to_upload = [
                    os.path.join(parent_dir, fn)
                    for fn in os.listdir(parent_dir)
                    if os.path.isfile(os.path.join(parent_dir, fn))
                ]
                if files_to_upload:
                    scp_cmd = scp_prefix + files_to_upload + [f"{host}:{remote_models_dir}/"]
                    run_cmd(scp_cmd, "upload PDK model files to Machine B", print_on_error=False)

            filename = os.path.basename(clean_path)
            new_remote_path = f"{remote_models_dir}/{filename}"
            path_replacements[clean_path] = new_remote_path
            path_replacements[clean_path.replace("\\", "/")] = new_remote_path
            path_replacements[clean_path.replace("/", "\\")] = new_remote_path

    if path_replacements:
        new_content = content
        for old_p, new_p in path_replacements.items():
            new_content = new_content.replace(old_p, new_p)

        if new_content != content:
            temp_sp = tempfile.NamedTemporaryFile("w", delete=False, suffix=".sp", encoding="utf-8")
            temp_sp.write(new_content)
            temp_sp.close()
            return temp_sp.name, True

    return input_sp, False


def deploy_osdi_libraries(local_binary, host, scp_prefix, ssh_prefix, remote_dir):
    """
    Find local .osdi Verilog-A model files and deploy them to Machine B.
    """
    osdi_search_dirs = [
        r"C:\EDA\GSPICE\osdi",
        r"C:\EDA\Tools\vacask_0.3.4.rc1\vacask_0.3.4.rc1_windows-x86_64\lib\mod",
    ]
    if local_binary:
        bin_dir = os.path.dirname(os.path.abspath(local_binary))
        osdi_search_dirs.insert(0, os.path.join(bin_dir, "..", "osdi"))
        osdi_search_dirs.insert(0, os.path.join(bin_dir, "osdi"))

    osdi_files = []
    seen_names = set()
    for d in osdi_search_dirs:
        if os.path.isdir(d):
            for fn in os.listdir(d):
                if fn.endswith(".osdi") and fn not in seen_names:
                    fp = os.path.join(d, fn)
                    if os.path.isfile(fp):
                        osdi_files.append(fp)
                        seen_names.add(fn)

    if not osdi_files:
        return

    print(f"Deploying {len(osdi_files)} OSDI Verilog-A model libraries to Machine B ...")
    remote_osdi_dir = f"{remote_dir}/osdi"
    remote_models_osdi_dir = f"{remote_dir}/models/osdi"
    remote_models_dir = f"{remote_dir}/models"

    create_remote_dir(ssh_prefix, host, remote_osdi_dir)
    create_remote_dir(ssh_prefix, host, remote_models_osdi_dir)

    run_cmd(scp_prefix + osdi_files + [f"{host}:{remote_osdi_dir}/"], "deploy OSDI models (osdi)", print_on_error=False)
    run_cmd(scp_prefix + osdi_files + [f"{host}:{remote_models_osdi_dir}/"], "deploy OSDI models (models/osdi)", print_on_error=False)
    run_cmd(scp_prefix + osdi_files + [f"{host}:{remote_models_dir}/"], "deploy OSDI models (models)", print_on_error=False)
    run_cmd(scp_prefix + osdi_files + [f"{host}:{remote_dir}/"], "deploy OSDI models (root)", print_on_error=False)


def deploy_binary_runtime_dependencies(local_binary, host, scp_prefix, remote_dir):
    """
    Upload DLLs next to the selected local binary. This is required for the
    SuiteSparse/KLU vcpkg build, where klu.dll and its dependencies sit beside
    gspice.exe after CMake's dependency deployment step.
    """
    if not local_binary:
        return
    bin_dir = os.path.dirname(os.path.abspath(local_binary))
    if not os.path.isdir(bin_dir):
        return
    dlls = [
        os.path.join(bin_dir, fn)
        for fn in os.listdir(bin_dir)
        if fn.lower().endswith(".dll") and os.path.isfile(os.path.join(bin_dir, fn))
    ]
    if not dlls:
        return
    print(f"Deploying {len(dlls)} runtime DLL(s) beside GSPICE on Machine B ...")
    run_cmd(scp_prefix + dlls + [f"{host}:{remote_dir}/"], "deploy binary runtime DLLs", print_on_error=False)


def run_remote_simulation_with_periodic_sync(ssh_prefix, scp_prefix, host, cmd, remote_raw, local_out_path, sync_pct_step=2.0):
    """
    Run remote simulation on Machine B, stream stdout/stderr live,
    and periodically dump partial .raw results to Machine A every 2% of progress or 4s.
    """
    import threading
    import time

    full_cmd = ssh_prefix + [host] + cmd
    print(f"Running simulation on Machine B: {' '.join(cmd)}")

    try:
        proc = subprocess.Popen(full_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
    except FileNotFoundError:
        print(f"ERROR: '{ssh_prefix[0]}' not found. Is SSH installed?", file=sys.stderr)
        sys.exit(5)

    last_synced_pct = 0.0
    last_sync_time = time.time()

    def sync_partial_results():
        nonlocal last_sync_time
        if not local_out_path or not remote_raw:
            return
        r = subprocess.run(scp_prefix + [f"{host}:{remote_raw}", local_out_path], capture_output=True, text=True, timeout=30)
        if r.returncode == 0:
            last_sync_time = time.time()
            print(f"[Periodic Sync] Synced partial raw results to Machine A: '{local_out_path}'")

    def bg_sync_loop():
        while proc.poll() is None:
            time.sleep(4.0)
            if proc.poll() is None and (time.time() - last_sync_time) >= 4.0:
                sync_partial_results()

    sync_thread = threading.Thread(target=bg_sync_loop, daemon=True)
    sync_thread.start()

    while True:
        line = proc.stdout.readline()
        if not line and proc.poll() is not None:
            break
        if line:
            sys.stdout.write(line)
            sys.stdout.flush()
            m = re.search(r"(\d+(?:\.\d+)?)%", line)
            if m:
                current_pct = float(m.group(1))
                if (current_pct - last_synced_pct) >= sync_pct_step:
                    last_synced_pct = current_pct
                    sync_partial_results()

    stderr_out = proc.stderr.read()
    if stderr_out:
        sys.stderr.write(stderr_out)
        sys.stderr.flush()

    if local_out_path and proc.returncode == 0:
        sync_partial_results()

    return proc.returncode


def main():
    parser = argparse.ArgumentParser(description="GSPICE remote simulation via SSH")
    parser.add_argument("input", help="Local .sp netlist file")
    parser.add_argument("--output", default="", help="Local output .raw file")
    parser.add_argument("--host", required=True, help="Remote hostname or IP")
    parser.add_argument("--user", required=True, help="SSH username")
    parser.add_argument("--key", default="", help="SSH private key path")
    parser.add_argument("--remote-gspice", default="gspice",
                        help="Path to gspice executable on remote (default: gspice)")
    parser.add_argument("--save", choices=["ALL", "SELECTED", "NONE", "all", "selected", "none"],
                        default="", help="Forward waveform save mode to remote gspice")
    parser.add_argument("--adaptive-maxstep", action="store_true",
                        help="Forward --adaptive-maxstep to remote gspice")
    parser.add_argument("--local-binary", default="", help="Local GSPICE binary path to deploy to remote host")
    parser.add_argument("--deploy-binary", action="store_true", help="Upload local GSPICE binary to remote host")
    parser.add_argument("--keep-remote", action="store_true", default=True, help="Retain results on remote machine B")
    parser.add_argument("--no-keep-remote", action="store_false", dest="keep_remote", help="Clean up remote directory after download")
    parser.add_argument("--remote-dir", default="", help="Custom working directory path on remote machine B")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: input file '{args.input}' not found", file=sys.stderr)
        sys.exit(2)

    host = f"{args.user}@{args.host}"
    scp_prefix = ["scp"]
    ssh_prefix = ["ssh"]
    if args.key:
        scp_prefix += ["-i", args.key]
        ssh_prefix += ["-i", args.key]

    if args.remote_dir:
        remote_dir = args.remote_dir.replace("\\", "/")
    else:
        remote_dir = f"/tmp/gspice_ssh_{os.getpid()}"

    remote_sp = f"{remote_dir}/{os.path.basename(args.input)}"
    remote_raw = remote_sp + ".raw"

    # 1) Create remote temp/work dir
    print(f"Connecting to {host} ...")
    if not create_remote_dir(ssh_prefix, host, remote_dir):
        print("ERROR: could not create remote working directory on Machine B", file=sys.stderr)
        sys.exit(3)

    # 2) Optionally deploy local GSPICE binary to Machine B (No pre-installation needed on Machine B!)
    remote_gspice_bin = args.remote_gspice
    if args.deploy_binary or args.local_binary:
        bin_path = args.local_binary if args.local_binary else r"C:\EDA\GSPICE\build\Release\gspice.exe"
        if os.path.isfile(bin_path):
            bin_name = os.path.basename(bin_path)
            remote_gspice_bin = f"{remote_dir}/{bin_name}"
            print(f"Deploying local GSPICE binary '{bin_path}' to Machine B ({remote_gspice_bin}) ...")
            r = run_cmd(scp_prefix + [bin_path, f"{host}:{remote_gspice_bin}"], "deploy binary to Machine B")
            if r.returncode == 0:
                deploy_binary_runtime_dependencies(bin_path, host, scp_prefix, remote_dir)
                run_cmd(ssh_prefix + [host, f"chmod +x {remote_gspice_bin}"], "chmod remote binary", print_on_error=False)
            else:
                print("WARNING: Could not upload local GSPICE binary; falling back to remote_gspice path", file=sys.stderr)
                remote_gspice_bin = args.remote_gspice

    # Deploy OSDI model libraries (psp103va.osdi, etc.) to Machine B
    deploy_osdi_libraries(args.local_binary, host, scp_prefix, ssh_prefix, remote_dir)

    # 3) Process netlist dependencies (PDK models) and upload netlist
    sp_to_upload, is_temp_sp = process_and_deploy_netlist_dependencies(args.input, host, scp_prefix, ssh_prefix, remote_dir)
    print(f"Uploading '{args.input}' to Machine B ...")
    r = run_cmd(scp_prefix + [sp_to_upload, f"{host}:{remote_sp}"], "upload netlist")
    if is_temp_sp and os.path.isfile(sp_to_upload):
        try: os.unlink(sp_to_upload)
        except OSError: pass
    if r.returncode != 0:
        if not args.keep_remote:
            cleanup_remote_dir(ssh_prefix, host, remote_dir)
        print("ERROR: upload failed", file=sys.stderr)
        sys.exit(3)

    # 4) Run simulation on remote Machine B with periodic 2% data dumping to Machine A
    cmd = [remote_gspice_bin, remote_sp]
    if args.output:
        cmd += ["--output", remote_raw]
    if args.save:
        cmd += ["--save", args.save.lower()]
    if args.adaptive_maxstep:
        cmd += ["--adaptive-maxstep"]
    out_path = args.output if args.output else os.path.basename(args.input) + ".raw"

    returncode = run_remote_simulation_with_periodic_sync(
        ssh_prefix, scp_prefix, host, cmd, remote_raw, out_path, sync_pct_step=2.0
    )
    if returncode != 0:
        if not args.keep_remote:
            cleanup_remote_dir(ssh_prefix, host, remote_dir)
        print(f"ERROR: remote simulation exited with code {returncode}", file=sys.stderr)
        sys.exit(returncode)

    # 5) Check if .raw was produced on Machine B
    if not check_remote_file_exists(ssh_prefix, host, remote_raw):
        if not args.keep_remote:
            cleanup_remote_dir(ssh_prefix, host, remote_dir)
        print("ERROR: remote simulation did not produce a .raw file", file=sys.stderr)
        sys.exit(4)

    # 6) Download result copy back to Machine A for GUI viewer
    out_path = args.output if args.output else os.path.basename(args.input) + ".raw"
    print(f"Downloading result to local Machine A: '{out_path}' ...")
    r = run_cmd(scp_prefix + [f"{host}:{remote_raw}", out_path], "download result")
    if r.returncode != 0:
        if not args.keep_remote:
            cleanup_remote_dir(ssh_prefix, host, remote_dir)
        print("ERROR: download failed", file=sys.stderr)
        sys.exit(3)

    # 7) Retain or cleanup remote Machine B results
    if args.keep_remote:
        print(f"Results retained on Machine B at: '{remote_dir}'")
    else:
        cleanup_remote_dir(ssh_prefix, host, remote_dir)

    print(f"Done. Local results: '{out_path}'")
    sys.exit(0)


if __name__ == "__main__":
    main()
