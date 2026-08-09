#!/usr/bin/env python3
"""Regenerate the PDF guide: python tools/gen_pdf.py"""

from fpdf import FPDF
import os

PDF_PATH = os.path.join(os.path.dirname(__file__), "..", "GSPICE_SSH_GUIDE.pdf")

doc = FPDF()
doc.set_auto_page_break(auto=True, margin=20)
doc.add_page()

doc.set_font("Helvetica", "B", 20)
doc.cell(0, 12, "GSPICE Remote Simulation via SSH", align="C", new_x="LMARGIN", new_y="NEXT")
doc.ln(6)
doc.set_font("Helvetica", "", 10)
doc.cell(0, 6, "Run GSPICE on a remote machine; results land on your local machine.", align="C", new_x="LMARGIN", new_y="NEXT")
doc.ln(8)

def heading(n, text):
    sizes = {1: 16, 2: 14, 3: 11}
    doc.set_font("Helvetica", "B", sizes.get(n, 14))
    doc.cell(0, 10, text, new_x="LMARGIN", new_y="NEXT")
    doc.ln(1)

def body(text):
    doc.set_font("Helvetica", "", 10)
    doc.multi_cell(0, 5, text)
    doc.ln(3)

def code(text):
    doc.set_font("Courier", "", 9)
    doc.multi_cell(0, 4.5, text)
    doc.ln(3)

# ---- 1 ----
heading(1, "1. How It Works")
body(
"1. Your local gspice calls tools/gspice_ssh.py\n"
"2. The script SCPs the .sp netlist to the remote\n"
"3. It SSHes in and runs gspice (just the binary, no install)\n"
"4. It SCPs the .raw result file back to your local machine\n"
"5. Remote temp files are cleaned up automatically\n\n"
"No GSPICE installation needed on the remote. Copy the binary once and you are done."
)

# ---- 2 ----
heading(1, "2. Quick Reference")
code(
"gspice test.sp --sim-env ssh              \\\n"
"    --ssh-host 192.168.1.100              \\\n"
"    --ssh-user alice                      \\\n"
"    --remote-gspice /home/alice/gspice    \\\n"
"    -o results.raw"
)
body("Omit any --ssh option and gspice will prompt interactively.")

# ---- 3 ----
heading(1, "3. Install OpenSSH")
heading(3, "On Windows 10/11 (local or remote)")
body(
"Settings -> Apps -> Optional Features -> Add a feature\n"
"Search for 'OpenSSH Server' and 'OpenSSH Client', install both.\n\n"
"Or via PowerShell (as Administrator):")
code(
"Add-WindowsCapability -Online -Name OpenSSH.Server~~~~0.0.1.0\n"
"Add-WindowsCapability -Online -Name OpenSSH.Client~~~~0.0.1.0")
heading(3, "On Linux (remote)")
code(
"sudo apt update && sudo apt install openssh-server\n"
"sudo systemctl enable ssh && sudo systemctl start ssh")

# ---- 4 ----
heading(1, "4. Password-less SSH Keys (Recommended)")
heading(3, "Generate key pair")
code('ssh-keygen -t ed25519 -f "$HOME/.ssh/gspice_remote"')
heading(3, "Copy public key to remote")
body("Windows PowerShell:")
code(
'type $env:USERPROFILE\\.ssh\\gspice_remote.pub | '
'ssh user@remote-host "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"')
body("Linux:")
code("ssh-copy-id -i ~/.ssh/gspice_remote.pub user@remote-host")
heading(3, "Test the connection")
code("ssh -i ~/.ssh/gspice_remote user@remote-host")
body("Should log in without a password.")

# ---- 5 ----
heading(1, "5. Copy GSPICE Binary to Remote")
code(
'scp -i ~/.ssh/gspice_remote C:\\EDA\\GSPICE\\build\\Release\\gspice.exe '
'user@remote-host:/home/alice/gspice')
body("If remote is Linux, make executable:")
code('ssh -i ~/.ssh/gspice_remote user@remote-host "chmod +x /home/alice/gspice"')

# ---- 6 ----
heading(1, "6. Firewall (Windows Remote Only)")
body("On the remote Windows machine, run PowerShell as Administrator:")
code(
"New-NetFirewallRule -DisplayName 'OpenSSH Server' "
"-Direction Inbound -Protocol TCP -LocalPort 22 -Action Allow")

# ---- 7 ----
heading(1, "7. Running a Remote Simulation")
heading(3, "All options on command line")
code(
"gspice test.sp --sim-env ssh                \\\n"
"    --ssh-host 192.168.1.100                \\\n"
"    --ssh-user alice                        \\\n"
'    --ssh-key \"%USERPROFILE%\\.ssh\\gspice_remote\" \\\n'
"    --remote-gspice /home/alice/gspice      \\\n"
"    -o results.raw")
heading(3, "Interactive prompt")
code(
"gspice test.sp --sim-env ssh\n"
"Remote host: 192.168.1.100\n"
"SSH user: alice\n"
"Remote gspice path: /home/alice/gspice")

# ---- 8 ----
doc.add_page()
heading(1, "8. Troubleshooting")

problems = [
("'ssh' is not recognized",
"OpenSSH not installed.\nFix: Settings -> Apps -> Optional Features -> install "
"'OpenSSH Client'. Or install Git for Windows / PuTTY."),

("Connection refused",
"SSH server not running or firewall blocking port 22.\n"
"- Remote: sudo systemctl status ssh\n"
"- Windows remote: services.msc -> start 'OpenSSH SSH Server'\n"
"- Firewall: see Section 6\n"
"- Test: ssh -v user@remote-host"),

("Permission denied (publickey)",
"SSH key misconfigured.\n"
"- Verify public key is in ~/.ssh/authorized_keys on remote\n"
"- Verify --ssh-key path is correct\n"
"- Remote: chmod 600 ~/.ssh/authorized_keys\n"
"- Test: ssh -i ~/.ssh/gspice_remote user@remote-host"),

("scp: file not found / upload failed",
"Source .sp missing or remote disk/temp issues.\n"
"- Check input .sp file exists\n"
"- Remote: df -h (disk space)\n"
"- Remote: touch /tmp/test_write (write permission)"),

("Remote gspice not found",
"--remote-gspice path is wrong.\n"
"- ssh user@remote-host ls -l /path/to/gspice\n"
"- Linux: chmod +x /path/to/gspice"),

("No .raw file returned",
"Simulation ran but no --output was given, or output path differs.\n"
"- Always use -o output.raw with --sim-env ssh\n"
"- Default: <input>.sp.raw"),

("SSH hangs before connecting",
"DNS or IPv6 timeout.\n"
"- Use IP address instead of hostname\n"
"- Add 'AddressFamily inet' to ~/.ssh/config\n"
"- ssh -4 user@host"),

("Windows: 'gspice' not found",
"gspice not in PATH.\n"
"Use full path or add to PATH:\n"
"$env:PATH += \";C:\\EDA\\GSPICE\\build\\Release\""),
]

for title, body_text in problems:
    heading(3, title)
    body(body_text)

# ---- 9 ----
heading(1, "9. Free Tools for Windows")
tools = [
("OpenSSH (built-in)", "ssh, scp, ssh-keygen -- comes with Windows 10/11"),
("PuTTY", "putty.exe, pscp.exe, pageant.exe -- https://www.putty.org"),
("WinSCP", "Graphical SCP/SFTP -- https://winscp.net"),
("Git for Windows", "Bundles OpenSSH + bash -- https://git-scm.com"),
]
for name, desc in tools:
    heading(3, f"- {name}")
    body(desc)

# ---- 10 ----
heading(1, "10. Files Changed")
code(
"tools/gspice_ssh.py    - SSH orchestration script (new)\n"
"src/core/main.cpp      - Added --sim-env, --ssh-host, --ssh-user,\n"
"                         --ssh-key, --remote-gspice options")

doc.ln(6)
doc.set_font("Helvetica", "I", 9)
doc.cell(0, 5, "Generated 2026-07-27", align="C")

doc.output(PDF_PATH)
print(f"PDF: {PDF_PATH}")
