# GSPICE OSDI / PSP Bring-Up

GSPICE has first-pass OSDI netlist plumbing for OpenVAF-style compact models:

- `.OSDI` and `.PRE_OSDI` loader directives
- `.MODEL` card storage
- `N...` compact-model instance parsing
- automatic `.osdi` discovery through `GSPICE_OSDI_DIR`, `NGSPICE_OSDI_DIR`, the deck folder, and local `osdi` folders

The IHP SG13G2 PSP Verilog-A sources are expected under:

```text
C:\EDA\LumenCircuitStudio\external\ihp_pdk\ihp-sg13g2\libs.tech\verilog-a\psp103
```

Compile them with:

```powershell
.\tools\build_ihp_psp_osdi.ps1 -OpenVafPath C:\path\to\openvaf.exe
```

Then run GSPICE with:

```powershell
$env:GSPICE_OSDI_DIR = "C:\EDA\GSPICE\osdi"
.\build\Release\gspice.exe path\to\deck.sp
```

Current limitation: GSPICE still needs full OpenVAF ABI verification, model/instance parameter binding, and charge/capacitance Jacobian handling before PSP simulation can be called production accurate.
