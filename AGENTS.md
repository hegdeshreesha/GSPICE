# GSPICE Agent Guide

GSPICE is a C++17 academic-beta SPICE-like simulator. Correctness, diagnostics,
and honest capability reporting matter more than feature count.

## Commands

Run from an x64 Native Tools environment when using MSVC/NMake:

```powershell
cd C:\EDA\GSPICE-clean-accuracy-fixes
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Focused examples:

```powershell
ctest --test-dir build -C Release -R "smoke_tran|regression_rc_step" --output-on-failure
python tools\validate_transient.py --gspice build\Release\gspice.exe --source .
build\Release\gspice.exe --capabilities
```

Release/vcpkg preset:

```powershell
$env:VCPKG_ROOT = "C:\path\to\vcpkg"
cmake --preset windows-vcpkg-release
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
```

## Repo Map

- `src/core`: CLI, parser, and simulator implementation.
- `include`: public/internal headers.
- `tests`: C++ tests, deck regressions, and PowerShell harnesses.
- `tools`: validation, oracle, build, and GMC utilities.
- `docs`: architecture, limitations, roadmap, and release status.

## Development Rules

- Add or adjust a regression test for every behavior change.
- Unsupported active syntax, model behavior, or analyses must fail loudly; do
  not silently substitute operating point or primitive MOS behavior.
- Keep clean-room boundaries: do not copy source, pseudocode, comments, tests,
  internal names, or structure from incompatible simulator projects.
- If a capability becomes more or less mature, update `--capabilities` output
  and `docs/LIMITATIONS.md` or related docs.
- Prefer focused ctest filters while iterating; run the full suite before
  calling broad simulator changes done.

## Verification Ladder

1. Build the touched target.
2. Run the narrowest relevant `ctest -R` filter.
3. Run validation tools when numerical behavior changes.
4. Run full `ctest --test-dir build -C Release --output-on-failure` for parser,
   solver, model, capability, or shared infrastructure changes.

