# VACASK Oracle Validation

GSPICE remains the production simulator.  VACASK is used only as an external
oracle for compact-model parity checks.

## Run Without VACASK

The harness can list the supported cases without a VACASK installation:

```powershell
python tools/vacask_parity.py --list-cases
```

## Produce A GSPICE CSV

This is useful when bringing up a VACASK converter/runner because it shows the
canonical column names expected by the comparator.

```powershell
python tools/vacask_parity.py `
  --gspice build-klu/Release/gspice.exe `
  --case ihp_lv_dc `
  --emit-gspice-csv build-klu/ihp_lv_dc_gspice.csv
```

## Compare Against A VACASK CSV

The VACASK CSV must have a header row.  Columns with matching names are compared
using the case's preferred signals, for example `V(drain)` for DC and
`V(drain).real` / `V(drain).imag` for AC.

```powershell
python tools/vacask_parity.py `
  --gspice build-klu/Release/gspice.exe `
  --case ihp_lv_dc `
  --vacask-csv path/to/ihp_lv_dc_vacask.csv
```

## Run VACASK Through A Command Template

Set `VACASK_COMMAND` to the local converter/runner.  The template receives:

- `{deck}`: the GSPICE test deck
- `{csv}`: the CSV file VACASK must write
- `{case}`: the validation case name
- `{workdir}`: a temporary work directory

Example:

```powershell
$env:VACASK_COMMAND = "python C:/tools/run_vacask_ihp.py --spice {deck} --csv {csv}"
ctest --test-dir build-klu -C Release -R validation_vacask --output-on-failure
```

Enable the optional CTest lane at configure time:

```powershell
cmake -S . -B build-klu -DGSPICE_ENABLE_VACASK_VALIDATION=ON
```
