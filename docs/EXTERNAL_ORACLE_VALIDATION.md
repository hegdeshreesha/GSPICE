# External Oracle Validation

GSPICE remains the production simulator. The external oracle is used only as an external
oracle for compact-model parity checks.

## Run Without An External Oracle

The harness can list the supported cases without an external simulator installation:

```powershell
python tools/external_oracle_parity.py --list-cases
```

## Produce A GSPICE CSV

This is useful when bringing up a converter/runner because it shows the
canonical column names expected by the comparator.

```powershell
python tools/external_oracle_parity.py `
  --gspice build-klu/Release/gspice.exe `
  --case ihp_lv_dc `
  --emit-gspice-csv build-klu/ihp_lv_dc_gspice.csv
```

## Compare Against An External CSV

The external CSV must have a header row. Columns with matching names are compared
using the case's preferred signals, for example `V(drain)` for DC and
`V(drain).real` / `V(drain).imag` for AC.

```powershell
python tools/external_oracle_parity.py `
  --gspice build-klu/Release/gspice.exe `
  --case ihp_lv_dc `
  --external_oracle-csv path/to/ihp_lv_dc_external_oracle.csv
```

## Run An External Oracle Through A Command Template

Set `EXTERNAL_ORACLE_COMMAND` to the local converter/runner. The template receives:

- `{deck}`: the GSPICE test deck
- `{csv}`: the CSV file the external tool must write
- `{case}`: the validation case name
- `{workdir}`: a temporary work directory

Example:

```powershell
$env:EXTERNAL_ORACLE_COMMAND = "python C:/tools/run_external_oracle_ihp.py --spice {deck} --csv {csv}"
ctest --test-dir build-klu -C Release -R validation_external_oracle --output-on-failure
```

Enable the optional CTest lane at configure time:

```powershell
cmake -S . -B build-klu -DGSPICE_ENABLE_EXTERNAL_ORACLE_VALIDATION=ON
```
