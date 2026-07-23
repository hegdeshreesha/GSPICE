param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Deck,
    [Parameter(Mandatory = $true)][string]$Output
)

$ErrorActionPreference = "Stop"
if (Test-Path -LiteralPath $Output) {
    Remove-Item -LiteralPath $Output -Force
}

& $Exe --threads 1 --format csv --output $Output $Deck | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "GSPICE failed with exit code $LASTEXITCODE"
}
if (-not (Test-Path -LiteralPath $Output)) {
    throw "CSV output was not created"
}

$lines = Get-Content -LiteralPath $Output
if ($lines.Count -lt 3) {
    throw "CSV output has too few rows: $($lines.Count)"
}
if ($lines[0] -notmatch '^time,' -or $lines[0] -notmatch '"V\(out\)"') {
    throw "Unexpected CSV header: $($lines[0])"
}
$columns = $lines[1].Split(',')
if ($columns.Count -lt 2) {
    throw "CSV data row has too few columns: $($lines[1])"
}
$parsed = 0.0
if (-not [double]::TryParse($columns[0], [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    throw "CSV time value is not numeric: $($columns[0])"
}
