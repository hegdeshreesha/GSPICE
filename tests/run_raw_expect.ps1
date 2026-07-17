param(
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Deck,
    [Parameter(Mandatory=$true)][string]$Signal,
    [Parameter(Mandatory=$true)][double]$Expected,
    [Parameter(Mandatory=$true)][double]$Tolerance
)

$rawPath = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_raw_expect_{0}.raw" -f ([System.Guid]::NewGuid().ToString("N")))
try {
    & $Exe --threads 1 -o $rawPath $Deck
    if ($LASTEXITCODE -ne 0) {
        throw "GSPICE exited with code $LASTEXITCODE"
    }
    if (-not (Test-Path $rawPath)) {
        throw "Expected RAW output was not created: $rawPath"
    }

    $variables = @{}
    $valuesStarted = $false
    $lastRow = $null
    foreach ($line in Get-Content $rawPath) {
        if ($line -match '^\s*(\d+)\s+(\S+)\s+\S+\s*$' -and -not $valuesStarted) {
            $variables[$matches[2]] = [int]$matches[1]
            continue
        }
        if ($line -match '^Values:') {
            $valuesStarted = $true
            continue
        }
        if ($valuesStarted -and $line.Trim().Length -gt 0) {
            $lastRow = $line.Trim() -split '\s+'
        }
    }

    if (-not $variables.ContainsKey($Signal)) {
        throw "Signal '$Signal' not found in RAW variables. Available: $($variables.Keys -join ', ')"
    }
    if ($null -eq $lastRow) {
        throw "RAW file contained no value rows"
    }
    $index = $variables[$Signal]
    if ($index -ge $lastRow.Count) {
        throw "Signal index $index for '$Signal' is outside final row with $($lastRow.Count) value(s)"
    }
    $actual = [double]$lastRow[$index]
    $error = [Math]::Abs($actual - $Expected)
    Write-Host ("RAW check {0}: actual={1:E12} expected={2:E12} tolerance={3:E12}" -f $Signal, $actual, $Expected, $Tolerance)
    if ($error -gt $Tolerance) {
        throw "RAW numeric check failed for '$Signal': error=$error"
    }
}
finally {
    if (Test-Path $rawPath) {
        Remove-Item -LiteralPath $rawPath -Force
    }
}
