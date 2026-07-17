param(
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Deck,
    [Parameter(Mandatory=$true)][int]$ExpectedVariables,
    [Parameter(Mandatory=$true)][string[]]$ExpectedSignals,
    [string[]]$AbsentSignals = @()
)

$rawPath = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_raw_header_{0}.raw" -f ([System.Guid]::NewGuid().ToString("N")))
try {
    & $Exe --threads 1 -o $rawPath $Deck
    if ($LASTEXITCODE -ne 0) {
        throw "GSPICE exited with code $LASTEXITCODE"
    }
    if (-not (Test-Path $rawPath)) {
        throw "Expected RAW output was not created: $rawPath"
    }

    $variables = New-Object System.Collections.Generic.List[string]
    $count = $null
    foreach ($line in Get-Content $rawPath) {
        if ($line -match '^No\.\s+Variables:\s+(\d+)\s*$') {
            $count = [int]$matches[1]
            continue
        }
        if ($line -match '^Values:') {
            break
        }
        if ($line -match '^\s*\d+\s+(\S+)\s+\S+\s*$') {
            $variables.Add($matches[1])
        }
    }

    if ($null -eq $count) {
        throw "RAW header did not contain No. Variables"
    }
    if ($count -ne $ExpectedVariables) {
        throw "Expected $ExpectedVariables RAW variables, found $count. Variables: $($variables -join ', ')"
    }
    foreach ($signal in $ExpectedSignals) {
        if (-not $variables.Contains($signal)) {
            throw "Expected signal '$signal' not found. Variables: $($variables -join ', ')"
        }
    }
    foreach ($signal in $AbsentSignals) {
        if ($variables.Contains($signal)) {
            throw "Signal '$signal' should not have been written. Variables: $($variables -join ', ')"
        }
    }
    Write-Host ("RAW header check passed: count={0}; variables={1}" -f $count, ($variables -join ', '))
}
finally {
    if (Test-Path $rawPath) {
        Remove-Item -LiteralPath $rawPath -Force
    }
}
