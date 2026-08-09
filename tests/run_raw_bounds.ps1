param(
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Deck,
    [Parameter(Mandatory=$true)][string]$Signal,
    [Parameter(Mandatory=$true)][double]$Min,
    [Parameter(Mandatory=$true)][double]$Max
)

$rawPath = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_raw_bounds_{0}.raw" -f ([System.Guid]::NewGuid().ToString("N")))
$output = & $Exe $Deck -o $rawPath --format raw --save all 2>&1
if ($LASTEXITCODE -ne 0) {
    $output | Write-Host
    throw "gspice exited with code $LASTEXITCODE"
}
if (!(Test-Path -LiteralPath $rawPath)) {
    throw "RAW output was not created: $rawPath"
}

$lines = Get-Content -LiteralPath $rawPath
$variables = @()
$mode = ""
foreach ($line in $lines) {
    $trim = $line.Trim()
    if ($trim -eq "Variables:") { $mode = "vars"; continue }
    if ($trim -eq "Values:") { break }
    if ($mode -eq "vars" -and $trim.Length -gt 0) {
        $parts = $trim -split "\s+"
        if ($parts.Count -ge 3 -and $parts[0] -match "^\d+$") {
            $variables += $parts[1]
        }
    }
}

$signalIndex = [Array]::IndexOf($variables, $Signal)
if ($signalIndex -lt 0) {
    throw "Signal '$Signal' not found in RAW variables: $($variables -join ', ')"
}

$valuesMode = $false
$seen = 0
$actualMin = [double]::PositiveInfinity
$actualMax = [double]::NegativeInfinity
foreach ($line in $lines) {
    $trim = $line.Trim()
    if ($trim -eq "Values:") { $valuesMode = $true; continue }
    if (!$valuesMode -or $trim.Length -eq 0) { continue }
    $parts = $trim -split "\s+"
    if ($parts.Count -eq $variables.Count + 1 -and $parts[0] -match "^\d+$") {
        $value = [double]$parts[$signalIndex + 1]
    } elseif ($parts.Count -eq $variables.Count) {
        $value = [double]$parts[$signalIndex]
    } else {
        continue
    }
    $actualMin = [Math]::Min($actualMin, $value)
    $actualMax = [Math]::Max($actualMax, $value)
    $seen += 1
}

Remove-Item -LiteralPath $rawPath -Force -ErrorAction SilentlyContinue
if ($seen -eq 0) {
    throw "No RAW values found for $Signal"
}
if ($actualMin -lt $Min -or $actualMax -gt $Max) {
    throw "$Signal out of bounds: min=$actualMin max=$actualMax expected [$Min, $Max]"
}
Write-Host "$Signal bounds OK: min=$actualMin max=$actualMax points=$seen"
