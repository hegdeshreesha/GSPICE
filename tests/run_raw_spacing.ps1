param(
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Deck,
    [Parameter(Mandatory=$true)][double]$MaxMinDt,
    [Parameter(Mandatory=$true)][int]$MinPoints
)

$rawPath = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_raw_spacing_{0}.raw" -f ([System.Guid]::NewGuid().ToString("N")))
$output = & $Exe $Deck -o $rawPath --format raw --save all 2>&1
if ($LASTEXITCODE -ne 0) {
    $output | Write-Host
    throw "gspice exited with code $LASTEXITCODE"
}
if (!(Test-Path -LiteralPath $rawPath)) {
    throw "RAW output was not created: $rawPath"
}

$times = New-Object System.Collections.Generic.List[double]
$valuesMode = $false
foreach ($line in Get-Content -LiteralPath $rawPath) {
    $trim = $line.Trim()
    if ($trim -eq "Values:") { $valuesMode = $true; continue }
    if (!$valuesMode -or $trim.Length -eq 0) { continue }
    $parts = $trim -split "\s+"
    if ($parts.Count -lt 2) { continue }
    if ($parts[0] -match "^\d+$") {
        $times.Add([double]$parts[1])
    } else {
        $times.Add([double]$parts[0])
    }
}

Remove-Item -LiteralPath $rawPath -Force -ErrorAction SilentlyContinue
if ($times.Count -lt $MinPoints) {
    throw "RAW point count $($times.Count) below required $MinPoints"
}

$minDt = [double]::PositiveInfinity
for ($i = 1; $i -lt $times.Count; $i++) {
    $dt = $times[$i] - $times[$i - 1]
    if ($dt -le 0.0) {
        throw "RAW times are not strictly increasing at index ${i}: $($times[$i - 1]) then $($times[$i])"
    }
    $minDt = [Math]::Min($minDt, $dt)
}
if ($minDt -gt $MaxMinDt) {
    throw "RAW minimum spacing $minDt exceeds required $MaxMinDt"
}
Write-Host "RAW spacing OK: points=$($times.Count) min_dt=$minDt"
