param(
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Deck,
    [string]$ExpectedRegex = "Simulation Completed Successfully"
)

$output = & $Exe --threads 1 $Deck 2>&1
$exitCode = $LASTEXITCODE
$text = $output -join "`n"
Write-Host $text

if ($exitCode -ne 0) {
    Write-Error "GSPICE exited with code $exitCode."
    exit 1
}

if ($text -match '(?im)\bwarning\b|WARNING:') {
    Write-Error "Unexpected warning in output."
    exit 1
}

if ($text -notmatch $ExpectedRegex) {
    Write-Error "Expected output to match '$ExpectedRegex'."
    exit 1
}

exit 0
