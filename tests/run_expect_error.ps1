param(
    [string]$Exe,
    [string]$Deck,
    [string]$ExpectedRegex
)

$output = & $Exe --threads 1 $Deck 2>&1
$exitCode = $LASTEXITCODE
$text = $output -join "`n"
Write-Host $text

if ($exitCode -eq 0) {
    Write-Error "Expected GSPICE to fail, but it exited with code 0."
    exit 1
}

if ($text -notmatch $ExpectedRegex) {
    Write-Error "Expected output to match '$ExpectedRegex'."
    exit 1
}

exit 0
