param(
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Deck
)

$outPath = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_sp_{0}.s2p" -f ([System.Guid]::NewGuid().ToString("N")))
try {
    & $Exe --threads 1 -o $outPath $Deck
    if ($LASTEXITCODE -ne 0) {
        throw "GSPICE exited with code $LASTEXITCODE"
    }
    if (-not (Test-Path $outPath)) {
        throw "Expected Touchstone output was not created: $outPath"
    }
    $content = Get-Content $outPath
    if (-not ($content | Where-Object { $_ -match '^#\s+Hz\s+S\s+MA\s+R\s+' })) {
        throw "Touchstone option line missing"
    }
    $data = $content | Where-Object { $_ -match '^\s*[+-]?\d' } | Select-Object -First 1
    if (-not $data) {
        throw "Touchstone data row missing"
    }
    $values = $data -split '\s+' | Where-Object { $_ -ne '' }
    if ($values.Count -ne 9) {
        throw "Expected 9 numeric values for 2-port MA row, found $($values.Count): $data"
    }
    Write-Host "Touchstone check passed: $data"
}
finally {
    if (Test-Path $outPath) {
        Remove-Item -LiteralPath $outPath -Force
    }
}
