param(
    [Parameter(Mandatory=$true)][string]$Exe
)

function Write-Deck([string]$Path, [string[]]$Lines) {
    Set-Content -LiteralPath $Path -Value $Lines -Encoding ASCII
}

function Read-Raw([string]$Path) {
    $vars = New-Object System.Collections.Generic.List[string]
    $rows = New-Object System.Collections.Generic.List[object]
    $inVars = $false
    $inValues = $false
    foreach ($line in Get-Content -LiteralPath $Path) {
        $s = $line.Trim()
        if ($s -eq "Variables:") { $inVars = $true; $inValues = $false; continue }
        if ($s -eq "Values:") { $inVars = $false; $inValues = $true; continue }
        if ($inVars -and $s -match '^\d+\s+(\S+)\s+\S+') {
            $vars.Add($matches[1])
            continue
        }
        if ($inValues -and $s) {
            $rows.Add(@($s -split '\s+' | ForEach-Object { [double]$_ }))
        }
    }
    if ($vars.Count -eq 0 -or $rows.Count -eq 0) {
        throw "Could not parse RAW file: $Path"
    }
    return @{ Vars = $vars; Rows = $rows }
}

function Value-At([hashtable]$Raw, [int]$Row, [string]$Name) {
    $idx = $Raw.Vars.IndexOf($Name)
    if ($idx -lt 0) { throw "Signal '$Name' not found. Variables: $($Raw.Vars -join ', ')" }
    return [double]$Raw.Rows[$Row][$idx]
}

function Assert-Close([string]$Label, [double]$Actual, [double]$Expected, [double]$Tol) {
    if ([double]::IsNaN($Actual) -or [double]::IsInfinity($Actual)) {
        throw "$Label is not finite: $Actual"
    }
    if ([math]::Abs($Actual - $Expected) -gt $Tol) {
        throw "$Label mismatch: actual=$Actual expected=$Expected tolerance=$Tol"
    }
}

$root = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_periodic_equiv_{0}" -f ([System.Guid]::NewGuid().ToString("N")))
New-Item -ItemType Directory -Path $root | Out-Null
try {
    $acDeck = Join-Path $root "linear_ac.sp"
    $pacDeck = Join-Path $root "linear_psspac.sp"
    $noiseDeck = Join-Path $root "linear_noise.sp"
    $pnoiseDeck = Join-Path $root "linear_pnoise.sp"
    $pssDeck = Join-Path $root "linear_pss.sp"

    Write-Deck $acDeck @(
        "* Linear AC oracle",
        "V1 in 0 DC 0 AC 1",
        "R1 in out 1k",
        "C1 out 0 1n",
        ".AC DEC 3 10 10k",
        ".END"
    )
    Write-Deck $pacDeck @(
        "* PSSPAC must reduce to AC for this static linear circuit",
        "V1 in 0 DC 0 AC 1",
        "R1 in out 1k",
        "C1 out 0 1n",
        ".PSSPAC 1k DEC 3 10 10k SIDEBANDS=1",
        ".END"
    )
    Write-Deck $noiseDeck @(
        "* Linear noise oracle",
        "V1 in 0 DC 0 AC 0",
        "R1 in out 1k",
        "R2 out 0 1k",
        ".NOISE V(out) V1 DEC 1 1 10",
        ".END"
    )
    Write-Deck $pnoiseDeck @(
        "* PNOISE must reduce to NOISE for this static linear circuit",
        "V1 in 0 DC 0 AC 0",
        "R1 in out 1k",
        "R2 out 0 1k",
        ".PNOISE V(out) V1 DEC 1 1 10 FUND=1k SIDEBANDS=1",
        ".END"
    )
    Write-Deck $pssDeck @(
        "* PSS static divider sanity",
        "V1 in 0 DC 1",
        "R1 in out 1k",
        "R2 out 0 1k",
        ".PSS 1k 1",
        ".END"
    )

    $acRaw = Join-Path $root "ac.raw"
    $pacRaw = Join-Path $root "psspac.raw"
    $noiseRaw = Join-Path $root "noise.raw"
    $pnoiseRaw = Join-Path $root "pnoise.raw"
    $pssRaw = Join-Path $root "pss.raw"

    & $Exe --threads 1 -o $acRaw $acDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "AC oracle failed with code $LASTEXITCODE" }
    & $Exe --threads 1 -o $pacRaw $pacDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "PSSPAC run failed with code $LASTEXITCODE" }
    & $Exe --threads 1 -o $noiseRaw $noiseDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "NOISE oracle failed with code $LASTEXITCODE" }
    & $Exe --threads 1 -o $pnoiseRaw $pnoiseDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "PNOISE run failed with code $LASTEXITCODE" }
    & $Exe --threads 1 -o $pssRaw $pssDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "PSS run failed with code $LASTEXITCODE" }

    $ac = Read-Raw $acRaw
    $pac = Read-Raw $pacRaw
    if ($ac.Rows.Count -ne $pac.Rows.Count) {
        throw "AC/PSSPAC point-count mismatch: $($ac.Rows.Count) vs $($pac.Rows.Count)"
    }
    for ($i = 0; $i -lt $ac.Rows.Count; $i++) {
        Assert-Close "PSSPAC frequency[$i]" (Value-At $pac $i "frequency") (Value-At $ac $i "frequency") 1e-9
        Assert-Close "PSSPAC V(out)[$i]" (Value-At $pac $i "V(out)") (Value-At $ac $i "V(out)") 1e-9
        Assert-Close "PSSPAC phase(V(out))[$i]" (Value-At $pac $i "phase(V(out))") (Value-At $ac $i "phase(V(out))") 1e-7
    }

    $noise = Read-Raw $noiseRaw
    $pnoise = Read-Raw $pnoiseRaw
    if ($noise.Rows.Count -ne $pnoise.Rows.Count) {
        throw "NOISE/PNOISE point-count mismatch: $($noise.Rows.Count) vs $($pnoise.Rows.Count)"
    }
    for ($i = 0; $i -lt $noise.Rows.Count; $i++) {
        Assert-Close "PNOISE frequency[$i]" (Value-At $pnoise $i "frequency") (Value-At $noise $i "frequency") 1e-12
        Assert-Close "PNOISE PSD[$i]" (Value-At $pnoise $i "pnoise_psd(V^2/Hz)") (Value-At $noise $i "onoise_psd(V^2/Hz)") 1e-24
    }

    $pss = Read-Raw $pssRaw
    Assert-Close "PSS V(out)" (Value-At $pss 0 "V(out)") 0.5 1e-9
    Assert-Close "PSS residual" (Value-At $pss 0 "PSS_residual") 0.0 1e-9

    Write-Host "Periodic equivalence checks passed."
}
finally {
    if (Test-Path -LiteralPath $root) {
        Remove-Item -LiteralPath $root -Recurse -Force
    }
}
