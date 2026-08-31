param(
    [Parameter(Mandatory=$true)][string]$GspiceExe,
    [Parameter(Mandatory=$true)][string]$NgspiceExe
)

function Write-Deck([string]$Path, [string[]]$Lines) {
    Set-Content -LiteralPath $Path -Value $Lines -Encoding ASCII
}

function Read-Raw([string]$Path) {
    $vars = New-Object System.Collections.Generic.List[string]
    $rows = New-Object System.Collections.Generic.List[object]
    $inVars = $false
    $inValues = $false
    $current = @()
    foreach ($line in Get-Content -LiteralPath $Path) {
        $s = $line.Trim()
        if ($s -eq "Variables:") { $inVars = $true; $inValues = $false; continue }
        if ($s -eq "Values:") { $inVars = $false; $inValues = $true; continue }
        if ($inVars -and $s -match '^\d+\s+(\S+)\s+\S+') {
            $vars.Add($matches[1])
            continue
        }
        if (-not $inValues -or -not $s) { continue }
        $parts = @($s -split '\s+')
        if ($parts.Count -ge 2 -and $parts[0] -match '^\d+$') {
            if ($current.Count -eq $vars.Count) { $rows.Add($current) }
            $current = @($parts[1..($parts.Count - 1)] | ForEach-Object { Convert-RawNumber $_ })
        } elseif ($parts.Count -eq $vars.Count) {
            if ($current.Count -eq $vars.Count) { $rows.Add($current) }
            $current = @($parts | ForEach-Object { Convert-RawNumber $_ })
        } elseif ($parts.Count -eq 1) {
            if ($current.Count -eq $vars.Count) { $rows.Add($current); $current = @() }
            $current += Convert-RawNumber $parts[0]
        }
    }
    if ($current.Count -eq $vars.Count) { $rows.Add($current) }
    if ($vars.Count -eq 0 -or $rows.Count -eq 0) { throw "Could not parse RAW file: $Path" }
    return @{ Vars = $vars; Rows = $rows }
}

function Convert-RawNumber([string]$Token) {
    if ($Token -match '^([^,]+),([^,]+)$') {
        return [pscustomobject]@{ Re = [double]$matches[1]; Im = [double]$matches[2] }
    }
    return [double]$Token
}

function Value-At([hashtable]$Raw, [int]$Row, [string]$Name) {
    $idx = $Raw.Vars.IndexOf($Name)
    if ($idx -lt 0) { throw "Signal '$Name' not found. Variables: $($Raw.Vars -join ', ')" }
    return $Raw.Rows[$Row][$idx]
}

function Mag([object]$Value) {
    if ($Value -is [double]) { return [math]::Abs($Value) }
    return [math]::Sqrt($Value.Re * $Value.Re + $Value.Im * $Value.Im)
}

function PhaseDeg([object]$Value) {
    if ($Value -is [double]) { return $(if ($Value -lt 0.0) { 180.0 } else { 0.0 }) }
    return [math]::Atan2($Value.Im, $Value.Re) * 180.0 / [math]::PI
}

function Assert-Close([string]$Label, [double]$Actual, [double]$Expected, [double]$Tol) {
    if ([double]::IsNaN($Actual) -or [double]::IsInfinity($Actual)) {
        throw "$Label is not finite: $Actual"
    }
    if ([math]::Abs($Actual - $Expected) -gt $Tol) {
        throw "$Label mismatch: actual=$Actual expected=$Expected tolerance=$Tol"
    }
}

$root = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_ngspice_oracle_{0}" -f ([System.Guid]::NewGuid().ToString("N")))
New-Item -ItemType Directory -Path $root | Out-Null
try {
    $gPacDeck = Join-Path $root "gspice_psspac.sp"
    $gNoiseDeck = Join-Path $root "gspice_pnoise.sp"
    $gPssDeck = Join-Path $root "gspice_pss.sp"
    $nAcDeck = Join-Path $root "ng_ac.sp"
    $nNoiseDeck = Join-Path $root "ng_noise.sp"
    $nOpDeck = Join-Path $root "ng_op.sp"

    Write-Deck $gPacDeck @(
        "* GSPICE PSSPAC linear reduction",
        "V1 in 0 DC 0 AC 1",
        "R1 in out 1k",
        "C1 out 0 1n",
        ".PSSPAC 1k DEC 3 10 10k SIDEBANDS=1",
        ".END"
    )
    Write-Deck $nAcDeck @(
        "* NGspice AC oracle",
        "V1 in 0 DC 0 AC 1",
        "R1 in out 1k",
        "C1 out 0 1n",
        ".control",
        "set filetype=ascii",
        "ac dec 3 10 10k",
        "write ng_ac.raw",
        "quit",
        ".endc",
        ".end"
    )
    Write-Deck $gNoiseDeck @(
        "* GSPICE PNOISE linear reduction",
        ".OPTIONS TEMP=25",
        "V1 in 0 DC 0 AC 0",
        "R1 in out 1k",
        "R2 out 0 1k",
        ".PNOISE V(out) V1 DEC 1 1 10 FUND=1k SIDEBANDS=1",
        ".END"
    )
    Write-Deck $nNoiseDeck @(
        "* NGspice NOISE oracle",
        ".OPTIONS TEMP=25",
        "V1 in 0 DC 0 AC 0",
        "R1 in out 1k",
        "R2 out 0 1k",
        ".control",
        "set filetype=ascii",
        "noise v(out) V1 dec 1 1 10",
        "setplot noise1",
        "write ng_noise.raw",
        "quit",
        ".endc",
        ".end"
    )
    Write-Deck $gPssDeck @(
        "* GSPICE PSS static divider",
        "V1 in 0 DC 1",
        "R1 in out 1k",
        "R2 out 0 1k",
        ".PSS 1k 1",
        ".END"
    )
    Write-Deck $nOpDeck @(
        "* NGspice OP oracle",
        "V1 in 0 DC 1",
        "R1 in out 1k",
        "R2 out 0 1k",
        ".control",
        "set filetype=ascii",
        "op",
        "write ng_op.raw",
        "quit",
        ".endc",
        ".end"
    )

    $gPacRaw = Join-Path $root "gspice_psspac.raw"
    $gNoiseRaw = Join-Path $root "gspice_pnoise.raw"
    $gPssRaw = Join-Path $root "gspice_pss.raw"
    & $GspiceExe --threads 1 -o $gPacRaw $gPacDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "GSPICE PSSPAC failed with code $LASTEXITCODE" }
    & $GspiceExe --threads 1 -o $gNoiseRaw $gNoiseDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "GSPICE PNOISE failed with code $LASTEXITCODE" }
    & $GspiceExe --threads 1 -o $gPssRaw $gPssDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "GSPICE PSS failed with code $LASTEXITCODE" }

    Push-Location $root
    & $NgspiceExe -b $nAcDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "NGspice AC failed with code $LASTEXITCODE" }
    & $NgspiceExe -b $nNoiseDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "NGspice NOISE failed with code $LASTEXITCODE" }
    & $NgspiceExe -b $nOpDeck | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "NGspice OP failed with code $LASTEXITCODE" }
    Pop-Location

    $gPac = Read-Raw $gPacRaw
    $nAc = Read-Raw (Join-Path $root "ng_ac.raw")
    if ($gPac.Rows.Count -ne $nAc.Rows.Count) {
        throw "PSSPAC/NGspice AC point-count mismatch: $($gPac.Rows.Count) vs $($nAc.Rows.Count)"
    }
    for ($i = 0; $i -lt $gPac.Rows.Count; $i++) {
        $nOut = Value-At $nAc $i "v(out)"
        Assert-Close "PSSPAC/NGspice frequency[$i]" (Value-At $gPac $i "frequency") (Mag (Value-At $nAc $i "frequency")) 1e-9
        Assert-Close "PSSPAC/NGspice V(out)[$i]" (Value-At $gPac $i "V(out)") (Mag $nOut) 1e-8
        Assert-Close "PSSPAC/NGspice phase(V(out))[$i]" (Value-At $gPac $i "phase(V(out))") (PhaseDeg $nOut) 1e-5
    }

    $gNoise = Read-Raw $gNoiseRaw
    $nNoise = Read-Raw (Join-Path $root "ng_noise.raw")
    if ($gNoise.Rows.Count -ne $nNoise.Rows.Count) {
        throw "PNOISE/NGspice NOISE point-count mismatch: $($gNoise.Rows.Count) vs $($nNoise.Rows.Count)"
    }
    for ($i = 0; $i -lt $gNoise.Rows.Count; $i++) {
        $ngOnoise = [double](Value-At $nNoise $i "onoise_spectrum")
        Assert-Close "PNOISE/NGspice frequency[$i]" (Value-At $gNoise $i "frequency") (Value-At $nNoise $i "frequency") 1e-12
        Assert-Close "PNOISE/NGspice sqrt[$i]" (Value-At $gNoise $i "pnoise_sqrt(V/rtHz)") $ngOnoise 1e-15
        Assert-Close "PNOISE/NGspice PSD[$i]" (Value-At $gNoise $i "pnoise_psd(V^2/Hz)") ($ngOnoise * $ngOnoise) 1e-23
    }

    $gPss = Read-Raw $gPssRaw
    $nOp = Read-Raw (Join-Path $root "ng_op.raw")
    Assert-Close "PSS/NGspice OP V(in)" (Value-At $gPss 0 "V(in)") (Value-At $nOp 0 "v(in)") 1e-9
    Assert-Close "PSS/NGspice OP V(out)" (Value-At $gPss 0 "V(out)") (Value-At $nOp 0 "v(out)") 1e-9

    Write-Host "NGspice periodic oracle checks passed."
}
finally {
    Pop-Location -ErrorAction SilentlyContinue
    if (Test-Path -LiteralPath $root) {
        Remove-Item -LiteralPath $root -Recurse -Force
    }
}
