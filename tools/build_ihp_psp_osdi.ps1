param(
    [string]$OpenVafPath = "",
    [string]$PdkRoot = "C:\EDA\LumenCircuitStudio\external\ihp_pdk\ihp-sg13g2",
    [string]$OutputDir = "C:\EDA\GSPICE\osdi"
)

$ErrorActionPreference = "Stop"

function Resolve-OpenVaf {
    param([string]$RequestedPath)

    if ($RequestedPath -and (Test-Path -LiteralPath $RequestedPath)) {
        return (Resolve-Path -LiteralPath $RequestedPath).Path
    }

    if ($env:OPENVAF -and (Test-Path -LiteralPath $env:OPENVAF)) {
        return (Resolve-Path -LiteralPath $env:OPENVAF).Path
    }

    $cmd = Get-Command openvaf -ErrorAction SilentlyContinue
    if ($cmd) {
        return $cmd.Source
    }

    $cmdReloaded = Get-Command openvaf-r -ErrorAction SilentlyContinue
    if ($cmdReloaded) {
        return $cmdReloaded.Source
    }

    $defaultReloaded = "C:\EDA\Tools\OpenVAF\openvaf-r.exe"
    if (Test-Path -LiteralPath $defaultReloaded) {
        return (Resolve-Path -LiteralPath $defaultReloaded).Path
    }

    throw "OpenVAF was not found. Pass -OpenVafPath, set OPENVAF, or add openvaf.exe/openvaf-r.exe to PATH."
}

$openvaf = Resolve-OpenVaf $OpenVafPath
$pspDir = Join-Path $PdkRoot "libs.tech\verilog-a\psp103"

if (!(Test-Path -LiteralPath $pspDir)) {
    throw "IHP PSP Verilog-A directory not found: $pspDir"
}

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
$resolvedOutput = (Resolve-Path -LiteralPath $OutputDir).Path

$models = @(
    @{ Source = "psp103.va"; Output = "psp103.osdi"; Alias = "psp103va.osdi" },
    @{ Source = "psp103_nqs.va"; Output = "psp103_nqs.osdi"; Alias = "pspnqs103va.osdi" }
)

foreach ($model in $models) {
    $source = Join-Path $pspDir $model.Source
    $output = Join-Path $resolvedOutput $model.Output
    $alias = Join-Path $resolvedOutput $model.Alias
    if (!(Test-Path -LiteralPath $source)) {
        throw "Missing Verilog-A source: $source"
    }

    Write-Host "Compiling $source with $openvaf"
    & $openvaf -I $pspDir -o $output $source
    if ($LASTEXITCODE -ne 0) {
        throw "OpenVAF failed for $source with exit code $LASTEXITCODE"
    }

    if (!(Test-Path -LiteralPath $output)) {
        throw "OpenVAF did not produce expected output: $output"
    }

    Copy-Item -LiteralPath $output -Destination $alias -Force
}

Write-Host ""
Write-Host "IHP PSP OSDI libraries copied to: $resolvedOutput"
Write-Host "Set this before running GSPICE if the deck does not contain explicit .OSDI lines:"
Write-Host "`$env:GSPICE_OSDI_DIR = `"$resolvedOutput`""
