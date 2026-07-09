param(
    [string]$GspiceExe = ""
)

$ErrorActionPreference = "Stop"
$repo = Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")
if ([string]::IsNullOrWhiteSpace($GspiceExe)) {
    $GspiceExe = Join-Path $repo "build\Release\gspice.exe"
}

if (!(Test-Path -LiteralPath $GspiceExe)) {
    throw "gspice.exe was not found at $GspiceExe"
}

$folder = Split-Path -Parent $GspiceExe
$requiredDlls = @("klu.dll", "amd.dll", "colamd.dll", "btf.dll", "suitesparseconfig.dll")
foreach ($dll in $requiredDlls) {
    $path = Join-Path $folder $dll
    if (!(Test-Path -LiteralPath $path)) {
        throw "Missing required KLU runtime DLL: $path"
    }
}

$deck = Join-Path $repo "tests\decks\solver_klu.sp"
$output = & $GspiceExe --threads 16 $deck 2>&1
$text = $output -join "`n"
Write-Host $text

if ($LASTEXITCODE -ne 0) {
    throw "GSPICE KLU smoke deck failed with exit code $LASTEXITCODE"
}
if ($text -notmatch "ext_klu_calls=[1-9]") {
    throw "GSPICE ran, but external KLU was not used."
}

Write-Host "GSPICE environment OK: SuiteSparse/KLU is active."
