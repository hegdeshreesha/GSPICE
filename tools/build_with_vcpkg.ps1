param(
    [string]$VcpkgRoot = $env:VCPKG_ROOT,
    [string]$Triplet = "x64-windows",
    [string]$BuildDir = "build-vcpkg",
    [switch]$SkipInstall
)

$ErrorActionPreference = "Stop"

function Find-CMake {
    $cmd = Get-Command cmake.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    $fallback = "C:\Program Files\CMake\bin\cmake.exe"
    if (Test-Path -LiteralPath $fallback) { return $fallback }
    throw "cmake.exe was not found. Install CMake or add it to PATH."
}

function Find-VcVars64 {
    $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
    if (Test-Path -LiteralPath $vswhere) {
        $install = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
        if ($install) {
            $candidate = Join-Path $install "VC\Auxiliary\Build\vcvars64.bat"
            if (Test-Path -LiteralPath $candidate) { return $candidate }
        }
    }

    $candidates = @(
        "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
    )
    foreach ($candidate in $candidates) {
        if (Test-Path -LiteralPath $candidate) { return $candidate }
    }
    throw "Visual Studio C++ Build Tools vcvars64.bat was not found."
}

function Find-Ninja([string]$Root) {
    $cmd = Get-Command ninja.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    $candidate = Get-ChildItem -Path (Join-Path $Root "downloads\tools") -Recurse -Filter ninja.exe -ErrorAction SilentlyContinue |
        Select-Object -First 1 -ExpandProperty FullName
    if ($candidate) { return $candidate }
    throw "ninja.exe was not found. Run vcpkg once, or install Ninja and add it to PATH."
}

if ([string]::IsNullOrWhiteSpace($VcpkgRoot)) {
    $VcpkgRoot = "C:\EDA\Tools\vcpkg"
}

$repo = Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")
$vcpkg = Join-Path $VcpkgRoot "vcpkg.exe"
if (!(Test-Path -LiteralPath $vcpkg)) {
    throw "vcpkg.exe was not found at $vcpkg. Clone and bootstrap vcpkg first, or set VCPKG_ROOT."
}

if (!$SkipInstall) {
    & $vcpkg install "suitesparse-klu:$Triplet"
}

$cmake = Find-CMake
$vcvars = Find-VcVars64
$ninja = Find-Ninja $VcpkgRoot
$toolchain = Join-Path $VcpkgRoot "scripts\buildsystems\vcpkg.cmake"
$buildPath = Join-Path $repo $BuildDir
$deployPath = Join-Path $repo "build\Release"

$configure = "call `"$vcvars`" && `"$cmake`" -S `"$repo`" -B `"$buildPath`" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=`"$ninja`" -DCMAKE_TOOLCHAIN_FILE=`"$toolchain`" -DVCPKG_TARGET_TRIPLET=$Triplet"
$build = "call `"$vcvars`" && `"$cmake`" --build `"$buildPath`" --config Release"

& cmd.exe /c $configure
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& cmd.exe /c $build
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

New-Item -ItemType Directory -Force -Path $deployPath | Out-Null
foreach ($name in @("gspice.exe", "klu.dll", "amd.dll", "colamd.dll", "btf.dll", "suitesparseconfig.dll")) {
    $source = Join-Path $buildPath $name
    if (Test-Path -LiteralPath $source) {
        Copy-Item -LiteralPath $source -Destination (Join-Path $deployPath $name) -Force
    }
}

Write-Host "GSPICE KLU build deployed to $deployPath"
