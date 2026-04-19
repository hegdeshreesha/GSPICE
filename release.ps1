# GSPICE Release Script
# Automates the creation of a standalone distribution package.

$Version = "0.4.0"
$DistDir = "dist/gspice_v$Version`_win64"

Write-Host "--- Packaging GSPICE v$Version ---" -ForegroundColor Cyan

# 1. Clean and Build in Release mode
Remove-Item -Path build -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path build
cd build
& "C:\Program Files\CMake\bin\cmake.exe" ..
& "C:\Program Files\CMake\bin\cmake.exe" --build . --config Release
cd ..

# 2. Create Distribution Folder
Remove-Item -Path $DistDir -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $DistDir -Force
New-Item -ItemType Directory -Path "$DistDir/examples" -Force
New-Item -ItemType Directory -Path "$DistDir/include" -Force

# 3. Copy Executable and Assets
Copy-Item "build/Release/gspice.exe" "$DistDir/"
Copy-Item "README.md" "$DistDir/"
Copy-Item "test_*.sp" "$DistDir/examples/"
Copy-Item "include/osdi.h" "$DistDir/include/"

# 4. Zip the release
$ZipFile = "dist/gspice_v$Version`_win64.zip"
Remove-Item $ZipFile -ErrorAction SilentlyContinue
Compress-Archive -Path "$DistDir/*" -DestinationPath $ZipFile

Write-Host "`nRelease Package Created: $ZipFile" -ForegroundColor Green
Write-Host "Standalone distribution ready for other computers."
