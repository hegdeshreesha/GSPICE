param(
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Deck,
    [Parameter(Mandatory=$true)][string]$Reference,
    [Parameter(Mandatory=$false)][double]$MaxRelError = 0.01
)

$ref = @()
foreach ($line in Get-Content $Reference) {
    $t = ($line -split '\s+') | Where-Object { $_ -ne '' }
    if ($t.Count -ge 2) {
        $ref += ,@([double]$t[0], [double]$t[$t.Count - 1])
    }
}
if ($ref.Count -lt 10) {
    throw "Reference file parsed only $($ref.Count) rows: $Reference"
}

$stdout = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_pmos_idvg_{0}.txt" -f ([System.Guid]::NewGuid().ToString("N")))
$stderr = Join-Path ([System.IO.Path]::GetTempPath()) ("gspice_pmos_idvg_{0}.err" -f ([System.Guid]::NewGuid().ToString("N")))
try {
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = $Exe
    $psi.Arguments = "--threads 1 `"$Deck`""
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $proc = [System.Diagnostics.Process]::Start($psi)
    $out = $proc.StandardOutput.ReadToEnd()
    $err = $proc.StandardError.ReadToEnd()
    $proc.WaitForExit()
    if ($proc.ExitCode -ne 0) {
        throw "GSPICE exited with code $($proc.ExitCode): $err"
    }

    $ours = @()
    foreach ($line in ($out -split "`r?`n")) {
        if ($line -match '^\s*([-0-9.eE+]+)\s*\|') {
            $t = ($line -split '\s+') | Where-Object { $_ -ne '' }
            if ($t.Count -ge 2 -and $t[0] -notmatch 'sweep') {
                $vg = [double]$t[0]
                if ($vg -ge -1.501 -and $vg -le 0.001) {
                    $ours += ,@($vg, [double]$t[$t.Count - 1])
                }
            }
        }
    }
    if ($ours.Count -lt 10) {
        throw "Parsed only $($ours.Count) sweep points from console output"
    }

    if ($ours.Count -ne $ref.Count) {
        throw "Point-count mismatch: ours=$($ours.Count) vs reference=$($ref.Count)"
    }

    $worst = 0.0
    $worstAt = -1.0
    for ($i = 0; $i -lt $ref.Count; $i++) {
        $rg = $ref[$i][0]
        $refId = $ref[$i][1]
        $ourId = $ours[$i][1]
        $err = [Math]::Abs(($ourId - $refId) / [Math]::Max([Math]::Abs($refId), 1e-15))
        if ($err -gt $worst) {
            $worst = $err
            $worstAt = $rg
        }
    }
    Write-Host ("PSP103 PMOS IdVg parity: points={0} max_rel_err={1:P3} at vg={2:F3} (limit={3:P2})" -f $ours.Count, $worst, $worstAt, $MaxRelError)
    if ($worst -gt $MaxRelError) {
        throw "PSP103 PMOS IdVg parity FAILED: max_rel_err=$worst exceeds $MaxRelError at Vg=$worstAt"
    }
}
finally {
    if (Test-Path -LiteralPath $stdout) { Remove-Item -LiteralPath $stdout -Force }
    if (Test-Path -LiteralPath $stderr) { Remove-Item -LiteralPath $stderr -Force }
}
