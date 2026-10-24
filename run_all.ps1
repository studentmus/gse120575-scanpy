# Run the whole pipeline inside the conda env (conda run sets up DLL paths on Windows).
# no ErrorActionPreference=Stop: PS 5.1 would treat Python warnings on stderr as errors
$env:PYTHONIOENCODING = "utf-8"
$conda = "$env:USERPROFILE\miniforge3\Scripts\conda.exe"
foreach ($s in "01_load", "02_qc", "03_cluster", "04_annotate", "05_proportions") {
    if (-not (Test-Path "scripts\$s.py")) { continue }
    Write-Host "== $s"
    & $conda run -n scanpy --no-capture-output python -u "scripts\$s.py"
    if ($LASTEXITCODE -ne 0) { throw "$s failed" }
}
