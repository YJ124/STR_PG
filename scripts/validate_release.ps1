$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location $ProjectRoot
python -m compileall -q src
python -m pytest -q
python -m strpg.cli --help | Out-Null
python -m strpg.cli genotype --help | Out-Null
PowerShell -ExecutionPolicy Bypass -File scripts\run_demo.ps1
PowerShell -ExecutionPolicy Bypass -File scripts\run_phmm_demo.ps1
