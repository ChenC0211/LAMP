param([string]$Preset = 'local', [switch]$Test, [int]$Jobs = 1, [string[]]$Target = @())
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$cmakeCommand = Get-Command cmake -ErrorAction SilentlyContinue
if ($cmakeCommand) {
    $cmakeExe = $cmakeCommand.Source
} else {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (-not (Test-Path -LiteralPath $vswhere)) { throw 'CMake or Visual Studio 2022 C++ tools are required.' }
    $vsPath = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $vsPath) { throw 'Visual Studio C++ tools were not found.' }
    $cmakeExe = Join-Path $vsPath 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
}
Push-Location $projectRoot
try {
    & $cmakeExe --preset $Preset
    if ($LASTEXITCODE) { throw "CMake configure failed ($LASTEXITCODE)." }
    # Normal development builds the application. Explicit -Test builds all
    # registered check executables before CTest; it never runs by default.
    if (-not $Target.Count) {
        $Target = if ($Preset -like '*-benchmarks') { @('lamp_benchmarks') } else { @('LAMP') }
    }
    $targetArgs = if ($Test) { @() } else { @('--target') + $Target }
    & $cmakeExe --build --preset $Preset @targetArgs --parallel $Jobs -- /nr:false
    if ($LASTEXITCODE) { throw "CMake build failed ($LASTEXITCODE)." }
    if ($Test) {
        & (Join-Path (Split-Path $cmakeExe) 'ctest.exe') --preset $Preset
        if ($LASTEXITCODE) { throw "CTest failed ($LASTEXITCODE)." }
    }
} finally {
    Pop-Location
}
