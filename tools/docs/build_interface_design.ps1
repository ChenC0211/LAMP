param([string]$XeLaTeX = 'xelatex')

$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../..'))
$sourceDir = Join-Path $projectRoot 'docs/history_todo'
$outputDir = Join-Path $projectRoot ('runs/checks/documents/'+(Get-Date -Format 'yyyyMMdd_HHmmss_fff')+'_interface_topology_design')
$document = 'QUANTUM_Interface_Topology_Benchmark_Design'
$compiler = (Get-Command $XeLaTeX -ErrorAction Stop).Source
New-Item -ItemType Directory -Path $outputDir -Force | Out-Null

Push-Location $sourceDir
try {
    $tocPath = Join-Path $outputDir "$document.toc"
    for ($pass = 1; $pass -le 4; $pass++) {
        $previousToc = if (Test-Path -LiteralPath $tocPath) {
            Get-Content -LiteralPath $tocPath -Raw
        } else { '' }
        $arguments = @(
            '--disable-installer', '-interaction=nonstopmode', '-halt-on-error',
            "-output-directory=$outputDir", "$document.tex"
        )
        $trace = Join-Path $outputDir "$document.pass$pass.txt"
        # MiKTeX may emit nonfatal user-log warnings on stderr.
        $ErrorActionPreference = 'Continue'
        try {
            & $compiler @arguments *> $trace
            $compilerExit = $LASTEXITCODE
        }
        finally { $ErrorActionPreference = 'Stop' }
        if ($compilerExit -ne 0) { throw "XeLaTeX failed. See $trace" }
        $logPath = Join-Path $outputDir "$document.log"
        $log = Get-Content -LiteralPath $logPath -Raw
        $currentToc = Get-Content -LiteralPath $tocPath -Raw
        $needsRerun = ($previousToc -cne $currentToc) -or
            ($log -match 'Rerun to|Rerun LaTeX|Label\(s\) may have changed')
        if ($pass -ge 2 -and -not $needsRerun) { break }
        if ($pass -eq 4 -and $needsRerun) { throw "References did not settle. See $logPath" }
    }
    if ($log -match 'undefined references|Citation .* undefined|Missing character:') {
        throw "Unresolved references or missing glyphs. See $logPath"
    }
    Write-Host "Built in $outputDir/$document.pdf"
}
finally { Pop-Location }
