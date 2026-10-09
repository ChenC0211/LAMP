param([string]$XeLaTeX = 'xelatex')

$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../..'))
$sourceDir = Join-Path $projectRoot 'docs/manuals'
$outputDir = Join-Path $projectRoot ('runs/checks/documents/'+(Get-Date -Format 'yyyyMMdd_HHmmss_fff')+'_manuals')
$compiler = (Get-Command $XeLaTeX -ErrorAction Stop).Source
$documents = @(
    'QUANTUM_Architecture_Models',
    'QUANTUM_User_Guide',
    'QUANTUM_Dictionary_Guide'
)
New-Item -ItemType Directory -Path $outputDir -Force | Out-Null

Push-Location $sourceDir
try {
    foreach ($document in $documents) {
        $tocPath = Join-Path $outputDir "$document.toc"
        for ($pass = 1; $pass -le 4; $pass++) {
            $previousToc = if (Test-Path -LiteralPath $tocPath) {
                Get-Content -LiteralPath $tocPath -Raw
            } else { '' }
            $trace = Join-Path $outputDir "$document.pass$pass.txt"
            $arguments = @(
                '--disable-installer', '-interaction=nonstopmode', '-halt-on-error',
                "-output-directory=$outputDir", "$document.tex"
            )
            # Native stderr may contain harmless MiKTeX user-log warnings.
            # The process exit code and TeX log decide whether compilation passed.
            $ErrorActionPreference = 'Continue'
            try {
                & $compiler @arguments *> $trace
                $compilerExit = $LASTEXITCODE
            }
            finally { $ErrorActionPreference = 'Stop' }
            if ($compilerExit -ne 0) {
                throw "XeLaTeX failed for $document (pass $pass). See $trace"
            }
            $logPath = Join-Path $outputDir "$document.log"
            $log = Get-Content -LiteralPath $logPath -Raw
            $currentToc = Get-Content -LiteralPath $tocPath -Raw
            $needsRerun = ($previousToc -cne $currentToc) -or
                ($log -match 'Rerun to|Rerun LaTeX|Label\(s\) may have changed')
            if ($pass -ge 2 -and -not $needsRerun) {
                break
            }
            if ($pass -eq 4 -and $needsRerun) {
                throw "References did not settle for $document. See $logPath"
            }
        }
        Write-Host "Built in $outputDir/manuals/$document.pdf"
    }
}
finally { Pop-Location }
