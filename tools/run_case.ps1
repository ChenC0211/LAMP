<#
.SYNOPSIS
Run a named configuration from cases/ in a fresh runs/<case>/<timestamp>/.
.EXAMPLE
.\tools\run_case.cmd hermite_gals_vortex_xy
.EXAMPLE
.\tools\run_case.cmd -List
#>
[CmdletBinding()]
param(
    [Parameter(Position=0)][string]$Case,
    [switch]$List,
    [ValidatePattern('^[A-Za-z0-9_-]*$')][string]$Label='',
    [ValidateRange(1,256)][int]$Threads=2,
    [string]$Program='',
    [switch]$PrepareOnly,
    [switch]$CheckOnly,
    [switch]$SkipCheck
)
$ErrorActionPreference='Stop'
$project=Split-Path -Parent $PSScriptRoot
$casesRoot=[IO.Path]::GetFullPath((Join-Path $project 'cases'))
$separator=[IO.Path]::DirectorySeparatorChar
if($List){
    Get-ChildItem -LiteralPath $casesRoot -Filter run.dict -File -Recurse |
        Where-Object { $_.Directory.Name -eq 'system' -and $_.FullName -notmatch '[\\/]output[\\/]' } |
        ForEach-Object { $_.Directory.Parent.FullName.Substring($casesRoot.Length+1).Replace('\','/') } |
        Sort-Object
    Get-ChildItem -LiteralPath $casesRoot -Filter benchmark.json -File -Recurse |
        ForEach-Object { $_.Directory.FullName.Substring($casesRoot.Length+1).Replace('\','/') } |
        Sort-Object
    return
}
if(-not $Case){throw 'Usage: tools\run_case.cmd <case-name> [-CheckOnly|-PrepareOnly]; use -List to show names.'}
if($PrepareOnly -and $CheckOnly){throw 'Use either -PrepareOnly or -CheckOnly.'}
if($SkipCheck -and $CheckOnly){throw '-SkipCheck cannot be combined with -CheckOnly.'}
$candidate=if(Test-Path -LiteralPath $Case){$Case}else{Join-Path $casesRoot $Case}
$casePath=(Resolve-Path -LiteralPath $candidate).Path.TrimEnd($separator)
if(-not $casePath.StartsWith($casesRoot+$separator,[StringComparison]::OrdinalIgnoreCase)){throw 'Select a configuration inside cases/.'}
$caseName=$casePath.Substring($casesRoot.Length+1).Replace('\','/')
if($caseName -match '(^|/)(output|runtime)(/|$)'){throw 'Select an original configuration, not a saved run.'}
for($cursor=Get-Item -LiteralPath $casePath; $cursor.FullName -ne $casesRoot; $cursor=$cursor.Parent){
    if($cursor.Attributes -band [IO.FileAttributes]::ReparsePoint){throw 'Case directory links are not supported.'}
}
$benchmarkPath=Join-Path $casePath 'benchmark.json'
$benchmark=$null
if(Test-Path -LiteralPath $benchmarkPath -PathType Leaf){
    if((Get-Item -LiteralPath $benchmarkPath).Attributes -band [IO.FileAttributes]::ReparsePoint){throw 'Benchmark descriptor links are not supported.'}
    $benchmark=Get-Content -LiteralPath $benchmarkPath -Raw -Encoding UTF8 | ConvertFrom-Json
    if($benchmark.program -notmatch '^LAMP_[A-Za-z0-9_]+_benchmark\.exe$' -or
        $benchmark.arguments -isnot [Array]){throw 'Benchmark requires a LAMP_<name>_benchmark.exe filename and an arguments array.'}
    foreach($argument in $benchmark.arguments){
        if($argument -isnot [string] -or $argument.Contains([string][char]0)){throw 'Benchmark arguments must be strings without NUL characters.'}
    }
}elseif(-not (Test-Path -LiteralPath (Join-Path $casePath 'system/run.dict') -PathType Leaf)){throw 'This case has no system/run.dict or benchmark.json; use -List.'}
if(-not $Program){$Program=if($benchmark){Join-Path $project ('bin/benchmarks/'+$benchmark.program)}else{Join-Path $project 'bin/LAMP.exe'}}
if(-not (Test-Path -LiteralPath $Program)){throw "Build the program first: tools\build.cmd; benchmark drivers: tools\build.cmd -Preset local-benchmarks"}
$Program=(Resolve-Path -LiteralPath $Program).Path
$stamp=Get-Date -Format 'yyyyMMdd_HHmmss_fff'
if($Label){$stamp+='_'+$Label}
$runRoot=Join-Path $project ('runs/'+$caseName.Replace('/','__'))
$run=Join-Path $runRoot ($stamp+'_'+[Guid]::NewGuid().ToString('N').Substring(0,8))
New-Item -ItemType Directory -Path $run | Out-Null
$inputFolders=@('system','constant','input','data','reference')
if($benchmark){Copy-Item -LiteralPath $benchmarkPath -Destination (Join-Path $run 'benchmark.json')}
foreach($folder in $inputFolders){
    $source=Join-Path $casePath $folder
    if(Test-Path -LiteralPath $source){
        $linked=@(Get-Item -LiteralPath $source; Get-ChildItem -LiteralPath $source -Recurse -Force) |
            Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint }
        if($linked){throw "Input links are not supported: $source"}
        Copy-Item -LiteralPath $source -Destination $run -Recurse
    }
}
$runtime=Join-Path $run 'runtime'
New-Item -ItemType Directory -Path $runtime | Out-Null
$exe=Join-Path $runtime ([IO.Path]::GetFileName($Program))
Copy-Item -LiteralPath $Program -Destination $exe
Get-ChildItem -LiteralPath (Split-Path -Parent $Program) -Filter '*.dll' -File | ForEach-Object {Copy-Item -LiteralPath $_.FullName -Destination $runtime}
$utf8=[Text.UTF8Encoding]::new($false)
$record=[ordered]@{
    case_name=$caseName;case=$casePath;run=$run;created=(Get-Date).ToString('o');status='prepared'
    check_only=[bool]$CheckOnly;skip_check=[bool]$SkipCheck;threads=$Threads
    kind=$(if($benchmark){'benchmark'}else{'case'});arguments=$(if($benchmark){@($benchmark.arguments)}else{@('--case',$run)})
    program=$exe;program_sha256=(Get-FileHash -LiteralPath $exe -Algorithm SHA256).Hash
    inputs=@(foreach($folder in $inputFolders){
        $path=Join-Path $run $folder
        if(Test-Path -LiteralPath $path){Get-ChildItem -LiteralPath $path -File -Recurse | Sort-Object FullName | ForEach-Object {
            [ordered]@{path=$_.FullName.Substring($run.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash}
        }}
    })
}
if($benchmark){$record.inputs+=@([ordered]@{path='benchmark.json';sha256=(Get-FileHash -LiteralPath (Join-Path $run 'benchmark.json') -Algorithm SHA256).Hash})}
function Save-Record {[IO.File]::WriteAllText((Join-Path $run 'run.json'),($record|ConvertTo-Json -Depth 6),$utf8)}
function Invoke-Program([string]$Log,[switch]$InputCheck){
    $arguments=@('--case',$run)
    if($benchmark){$arguments=@($benchmark.arguments)}
    if($InputCheck){$arguments+='--check-case'}
    $previousPreference=$ErrorActionPreference
    try{
        $ErrorActionPreference='Continue'
        $global:LASTEXITCODE=$null
        & $exe @arguments *> $Log
        if($null -eq $global:LASTEXITCODE){throw 'The program did not return an exit code.'}
        return $global:LASTEXITCODE
    }finally{$ErrorActionPreference=$previousPreference}
}
Save-Record
[IO.File]::WriteAllText((Join-Path $runRoot 'latest.txt'),$run,$utf8)
Write-Output "Run directory: $run"
if($PrepareOnly){Write-Output 'Prepared inputs and program; no physical steps.';return}
$saved=@{}
foreach($key in @('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OMP_DYNAMIC','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')){$saved[$key]=[Environment]::GetEnvironmentVariable($key,'Process')}
$env:OMP_NUM_THREADS=[string]$Threads;$env:OMP_THREAD_LIMIT=[string]$Threads;$env:OMP_DYNAMIC='FALSE';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
$watch=[Diagnostics.Stopwatch]::StartNew()
try{
    if($benchmark -and $CheckOnly){$record.status='input_checked';$record.exit_code=0;Write-Output 'Benchmark descriptor and runtime snapshot checked; no physical steps.';return}
    if(-not $SkipCheck -and -not $benchmark){
        $record.status='checking';Save-Record
        $log=Join-Path $run 'check.log'
        $code=Invoke-Program $log -InputCheck
        $record.check_exit_code=$code
        if($code -ne 0){
            $record.status='input_failed';$record.exit_code=$code
            Get-Content -LiteralPath $log -Tail 12
            Write-Output "Input check failed ($code): $log"
            exit $code
        }
        if($CheckOnly){$record.status='input_checked';$record.exit_code=0;Write-Output "Input check OK: $log";return}
    }
    $record.status='running';Save-Record
    $log=Join-Path $run 'run.log'
    Write-Output "Running; log: $log"
    Push-Location $run
    try {$code=Invoke-Program $log} finally {Pop-Location}
    $record.status=if($code -eq 0){'completed'}else{'failed'}
    $record.exit_code=$code
    if($code -ne 0){Get-Content -LiteralPath $log -Tail 12;Write-Output "Run failed ($code): $log";exit $code}
    if($benchmark){Write-Output "Completed benchmark; results: $log"}else{Write-Output "Completed: $(Join-Path $run 'output')"}
    $pvd=Join-Path $run 'output/fields/solution.pvd'
    if(Test-Path -LiteralPath $pvd){Write-Output "ParaView: $pvd"}
}catch{
    $record.status='launcher_failed';$record.launcher_error=$_.Exception.Message
    throw
}finally{
    $watch.Stop();$record.seconds=$watch.Elapsed.TotalSeconds;$record.finished=(Get-Date).ToString('o');Save-Record
    foreach($key in $saved.Keys){[Environment]::SetEnvironmentVariable($key,$saved[$key],'Process')}
}
