<#
.SYNOPSIS
Run a named configuration from cases/ in a fresh runs/<case>/<timestamp>/.
.EXAMPLE
.\tools\run_case.cmd vof_a_liquid_gas_drop_deformation_amr -Resume
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
    [string]$Restart='',
    [switch]$Resume,
    [string]$EndTime='',
    [Alias('h')][switch]$Help,
    [switch]$PrepareOnly,
    [switch]$CheckOnly,
    [switch]$SkipCheck
)
$ErrorActionPreference='Stop'
$threadsExplicit=$PSBoundParameters.ContainsKey('Threads')
$endTimeExplicit=$PSBoundParameters.ContainsKey('EndTime')
if($Help){
    Write-Output @'
Usage: tools\run_case.cmd <case-name> [options]
  -List                   List available cases.
  -Resume                 Choose a saved run interactively, or paste its path.
  -Restart <path>          Resume an old run or explicit checkpoint directory.
  -EndTime <code-time>     Extend the end time with -Resume or -Restart.
  -PrepareOnly            Copy inputs/runtime without checking or advancing.
  -CheckOnly              Validate inputs without physical steps.
  -Program <exe>          Select an executable; default: bin/LAMP.exe.
  -Threads <N>            Override Parallel.threads.
  -Help                   Show this help.

run.dict: Checkpoint { everySteps 200; } // accepted steps; 0 disables saving
Restart supports Hydro/A VOF checkpoints. Old runs without checkpoints and B1
cannot use this restart entry. A restart creates a new run using frozen inputs.
'@
    return
}
# Hash directly through .NET so the launcher also works in PowerShell hosts
# where the optional Get-FileHash command is unavailable.
function Get-CaseSha256([string]$Path){
    $stream=[IO.File]::OpenRead($Path)
    $algorithm=[Security.Cryptography.SHA256]::Create()
    try{return [BitConverter]::ToString($algorithm.ComputeHash($stream)).Replace('-','')}
    finally{$algorithm.Dispose();$stream.Dispose()}
}
function Assert-NoRestartLinks([string]$Path){
    $cursor=Get-Item -LiteralPath $Path -Force
    while($null -ne $cursor){
        if($cursor.Attributes -band [IO.FileAttributes]::ReparsePoint){throw "Restart directory links are not supported: $($cursor.FullName)"}
        $cursor=if($cursor.PSIsContainer){$cursor.Parent}else{$cursor.Directory}
    }
}
function Resolve-RestartSource([string]$Requested,[string]$ExpectedCase,[string]$Project){
    $entry=Get-Item -LiteralPath $Requested -Force
    if(-not $entry.PSIsContainer){throw '-Restart requires an old run or checkpoint directory.'}
    $sourcePath=$entry.FullName.TrimEnd([IO.Path]::DirectorySeparatorChar)
    $runs=[IO.Path]::GetFullPath((Join-Path $Project 'runs')).TrimEnd([IO.Path]::DirectorySeparatorChar)
    if(-not $sourcePath.StartsWith($runs+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)){
        throw '-Restart requires a saved run inside this project runs/ directory.'
    }
    Assert-NoRestartLinks $sourcePath
    $explicitCheckpoint=$entry.Name -cmatch '\Acheckpoint_[0-9]+\z'
    $parentRun=$entry
    if($explicitCheckpoint){
        while($null -ne $parentRun -and $parentRun.FullName.StartsWith($runs+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase) -and
            -not (Test-Path -LiteralPath (Join-Path $parentRun.FullName 'run.json') -PathType Leaf)){
            $parentRun=$parentRun.Parent
        }
    }
    if($null -eq $parentRun -or -not $parentRun.FullName.StartsWith($runs+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase) -or
        -not (Test-Path -LiteralPath (Join-Path $parentRun.FullName 'run.json') -PathType Leaf)){
        throw 'Restart source has no owning saved run.json.'
    }
    Assert-NoRestartLinks (Join-Path $parentRun.FullName 'run.json')
    $old=Get-Content -LiteralPath (Join-Path $parentRun.FullName 'run.json') -Raw -Encoding UTF8 | ConvertFrom-Json
    if(-not $old.case_name -or -not ([string]$old.case_name).Equals($ExpectedCase,[StringComparison]::OrdinalIgnoreCase)){
        throw 'Restart source case_name differs from the selected case.'
    }
    $checkpointRoot=Join-Path $parentRun.FullName 'output/checkpoints'
    if($explicitCheckpoint){
        if(-not $entry.Parent.FullName.Equals([IO.Path]::GetFullPath($checkpointRoot),[StringComparison]::OrdinalIgnoreCase)){
            throw 'Explicit restart checkpoint must belong to its run output/checkpoints/ directory.'
        }
        $checkpoint=$sourcePath
    }else{
        $latest=Join-Path $checkpointRoot 'latest.txt'
        if(-not (Test-Path -LiteralPath $latest -PathType Leaf)){throw 'Saved run has no completed checkpoint latest.txt.'}
        Assert-NoRestartLinks $latest
        $name=([IO.File]::ReadAllText($latest,[Text.Encoding]::UTF8)).Trim()
        if($name -cnotmatch '\Acheckpoint_[0-9]+\z'){throw 'Checkpoint latest.txt must contain only a checkpoint_<step> directory name.'}
        $checkpoint=Join-Path $checkpointRoot $name
        if(-not (Test-Path -LiteralPath $checkpoint -PathType Container)){throw 'Newest checkpoint directory is missing; no automatic rollback is performed.'}
        Assert-NoRestartLinks $checkpoint
    }
    $manifest=Join-Path $checkpoint 'manifest.txt'
    if(-not (Test-Path -LiteralPath $manifest -PathType Leaf)){throw 'Newest selected checkpoint manifest.txt is missing; no automatic rollback is performed.'}
    Assert-NoRestartLinks $manifest
    if(-not (Test-Path -LiteralPath (Join-Path $parentRun.FullName 'system/run.dict') -PathType Leaf)){
        throw 'Restart source has no frozen system/run.dict.'
    }
    return [ordered]@{source_run=$parentRun.FullName;source_case_name=[string]$old.case_name;
        requested=$sourcePath;checkpoint_directory=[IO.Path]::GetFullPath($checkpoint);
        checkpoint_manifest_sha256=(Get-CaseSha256 $manifest);explicit_checkpoint=[bool]$explicitCheckpoint}
}
function Read-RestartSource([string]$CaseName,[string]$Project){
    $root=Join-Path $Project ('runs/'+$CaseName.Replace('/','__'))
    $choices=@()
    if(Test-Path -LiteralPath $root -PathType Container){
        $choices=@(Get-ChildItem -LiteralPath $root -Directory |
            Where-Object { Test-Path -LiteralPath (Join-Path $_.FullName 'output/checkpoints/latest.txt') -PathType Leaf } |
            Sort-Object Name -Descending | Select-Object -First 10)
    }
    $previousInput=[Console]::InputEncoding;$previousOutput=[Console]::OutputEncoding
    $previousPipelineEncoding=$OutputEncoding;$promptEncoding=[Text.UTF8Encoding]::new($false)
    try{
        [Console]::InputEncoding=$promptEncoding;[Console]::OutputEncoding=$promptEncoding;$OutputEncoding=$promptEncoding
        [Console]::WriteLine('[RESTART] Saved runs with a checkpoint index for '+$CaseName+':')
        if($choices.Count -eq 0){[Console]::WriteLine('  No indexed runs found. You may paste an explicit saved checkpoint path.')}
        for($i=0;$i -lt $choices.Count;$i++){[Console]::WriteLine('  ['+($i+1)+'] '+$choices[$i].FullName)}
        [Console]::WriteLine('The selected checkpoint will be validated. No automatic fallback is performed.')
        [Console]::Write('Enter a number or saved run/checkpoint path (q cancels): ')
        $answer=[Console]::ReadLine()
    }finally{
        [Console]::InputEncoding=$previousInput;[Console]::OutputEncoding=$previousOutput;$OutputEncoding=$previousPipelineEncoding
    }
    if($null -eq $answer -or [string]::IsNullOrWhiteSpace($answer)){throw 'No restart source selected; use -Restart <path> for noninteractive execution.'}
    $answer=$answer.Trim()
    if($answer -ieq 'q'){return $null}
    if($answer -match '\A[0-9]+\z'){
        $selection=0
        if(-not [int]::TryParse($answer,[ref]$selection) -or $selection -lt 1 -or $selection -gt $choices.Count){
            throw 'Restart selection is outside the displayed range.'
        }
        return $choices[$selection-1].FullName
    }
    if($answer.Length -ge 2 -and (($answer[0] -eq '"' -and $answer[$answer.Length-1] -eq '"') -or
        ($answer[0] -eq "'" -and $answer[$answer.Length-1] -eq "'"))){$answer=$answer.Substring(1,$answer.Length-2)}
    return $answer
}
function Get-CaseDictionaryTokens([string]$Text){
    # The formal dictionary lexer has single-line quotes and # / // comments.
    # Track actual top-level brace tokens so an old Restart block can be
    # replaced without matching comments or a quoted path/formula.
    $pattern='#[^\r\n]*|//[^\r\n]*|"[^"\r\n]*"|''[^''\r\n]*''|[{}]|[^\s{}();,]+'
    return @([regex]::Matches($Text,$pattern) | Where-Object { -not ($_.Value.StartsWith('#') -or $_.Value.StartsWith('//')) })
}
function Set-RestartDirectory([string]$Text,[string]$Directory){
    if($Directory.IndexOfAny([char[]]@('"',"`r","`n",[char]0)) -ge 0){throw 'Restart path cannot contain quotes, newlines or NUL.'}
    $tokens=@(Get-CaseDictionaryTokens $Text)
    $depth=0;$begin=-1;$end=-1
    for($i=0;$i -lt $tokens.Count;$i++){
        $token=$tokens[$i]
        if($depth -eq 0 -and $token.Value -ceq 'Restart'){
            if($begin -ge 0){throw 'Frozen run.dict contains duplicate top-level Restart entries.'}
            if($i+1 -ge $tokens.Count -or $tokens[$i+1].Value -ne '{'){throw 'Frozen Restart must be a block.'}
            $begin=$token.Index
        }
        if($token.Value -eq '{'){$depth++}
        elseif($token.Value -eq '}'){
            $depth--;if($depth -lt 0){throw 'Frozen run.dict has unmatched braces.'}
            if($depth -eq 0 -and $begin -ge 0 -and $end -lt 0){$end=$token.Index+$token.Length}
        }
    }
    if($depth -ne 0 -or ($begin -ge 0 -and $end -lt 0)){throw 'Frozen run.dict has an incomplete block.'}
    $block='Restart { directory "'+$Directory.Replace('\','/')+'"; }'
    if($begin -ge 0){return $Text.Substring(0,$begin)+$block+$Text.Substring($end)}
    return $Text+"`n"+$block+"`n"
}
function Set-RestartEndTime([string]$Text,[double]$Value){
    $tokens=@(Get-CaseDictionaryTokens $Text)
    $depth=0;$timeStart=-1;$inTime=$false;$valueToken=$null
    for($i=0;$i -lt $tokens.Count;$i++){
        $token=$tokens[$i]
        if($depth -eq 0 -and $token.Value -ceq 'Time'){
            if($timeStart -ge 0){throw 'Frozen run.dict contains duplicate top-level Time entries.'}
            if($i+1 -ge $tokens.Count -or $tokens[$i+1].Value -ne '{'){throw 'Frozen Time must be a block.'}
            $timeStart=$i;$inTime=$true
        }
        if($depth -eq 1 -and $inTime -and $token.Value -ceq 'endTime'){
            if($null -ne $valueToken -or $i+1 -ge $tokens.Count){throw 'Frozen Time must contain exactly one endTime value.'}
            $valueToken=$tokens[$i+1]
        }
        if($token.Value -eq '{'){$depth++}
        elseif($token.Value -eq '}'){$depth--;if($depth -lt 0){throw 'Frozen run.dict has unmatched braces.'};if($depth -eq 0){$inTime=$false}}
    }
    if($depth -ne 0 -or $timeStart -lt 0 -or $null -eq $valueToken){throw 'Frozen run.dict has no complete Time.endTime value.'}
    $previous=0.0;$raw=$valueToken.Value
    if($raw.Length -ge 2 -and ($raw[0] -eq '"' -or $raw[0] -eq "'")){$raw=$raw.Substring(1,$raw.Length-2)}
    if(-not [double]::TryParse($raw,[Globalization.NumberStyles]::Float,[Globalization.CultureInfo]::InvariantCulture,[ref]$previous) -or
        [double]::IsNaN($previous) -or [double]::IsInfinity($previous) -or $previous -le 0){throw 'Frozen Time.endTime is not a finite positive code time.'}
    $replacement=$Value.ToString('R',[Globalization.CultureInfo]::InvariantCulture)
    return [pscustomobject]@{Text=$Text.Substring(0,$valueToken.Index)+$replacement+$Text.Substring($valueToken.Index+$valueToken.Length);
        OldValue=$previous;NewValue=$Value}
}
$project=Split-Path -Parent $PSScriptRoot
$casesRoot=[IO.Path]::GetFullPath((Join-Path $project 'cases'))
$separator=[IO.Path]::DirectorySeparatorChar
if($List){
    Get-ChildItem -LiteralPath $casesRoot -Filter run.dict -File -Recurse |
        Where-Object { $_.Directory.Name -eq 'system' -and $_.FullName -notmatch '[\\/]output[\\/]' -and
            -not (Test-Path -LiteralPath (Join-Path $_.Directory.Parent.FullName 'retired.json')) } |
        ForEach-Object { $_.Directory.Parent.FullName.Substring($casesRoot.Length+1).Replace('\','/') } |
        Sort-Object
    Get-ChildItem -LiteralPath $casesRoot -Filter benchmark.json -File -Recurse |
        Where-Object { -not (Test-Path -LiteralPath (Join-Path $_.Directory.FullName 'retired.json')) } |
        ForEach-Object { $_.Directory.FullName.Substring($casesRoot.Length+1).Replace('\','/') } |
        Sort-Object
    return
}
if(-not $Case){throw 'Usage: tools\run_case.cmd <case-name> [-Resume | -Restart <old-run-or-checkpoint>] [-EndTime <code-time>]; use -Help or -List.'}
if($Resume -and $PSBoundParameters.ContainsKey('Restart')){throw 'Use either -Resume or -Restart <path>.'}
if($PSBoundParameters.ContainsKey('Restart') -and [string]::IsNullOrWhiteSpace($Restart)){throw '-Restart requires a nonempty saved run or checkpoint directory.'}
$restartEndTime=0.0
if($endTimeExplicit){
    if(-not $Restart -and -not $Resume){throw '-EndTime is only supported together with -Restart or -Resume.'}
    if(-not [double]::TryParse($EndTime,[Globalization.NumberStyles]::Float,[Globalization.CultureInfo]::InvariantCulture,[ref]$restartEndTime) -or
        [double]::IsNaN($restartEndTime) -or [double]::IsInfinity($restartEndTime) -or $restartEndTime -le 0){
        throw '-EndTime requires a finite positive code time.'
    }
}
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
$retirementPath=Join-Path $casePath 'retired.json'
if(Test-Path -LiteralPath $retirementPath -PathType Leaf){
    $retirement=Get-Content -LiteralPath $retirementPath -Raw -Encoding UTF8 | ConvertFrom-Json
    throw "This case is retired: $($retirement.reason)"
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
$restartSource=$null;$inputCasePath=$casePath
if($Resume){
    if($benchmark){throw '-Resume is supported only by formal Hydro/A VOF cases, not benchmark drivers.'}
    $selectedRestart=Read-RestartSource $caseName $project
    if($null -eq $selectedRestart){Write-Output 'Restart cancelled; no run directory created.';return}
    $Restart=[string]$selectedRestart
}
if($Restart){
    if($benchmark){throw '-Restart is supported only by formal Hydro/A VOF cases, not benchmark drivers.'}
    $restartSource=Resolve-RestartSource $Restart $caseName $project
    $inputCasePath=$restartSource.source_run
}
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
    $source=Join-Path $inputCasePath $folder
    if(Test-Path -LiteralPath $source){
        $linked=@(Get-Item -LiteralPath $source; Get-ChildItem -LiteralPath $source -Recurse -Force) |
            Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint }
        if($linked){throw "Input links are not supported: $source"}
        Copy-Item -LiteralPath $source -Destination $run -Recurse
    }
}
$utf8=[Text.UTF8Encoding]::new($false)
if($restartSource){
    $runDictionary=Join-Path $run 'system/run.dict'
    $frozen=[IO.File]::ReadAllText($runDictionary,[Text.Encoding]::UTF8)
    $frozen=Set-RestartDirectory $frozen $restartSource.checkpoint_directory
    if($endTimeExplicit){
        $changed=Set-RestartEndTime $frozen $restartEndTime;$frozen=$changed.Text
        $restartSource['end_time_override']=[ordered]@{old=$changed.OldValue;new=$changed.NewValue;units='code_time'}
    }
    [IO.File]::WriteAllText($runDictionary,$frozen,$utf8)
}
$runtime=Join-Path $run 'runtime'
New-Item -ItemType Directory -Path $runtime | Out-Null
$exe=Join-Path $runtime ([IO.Path]::GetFileName($Program))
Copy-Item -LiteralPath $Program -Destination $exe
Get-ChildItem -LiteralPath (Split-Path -Parent $Program) -Filter '*.dll' -File | ForEach-Object {Copy-Item -LiteralPath $_.FullName -Destination $runtime}
$record=[ordered]@{
    case_name=$caseName;case=$casePath;run=$run;created=(Get-Date).ToString('o');status='prepared'
    check_only=[bool]$CheckOnly;skip_check=[bool]$SkipCheck
    requested_threads=$(if($threadsExplicit){$Threads}else{$null});threads=$(if($benchmark -or $threadsExplicit){$Threads}else{$null})
    threads_source=$(if($threadsExplicit){'command_line'}elseif($benchmark){'default'}else{'pending_runtime'})
    kind=$(if($benchmark){'benchmark'}else{'case'});arguments=$(if($benchmark){@($benchmark.arguments)}elseif($threadsExplicit){@('--case',$run,'--threads',[string]$Threads)}else{@('--case',$run)})
    program=$exe;program_sha256=(Get-CaseSha256 $exe)
    restart=$restartSource;input_source=$inputCasePath
    inputs=@(foreach($folder in $inputFolders){
        $path=Join-Path $run $folder
        if(Test-Path -LiteralPath $path){Get-ChildItem -LiteralPath $path -File -Recurse | Sort-Object FullName | ForEach-Object {
            [ordered]@{path=$_.FullName.Substring($run.Length+1).Replace('\','/');sha256=(Get-CaseSha256 $_.FullName)}
        }}
    })
}
if($benchmark){$record.inputs+=@([ordered]@{path='benchmark.json';sha256=(Get-CaseSha256 (Join-Path $run 'benchmark.json'))})}
function Save-Record {[IO.File]::WriteAllText((Join-Path $run 'run.json'),($record|ConvertTo-Json -Depth 6),$utf8)}
function Invoke-Program([string]$Log,[switch]$InputCheck,[switch]$QuietConsole){
    $arguments=@('--case',$run)
    if($benchmark){$arguments=@($benchmark.arguments)}
    if(-not $benchmark -and $threadsExplicit){$arguments+=@('--threads',[string]$Threads)}
    if($InputCheck){$arguments+='--check-case'}
    $previousPreference=$ErrorActionPreference
    $previousConsoleInput=[Console]::InputEncoding
    $previousConsoleOutput=[Console]::OutputEncoding
    $previousOutputEncoding=$OutputEncoding
    $writer=[IO.StreamWriter]::new($Log,$false,$utf8)
    $writer.AutoFlush=$true
    try{
        # Native stderr is merged into the same live stream. Console writes do
        # not enter this function's return pipeline: its sole result is the exit code.
        # Windows PowerShell uses Console.OutputEncoding to decode native
        # pipes. Match LAMP's UTF-8 before it starts, including stderr/logs.
        [Console]::InputEncoding=$utf8
        [Console]::OutputEncoding=$utf8
        $OutputEncoding=$utf8
        $ErrorActionPreference='Continue'
        $global:LASTEXITCODE=$null
        & $exe @arguments 2>&1 | ForEach-Object {
            $line=[string]$_
            if($line -match '^\[Runtime\] threads=(\d+) requested_threads=(\d+) source=([A-Za-z_.]+)$'){
                $record.threads=[int]$Matches[1];$record.threads_source=$Matches[3]
            }
            $writer.WriteLine($line)
            if(-not $QuietConsole){[Console]::WriteLine($line)}
        }
        if($null -eq $global:LASTEXITCODE){throw 'The program did not return an exit code.'}
        return $global:LASTEXITCODE
    }finally{
        $writer.Dispose()
        $ErrorActionPreference=$previousPreference
        [Console]::InputEncoding=$previousConsoleInput
        [Console]::OutputEncoding=$previousConsoleOutput
        $OutputEncoding=$previousOutputEncoding
    }
}
Save-Record
[IO.File]::WriteAllText((Join-Path $runRoot 'latest.txt'),$run,$utf8)
Write-Output "Run directory: $run"
if($PrepareOnly){Write-Output 'Prepared inputs and program; no physical steps.';return}
$saved=@{}
foreach($key in @('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OMP_DYNAMIC','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')){$saved[$key]=[Environment]::GetEnvironmentVariable($key,'Process')}
if($benchmark){
    $env:OMP_NUM_THREADS=[string]$Threads;$env:OMP_THREAD_LIMIT=[string]$Threads
}else{
    # The formal Runtime resolves Parallel.threads or --threads before OpenMP
    # initializes. Inherited limits must not silently cap the resolved input.
    [Environment]::SetEnvironmentVariable('OMP_NUM_THREADS',$null,'Process')
    [Environment]::SetEnvironmentVariable('OMP_THREAD_LIMIT',$null,'Process')
}
$env:OMP_DYNAMIC='FALSE';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
$watch=[Diagnostics.Stopwatch]::StartNew()
try{
    if($benchmark -and $CheckOnly){$record.status='input_checked';$record.exit_code=0;Write-Output 'Benchmark descriptor and runtime snapshot checked; no physical steps.';return}
    if(-not $SkipCheck -and -not $benchmark){
        $record.status='checking';Save-Record
        $log=Join-Path $run 'check.log'
        Write-Output '[CHECK] Validating inputs...'
        # Automatic preflight keeps its full log without repeating the App
        # startup on screen. Explicit -CheckOnly remains directly observable.
        $code=Invoke-Program $log -InputCheck -QuietConsole:(!$CheckOnly)
        $record.check_exit_code=$code
        if($code -ne 0){
            $record.status='input_failed';$record.exit_code=$code
            if(-not $CheckOnly){
                foreach($line in [IO.File]::ReadLines($log,$utf8)){[Console]::WriteLine($line)}
            }
            Write-Output "Input check failed ($code): $log"
            exit $code
        }
        if($CheckOnly){$record.status='input_checked';$record.exit_code=0;Write-Output "Input check OK: $log";return}
        Write-Output '[CHECK] OK'
    }
    $record.status='running';Save-Record
    $log=Join-Path $run 'run.log'
    Write-Output "Running; log: $log"
    Push-Location $run
    try {$code=Invoke-Program $log} finally {Pop-Location}
    $record.status=if($code -eq 0){'completed'}else{'failed'}
    $record.exit_code=$code
    if($code -ne 0){Write-Output "Run failed ($code): $log";exit $code}
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
