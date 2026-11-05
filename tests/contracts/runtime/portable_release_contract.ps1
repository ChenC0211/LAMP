<#
.SYNOPSIS
Exercise selected-case release orchestration without building a solver.
.DESCRIPTION
The package consumer is an isolated fixture. A build trap ensures compilation
is never invoked. The real package acceptance
owns PE dependencies, ZIP contents and solver/restart verification.
#>
[CmdletBinding()]
param(
    [string]$Source='',
    [string]$WorkDirectory=''
)
$ErrorActionPreference='Stop'
if(-not $Source){$Source=Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $PSScriptRoot))}
$Source=(Get-Item -LiteralPath $Source).FullName
if(-not $WorkDirectory){$WorkDirectory=Join-Path $Source 'runs/checks/portable_release_contract'}
$WorkDirectory=[IO.Path]::GetFullPath($WorkDirectory)
if($WorkDirectory -notmatch '\AD:[\\/]'){throw 'Contract outputs must be on D:.'}
$work=Join-Path $WorkDirectory ((Get-Date -Format 'yyyyMMdd_HHmmss_fff')+'_'+[Guid]::NewGuid().ToString('N').Substring(0,8))
[void][IO.Directory]::CreateDirectory($work)
$utf8=[Text.UTF8Encoding]::new($false)
$powershell=Join-Path $env:SystemRoot 'System32/WindowsPowerShell/v1.0/powershell.exe'
$results=[Collections.Generic.List[object]]::new()
function Require([bool]$Condition,[string]$Message){if(-not $Condition){throw $Message}}
function Write-Fixture([string]$Path,[string]$Text){
    [void][IO.Directory]::CreateDirectory((Split-Path -Parent $Path))
    [IO.File]::WriteAllText($Path,$Text,$utf8)
}
function Get-Sha([string]$Path){return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()}
function New-Fixture([string]$Name){
    $root=Join-Path $work ($Name+' moved '+[char]0x6db2+[char]0x6ef4+' project')
    [void][IO.Directory]::CreateDirectory((Join-Path $root 'tools'))
    foreach($file in @('release_portable.ps1','release_portable.cmd')){
        Copy-Item -LiteralPath (Join-Path $Source ('tools/'+$file)) -Destination (Join-Path $root ('tools/'+$file))
    }
    Write-Fixture (Join-Path $root 'bin/LAMP.exe') 'existing compiled executable fixture'
    Write-Fixture (Join-Path $root 'tools/build.ps1') @'
$root=Split-Path -Parent $PSScriptRoot
$event=[ordered]@{consumer='build'}
[IO.File]::AppendAllText((Join-Path $root 'events.jsonl'),(($event|ConvertTo-Json -Compress)+[Environment]::NewLine))
throw 'Compilation is forbidden by this packaging contract.'
'@
    Write-Fixture (Join-Path $root 'tools/package_portable.ps1') @'
param([string]$Program,[string]$OutputDirectory,[string]$Version,[string[]]$Cases,
 [string]$OpenMpRuntime,[string]$OpenMpLicenseDirectory,[string]$VcRuntimeDirectory,
 [string]$VcRuntimeLicensePath,[string]$AdditionalLicenseDirectory)
$root=Split-Path -Parent $PSScriptRoot
$event=[ordered]@{consumer='package';program=$Program;output=$OutputDirectory;version=$Version;
 cases=@($Cases);cases_explicit=$PSBoundParameters.ContainsKey('Cases');
 omp=$OpenMpRuntime;omp_licenses=$OpenMpLicenseDirectory;crt=$VcRuntimeDirectory;
 crt_license=$VcRuntimeLicensePath;third_party=$AdditionalLicenseDirectory;
 program_text=[IO.File]::ReadAllText($Program)}
[IO.File]::AppendAllText((Join-Path $root 'events.jsonl'),(($event|ConvertTo-Json -Compress)+[Environment]::NewLine))
$modeFile=Join-Path $root 'package.mode'
$mode=if(Test-Path -LiteralPath $modeFile){[IO.File]::ReadAllText($modeFile)}else{'success'}
if($mode -eq 'exit9'){exit 9}
$destination=Join-Path $OutputDirectory ('LAMP-Portable-'+$Version)
[void][IO.Directory]::CreateDirectory($destination)
if($mode -ne 'missing-manifest'){
 [IO.File]::WriteAllText((Join-Path $destination 'portable_manifest.json'),'{"fixture":true}')
}
if($mode -ne 'missing-zip'){[IO.File]::WriteAllText(($destination+'.zip'),'fixture archive marker')}
Write-Output 'fixture package published'
'@
    foreach($case in @('alpha','beta','gamma','retired','bench','unsupported','z_electrospray')){
        $app=if($case -eq 'unsupported'){'LegacyApp'}elseif($case -eq 'gamma'){'Electric'}elseif($case -eq 'z_electrospray'){'Electrospray'}else{'Hydrodynamics'}
        Write-Fixture (Join-Path $root ('cases/'+$case+'/system/run.dict')) ('app '+$app+';')
    }
    Write-Fixture (Join-Path $root 'cases/retired/retired.json') '{}'
    Write-Fixture (Join-Path $root 'cases/bench/benchmark.json') '{}'
    Write-Fixture (Join-Path $root 'runs/history/do-not-package.txt') 'old numerical output fixture'
    Write-Fixture (Join-Path $root 'build/do-not-package.txt') 'compiler object fixture'
    $runtime=Join-Path $root '.portable/runtime'
    foreach($relative in @('bin/libomp140.x86_64.dll','bin/msvcp140.dll','bin/vcruntime140.dll','bin/vcruntime140_1.dll',
        'licenses/CRT/LICENSE-VCRUNTIME.rtf','licenses/OpenMP/LICENSE.TXT','licenses/ThirdParty/NOTICES.txt')){
        Write-Fixture (Join-Path $runtime $relative) ('runtime fixture: '+$relative)
    }
    $entries=@(Get-ChildItem -LiteralPath $runtime -File -Recurse|ForEach-Object {
        [ordered]@{path=$_.FullName.Substring($runtime.Length+1).Replace('\','/');sha256=(Get-Sha $_.FullName);bytes=$_.Length}
    })
    Write-Fixture (Join-Path $runtime 'runtime-cache.json') (([ordered]@{schema_version=1;files=$entries})|ConvertTo-Json -Depth 6)
    [void][IO.Directory]::CreateDirectory((Join-Path $root 'unrelated caller'))
    return [pscustomobject]@{Root=$root;Runtime=$runtime;Output=(Join-Path $root 'dist');Version=$Name}
}
function Quote-Argument([string]$Value){
    if($Value.Contains('"')){throw 'Fixture argument contains an unexpected quote.'}
    return '"'+$Value+'"'
}
function Invoke-Fixture($Fixture,[string[]]$Arguments=@(),[AllowNull()][string]$InputText=$null,[switch]$Cmd){
    $program=Join-Path $Fixture.Root 'bin/LAMP.exe'
    $originalSha=if(Test-Path -LiteralPath $program -PathType Leaf){Get-Sha $program}else{$null}
    $info=[Diagnostics.ProcessStartInfo]::new()
    if($Cmd){
        $info.FileName=Join-Path $env:SystemRoot 'System32/cmd.exe'
        $command=(Quote-Argument (Join-Path $Fixture.Root 'tools/release_portable.cmd'))
        if($Arguments.Count){$command+=' '+(($Arguments|ForEach-Object {Quote-Argument $_}) -join ' ')}
        $info.Arguments='/d /s /c "'+$command+'"'
    }else{
        $info.FileName=$powershell
        $all=@('-NoProfile','-ExecutionPolicy','Bypass','-File',(Join-Path $Fixture.Root 'tools/release_portable.ps1'))+$Arguments
        $info.Arguments=($all|ForEach-Object {Quote-Argument $_}) -join ' '
    }
    $info.WorkingDirectory=Join-Path $Fixture.Root 'unrelated caller'
    $info.UseShellExecute=$false;$info.CreateNoWindow=$true
    $info.RedirectStandardInput=$true;$info.RedirectStandardOutput=$true;$info.RedirectStandardError=$true
    $info.EnvironmentVariables['PATH']=$env:SystemRoot+';'+(Join-Path $env:SystemRoot 'System32')
    $temp=Join-Path $Fixture.Root 'temp';[void][IO.Directory]::CreateDirectory($temp)
    $info.EnvironmentVariables['TEMP']=$temp;$info.EnvironmentVariables['TMP']=$temp
    $process=[Diagnostics.Process]::new();$process.StartInfo=$info
    try{
        [void]$process.Start()
        $outTask=$process.StandardOutput.ReadToEndAsync();$errTask=$process.StandardError.ReadToEndAsync()
        if($null -ne $InputText){$process.StandardInput.Write($InputText)}
        $process.StandardInput.Close()
        if(-not $process.WaitForExit(30000)){$process.Kill();throw 'Release command did not terminate within 30 seconds.'}
        $result=[pscustomobject]@{exit=$process.ExitCode;stdout=$outTask.Result;stderr=$errTask.Result}
    }finally{$process.Dispose()}
    if($originalSha){Require ((Test-Path -LiteralPath $program -PathType Leaf) -and (Get-Sha $program) -eq $originalSha) 'Packaging modified the existing executable.'}
    Write-Fixture (Join-Path $Fixture.Root 'process-result.json') ($result|ConvertTo-Json -Depth 4)
    return $result
}
function Read-Events($Fixture){
    $path=Join-Path $Fixture.Root 'events.jsonl'
    return @(if(Test-Path -LiteralPath $path){[IO.File]::ReadAllLines($path)|Where-Object {$_}|ForEach-Object {$_|ConvertFrom-Json}})
}
function Assert-Cases($Event,[string[]]$Expected){
    Require (($Event.cases -join ',') -ceq ($Expected -join ',')) ('Unexpected cases: '+($Event.cases -join ','))
    Require $Event.cases_explicit 'Packager did not receive an explicit Cases selection.'
}
function Read-Record($Fixture){
    $records=@(Get-ChildItem -LiteralPath (Join-Path $Fixture.Root 'runs/packaging') -Filter 'release.json' -Recurse)
    Require ($records.Count -eq 1) 'Expected one durable release record.'
    return [IO.File]::ReadAllText($records[0].FullName)|ConvertFrom-Json
}
function Complete-Case([string]$Name,$Result){
    $results.Add([ordered]@{name=$Name;passed=$true;exit=$Result.exit})
    Write-Output ('PASS '+$Name)
}

try{
    $fixture=New-Fixture 'default_package'
    $result=Invoke-Fixture $fixture @('-Cases','gamma,alpha','-Version',$fixture.Version)
    Require ($result.exit -eq 0) ('Default release failed: '+$result.stderr)
    $events=@(Read-Events $fixture)
    Require ($events.Count -eq 1 -and $events[0].consumer -eq 'package') 'Default packaging invoked compilation.'
    Assert-Cases $events[0] @('gamma','alpha')
    Require ($events[0].program_text -eq 'existing compiled executable fixture') 'Packager did not consume the existing executable.'
    Require ($events[0].omp -eq (Join-Path $fixture.Runtime 'bin/libomp140.x86_64.dll')) 'Runtime cache was not delegated.'
    $record=Read-Record $fixture
    Require ($record.status -eq 'completed' -and $record.mode -eq 'existing_binary' -and ($record.cases -join ',') -eq 'gamma,alpha') 'Completed record lost existing-binary/selection provenance.'
    Require (@($record.PSObject.Properties.Name|Where-Object {$_ -in @('skip_build','preset','jobs','build_directory')}).Count -eq 0) 'Packaging record incorrectly reports build settings.'
    Require ($record.program_sha256 -eq (Get-Sha (Join-Path $fixture.Root 'bin/LAMP.exe'))) 'Release record executable hash differs from packaged input.'
    Complete-Case 'default packaging preserves existing executable and explicit selection' $result

    $fixture=New-Fixture 'electrospray_explicit'
    $result=Invoke-Fixture $fixture @('-Cases','z_electrospray','-Version',$fixture.Version)
    Require ($result.exit -eq 0) ('Explicit Electrospray release failed: '+$result.stderr)
    $events=@(Read-Events $fixture)
    Require ($events.Count -eq 1 -and $events[0].consumer -eq 'package') 'Electrospray selection invoked compilation or did not delegate exactly once.'
    Assert-Cases $events[0] @('z_electrospray')
    $record=Read-Record $fixture
    Require ($record.status -eq 'completed' -and $record.mode -eq 'existing_binary' -and ($record.cases -join ',') -eq 'z_electrospray') 'Completed release lost the explicit Electrospray selection.'
    Complete-Case 'Electrospray explicit selection packages only the chosen case' $result

    $fixture=New-Fixture 'path_overrides'
    $alternative=Join-Path $fixture.Root 'alternative runtime'
    Copy-Item -LiteralPath $fixture.Runtime -Destination $alternative -Recurse
    $output=Join-Path $fixture.Root 'chosen output'
    $result=Invoke-Fixture $fixture @('-Cases','beta','-RuntimeDirectory',$alternative,'-OutputDirectory',$output,'-Version','selected-only')
    Require ($result.exit -eq 0) ('Selected release failed: '+$result.stderr)
    $events=@(Read-Events $fixture)
    Require ($events.Count -eq 1 -and $events[0].consumer -eq 'package') 'Packaging invoked the build consumer.'
    Assert-Cases $events[0] @('beta')
    Require ($events[0].program_text -eq 'existing compiled executable fixture' -and $events[0].output -eq $output -and $events[0].version -eq 'selected-only') 'Explicit output/version or current executable was lost.'
    Require ($events[0].crt -eq (Join-Path $alternative 'bin')) 'Explicit runtime directory was ignored.'
    Require (Test-Path -LiteralPath (Join-Path $output 'LAMP-Portable-selected-only.zip')) 'Selected package was not published in requested output directory.'
    Complete-Case 'only selected case and override paths' $result

    $fixture=New-Fixture 'interactive_range'
    $result=Invoke-Fixture $fixture @('-Version',$fixture.Version) "3,1-2,1`r`n"
    Require ($result.exit -eq 0) ('Interactive release failed: '+$result.stderr)
    $events=@(Read-Events $fixture)
    Require ($events.Count -eq 1 -and $events[0].consumer -eq 'package') 'Interactive selection did not package without compilation.'
    Assert-Cases $events[0] @('gamma','alpha','beta')
    Require ($result.stdout -match '(?m)^\s*\[4\] z_electrospray\s*$') 'Electrospray was absent from the supported selection menu or changed the original case numbering.'
    Require ($result.stdout -notmatch '\[\d+\] (retired|bench|unsupported)') 'Unsupported/retired cases appeared in the selection menu.'
    Complete-Case 'interactive CSV ranges preserve order and deduplicate' $result

    foreach($cancel in @('q','eof','blank')){
        $fixture=New-Fixture ('cancel_'+$cancel)
        $taskInput=if($cancel -eq 'q'){"q`r`n"}elseif($cancel -eq 'blank'){"`r`n"}else{$null}
        $result=Invoke-Fixture $fixture @() $taskInput
        Require ($result.exit -eq 0 -and @(Read-Events $fixture).Count -eq 0) ('Cancellation started work: '+$cancel)
        Require (-not (Test-Path -LiteralPath (Join-Path $fixture.Root 'runs/packaging'))) 'Cancellation created a release transaction.'
        Complete-Case ('cancel '+$cancel+' before packaging') $result
    }

    foreach($bad in @('unknown','retired','bench','unsupported','0','2-9','3-1','')){
        $fixture=New-Fixture ('bad_case_'+[Guid]::NewGuid().ToString('N').Substring(0,6))
        $result=Invoke-Fixture $fixture @('-Cases',$bad,'-Version',$fixture.Version)
        Require ($result.exit -ne 0 -and @(Read-Events $fixture).Count -eq 0) ('Invalid explicit selection started work: '+$bad)
        Complete-Case ('reject explicit selection ['+$bad+']') $result
    }

    $fixture=New-Fixture 'missing_executable'
    Remove-Item -LiteralPath (Join-Path $fixture.Root 'bin/LAMP.exe')
    $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version)
    Require ($result.exit -ne 0 -and @(Read-Events $fixture).Count -eq 0) 'Missing executable invoked compilation or packaging.'
    Require ($result.stderr -match 'LAMP.exe.*missing') 'Missing executable did not produce a useful diagnostic.'
    Complete-Case 'missing executable refuses packaging' $result

    $fixture=New-Fixture 'no_build_environment'
    Remove-Item -LiteralPath (Join-Path $fixture.Root 'tools/build.ps1')
    $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version)
    $events=@(Read-Events $fixture)
    Require ($result.exit -eq 0 -and $events.Count -eq 1 -and $events[0].consumer -eq 'package') 'Packaging required build scripts, CMake presets or build cache.'
    Complete-Case 'no build script presets or cache required' $result

    $fixture=New-Fixture 'ignores_build_cache'
    Write-Fixture (Join-Path $fixture.Root 'CMakeUserPresets.json') 'invalid unused preset fixture'
    Write-Fixture (Join-Path $fixture.Root 'build/windows-local/CMakeCache.txt') "LAMP_RUNTIME_DIRECTORY:PATH=D:\unrelated-output`r`n"
    $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version)
    $events=@(Read-Events $fixture)
    Require ($result.exit -eq 0 -and $events.Count -eq 1 -and $events[0].consumer -eq 'package') 'Packaging incorrectly depended on CMake presets/cache.'
    Complete-Case 'existing executable independent of presets/cache' $result

    foreach($failure in @('corrupt','unlisted','missing','traversal')){
        $fixture=New-Fixture ('cache_'+$failure)
        $manifestPath=Join-Path $fixture.Runtime 'runtime-cache.json'
        if($failure -eq 'corrupt'){Write-Fixture (Join-Path $fixture.Runtime 'bin/vcruntime140.dll') 'tampered runtime'}
        elseif($failure -eq 'unlisted'){Write-Fixture (Join-Path $fixture.Runtime 'bin/unlisted.dll') 'unlisted runtime'}
        elseif($failure -eq 'missing'){
            $inventory=[IO.File]::ReadAllText($manifestPath)|ConvertFrom-Json
            $inventory.files=@($inventory.files|Where-Object {$_.path -ne 'bin/vcruntime140.dll'})
            Write-Fixture $manifestPath ($inventory|ConvertTo-Json -Depth 6)
        }else{
            $inventory=[IO.File]::ReadAllText($manifestPath)|ConvertFrom-Json
            $inventory.files[0].path='../escape.dll'
            Write-Fixture $manifestPath ($inventory|ConvertTo-Json -Depth 6)
        }
        $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version)
        Require ($result.exit -ne 0 -and @(Read-Events $fixture).Count -eq 0) ('Invalid runtime cache started packaging: '+$failure)
        Complete-Case ('cache '+$failure+' rejected before packaging') $result
    }

    foreach($failure in @('exit9','missing-zip','missing-manifest')){
        $fixture=New-Fixture ('package_'+$failure)
        Write-Fixture (Join-Path $fixture.Root 'package.mode') $failure
        $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version)
        Require ($result.exit -ne 0) ('Incomplete/failed packager reported success: '+$failure)
        Require ((Read-Record $fixture).status -eq 'failed') 'Packaging failure did not leave a failed transaction record.'
        Complete-Case ('packager '+$failure+' retained as failure') $result
    }

    $fixture=New-Fixture 'existing_release'
    Write-Fixture (Join-Path $fixture.Output ('LAMP-Portable-'+$fixture.Version+'.zip')) 'existing archive must survive'
    $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version)
    Require ($result.exit -ne 0 -and @(Read-Events $fixture).Count -eq 0) 'Existing version was overwritten.'
    Require ([IO.File]::ReadAllText((Join-Path $fixture.Output ('LAMP-Portable-'+$fixture.Version+'.zip'))) -eq 'existing archive must survive') 'Existing archive was modified.'
    Complete-Case 'existing version remains immutable' $result

    $fixture=New-Fixture 'cmd_relative'
    $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version) $null -Cmd
    Require ($result.exit -eq 0) ('CMD entry failed from another working directory: '+$result.stderr)
    $events=@(Read-Events $fixture);Require ($events.Count -eq 1) 'CMD entry did not delegate exactly once.'
    Assert-Cases $events[0] @('alpha')
    Complete-Case 'CMD relative path and clean PATH' $result

    $fixture=New-Fixture 'cmd_failure'
    Write-Fixture (Join-Path $fixture.Root 'package.mode') 'exit9'
    $result=Invoke-Fixture $fixture @('-Cases','alpha','-Version',$fixture.Version) $null -Cmd
    $events=@(Read-Events $fixture)
    Require ($result.exit -ne 0 -and $events.Count -eq 1 -and $events[0].consumer -eq 'package') 'CMD entry lost packaging failure or invoked compilation.'
    Complete-Case 'CMD preserves failure exit status' $result

    $fixture=New-Fixture 'help'
    $result=Invoke-Fixture $fixture @('-Help') $null -Cmd
    Require ($result.exit -eq 0 -and $result.stdout -match 'existing executable' -and $result.stdout -match 'never invokes' -and @(Read-Events $fixture).Count -eq 0) 'Help started work or omitted the existing binary workflow.'
    Complete-Case 'help performs no packaging' $result

    $summary=[ordered]@{passed=$true;checks=$results.Count;scope='orchestration fixtures only; real PE/ZIP/solver checks are separate';results=$results}
    Write-Fixture (Join-Path $work 'verification.json') ($summary|ConvertTo-Json -Depth 7)
    Write-Output ('Portable release contract passed '+$results.Count+' checks. Evidence: '+$work)
}catch{
    $summary=[ordered]@{passed=$false;checks=$results.Count;error=$_.Exception.Message;results=$results}
    Write-Fixture (Join-Path $work 'verification.json') ($summary|ConvertTo-Json -Depth 7)
    throw
}
