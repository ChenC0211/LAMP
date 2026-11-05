<#
.SYNOPSIS
Start a case from a movable LAMP portable distribution.
.EXAMPLE
.\Start.cmd -Case vof_a_drop_deformation_amr -Action resume -NoPause
#>
[CmdletBinding()]
param(
    [string]$Case='',
    [ValidateSet('','new','resume','check')][string]$Action='',
    [string]$Restart='',
    [string]$EndTime='',
    [ValidateRange(1,256)][int]$Threads,
    [switch]$NoPause,
    [Alias('h')][switch]$Help
)
$ErrorActionPreference='Stop'
Set-StrictMode -Version 2.0

function Get-PortableSha256([string]$Path){
    $stream=[IO.File]::OpenRead($Path)
    $algorithm=[Security.Cryptography.SHA256]::Create()
    try { return [BitConverter]::ToString($algorithm.ComputeHash($stream)).Replace('-','').ToLowerInvariant() }
    finally { $algorithm.Dispose(); $stream.Dispose() }
}
function Assert-PortableAssets([string]$Root){
    $manifestPath=Join-Path $Root 'portable_manifest.json'
    if(-not (Test-Path -LiteralPath $manifestPath -PathType Leaf)){
        throw 'portable_manifest.json is missing. Extract the complete portable ZIP before starting.'
    }
    $manifest=[IO.File]::ReadAllText($manifestPath,[Text.Encoding]::UTF8) | ConvertFrom-Json
    if($manifest.schema_version -ne 1 -or $null -eq $manifest.files){ throw 'Unsupported portable manifest schema.' }
    $seen=@{}
    $prefix=$Root.TrimEnd([IO.Path]::DirectorySeparatorChar)+[IO.Path]::DirectorySeparatorChar
    foreach($entry in $manifest.files){
        $relative=[string]$entry.path
        if(-not ($relative -match '\Abin/' -or $relative -in @('tools/run_case.ps1','tools/run_case.cmd'))){ continue }
        if($seen.ContainsKey($relative)){ throw "Duplicate executable asset in manifest: $relative" }
        $seen[$relative]=$true
        if([IO.Path]::IsPathRooted($relative) -or $relative.Contains('\') -or $relative -match '(^|/)\.\.(/|$)'){
            throw "Invalid executable asset path: $relative"
        }
        $path=[IO.Path]::GetFullPath((Join-Path $Root $relative))
        if(-not $path.StartsWith($prefix,[StringComparison]::OrdinalIgnoreCase)){
            throw "Executable asset escapes the portable directory: $relative"
        }
        if(-not (Test-Path -LiteralPath $path -PathType Leaf)){ throw "Required portable executable asset is missing: $relative" }
        if([string]$entry.sha256 -notmatch '\A[0-9a-fA-F]{64}\z' -or
            (Get-PortableSha256 $path) -ine [string]$entry.sha256){
            throw "Portable executable asset is damaged or changed: $relative. Extract the complete ZIP again."
        }
    }
    foreach($required in @('bin/LAMP.exe','tools/run_case.ps1','tools/run_case.cmd')){
        if(-not $seen.ContainsKey($required)){ throw "Portable manifest has no required executable asset: $required" }
    }
    # Reject unlisted runtime libraries too: the DLL loader can select these
    # even if the listed runtime files themselves have matching hashes.
    foreach($asset in Get-ChildItem -LiteralPath (Join-Path $Root 'bin') -File){
        if($asset.Extension -in @('.exe','.dll') -and -not $seen.ContainsKey('bin/'+$asset.Name)){
            throw "Portable runtime has an unlisted executable asset: bin/$($asset.Name)"
        }
    }
}
function Read-PortableSelection([string]$Prompt){
    [Console]::Write($Prompt)
    $answer=[Console]::ReadLine()
    if($null -eq $answer -or [string]::IsNullOrWhiteSpace($answer) -or $answer.Trim() -ieq 'q'){ return $null }
    return $answer.Trim()
}
function Get-PortableCases([string]$Root){
    $casesRoot=Join-Path $Root 'cases'
    if(-not (Test-Path -LiteralPath $casesRoot -PathType Container)){ throw 'The portable cases/ directory is missing.' }
    $prefix=$casesRoot.TrimEnd([IO.Path]::DirectorySeparatorChar)+[IO.Path]::DirectorySeparatorChar
    $names=@(Get-ChildItem -LiteralPath $casesRoot -Recurse -File -Filter 'run.dict' |
        Where-Object { $_.Directory.Name -eq 'system' -and
            -not (Test-Path -LiteralPath (Join-Path $_.Directory.Parent.FullName 'retired.json')) } |
        ForEach-Object { $_.Directory.Parent.FullName.Substring($prefix.Length).Replace('\','/') } |
        Sort-Object)
    if($names.Count -eq 0){ throw 'No runnable cases with system/run.dict were packaged.' }
    return $names
}

try {
    if($Help){
        Write-Output @'
Usage: Start.cmd [-Case <name>] [-Action new|resume|check] [options]
  No arguments          Choose a case and action interactively; q cancels.
  -Restart <path>        Old run/checkpoint inside this package's runs/.
  -EndTime <code-time>   Extend the saved end time with -Action resume.
  -Threads <1..256>      Override the case's requested thread count.
  -NoPause              Return immediately when Start.cmd completes.

Examples:
  Start.cmd -Case vof_a_drop_deformation_amr -Action new -NoPause
  Start.cmd -Case vof_a_drop_deformation_amr -Action resume
  Start.cmd -Case vof_a_drop_deformation_amr -Action check -NoPause

Edit cases/<name>/system/*.dict before a new run. A resume uses the old run's
frozen dictionaries. Every action creates a fresh run directory via tools/run_case.
Restart requires a completed Hydro/A VOF checkpoint; B1 restart is unsupported.
'@
        exit 0
    }
    $root=[IO.Path]::GetFullPath($PSScriptRoot)
    Assert-PortableAssets $root
    $cases=@(Get-PortableCases $root)
    $oldInput=[Console]::InputEncoding
    $oldOutput=[Console]::OutputEncoding
    try {
        [Console]::InputEncoding=[Text.UTF8Encoding]::new($false)
        [Console]::OutputEncoding=[Text.UTF8Encoding]::new($false)
        if(-not $Case){
            [Console]::WriteLine('LAMP portable - cases:')
            for($i=0;$i -lt $cases.Count;$i++){ [Console]::WriteLine('  ['+($i+1)+'] '+$cases[$i]) }
            $answer=Read-PortableSelection 'Select a case number or exact name (q cancels): '
            if($null -eq $answer){ [Console]::WriteLine('[PORTABLE] Cancelled; no run created.'); exit 0 }
            if($answer -match '\A[0-9]+\z'){
                $number=0
                if(-not [int]::TryParse($answer,[ref]$number) -or $number -lt 1 -or $number -gt $cases.Count){
                    throw 'Case selection is outside the displayed range.'
                }
                $Case=$cases[$number-1]
            }else{ $Case=$answer }
        }
        $matched=@($cases | Where-Object { $_ -ieq $Case })
        if($matched.Count -ne 1){ throw "Unknown packaged case: $Case" }
        $Case=$matched[0]
        if(-not $Action){
            [Console]::WriteLine('  [1] New run    [2] Resume from checkpoint    [3] Validate dictionaries    [q] Cancel')
            $answer=Read-PortableSelection 'Select action: '
            if($null -eq $answer){ [Console]::WriteLine('[PORTABLE] Cancelled; no run created.'); exit 0 }
            switch -Exact ($answer){
                '1' { $Action='new' }
                '2' { $Action='resume' }
                '3' { $Action='check' }
                default { throw 'Action selection must be 1, 2, 3 or q.' }
            }
        }
    }finally{
        [Console]::InputEncoding=$oldInput
        [Console]::OutputEncoding=$oldOutput
    }
    if($Action -ne 'resume' -and ($Restart -or $PSBoundParameters.ContainsKey('EndTime'))){
        throw '-Restart and -EndTime require -Action resume.'
    }
    $arguments=@{ Case=$Case; Program=(Join-Path $root 'bin/LAMP.exe') }
    if($PSBoundParameters.ContainsKey('Threads')){ $arguments.Threads=$Threads }
    switch($Action){
        'check' { $arguments.CheckOnly=$true }
        'resume' {
            if($Restart){
                $arguments.Restart=if([IO.Path]::IsPathRooted($Restart)){$Restart}else{Join-Path $root $Restart}
            }else{ $arguments.Resume=$true }
            if($PSBoundParameters.ContainsKey('EndTime')){ $arguments.EndTime=$EndTime }
        }
    }
    [Console]::WriteLine('[PORTABLE] '+$Action+' '+$Case)
    $global:LASTEXITCODE=0
    # The existing launcher owns its parameter/default rules; do not impose
    # this entry script's StrictMode on a child script.
    Set-StrictMode -Off
    & (Join-Path $root 'tools/run_case.ps1') @arguments
    exit $LASTEXITCODE
}catch{
    [Console]::Error.WriteLine('[PORTABLE] ERROR: '+$_.Exception.Message)
    exit 2
}
