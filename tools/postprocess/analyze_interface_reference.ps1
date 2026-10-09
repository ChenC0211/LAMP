<#
.SYNOPSIS
Verify a saved Interface run against its frozen analytic reference.
.EXAMPLE
.\tools\postprocess\analyze_interface_reference.ps1 -Run runs/<case>/<timestamp>
#>
[CmdletBinding()]
param([string]$Run='', [string]$Python='', [switch]$AcceptedPrefix, [switch]$SelfTest)
$ErrorActionPreference='Stop'
$project=Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
if(-not $Python){
    $candidate=Get-Command pvpython.exe -ErrorAction SilentlyContinue
    if($candidate){$Python=$candidate.Source}
    else{
        $roots=@((Join-Path $env:LOCALAPPDATA 'Programs'),$env:ProgramFiles)
        foreach($root in $roots){
            $found=Get-ChildItem -LiteralPath $root -Directory -Filter 'ParaView*' -ErrorAction SilentlyContinue |
                Sort-Object Name -Descending | ForEach-Object {Join-Path $_.FullName 'bin/pvpython.exe'} |
                Where-Object {Test-Path -LiteralPath $_ -PathType Leaf} | Select-Object -First 1
            if($found){$Python=$found;break}
        }
        if(-not $Python){
            $candidate=Get-Command python.exe -ErrorAction SilentlyContinue
            if($candidate -and $candidate.Source -notmatch 'WindowsApps'){$Python=$candidate.Source}
        }
    }
}
if(-not $Python){throw 'Provide -Python <python.exe> with NumPy/Matplotlib, or install ParaView with pvpython.'}
$Python=(Resolve-Path -LiteralPath $Python).Path
$arguments=@()
if([IO.Path]::GetFileName($Python) -eq 'pvpython.exe'){$arguments+='--disable-registry'}
$arguments+=Join-Path $PSScriptRoot 'analyze_interface_reference.py'
if($SelfTest){$arguments+='--self-test'}else{
    if(-not $Run){throw 'Provide -Run <timestamped run directory>, or -SelfTest.'}
    $arguments+=@('--run',(Resolve-Path -LiteralPath $Run).Path)
}
if($AcceptedPrefix){$arguments+='--accepted-prefix'}
$cache=Join-Path $project ('runs/checks/interface_reference/'+(Get-Date -Format 'yyyyMMdd_HHmmss_fff')+'/matplotlib')
New-Item -ItemType Directory -Force -Path $cache | Out-Null
$previous=[Environment]::GetEnvironmentVariable('MPLCONFIGDIR','Process')
try{$env:MPLCONFIGDIR=$cache; & $Python @arguments; $code=$LASTEXITCODE}
finally{[Environment]::SetEnvironmentVariable('MPLCONFIGDIR',$previous,'Process')}
exit $code
