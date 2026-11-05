<#
.SYNOPSIS
Create an app-local Windows x64 LAMP distribution without developer tools.
.EXAMPLE
tools\package_portable.cmd -OpenMpRuntime D:\runtime\libomp140.x86_64.dll -OpenMpLicenseDirectory D:\runtime\licenses -VcRuntimeDirectory D:\redist\x64\Microsoft.VC143.CRT -VcRuntimeLicensePath D:\redist\LICENSE.rtf -AdditionalLicenseDirectory D:\licenses
#>
[CmdletBinding()]
param(
    [string]$Program='',
    [string]$OutputDirectory='',
    [Parameter(Mandatory=$true)][string]$OpenMpRuntime,
    [Parameter(Mandatory=$true)][string]$OpenMpLicenseDirectory,
    [string]$VcRuntimeDirectory='',
    [Parameter(Mandatory=$true)][string]$VcRuntimeLicensePath,
    [Parameter(Mandatory=$true)][string]$AdditionalLicenseDirectory,
    [string]$ManifestTool='',
    [string[]]$Cases=@(),
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9._-]*$')][string]$Version=(Get-Date -Format 'yyyyMMdd-HHmmss')
)
$ErrorActionPreference='Stop'
. (Join-Path $PSScriptRoot 'portable/Read-PortablePe.ps1')
. (Join-Path $PSScriptRoot 'portable/Test-PortableInputs.ps1')
$project=Split-Path -Parent $PSScriptRoot
$utf8=[Text.UTF8Encoding]::new($false)
function Get-PortableSha([string]$Path){
    $stream=[IO.File]::OpenRead($Path);$hash=[Security.Cryptography.SHA256]::Create()
    try{return [BitConverter]::ToString($hash.ComputeHash($stream)).Replace('-','').ToLowerInvariant()}
    finally{$hash.Dispose();$stream.Dispose()}
}
function Assert-PortableNoLinks([string]$Path){
    $entry=Get-Item -LiteralPath $Path -Force
    $items=@($entry)
    if($entry.PSIsContainer){$items+=@(Get-ChildItem -LiteralPath $Path -Recurse -Force)}
    foreach($item in $items){if($item.Attributes -band [IO.FileAttributes]::ReparsePoint){throw "Package inputs cannot contain links: $($item.FullName)"}}
}
function Copy-PortableFile([string]$Source,[string]$Destination){
    Assert-PortableNoLinks $Source
    [void][IO.Directory]::CreateDirectory((Split-Path -Parent $Destination))
    Copy-Item -LiteralPath $Source -Destination $Destination
}
function Invoke-PortableManifestTool([string[]]$Arguments,[string]$Log){
    $previousPreference=$ErrorActionPreference;$ErrorActionPreference='Continue';$global:LASTEXITCODE=$null
    $previousTemp=$env:TEMP;$previousTmp=$env:TMP;$manifestTemp=Join-Path $staging 'temp'
    [void][IO.Directory]::CreateDirectory($manifestTemp);$env:TEMP=$manifestTemp;$env:TMP=$manifestTemp
    try{$lines=@(& $ManifestTool @Arguments 2>&1);$code=$LASTEXITCODE}
    finally{$ErrorActionPreference=$previousPreference;$env:TEMP=$previousTemp;$env:TMP=$previousTmp}
    [IO.File]::AppendAllLines($Log,[string[]]@($lines|ForEach-Object {[string]$_}),$utf8)
    if($null -eq $code -or $code -ne 0){throw "Windows manifest tool failed (exit $code); see $Log"}
}
if(-not $Program){$Program=Join-Path $project 'bin/LAMP.exe'}
$Program=(Get-Item -LiteralPath $Program).FullName
$OpenMpRuntime=(Get-Item -LiteralPath $OpenMpRuntime).FullName
$OpenMpLicenseDirectory=(Get-Item -LiteralPath $OpenMpLicenseDirectory).FullName
$VcRuntimeLicensePath=(Get-Item -LiteralPath $VcRuntimeLicensePath).FullName
$AdditionalLicenseDirectory=(Get-Item -LiteralPath $AdditionalLicenseDirectory).FullName
if(-not $ManifestTool){
    $sdk=Join-Path ${env:ProgramFiles(x86)} 'Windows Kits/10/bin'
    $tools=@(if(Test-Path -LiteralPath $sdk){Get-ChildItem -LiteralPath $sdk -Directory|Where-Object {$_.Name -match '\A[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+\z'}|
        Sort-Object {[version]$_.Name} -Descending|ForEach-Object {Join-Path $_.FullName 'x64/mt.exe'}|Where-Object {Test-Path -LiteralPath $_ -PathType Leaf}})
    if($tools.Count){$ManifestTool=$tools[0]}
}
if(-not $ManifestTool -or -not (Test-Path -LiteralPath $ManifestTool -PathType Leaf)){throw 'Packaging requires Windows SDK mt.exe to embed the UTF-8 process manifest. Supply -ManifestTool; the destination computer does not need the SDK.'}
$ManifestTool=(Get-Item -LiteralPath $ManifestTool).FullName
if(-not $OutputDirectory){$OutputDirectory=Join-Path $project 'dist'}
$OutputDirectory=[IO.Path]::GetFullPath($OutputDirectory)
if($OutputDirectory -notmatch '\AD:[\\/]'){throw 'Store package output on D: as required by the project storage policy.'}
if($OpenMpRuntime -match '(?i)[\\/]debug_nonredist[\\/]'){throw 'Visual Studio debug_nonredist OpenMP cannot be redistributed.'}
if((Get-PortableSha $OpenMpRuntime) -eq '9ea4ac3d55a425d5137a5590e3b7d2a2289b3bbff7941abe2b96fa148b75eef3'){
    throw 'The original Visual Studio nonredistributable OpenMP DLL cannot be packaged.'
}
$programSha=Get-PortableSha $Program
$programImage=Read-PortablePe $Program;Assert-PortableReleaseImports $programImage
if($programImage.IsDll){throw 'Program must be an executable, not a DLL.'}
$ompImage=Read-PortablePe $OpenMpRuntime;Assert-PortableReleaseImports $ompImage
if(-not $ompImage.IsDll){throw 'OpenMpRuntime must be a Windows x64 DLL.'}
$ompImports=@($programImage.Imports|Where-Object {$_ -match '\Alibomp.*\.dll\z'})
if($ompImports.Count -ne 1){throw 'Expected exactly one LLVM OpenMP import in LAMP.exe.'}
if(-not $VcRuntimeDirectory){
    $roots=@(${env:ProgramFiles},${env:ProgramFiles(x86)})|Where-Object {$_}
    $candidates=@(foreach($root in $roots){
        $vs=Join-Path $root 'Microsoft Visual Studio/2022'
        if(Test-Path -LiteralPath $vs){Get-ChildItem -LiteralPath $vs -Filter 'Microsoft.VC143.CRT' -Directory -Recurse -ErrorAction SilentlyContinue |
            Where-Object {$_.FullName -match '(?i)[\\/]Redist[\\/]MSVC[\\/][^\\/]+[\\/]x64[\\/]' -and $_.FullName -notmatch '(?i)debug'}}
    })
    $eligible=@(foreach($candidate in $candidates){
        $vc=Join-Path $candidate.FullName 'vcruntime140.dll'
        if(Test-Path -LiteralPath $vc){$image=Read-PortablePe $vc;if($image.Version -ge [version]'14.41'){[pscustomobject]@{Path=$candidate.FullName;Version=$image.Version}}}
    })
    $best=$eligible|Sort-Object Version -Descending|Select-Object -First 1
    if(-not $best){throw 'No formal VS Redist x64 CRT >=14.41 found. Supply -VcRuntimeDirectory from a current official redistributable.'}
    $VcRuntimeDirectory=$best.Path
}
$VcRuntimeDirectory=(Get-Item -LiteralPath $VcRuntimeDirectory).FullName
if($VcRuntimeDirectory -match '(?i)[\\/](?:debug_nonredist|debug)[\\/]'){throw 'Debug CRT directories cannot be packaged.'}
Assert-PortableNoLinks $VcRuntimeDirectory
$available=@{}
foreach($dll in Get-ChildItem -LiteralPath $VcRuntimeDirectory -Filter '*.dll' -File){$available[$dll.Name.ToLowerInvariant()]=$dll.FullName}
$available[$ompImports[0]]=$OpenMpRuntime
$closure=@{};$pending=[Collections.Generic.Queue[string]]::new();$pending.Enqueue($Program)
while($pending.Count){
    $path=$pending.Dequeue();$image=Read-PortablePe $path;Assert-PortableReleaseImports $image
    foreach($name in $image.Imports){
        if(Test-PortableSystemDll $name){continue}
        if($closure.ContainsKey($name)){continue}
        if(-not $available.ContainsKey($name)){throw "Unresolved non-system runtime dependency: $name required by $path"}
        $dependency=Read-PortablePe $available[$name];Assert-PortableReleaseImports $dependency
        if(-not $dependency.IsDll){throw "Runtime dependency is not a DLL: $name"}
        if($name -match '\A(?:msvcp|vcruntime|concrt|vccorlib)140.*\.dll\z' -and $dependency.Version -lt [version]'14.41'){
            throw "CRT must be >=14.41 for this MSVC 19.41 release: $name version $($dependency.Version)"
        }
        $closure[$name]=[pscustomobject]@{Path=$available[$name];Image=$dependency;Sha256=(Get-PortableSha $available[$name])};$pending.Enqueue($available[$name])
    }
}
$name='LAMP-Portable-'+$Version;$final=Join-Path $OutputDirectory $name;$zip=$final+'.zip'
if((Test-Path -LiteralPath $final) -or (Test-Path -LiteralPath $zip)){throw "Release already exists; choose a new -Version: $final"}
[void][IO.Directory]::CreateDirectory($OutputDirectory)
$staging=Join-Path $OutputDirectory ('.pending-'+[Guid]::NewGuid().ToString('N'))
[void][IO.Directory]::CreateDirectory($staging)
$package=Join-Path $staging $name
foreach($folder in @('bin','tools','cases','runs','licenses','validation')){[void][IO.Directory]::CreateDirectory((Join-Path $package $folder))}
Copy-PortableFile $Program (Join-Path $package 'bin/LAMP.exe')
if((Get-PortableSha (Join-Path $package 'bin/LAMP.exe')) -ne $programSha -or (Get-PortableSha $Program) -ne $programSha){throw 'Source executable changed during packaging; no release was published.'}
$deployedProgram=Join-Path $package 'bin/LAMP.exe'
$deploymentEvidence=Join-Path $package 'validation/deployment'
[void][IO.Directory]::CreateDirectory($deploymentEvidence)
$originalManifest=Join-Path $deploymentEvidence 'original.manifest'
$embeddedManifest=Join-Path $deploymentEvidence 'embedded.manifest'
$manifestFragment=Join-Path $deploymentEvidence 'LAMP.Portable.manifest'
$manifestLog=Join-Path $deploymentEvidence 'manifest-tool.log'
Copy-PortableFile (Join-Path $PSScriptRoot 'portable/LAMP.Portable.manifest') $manifestFragment
Invoke-PortableManifestTool @('-nologo',('-inputresource:'+$deployedProgram+';#1'),('-out:'+$originalManifest)) $manifestLog
Invoke-PortableManifestTool @('-nologo',('-inputresource:'+$deployedProgram+';#1'),'-manifest',$manifestFragment,('-outputresource:'+$deployedProgram+';#1')) $manifestLog
Invoke-PortableManifestTool @('-nologo',('-inputresource:'+$deployedProgram+';#1'),('-out:'+$embeddedManifest)) $manifestLog
$deployedImage=Read-PortablePe $deployedProgram
Assert-PortableResourceOnlyTransform $programImage $deployedImage
[xml]$originalXml=[IO.File]::ReadAllText($originalManifest)
[xml]$embeddedXml=[IO.File]::ReadAllText($embeddedManifest)
$originalTrust=$originalXml.SelectSingleNode('//*[local-name()="trustInfo"]')
$embeddedTrust=$embeddedXml.SelectSingleNode('//*[local-name()="trustInfo"]')
if($null -eq $originalTrust -or $null -eq $embeddedTrust -or $originalTrust.OuterXml -cne $embeddedTrust.OuterXml){throw 'UTF-8 deployment did not preserve original executable trustInfo.'}
$activeCodePage=$embeddedXml.SelectSingleNode('//*[local-name()="activeCodePage" and namespace-uri()="http://schemas.microsoft.com/SMI/2019/WindowsSettings"]')
if($null -eq $activeCodePage -or $activeCodePage.InnerText -cne 'UTF-8'){throw 'UTF-8 deployment manifest is missing its activeCodePage setting.'}
if((Get-PortableSha $Program) -ne $programSha){throw 'Source executable changed while creating the deployment manifest.'}
$packagedProgramSha=Get-PortableSha $deployedProgram
foreach($key in $closure.Keys){
    $target=Join-Path $package ('bin/'+$key);Copy-PortableFile $closure[$key].Path $target
    if((Get-PortableSha $target) -ne $closure[$key].Sha256 -or (Get-PortableSha $closure[$key].Path) -ne $closure[$key].Sha256){throw "Runtime changed during packaging: $key"}
}
foreach($file in @('run_case.cmd','run_case.ps1')){Copy-PortableFile (Join-Path $PSScriptRoot $file) (Join-Path $package ('tools/'+$file))}
foreach($file in @('start.cmd','start.ps1')){Copy-PortableFile (Join-Path $PSScriptRoot ('portable/'+$file)) (Join-Path $package $file)}
Copy-PortableFile (Join-Path $PSScriptRoot 'portable/README.txt') (Join-Path $package 'README.txt')
Assert-PortableNoLinks $OpenMpLicenseDirectory
if(-not @(Get-ChildItem -LiteralPath $OpenMpLicenseDirectory -File -Recurse).Count){throw 'OpenMP license directory is empty.'}
Copy-Item -LiteralPath $OpenMpLicenseDirectory -Destination (Join-Path $package 'licenses/OpenMP') -Recurse
Copy-PortableFile $VcRuntimeLicensePath (Join-Path $package ('licenses/Microsoft-Visual-Cpp-Runtime'+[IO.Path]::GetExtension($VcRuntimeLicensePath)))
Assert-PortableNoLinks $AdditionalLicenseDirectory
if(-not @(Get-ChildItem -LiteralPath $AdditionalLicenseDirectory -File -Recurse).Count){throw 'Additional license directory must contain static dependency licenses and notices.'}
Copy-Item -LiteralPath $AdditionalLicenseDirectory -Destination (Join-Path $package 'licenses/ThirdParty') -Recurse
if(Test-Path -LiteralPath (Join-Path $project 'LICENSE')){Copy-PortableFile (Join-Path $project 'LICENSE') (Join-Path $package 'licenses/LAMP-LICENSE')}
$casesRoot=Join-Path $project 'cases';$explicit=$Cases.Count -gt 0
if(-not $explicit){$Cases=@(Get-ChildItem -LiteralPath $casesRoot -Directory|Where-Object {
    -not (Test-Path -LiteralPath (Join-Path $_.FullName 'retired.json')) -and
    (Test-Path -LiteralPath (Join-Path $_.FullName 'system/run.dict'))}|Select-Object -ExpandProperty Name|Sort-Object)}
$included=[Collections.Generic.List[string]]::new();$excluded=[Collections.Generic.List[object]]::new()
$oldPath=$env:PATH;$oldTmp=$env:TMP;$oldTemp=$env:TEMP
$oldConsoleEncoding=[Console]::OutputEncoding;$oldPipelineEncoding=$OutputEncoding
$env:PATH=Join-Path $env:SystemRoot 'System32';$env:TEMP=Join-Path $staging 'temp';$env:TMP=$env:TEMP
[Console]::OutputEncoding=$utf8;$OutputEncoding=$utf8
[void][IO.Directory]::CreateDirectory($env:TEMP)
try{
    foreach($case in $Cases){
        if($case -notmatch '\A[A-Za-z0-9_-]+\z'){throw "Select a top-level formal case name: $case"}
        $source=Join-Path $casesRoot $case
        if(-not (Test-Path -LiteralPath (Join-Path $source 'system/run.dict')) -or
            (Test-Path -LiteralPath (Join-Path $source 'retired.json')) -or (Test-Path -LiteralPath (Join-Path $source 'benchmark.json'))){throw "Case is missing, retired, or a benchmark: $case"}
        $runText=[IO.File]::ReadAllText((Join-Path $source 'system/run.dict'))
        $app=([regex]::Match($runText,'\bapp\s+([^\s;{}]+)')).Groups[1].Value
        if($app -notin @('Electric','Hydrodynamics','VofElectrohydrodynamics','Electrospray')){
            if($explicit){throw "Case uses a retired/unsupported app: $case ($app)"}
            $excluded.Add([ordered]@{case=$case;reason="unsupported app $app"});continue
        }
        Assert-PortableNoLinks $source;Assert-PortableInputResources $source
        $destination=Join-Path $package ('cases/'+$case);[void][IO.Directory]::CreateDirectory($destination)
        foreach($folder in @('system','constant','input','data','reference')){
            $input=Join-Path $source $folder
            if(Test-Path -LiteralPath $input){Copy-Item -LiteralPath $input -Destination $destination -Recurse}
        }
        if(Test-Path -LiteralPath (Join-Path $source 'README.md')){Copy-PortableFile (Join-Path $source 'README.md') (Join-Path $destination 'README.md')}
        $log=Join-Path $package ('validation/'+$case+'.check.log')
        $nativeErrorPreference=$ErrorActionPreference;$ErrorActionPreference='Continue'
        $global:LASTEXITCODE=$null
        try{$lines=@(& (Join-Path $package 'bin/LAMP.exe') --case $destination --check-case 2>&1);$code=$LASTEXITCODE}
        finally{$ErrorActionPreference=$nativeErrorPreference}
        [IO.File]::WriteAllLines($log,[string[]]@($lines|ForEach-Object {[string]$_}),$utf8)
        if($null -eq $code -or $code -ne 0){
            if($explicit){throw "Packaged executable rejected case $case (exit $code); see $log"}
            $excluded.Add([ordered]@{case=$case;reason="input check exit $code"})
            # This path was created in this staging package, never a source.
            if(-not $destination.StartsWith($package+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)){throw 'Invalid staging case cleanup path.'}
            Remove-Item -LiteralPath $destination -Recurse -Force
            continue
        }
        $included.Add($case);Write-Output "[PACKAGE CHECK] $case OK"
    }
}finally{
    $env:PATH=$oldPath;$env:TEMP=$oldTemp;$env:TMP=$oldTmp
    [Console]::OutputEncoding=$oldConsoleEncoding;$OutputEncoding=$oldPipelineEncoding
}
if(-not $included.Count){throw "No formal cases passed; staging evidence retained at $staging"}
$files=@(Get-ChildItem -LiteralPath $package -File -Recurse|Sort-Object FullName|ForEach-Object {
    $entry=[ordered]@{path=$_.FullName.Substring($package.Length+1).Replace('\','/');bytes=$_.Length;sha256=(Get-PortableSha $_.FullName)}
    if($_.Extension -in @('.exe','.dll')){$info=Read-PortablePe $_.FullName;$entry['file_version']=$info.FileVersion;$entry['imports']=@($info.Imports)}
    $entry
})
$manifest=[ordered]@{schema_version=1;name=$name;created=(Get-Date).ToString('o');platform='Windows 10 version 1903 or later / Windows 11, x64';
    source_program=[ordered]@{path=$Program;sha256=$programSha};additional_license_source=$AdditionalLicenseDirectory;
    packaged_program=[ordered]@{path='bin/LAMP.exe';sha256=$packagedProgramSha};
    deployment_transform=[ordered]@{type='Embedded activeCodePage=UTF-8 application manifest';manifest_tool=$ManifestTool;
        tool_file_version=([Diagnostics.FileVersionInfo]::GetVersionInfo($ManifestTool)).FileVersion;
        fragment_sha256=(Get-PortableSha $manifestFragment);original_trust_info_preserved=$true;imports_unchanged=$true;non_resource_sections_unchanged=$true;
        non_resource_sections=@($programImage.Sections|Where-Object {$_.Name -cne $programImage.ResourceSection}|ForEach-Object {[ordered]@{name=$_.Name;sha256=$_.Sha256}})};
    runtime=[ordered]@{openmp_source=$OpenMpRuntime;openmp_sha256=(Get-PortableSha $OpenMpRuntime);vc_runtime_source=$VcRuntimeDirectory;
        minimum_vc_version='14.41';non_system_dependencies=@($closure.Keys|Sort-Object);system_imports='Windows system DLL/API-set allowlist'};
    cases=$included.ToArray();excluded_cases=$excluded.ToArray();files=$files;
    validation=[ordered]@{dependency_closure='passed';release_runtime_imports='passed';case_preflight='passed';path='System32 only';
        clean_windows_machine='Not performed by this packaging script; see release verification record.'}}
[IO.File]::WriteAllText((Join-Path $package 'portable_manifest.json'),($manifest|ConvertTo-Json -Depth 9),$utf8)
Add-Type -AssemblyName System.IO.Compression.FileSystem
$pendingZip=Join-Path $staging ($name+'.zip')
[IO.Compression.ZipFile]::CreateFromDirectory($package,$pendingZip,[IO.Compression.CompressionLevel]::Optimal,$true)
Move-Item -LiteralPath $package -Destination $final
Move-Item -LiteralPath $pendingZip -Destination $zip
Write-Output "Portable directory: $final"
Write-Output "Portable ZIP: $zip"
Write-Output ('SHA256: '+(Get-PortableSha $zip))
