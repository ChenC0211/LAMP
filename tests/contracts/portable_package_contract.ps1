[CmdletBinding()]
param([string]$WorkingDirectory='', [string]$ManifestTool='')
$ErrorActionPreference='Stop'
$project=Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
. (Join-Path $project 'tools/portable/Read-PortablePe.ps1')
. (Join-Path $project 'tools/portable/Test-PortableInputs.ps1')
if(-not $WorkingDirectory){$WorkingDirectory=Join-Path $project 'runs/checks/portable_release_20261019/packager/contract'}
[void][IO.Directory]::CreateDirectory($WorkingDirectory)
$script:checks=0
function Assert-Contract([bool]$Condition,[string]$Message){if(-not $Condition){throw $Message};$script:checks++}
function Assert-Rejected([scriptblock]$Action,[string]$Message){
    $rejected=$false
    try{& $Action|Out-Null}catch{$rejected=$true}
    Assert-Contract $rejected $Message
}
function New-ContractImage([string]$Name,[string]$Import='kernel32.dll',[switch]$Delay){
    $image=[byte[]]::new(2048)
    function Put16([int]$Offset,[uint16]$Value){[BitConverter]::GetBytes($Value).CopyTo($image,$Offset)}
    function Put32([int]$Offset,[uint32]$Value){[BitConverter]::GetBytes($Value).CopyTo($image,$Offset)}
    Put16 0 0x5a4d;Put32 0x3c 128;Put32 128 0x4550
    Put16 132 0x8664;Put16 134 1;Put16 148 240;Put16 150 0x2022
    Put16 152 0x20b;Put32 (152+60) 512;Put32 (152+108) 16
    $index=if($Delay){13}else{1};$stride=if($Delay){32}else{20}
    Put32 (152+112+$index*8) 4096;Put32 (152+116+$index*8) (2*$stride)
    # Section maps 0x1000..0x15ff to byte 0x200..0x7ff.
    Put32 (392+8) 1536;Put32 (392+12) 4096;Put32 (392+16) 1536;Put32 (392+20) 512
    if($Delay){Put32 512 1;Put32 516 4352}else{Put32 524 4352}
    [Text.Encoding]::ASCII.GetBytes($Import+[char]0).CopyTo($image,768)
    $path=Join-Path $WorkingDirectory $Name;[IO.File]::WriteAllBytes($path,$image)
    return $path
}
$normal=New-ContractImage 'normal.dll'
$read=Read-PortablePe $normal
Assert-Contract ($read.Machine -eq 'AMD64' -and $read.IsDll -and $read.Imports.Count -eq 1 -and $read.Imports[0] -eq 'kernel32.dll') 'Native x64 import was not read.'
$delayed=Read-PortablePe (New-ContractImage 'delay.dll' 'vcruntime140.dll' -Delay)
Assert-Contract ($delayed.Imports.Count -eq 1 -and $delayed.Imports[0] -eq 'vcruntime140.dll') 'Delay dependency was not included.'
Assert-Contract (Test-PortableSystemDll 'kernel32.dll') 'System DLL rejected.'
Assert-Contract (Test-PortableSystemDll 'api-ms-win-crt-runtime-l1-1-0.dll') 'API-set rejected.'
Assert-Contract (-not (Test-PortableSystemDll 'msvcp140.dll')) 'App-local CRT was incorrectly treated as a system DLL.'
Assert-Contract (-not (Test-PortableSystemDll 'surprise.dll')) 'Unrecognized dependency was incorrectly exempted.'
foreach($name in @('vcruntime140d.dll','msvcp140d_atomic_wait.dll','ucrtbased.dll')){
    $debug=Read-PortablePe (New-ContractImage ($name+'.fixture') $name)
    Assert-Rejected {Assert-PortableReleaseImports $debug} "Debug runtime was accepted: $name"
}
$bad=[IO.File]::ReadAllBytes($normal);$bad[132]=0x4c;$bad[133]=1
$path=Join-Path $WorkingDirectory 'wrong-machine.dll';[IO.File]::WriteAllBytes($path,$bad)
Assert-Rejected {Read-PortablePe $path} 'An x86 image was accepted.'
$bad=[IO.File]::ReadAllBytes($normal);$bad=$bad[0..530]
$path=Join-Path $WorkingDirectory 'truncated.dll';[IO.File]::WriteAllBytes($path,$bad)
Assert-Rejected {Read-PortablePe $path} 'Truncated imports were accepted.'
$bad=[IO.File]::ReadAllBytes($normal);[BitConverter]::GetBytes([uint32]999999).CopyTo($bad,524)
$path=Join-Path $WorkingDirectory 'unmapped-name.dll';[IO.File]::WriteAllBytes($path,$bad)
Assert-Rejected {Read-PortablePe $path} 'An unmapped import name RVA was accepted.'
$bad=[IO.File]::ReadAllBytes($normal);[BitConverter]::GetBytes([uint32]20).CopyTo($bad,(152+116+8))
$path=Join-Path $WorkingDirectory 'unterminated.dll';[IO.File]::WriteAllBytes($path,$bad)
Assert-Rejected {Read-PortablePe $path} 'An unterminated imports table was accepted.'
foreach($case in @('electric_charged_circle_xy','vof_a_drop_deformation_amr')){
    Assert-PortableInputResources (Join-Path $project ('cases/'+$case));$script:checks++
}
$resourceCase=Join-Path $WorkingDirectory 'resource-case'
[void][IO.Directory]::CreateDirectory((Join-Path $resourceCase 'system'))
[void][IO.Directory]::CreateDirectory((Join-Path $resourceCase 'input'))
$resourceDictionary=Join-Path $resourceCase 'system/run.dict';$encoding=[Text.UTF8Encoding]::new($false)
[IO.File]::WriteAllText((Join-Path $resourceCase 'input/shape.stl'),'solid shape',$encoding)
[IO.File]::WriteAllText($resourceDictionary,'Output { Fields { profile physical; } } Geometry { meshFile "input/shape.stl"; }',$encoding)
Assert-PortableInputResources $resourceCase;$script:checks++
foreach($badResource in @('D:/outside/shape.stl','../shape.stl','input/missing.stl','geometry/shape.stl')){
    [IO.File]::WriteAllText($resourceDictionary,('Geometry { meshFile "'+$badResource+'"; }'),$encoding)
    Assert-Rejected {Assert-PortableInputResources $resourceCase} "Unsafe/missing case resource accepted: $badResource"
}
[IO.File]::WriteAllText($resourceDictionary,'// file "D:/outside/shape.stl";'+"`n"+'Initial { phi "(x-.5)^2+(y-.5)^2-.1^2"; } profile physical;',$encoding)
Assert-PortableInputResources $resourceCase;$script:checks++
foreach($file in @('tools/package_portable.ps1','tools/portable/Read-PortablePe.ps1','tools/portable/Test-PortableInputs.ps1')){
    $tokens=$null;$errors=$null;[void][Management.Automation.Language.Parser]::ParseFile((Join-Path $project $file),[ref]$tokens,[ref]$errors)
    Assert-Contract ($errors.Count -eq 0) "PowerShell syntax errors in $file"
}
if($ManifestTool){
    $source=Join-Path $project 'bin/LAMP.exe';$copy=Join-Path $WorkingDirectory 'utf8-deployed.exe'
    $sourceHash=(Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash
    Copy-Item -LiteralPath $source -Destination $copy
    $originalManifest=Join-Path $WorkingDirectory 'original.manifest';$mergedManifest=Join-Path $WorkingDirectory 'merged.manifest'
    & $ManifestTool -nologo ('-inputresource:'+$copy+';#1') ('-out:'+$originalManifest)|Out-Null
    if($LASTEXITCODE -ne 0){throw 'Failed to extract the contract source manifest.'}
    & $ManifestTool -nologo ('-inputresource:'+$copy+';#1') -manifest (Join-Path $project 'tools/portable/LAMP.Portable.manifest') ('-outputresource:'+$copy+';#1')|Out-Null
    if($LASTEXITCODE -ne 0){throw 'Failed to merge the contract UTF-8 manifest.'}
    & $ManifestTool -nologo ('-inputresource:'+$copy+';#1') ('-out:'+$mergedManifest)|Out-Null
    if($LASTEXITCODE -ne 0){throw 'Failed to extract the contract merged manifest.'}
    $before=Read-PortablePe $source;$after=Read-PortablePe $copy
    Assert-PortableResourceOnlyTransform $before $after;$script:checks++
    Assert-Contract ((Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash -ceq $sourceHash) 'Manifest deployment changed the original source executable.'
    [xml]$original=[IO.File]::ReadAllText($originalManifest);[xml]$merged=[IO.File]::ReadAllText($mergedManifest)
    Assert-Contract ($original.SelectSingleNode('//*[local-name()="trustInfo"]').OuterXml -ceq $merged.SelectSingleNode('//*[local-name()="trustInfo"]').OuterXml) 'Original trustInfo was not preserved.'
    Assert-Contract ($merged.SelectSingleNode('//*[local-name()="activeCodePage"]').InnerText -ceq 'UTF-8') 'UTF-8 process code page was not embedded.'
    $bad=[IO.File]::ReadAllBytes($copy);$codeSection=$after.Sections|Where-Object {$_.Name -ceq '.text'}|Select-Object -First 1
    $bad[$codeSection.Raw]=$bad[$codeSection.Raw] -bxor 1
    $tampered=Join-Path $WorkingDirectory 'changed-code.exe';[IO.File]::WriteAllBytes($tampered,$bad)
    Assert-Rejected {Assert-PortableResourceOnlyTransform $before (Read-PortablePe $tampered)} 'A deployment transform that changed code was accepted.'
    $after.Imports+=@('unexpected.dll')
    Assert-Rejected {Assert-PortableResourceOnlyTransform $before $after} 'A deployment transform that changed imports was accepted.'
}
Write-Output "Portable package contract PASS: $script:checks checks"
