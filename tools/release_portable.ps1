<#
.SYNOPSIS
Package the existing compiled application and selected cases.
.EXAMPLE
.\tools\release_portable.cmd
.EXAMPLE
.\tools\release_portable.cmd -Cases "hydro_static_drop,vof_a_charged_drop"
#>
[CmdletBinding()]
param(
    [string]$Cases='',
    [string]$RuntimeDirectory='',
    [string]$OutputDirectory='',
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9._-]*$')][string]$Version,
    [Alias('h')][switch]$Help
)
$ErrorActionPreference='Stop'
$project=Split-Path -Parent $PSScriptRoot
$utf8=[Text.UTF8Encoding]::new($false)
$record=$null;$recordPath=$null
$savedTemp=$env:TEMP;$savedTmp=$env:TMP
function Get-ReleaseSha([string]$Path){
    $stream=[IO.File]::OpenRead($Path);$hash=[Security.Cryptography.SHA256]::Create()
    try{return [BitConverter]::ToString($hash.ComputeHash($stream)).Replace('-','').ToLowerInvariant()}
    finally{$hash.Dispose();$stream.Dispose()}
}
function Save-ReleaseRecord {
    if($recordPath){[IO.File]::WriteAllText($recordPath,($record|ConvertTo-Json -Depth 7),$utf8)}
}
function Get-ReleaseCases {
    $root=Join-Path $project 'cases'
    return @(Get-ChildItem -LiteralPath $root -Directory|Where-Object {
        $run=Join-Path $_.FullName 'system/run.dict'
        if(-not (Test-Path -LiteralPath $run -PathType Leaf) -or
            (Test-Path -LiteralPath (Join-Path $_.FullName 'retired.json')) -or
            (Test-Path -LiteralPath (Join-Path $_.FullName 'benchmark.json'))){return $false}
        $text=[IO.File]::ReadAllText($run)
        # Exclude comments before identifying the formal application.
        $text=[regex]::Replace($text,'#[^\r\n]*|//[^\r\n]*','')
        $app=([regex]::Match($text,'\bapp\s+([^\s;{}]+)')).Groups[1].Value
        return $app -in @('Electric','Hydrodynamics','VofElectrohydrodynamics','Electrospray')
    }|Select-Object -ExpandProperty Name|Sort-Object)
}
function Resolve-ReleaseCases([string]$Expression,[string[]]$Available){
    if([string]::IsNullOrWhiteSpace($Expression)){throw 'Select at least one case; -Cases cannot be empty.'}
    $selected=[Collections.Generic.List[string]]::new()
    foreach($token in @($Expression.Trim() -split '[,;\s]+')){
        $indices=@()
        if($token -match '\A[0-9]+-[0-9]+\z'){
            $ends=$token.Split('-');$first=0;$last=0
            if(-not [int]::TryParse($ends[0],[ref]$first) -or -not [int]::TryParse($ends[1],[ref]$last) -or
                $first -lt 1 -or $last -lt $first -or $last -gt $Available.Count){throw "Invalid case range: $token"}
            $indices=@($first..$last)
        }elseif($token -match '\A[0-9]+\z'){
            $number=0
            if(-not [int]::TryParse($token,[ref]$number) -or $number -lt 1 -or $number -gt $Available.Count){throw "Case number outside the displayed range: $token"}
            $indices=@($number)
        }else{
            $matches=@($Available|Where-Object {$_ -ieq $token})
            if($matches.Count -ne 1){throw "Unknown or unsupported formal case: $token"}
            if($selected -notcontains $matches[0]){$selected.Add($matches[0])}
        }
        foreach($number in $indices){
            $name=$Available[$number-1]
            if($selected -notcontains $name){$selected.Add($name)}
        }
    }
    if(-not $selected.Count){throw 'Select at least one case.'}
    return $selected.ToArray()
}
function Assert-ReleaseRuntime([string]$Root){
    $manifest=Join-Path $Root 'runtime-cache.json'
    if(-not (Test-Path -LiteralPath $manifest -PathType Leaf)){
        throw "Portable runtime cache is missing: $manifest. Restore the complete .portable/runtime directory or supply -RuntimeDirectory."
    }
    $inventory=[IO.File]::ReadAllText($manifest,[Text.Encoding]::UTF8)|ConvertFrom-Json
    if($inventory.schema_version -ne 1 -or -not $inventory.files){throw 'Unsupported or empty portable runtime cache manifest.'}
    $prefix=$Root.TrimEnd('\','/')+[IO.Path]::DirectorySeparatorChar
    $seen=@{}
    foreach($entry in $inventory.files){
        $relative=[string]$entry.path
        if(-not $relative -or [IO.Path]::IsPathRooted($relative) -or $relative.Contains('\') -or $relative -match '(^|/)\.\.(/|$)' -or
            $seen.ContainsKey($relative)){throw "Invalid/duplicate runtime cache entry: $relative"}
        $seen[$relative]=$true
        $file=[IO.Path]::GetFullPath((Join-Path $Root $relative))
        if(-not $file.StartsWith($prefix,[StringComparison]::OrdinalIgnoreCase) -or -not (Test-Path -LiteralPath $file -PathType Leaf)){
            throw "Runtime cache file is missing or outside its directory: $relative"
        }
        if([string]$entry.sha256 -notmatch '\A[0-9a-fA-F]{64}\z' -or (Get-ReleaseSha $file) -ine [string]$entry.sha256){
            throw "Runtime cache file failed SHA256 verification: $relative"
        }
    }
    foreach($required in @('bin/libomp140.x86_64.dll','bin/msvcp140.dll','bin/vcruntime140.dll','bin/vcruntime140_1.dll','licenses/CRT/LICENSE-VCRUNTIME.rtf')){
        if(-not $seen.ContainsKey($required)){throw "Runtime cache manifest is missing required file: $required"}
    }
    foreach($folder in @('licenses/OpenMP/','licenses/ThirdParty/')){
        if(-not @($seen.Keys|Where-Object {$_.StartsWith($folder,[StringComparison]::OrdinalIgnoreCase)}).Count){
            throw "Runtime cache manifest is missing licenses in $folder"
        }
    }
    # Each cache member must be inventoried, including all bundled notices.
    foreach($file in Get-ChildItem -LiteralPath $Root -File -Recurse){
        $relative=$file.FullName.Substring($prefix.Length).Replace('\','/')
        if($relative -ne 'runtime-cache.json' -and -not $seen.ContainsKey($relative)){throw "Unlisted file in portable runtime cache: $relative"}
    }
}
try{
    if($Help){
        Write-Output @'
Usage: tools\release_portable.cmd [options]
  No arguments           Choose cases interactively, then package bin/LAMP.exe.
  -Cases "case1,case2"    Package only these formal cases; numbers/ranges also work.
  -RuntimeDirectory dir  Override the reusable .portable/runtime cache.
  -OutputDirectory dir   Override the default dist/ destination (on D:).
  -Version name          Set a new version; existing releases are never overwritten.

Examples:
  tools\release_portable.cmd
  tools\release_portable.cmd -Cases "hydro_static_drop,vof_a_charged_drop"
  tools\release_portable.cmd -Cases hydro_static_drop

Packaging uses the existing executable and Windows SDK mt.exe; it never invokes
a C++ compiler or CMake. The generated
portable ZIP contains only runtime files, licenses and the selected case inputs.
No cases are selected automatically. q or EOF cancels before packaging.
'@
        exit 0
    }
    $available=@(Get-ReleaseCases)
    if(-not $available.Count){throw 'No supported formal cases found in cases/.'}
    if(-not $PSBoundParameters.ContainsKey('Cases')){
        [Console]::WriteLine('Select cases for the portable package:')
        for($i=0;$i -lt $available.Count;$i++){[Console]::WriteLine('  ['+($i+1)+'] '+$available[$i])}
        [Console]::Write('Numbers/names separated by commas (example: 1,3-5); q cancels: ')
        $answer=[Console]::ReadLine()
        if($null -eq $answer -or [string]::IsNullOrWhiteSpace($answer) -or $answer.Trim() -ieq 'q'){
            Write-Output '[RELEASE] Cancelled; no package started.';exit 0
        }
        $Cases=$answer
    }
    $selected=@(Resolve-ReleaseCases $Cases $available)
    $program=Join-Path $project 'bin/LAMP.exe'
    if(-not (Test-Path -LiteralPath $program -PathType Leaf)){throw 'bin/LAMP.exe is missing; compile the application separately before packaging.'}
    if(-not $RuntimeDirectory){$RuntimeDirectory=Join-Path $project '.portable/runtime'}
    $RuntimeDirectory=(Get-Item -LiteralPath $RuntimeDirectory).FullName
    Assert-ReleaseRuntime $RuntimeDirectory
    if(-not $OutputDirectory){$OutputDirectory=Join-Path $project 'dist'}
    $OutputDirectory=[IO.Path]::GetFullPath($OutputDirectory)
    if($OutputDirectory -notmatch '\AD:[\\/]'){throw 'Store package output on D: as required by the project storage policy.'}
    if(-not $Version){$Version=(Get-Date -Format 'yyyyMMdd-HHmmss-fff')+'-'+[Guid]::NewGuid().ToString('N').Substring(0,8)+'-x64'}
    $destination=Join-Path $OutputDirectory ('LAMP-Portable-'+$Version)
    if((Test-Path -LiteralPath $destination) -or (Test-Path -LiteralPath ($destination+'.zip'))){throw "Release already exists: $destination"}
    $work=Join-Path $project ('runs/packaging/'+$Version+'-'+[Guid]::NewGuid().ToString('N').Substring(0,8))
    [void][IO.Directory]::CreateDirectory((Join-Path $work 'temp'))
    $env:TEMP=Join-Path $work 'temp';$env:TMP=$env:TEMP
    $recordPath=Join-Path $work 'release.json'
    $record=[ordered]@{status='preparing';created=(Get-Date).ToString('o');cases=$selected;mode='existing_binary';
        runtime_cache=$RuntimeDirectory;output_directory=$OutputDirectory;version=$Version}
    Save-ReleaseRecord
    Write-Output ('[RELEASE] Selected cases: '+($selected -join ', '))
    Write-Output '[RELEASE] Using the existing compiled bin/LAMP.exe.'
    $record['program_sha256']=Get-ReleaseSha $program
    $record.status='packaging';Save-ReleaseRecord
    $packageOptions=@{Program=$program;OutputDirectory=$OutputDirectory;Version=$Version;Cases=[string[]]$selected;
        OpenMpRuntime=(Join-Path $RuntimeDirectory 'bin/libomp140.x86_64.dll');
        OpenMpLicenseDirectory=(Join-Path $RuntimeDirectory 'licenses/OpenMP');
        VcRuntimeDirectory=(Join-Path $RuntimeDirectory 'bin');
        VcRuntimeLicensePath=(Join-Path $RuntimeDirectory 'licenses/CRT/LICENSE-VCRUNTIME.rtf');
        AdditionalLicenseDirectory=(Join-Path $RuntimeDirectory 'licenses/ThirdParty')}
    $global:LASTEXITCODE=0
    & (Join-Path $PSScriptRoot 'package_portable.ps1') @packageOptions | Tee-Object -FilePath (Join-Path $work 'package.log')
    if($LASTEXITCODE -ne 0){throw "Packaging failed ($LASTEXITCODE)."}
    $zip=$destination+'.zip'
    if(-not (Test-Path -LiteralPath $zip -PathType Leaf) -or -not (Test-Path -LiteralPath (Join-Path $destination 'portable_manifest.json'))){
        throw 'Packager returned without publishing both the ZIP and portable manifest.'
    }
    $record.status='completed';$record['zip']=$zip;$record['zip_sha256']=Get-ReleaseSha $zip
    $record['finished']=(Get-Date).ToString('o');Save-ReleaseRecord
    Write-Output ('[RELEASE] Ready: '+$zip)
    exit 0
}catch{
    if($record){$record.status='failed';$record['error']=$_.Exception.Message;Save-ReleaseRecord}
    [Console]::Error.WriteLine('[RELEASE] ERROR: '+$_.Exception.Message)
    exit 1
}finally{
    $env:TEMP=$savedTemp;$env:TMP=$savedTmp
}
