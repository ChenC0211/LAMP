# Pure PowerShell PE import reader. Packaging needs no Visual Studio/dumpbin.
function Read-PortablePe([string]$Path){
    $bytes=[IO.File]::ReadAllBytes($Path)
    function Assert-PeRange([long]$Offset,[long]$Size){
        if($Offset -lt 0 -or $Size -lt 0 -or $Offset -gt $bytes.LongLength-$Size){throw "Truncated PE structure: $Path"}
    }
    function Read-PeU16([long]$Offset){Assert-PeRange $Offset 2;return [BitConverter]::ToUInt16($bytes,[int]$Offset)}
    function Read-PeU32([long]$Offset){Assert-PeRange $Offset 4;return [BitConverter]::ToUInt32($bytes,[int]$Offset)}
    Assert-PeRange 0 64
    if((Read-PeU16 0) -ne 0x5a4d){throw "Not a Windows PE image: $Path"}
    $pe=[long](Read-PeU32 0x3c)
    Assert-PeRange $pe 24
    if((Read-PeU32 $pe) -ne 0x4550){throw "Invalid PE signature: $Path"}
    $machine=Read-PeU16 ($pe+4);$sectionCount=Read-PeU16 ($pe+6)
    $optionalSize=Read-PeU16 ($pe+20);$characteristics=Read-PeU16 ($pe+22);$optional=$pe+24
    Assert-PeRange $optional $optionalSize
    if($machine -ne 0x8664 -or $optionalSize -lt 112 -or (Read-PeU16 $optional) -ne 0x20b){throw "Portable release requires a Windows x64 PE32+ image: $Path"}
    $directoryCount=Read-PeU32 ($optional+108)
    if($directoryCount -gt 16 -or $optionalSize -lt 112+8*$directoryCount){throw "Invalid PE data directory: $Path"}
    $headerSize=Read-PeU32 ($optional+60);$sectionStart=$optional+$optionalSize
    if($sectionCount -lt 1 -or $sectionCount -gt 96){throw "Invalid PE section count: $Path"}
    Assert-PeRange $sectionStart (40*$sectionCount)
    $sections=@(for($i=0;$i -lt $sectionCount;$i++){
        $s=$sectionStart+40*$i
        $section=[pscustomobject]@{Name=[Text.Encoding]::ASCII.GetString($bytes,[int]$s,8).TrimEnd([char]0);
            Virtual=[long](Read-PeU32 ($s+12));VirtualSize=[long](Read-PeU32 ($s+8));RawSize=[long](Read-PeU32 ($s+16));Raw=[long](Read-PeU32 ($s+20));
            Characteristics=[long](Read-PeU32 ($s+36));Sha256=''}
        Assert-PeRange $section.Raw $section.RawSize
        $hash=[Security.Cryptography.SHA256]::Create()
        try{$section.Sha256=[BitConverter]::ToString($hash.ComputeHash($bytes,[int]$section.Raw,[int]$section.RawSize)).Replace('-','').ToLowerInvariant()}
        finally{$hash.Dispose()}
        $section
    })
    function Resolve-PeRva([long]$Rva,[long]$Size){
        if($Rva -lt $headerSize){Assert-PeRange $Rva $Size;return $Rva}
        foreach($s in $sections){
            if($Rva -ge $s.Virtual -and $Rva-$s.Virtual -lt [Math]::Max($s.VirtualSize,$s.RawSize)){
                $delta=$Rva-$s.Virtual
                if($Size -gt $s.RawSize-$delta){throw "PE RVA points outside file-backed section: $Path"}
                $offset=$s.Raw+$delta;Assert-PeRange $offset $Size;return $offset
            }
        }
        throw "Unmapped PE RVA: $Path"
    }
    function Read-PeDllName([long]$Rva){
        $offset=Resolve-PeRva $Rva 1;$name=[Text.StringBuilder]::new()
        for($j=0;$j -lt 260;$j++){
            # Resolve each byte so an unterminated name cannot cross raw data.
            $at=Resolve-PeRva ($Rva+$j) 1;$value=$bytes[$at]
            if($value -eq 0){
                $result=$name.ToString()
                if($result -notmatch '\A[A-Za-z0-9_.-]+\.dll\z'){throw "Invalid PE DLL import name '$result': $Path"}
                return $result.ToLowerInvariant()
            }
            if($value -gt 127 -or $value -lt 33){throw "Invalid PE DLL import bytes: $Path"}
            [void]$name.Append([char]$value)
        }
        throw "Unterminated PE import name: $Path"
    }
    $imports=[Collections.Generic.List[string]]::new()
    foreach($entry in @(@{Index=1;Stride=20;NameOffset=12},@{Index=13;Stride=32;NameOffset=4})){
        if($directoryCount -le $entry.Index){continue}
        $directory=$optional+112+8*$entry.Index
        $rva=[long](Read-PeU32 $directory);$size=[long](Read-PeU32 ($directory+4))
        if($rva -eq 0 -and $size -eq 0){continue}
        if($rva -eq 0 -or $size -lt $entry.Stride){throw "Invalid PE imports directory: $Path"}
        $terminated=$false
        for($relative=0;$relative -le $size-$entry.Stride;$relative+=$entry.Stride){
            $descriptor=Resolve-PeRva ($rva+$relative) $entry.Stride
            $allZero=$true
            for($j=0;$j -lt $entry.Stride;$j+=4){if((Read-PeU32 ($descriptor+$j)) -ne 0){$allZero=$false}}
            if($allZero){$terminated=$true;break}
            if($entry.Index -eq 13 -and (Read-PeU32 $descriptor) -ne 1){throw "Unsupported legacy VA delay imports in x64 image: $Path"}
            $imports.Add((Read-PeDllName (Read-PeU32 ($descriptor+$entry.NameOffset))))
        }
        if(-not $terminated){throw "Unterminated PE imports directory: $Path"}
    }
    $version=[Diagnostics.FileVersionInfo]::GetVersionInfo($Path)
    $resourceSection=''
    if($directoryCount -gt 2){
        $resourceRva=[long](Read-PeU32 ($optional+112+16))
        if($resourceRva){foreach($section in $sections){if($resourceRva -ge $section.Virtual -and $resourceRva-$section.Virtual -lt [Math]::Max($section.VirtualSize,$section.RawSize)){$resourceSection=$section.Name}}}
    }
    [pscustomobject]@{Path=[IO.Path]::GetFullPath($Path);Machine='AMD64';IsDll=[bool]($characteristics -band 0x2000);
        Sections=$sections;ResourceSection=$resourceSection;
        Imports=@($imports|Sort-Object -Unique);FileVersion=$version.FileVersion;ProductVersion=$version.ProductVersion;
        Version=[version]::new([Math]::Max(0,$version.FileMajorPart),[Math]::Max(0,$version.FileMinorPart),[Math]::Max(0,$version.FileBuildPart),[Math]::Max(0,$version.FilePrivatePart))}
}
function Assert-PortableResourceOnlyTransform($Original,$Transformed){
    if(-not $Original.ResourceSection -or $Original.ResourceSection -ne $Transformed.ResourceSection){throw 'Manifest deployment changed or lost the resource section identity.'}
    if(($Original.Imports -join '|') -cne ($Transformed.Imports -join '|')){throw 'Manifest deployment changed executable imports.'}
    $before=@($Original.Sections|Where-Object {$_.Name -cne $Original.ResourceSection})
    $after=@($Transformed.Sections|Where-Object {$_.Name -cne $Transformed.ResourceSection})
    if($before.Count -ne $after.Count){throw 'Manifest deployment changed non-resource section count.'}
    for($i=0;$i -lt $before.Count;$i++){
        foreach($property in @('Name','Virtual','VirtualSize','RawSize','Characteristics','Sha256')){
            if($before[$i].$property -cne $after[$i].$property){throw "Manifest deployment changed non-resource section $($before[$i].Name): $property"}
        }
    }
}
function Test-PortableSystemDll([string]$Name){
    if($Name -match '\A(?:api-ms-win-|ext-ms-win-)[a-z0-9_.-]+\.dll\z'){return $true}
    return $Name -in @('kernel32.dll','kernelbase.dll','ntdll.dll','advapi32.dll','user32.dll','gdi32.dll','gdi32full.dll',
        'shell32.dll','shlwapi.dll','ole32.dll','oleaut32.dll','combase.dll','comdlg32.dll','rpcrt4.dll','ws2_32.dll',
        'bcrypt.dll','bcryptprimitives.dll','crypt32.dll','secur32.dll','version.dll','winmm.dll','imm32.dll','setupapi.dll',
        'cfgmgr32.dll','userenv.dll','psapi.dll','dbghelp.dll','powrprof.dll','normaliz.dll','msvcrt.dll','ucrtbase.dll','wintrust.dll','imagehlp.dll')
}
function Assert-PortableReleaseImports($Image){
    foreach($name in $Image.Imports){
        if($name -match '\A(?:msvcp|msvcr|vcruntime|concrt|vccorlib)[0-9_]+d(?:_[a-z0-9_]+)?\.dll\z' -or $name -eq 'ucrtbased.dll'){
            throw "Debug CRT is not supported in a portable release: $name in $($Image.Path)"
        }
    }
}
