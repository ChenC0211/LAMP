# The run launcher freezes these input directories for each fresh run.
function Assert-PortableInputResources([string]$CaseRoot){
    $allowed=@('system','constant','input','data','reference')
    foreach($dictionary in Get-ChildItem -LiteralPath (Join-Path $CaseRoot 'system') -Filter '*.dict' -File -Recurse){
        $text=[IO.File]::ReadAllText($dictionary.FullName,[Text.Encoding]::UTF8)
        $tokens=@([regex]::Matches($text,'#[^\r\n]*|//[^\r\n]*|"[^"\r\n]*"|''[^''\r\n]*''|[^\s{}();,]+') |
            Where-Object {-not ($_.Value.StartsWith('#') -or $_.Value.StartsWith('//'))})
        for($i=0;$i -lt $tokens.Count;$i++){
            $value=$tokens[$i].Value.Trim('"',"'")
            if($value -match '\A(?:[A-Za-z]:[\\/]|[\\/]{2}|/)' -or $value -match '\A(?:file|https?)://'){
                throw "Portable inputs cannot reference absolute/external paths: $($dictionary.FullName): $value"
            }
            $key=if($i -gt 0){$tokens[$i-1].Value}else{''}
            $pathKey=$key -match '(?i)\A(?:file|filename|path|directory|table|meshfile|geometryfile|datafile|sourcefile)\z'
            $resource=$value -match '(?i)\.(?:stl|obj|off|ply|vtk|vtu|msh|mesh|csv|dat|txt|json|bin|h5|hdf5|npy|dict)\z'
            $relative=$value -match '\A\.\.?[\\/]'
            if(-not ($pathKey -or $resource -or $relative)){continue}
            if($value -match '\A[+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?\z'){continue}
            $resolved=[IO.Path]::GetFullPath((Join-Path $CaseRoot $value))
            $prefix=$CaseRoot.TrimEnd('\','/')+[IO.Path]::DirectorySeparatorChar
            if(-not $resolved.StartsWith($prefix,[StringComparison]::OrdinalIgnoreCase)){throw "Input resource escapes its case: $value"}
            $relativeName=$resolved.Substring($prefix.Length).Replace('\','/')
            $folder=$relativeName.Split('/')[0]
            if($folder -notin $allowed){throw "Resource must be in system/constant/input/data/reference so the launcher freezes it: $value"}
            if(-not (Test-Path -LiteralPath $resolved)){throw "Missing case resource: $value in $($dictionary.FullName)"}
        }
    }
}
