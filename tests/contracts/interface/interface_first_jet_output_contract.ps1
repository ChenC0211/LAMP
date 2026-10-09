param([Parameter(Mandatory=$true)][string]$Program,[Parameter(Mandatory=$true)][string]$RunRoot)
$ErrorActionPreference='Stop'
$culture=[Globalization.CultureInfo]::InvariantCulture
function Assert([bool]$ok,[string]$why){if(!$ok){throw $why}}
function Number([string]$s){[double]::Parse($s,$culture)}
function Read-Log([string]$p){(Get-Content -Raw -LiteralPath $p).Replace("`r`n","`n")}
function Write-Text([string]$p,[string]$s){[IO.File]::WriteAllText($p,$s,[Text.UTF8Encoding]::new($false))}
function Case([string]$name,[hashtable]$body){
 $path=Join-Path $RunRoot $name
 New-Item -ItemType Directory -Path "$path/system" -Force|Out-Null
 foreach($entry in $body.GetEnumerator()){Write-Text "$path/system/$($entry.Key)" $entry.Value}
 return $path
}
function Fixture([string]$name){
 $map=@{staticCircleXY='static_xy';staticSphereRZ='static_rz';translationCircleXY='translation_xy';translationSphereRZ='translation_rz';solidRotationXY='rotation_xy';polynomialVortexXY='vortex_xy'}
 $folder=Join-Path $PSScriptRoot ('../../fixtures/interface_transport/'+$map[$name]+'/system')
 $body=@{};foreach($file in Get-ChildItem -LiteralPath $folder -Filter '*.dict' -File){$body[$file.Name]=[IO.File]::ReadAllText($file.FullName)}
 return $body
}
function Invoke-Case([string]$path,[bool]$success=$true){
 $ErrorActionPreference='Continue'
 & $Program --case $path *> "$path/run.log"
 $rc=$LASTEXITCODE
 $ErrorActionPreference='Stop'
 Write-Text "$path/exit.txt" "$rc`n"
 Assert (($rc -eq 0) -eq $success) "Unexpected run status $rc at $path"
 return Read-Log "$path/run.log"
}
function Array-Data($piece,[string]$xpath){
 $node=$piece.SelectSingleNode($xpath);Assert ($null -ne $node) "Missing VTK array $xpath"
 return @($node.InnerText.Trim() -split '\s+'|Where-Object {$_ -ne ''}|ForEach-Object {Number $_})
}
function Numeric-Output([string]$out){
 return @([regex]::Matches($out,'(?m)^(\[Step [^\n]+|FIRST_JET_APP_(?:RESULT|REGRID),[^\n]+)$')|ForEach-Object {$_.Value}) -join "`n"
}
$framesChecked=0;$vtkFiles=0;$pointSamples=0;$frequencyCases=0
foreach($problem in @('staticCircleXY','staticSphereRZ','translationCircleXY','translationSphereRZ','solidRotationXY','polynomialVortexXY')){
 $radial=$problem -like '*RZ';$static=$problem -like 'static*'
 $dt=if($static){'.005'}elseif($problem -eq 'polynomialVortexXY'){'.001'}else{'.01'}
 $mesh=if($radial){'amr'}else{'uniform'}
 $base=Fixture $problem
 $base['mesh.dict']=$base['mesh.dict'].Replace('initialMesh uniform;',"initialMesh $mesh; Refinement { min (.25 .25); max (.75 .75); level 1; } Regrid { atSteps (1 2 3); }")
 $baseline=$null;$curvatureBaseline=$null
 foreach($frequency in @(0,1,3)){
  $body=$base.Clone()
  if($frequency -gt 0){$body['run.dict']=$body['run.dict'].Replace('enabled false;',"enabled true; writeEverySteps $frequency;")}
  $path=Case "$($problem)_$frequency" $body
  $out=Invoke-Case $path
  $numeric=Numeric-Output $out
  if($frequency -eq 0){$baseline=$numeric}else{Assert ($numeric -ceq $baseline) "Output schedule changed numeric result: $problem/$frequency"}
  ++$frequencyCases
  $meta="$path/output/meta";$series="$path/output/series"
  $manifest=Get-Content -Raw -Encoding UTF8 "$meta/manifest.json"|ConvertFrom-Json
  Assert ($manifest.exit_status -eq 'success' -and $manifest.accepted_steps -eq 4 -and $manifest.final_time -eq 4*(Number $dt)) 'Wrong successful manifest'
  Assert ($manifest.selected_methods.transport -eq 'FirstJet' -and $null -eq $manifest.inventory -and $null -eq $manifest.phase_flux_ledger) 'False inventory metadata'
  Assert ($manifest.git_commit -eq 'unavailable' -and $manifest.compiler -and $manifest.build_type -and $manifest.run_id -match '^[a-f0-9]{64}$' -and
   $manifest.start_time -match '^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$' -and $manifest.case_path) 'Missing run/build provenance'
  Assert (@($manifest.input_file_hashes.PSObject.Properties).Count -eq 5) 'Expected all five input hashes'
  foreach($prop in $manifest.input_file_hashes.PSObject.Properties){Assert ($prop.Value -eq 'sha256:'+(Get-FileHash -LiteralPath "$path/$($prop.Name)").Hash.ToLowerInvariant()) 'Input SHA-256 mismatch'}
  $trace=@(Import-Csv "$meta/configuration_trace.csv")
  Assert (@($trace|Where-Object {$_.used -ne 'true' -or !$_.consumer -or !$_.units -or !$_.file -or [int]$_.line -lt 1}).Count -eq 0) 'Unused/missing configuration provenance'
  Assert (@($trace|Where-Object {$_.field -eq 'Geometry.coordinates' -and $_.effective_value -eq $(if($radial){'rz'}else{'xy'})}).Count -eq 1) 'Coordinate provenance missing'
  Assert ((Read-Log "$meta/effective_config.dict").Contains('Initial { phi "') -and !(Read-Log "$meta/effective_config.dict").Contains('problem ')) 'Generic effective data absent'
  $rows=@(Import-Csv "$series/conservation.csv");Assert ($rows.Count -eq 5) 'Missing initial/accepted conservation row'
  foreach($r in $rows){
   Assert ([int]$r.curvature_quadrature_samples -gt 0 -and [int]$r.curvature_unavailable_samples -eq 0 -and
    (Number $r.curvature_measure) -gt 0 -and (Number $r.kappa_rms) -ge 0 -and
    (Number $r.kappa_rms) -le (Number $r.kappa_max)+1e-12) 'Missing/invalid same-publication curvature statistics'
   if($static -or $problem -like 'translation*'){
    $expected=if($radial){2/.27}else{1/.19}
    Assert ([Math]::Abs((Number $r.kappa_rms)-$expected) -lt 1e-10) 'XY/RZ analytic curvature RMS mismatch'
   }
  }
  $curvatureRows=($rows|ForEach-Object {"$($_.step),$($_.geometry_identity),$($_.kappa_rms),$($_.curvature_quadrature_samples),$($_.curvature_unavailable_samples),$($_.curvature_measure)"}) -join "`n"
  if($frequency -eq 0){$curvatureBaseline=$curvatureRows}else{Assert ($curvatureRows -ceq $curvatureBaseline) 'Output schedule changed same-publication curvature diagnostics'}
  $logs=@([regex]::Matches($out,'(?m)^\[Step (\d+)\] t=(\S+) dt=(\S+) cells=(\d+) geometry=(\d+) V0=(\S+) V=(\S+) dV/V=(\S+) kappaMax=(\S+) inventory=NA flux_ledger=NA$'))
  Assert ($logs.Count -eq 4) 'Missing stdout accepted rows'
  for($k=1;$k -le 4;++$k){$g=$logs[$k-1].Groups;$r=$rows[$k]
   foreach($pair in @(@('step',1),@('time',2),@('dt',3),@('geometry_identity',5),@('initial_phase_volume',6),@('phase_volume',7),@('relative_volume_drift',8),@('kappa_max',9))){
    Assert ($r.($pair[0]) -ceq $g[$pair[1]].Value) "Series/stdout difference: $($pair[0])"
   }
  }
  foreach($r in $rows){foreach($name in @('analytic_volume','absolute_volume_error','phase_volume_boundary_flux','phase_volume_closure','mass','mass_boundary_flux','mass_closure','divergence_norm','momentum_balance_norm','inventory','flux_ledger')){Assert ($r.$name -eq 'NA') "False physical field $name"}}
  $solver=@(Import-Csv "$series/solver.csv");Assert ($solver.Count -eq 5 -and @($solver|Where-Object {$_.status -ne 'accepted' -or $_.linear_iterations -ne 'NA'}).Count -eq 0) 'False solver diagnostics'
  $events=@(Import-Csv "$series/events.csv");Assert ($events.Count -eq 3 -and @($events|Where-Object {$_.event -ne 'regrid'}).Count -eq 0) 'Regrid events missing'
  $perf=@(Import-Csv "$series/performance.csv");Assert ($perf.Count -eq 52) 'Missing App/producer stage performance'
  foreach($p in $perf){Assert ((Number $p.wall_seconds) -ge 0 -and [int]$p.leaf_cells -gt 0 -and [int]$p.native_faces -gt 0 -and
   [int]$p.primary_vertices -gt 0 -and [int]$p.primary_edges -gt 0 -and [int]$p.centers -eq [int]$p.leaf_cells -and $p.nnz -eq 'NA') 'Invalid structural/performance data'}
  $phases=@('interface_input','interface_boundary','interface_vertices','interface_edges','interface_centers','interface_mixed_closure','interface_publication')
  for($step=1;$step -le 4;++$step){$byStage=@{};foreach($p in @($perf|Where-Object {[int]$_.step -eq $step})){
    Assert (!$byStage.ContainsKey($p.stage)) 'Duplicate timer';$byStage[$p.stage]=Number $p.wall_seconds
   }
   $sum=0.;foreach($phase in $phases){Assert ($byStage.ContainsKey($phase)) 'Missing producer phase';$sum+=$byStage[$phase]}
   Assert ($sum -le $byStage.interface_producer_total+1e-9 -and $byStage.interface_producer_total -le $byStage.interface_step+1e-9 -and
    $byStage.interface_step+$byStage.interface_mesh_transfer+$byStage.io -le $byStage.step_total+1e-9) 'Timer partition/double accounting error'
   $source=@($perf|Where-Object {[int]$_.step -eq $step-1 -and $_.stage -eq 'io'})[0]
   foreach($p in @($perf|Where-Object {[int]$_.step -eq $step -and ($_.stage -in $phases -or $_.stage -in @('interface_step','interface_producer_total'))})){
    foreach($name in @('leaf_cells','native_faces','cut_cells','interface_panels','primary_vertices','primary_edges','centers')){
     Assert ($p.$name -eq $source.$name) 'Producer timer incorrectly labeled with post-regrid geometry size'
    }
   }
  }
  $accuracy=Get-Content -Raw -Encoding UTF8 "$meta/accuracy.json"|ConvertFrom-Json
  $tokens=[regex]::Match($out,'(?m)^FIRST_JET_APP_RESULT,([^\n]+)').Groups[1].Value.Split(',')
  foreach($name in @('reference','position','gradient','curvature','absolute_volume_error')){Assert ($null -eq $accuracy.$name) 'False analytic reference metadata'}
  Assert ($accuracy.peak_relative_volume_drift -eq (Number $tokens[9]) -and $tokens[5] -eq 'NA' -and $tokens[11] -eq 'COMPLETED') 'Final diagnostics/completion mismatch'
  $frames=@(Import-Csv "$meta/frames.csv")
  if($frequency -eq 0){Assert ($frames.Count -eq 0 -and $manifest.accepted_frames -eq 0 -and !(Test-Path "$path/output/fields/solution.pvd")) 'Fields disabled but registered';continue}
  $expectedSteps=if($frequency -eq 1){@(0,1,2,3,4)}else{@(0,3,4)}
  Assert ($frames.Count -eq $expectedSteps.Count -and $manifest.accepted_frames -eq $frames.Count) 'Incorrect frame frequency'
  [xml]$pvd=Get-Content -Raw -Encoding UTF8 "$path/output/fields/solution.pvd"
  $datasets=@($pvd.VTKFile.Collection.DataSet);Assert ($datasets.Count -eq 2*$frames.Count) 'Missing PVD groups'
  for($j=0;$j -lt $frames.Count;++$j){$f=$frames[$j];$r=$rows[[int]$f.step]
   Assert ([int]$f.frame -eq $j -and [int]$f.step -eq $expectedSteps[$j] -and $f.time -ceq $r.time -and $f.geometry_identity -ceq $r.geometry_identity) 'Frame identity/time mismatch'
   foreach($group in @(0,1)){$ds=$datasets[2*$j+$group];Assert ($ds.timestep -ceq $f.time -and (Test-Path "$path/output/fields/$($ds.file)")) 'Invalid PVD reference'}
   [xml]$vtu=Get-Content -Raw -Encoding UTF8 "$path/output/$($f.volume)"
   [xml]$vtp=Get-Content -Raw -Encoding UTF8 "$path/output/$($f.interface)"
   foreach($x in @($vtu,$vtp)){foreach($pair in @(@('TimeValue',$f.time),@('geometry_identity',$f.geometry_identity),@('accepted_step',$f.step),@('mesh_revision',$f.mesh_revision))){
    Assert ($x.SelectSingleNode("//FieldData/DataArray[@Name='$($pair[0])']").InnerText -ceq $pair[1]) 'VTK FieldData mismatch'
   }}
   $v=$vtu.VTKFile.UnstructuredGrid.Piece;$cells=[int]$v.NumberOfCells;$vertices=[int]$v.NumberOfPoints
   $points=Array-Data $v "Points/DataArray";$conn=Array-Data $v "Cells/DataArray[@Name='connectivity']";$offsets=Array-Data $v "Cells/DataArray[@Name='offsets']";$types=Array-Data $v "Cells/DataArray[@Name='types']"
   Assert ($points.Count -eq 3*$vertices -and $conn.Count -eq 4*$cells -and $offsets.Count -eq $cells -and $types.Count -eq $cells) 'Broken native grid arrays'
   for($k=0;$k -lt $cells;++$k){Assert ($offsets[$k] -eq 4*($k+1) -and $types[$k] -eq 9) 'Bad quad connectivity'}
   Assert (@($conn|Where-Object {$_ -lt 0 -or $_ -ge $vertices -or $_ -ne [int]$_}).Count -eq 0) 'Bad native vertex id'
   $vol=Array-Data $v "CellData/DataArray[@Name='geometric_liquid_volume']";$cap=Array-Data $v "CellData/DataArray[@Name='cell_capacity']";$alpha=Array-Data $v "CellData/DataArray[@Name='phase_volume_fraction']"
   foreach($arr in @($vol,$cap,$alpha)){Assert ($arr.Count -eq $cells) 'Wrong cell data count'}
   $sum=0.;$capSum=0.;for($k=0;$k -lt $cells;++$k){Assert ($cap[$k] -gt 0 -and [Math]::Abs($alpha[$k]*$cap[$k]-$vol[$k]) -lt 1e-14) 'Inconsistent geometric volume fraction';$sum+=$vol[$k];$capSum+=$cap[$k]}
   Assert ([Math]::Abs($sum-(Number $r.phase_volume)) -lt 1e-13) 'VTU and publication sum differ'
   $expectedCapacity=if($radial){[Math]::PI}else{1.};Assert ([Math]::Abs($capSum-$expectedCapacity) -lt 1e-12) 'XY/RZ capacity wrong'
   foreach($name in @('level','cell_id','phi_center')){Assert ((Array-Data $v "CellData/DataArray[@Name='$name']").Count -eq $cells) 'Missing cell field'}
   $v=$vtp.VTKFile.PolyData.Piece;$np=[int]$v.NumberOfPoints;$nl=[int]$v.NumberOfLines
   $p=Array-Data $v "Points/DataArray";$n=Array-Data $v "PointData/DataArray[@Name='normal']";$kappa=Array-Data $v "PointData/DataArray[@Name='curvature']"
   $meridional=Array-Data $v "PointData/DataArray[@Name='meridional_curvature']";$azimuthal=Array-Data $v "PointData/DataArray[@Name='azimuthal_curvature']"
   $res=Array-Data $v "PointData/DataArray[@Name='level_set_residual']";$indices=Array-Data $v "Lines/DataArray[@Name='connectivity']";$ends=Array-Data $v "Lines/DataArray[@Name='offsets']"
   Assert ($np -eq 9*$nl -and $p.Count -eq 3*$np -and $n.Count -eq 3*$np -and $kappa.Count -eq $np -and $indices.Count -eq $np -and $ends.Count -eq $nl) 'Broken interface polyline arrays'
   for($k=0;$k -lt $np;++$k){Assert ([double]::IsNaN($kappa[$k]) -eq $false -and [double]::IsInfinity($kappa[$k]) -eq $false -and
    [Math]::Abs($n[3*$k]*$n[3*$k]+$n[3*$k+1]*$n[3*$k+1]-1) -lt 1e-12 -and $n[3*$k+2] -eq 0 -and
    [Math]::Abs($kappa[$k]-$meridional[$k]-$azimuthal[$k]) -lt 1e-12 -and $res[$k] -le 1e-9 -and $indices[$k] -eq $k) 'Invalid same-scalar geometry data'
    if($static -or $problem -like 'translation*'){
     $radius=if($radial){.27}else{.19};$cx=if($radial){0}else{.43};$cy=.51
     if($problem -like 'translation*'){if($radial){$cy+=.2*(Number $f.time)}else{$cx+=.2*(Number $f.time)}}
     $x=$p[3*$k]-$cx;$y=$p[3*$k+1]-$cy;$factor=if($radial){2}else{1}
     Assert ([Math]::Abs([Math]::Sqrt($x*$x+$y*$y)-$radius) -lt 1e-12 -and [Math]::Abs($kappa[$k]-$factor/$radius) -lt 1e-10 -and
      $n[3*$k]*$x+$n[3*$k+1]*$y -gt 0) 'Wrong interface position/curvature/orientation'
    }
   }
   for($k=0;$k -lt $nl;++$k){Assert ($ends[$k] -eq 9*($k+1)) 'Bad interface offsets'}
   ++$framesChecked;$vtkFiles+=2;$pointSamples+=$np
  }
 }
}
# New shape, non-unit domain and affine shear require no C++ enum.
# Five files are prepared directly; the test never adapts an old schema.
$ellipse=@{
 'run.dict'='Run { app Interface; } Time { units code; dt .01; endTime .04; } Output { Fields { enabled true; } }'
 'geometry.dict'='Geometry { coordinates xy; Domain { min (-.2 -.1); max (1.8 1.9); } Initial { phi "((x-cx)/ax)^2+((y-cy)/ay)^2-1"; } }'
 'mesh.dict'='Mesh { resolution 8; }'
 'fields.dict'='Fields { Constants { ax .32; ay .27; cx .7; cy .8; } Velocity { type constant; value (.07 0); } }'
 'interface.dict'='Interface { method hermiteGals; Boundary { type prescribed; phi "((x-cx-.07*t)/ax)^2+((y-cy)/ay)^2-1"; } }'
}
$constant=$null
foreach($provider in @('constant','polynomial','shear')){
 $body=$ellipse.Clone()
 if($provider -eq 'polynomial'){$body['fields.dict']=$body['fields.dict'].Replace('type constant; value (.07 0);','type polynomial; origin (.2 .3); scale (2 3); Coefficients { c00 { x ".07"; y "0"; } }')}
 if($provider -eq 'shear'){
  $body['fields.dict']=$body['fields.dict'].Replace('type constant; value (.07 0);','type polynomial; origin (.2 .3); scale (2 3); Coefficients { c00 { x ".055"; y "0"; } c01 { x ".09"; y "0"; } }')
  $body['interface.dict']=$body['interface.dict'].Replace('x-cx-.07*t','x-cx-(.07+.03*(y-cy))*t')
 }
 $path=Case "new_ellipse_$provider" $body;$null=Invoke-Case $path
 $frames=@(Import-Csv "$path/output/meta/frames.csv");$f=$frames[-1]
 [xml]$surface=Get-Content -Raw -Encoding UTF8 "$path/output/$($f.interface)"
 $piece=$surface.VTKFile.PolyData.Piece;$points=Array-Data $piece 'Points/DataArray'
 for($k=0;$k -lt $points.Count/3;++$k){
  $speed=if($provider -eq 'shear'){.07+.03*($points[3*$k+1]-.8)}else{.07}
  $x=($points[3*$k]-.7-$speed*.04)/.32;$y=($points[3*$k+1]-.8)/.27
  Assert ([Math]::Abs($x*$x+$y*$y-1) -lt 1e-11) 'Generic new ellipse failed independent translation reference'
 }
 $rows=@(Import-Csv "$path/output/series/conservation.csv");Assert ([Math]::Abs((Number $rows[-1].phase_volume)-[Math]::PI*.32*.27) -lt 1e-12) 'Generic ellipse area mismatch'
 if($provider -eq 'constant'){$constant=$points}elseif($provider -eq 'polynomial'){Assert ($points.Count -eq $constant.Count) 'Constant/polynomial topology mismatch';for($k=0;$k -lt $points.Count;++$k){Assert ([Math]::Abs($points[$k]-$constant[$k]) -lt 1e-12) 'Constant/polynomial provider mismatch'}}
 $replay=Join-Path $RunRoot "replay_$provider";New-Item -ItemType Directory -Path "$replay/system" -Force|Out-Null
 $effective=@(Get-ChildItem -LiteralPath "$path/output/meta/effective" -Filter '*.dict' -File)
 Assert ($effective.Count -eq 5) 'Effective configuration lost file ownership'
 foreach($file in $effective){Copy-Item -LiteralPath $file.FullName -Destination "$replay/system/$($file.Name)"}
 $ErrorActionPreference='Continue';& $Program --case $replay --check-case *> "$replay/check.log";$rc=$LASTEXITCODE;$ErrorActionPreference='Stop'
 Assert ($rc -eq 0 -and !(Test-Path "$replay/output")) 'Effective five-file configuration cannot be read again'
}
# Runtime numerical rejection retains initial accepted frame and time.
$base=Fixture 'translationCircleXY'
$base['run.dict']=$base['run.dict'].Replace('dt 0.01;','dt 1;').Replace('endTime 0.040000000000000001;','endTime 2;').Replace('enabled false;','enabled true;')
$path=Case 'cfl_rejected' $base
$out=Invoke-Case $path $false
$failure=Get-Content -Raw -Encoding UTF8 "$path/output/meta/failure.json"|ConvertFrom-Json
$manifest=Get-Content -Raw -Encoding UTF8 "$path/output/meta/manifest.json"|ConvertFrom-Json
$frames=@(Import-Csv "$path/output/meta/frames.csv")
Assert ($failure.accepted_step -eq 0 -and $failure.accepted_time -eq 0 -and $failure.attempted_step -eq 1 -and $failure.reason -match 'half a cell' -and
 $manifest.exit_status -eq 'failed' -and $manifest.accepted_steps -eq 0 -and $frames.Count -eq 1 -and $frames[0].step -eq '0' -and $out -notmatch '\[Step ') 'False accepted step on rejection'
# A fresh-root violation preserves the existing run byte for byte.
$existing=@(Get-ChildItem -LiteralPath "$path/output" -File -Recurse|ForEach-Object {[pscustomobject]@{path=$_.FullName;hash=(Get-FileHash $_.FullName).Hash}})
$out=Invoke-Case $path $false
Assert ($out -match 'fresh run directory') 'Existing output not rejected'
foreach($file in $existing){Assert ((Get-FileHash -LiteralPath $file.path).Hash -eq $file.hash) 'Existing run overwritten'}
# Inject a real filesystem write failure after the run has accepted its initial
# frame. This is an external test action, not a numerical/debug input option.
$body=Fixture 'staticCircleXY';$body['mesh.dict']=$body['mesh.dict'].Replace('resolution 8;','resolution 16;')
$body['run.dict']=$body['run.dict'].Replace('endTime 0.02;','endTime .04;').Replace('enabled false;','enabled true;')
$path=Case 'io_rejected' $body
$process=Start-Process -FilePath $Program -ArgumentList @('--case',('"'+$path+'"')) -PassThru -WindowStyle Hidden -RedirectStandardOutput "$path/run.log" -RedirectStandardError "$path/error.log"
$null=$process.Handle # Retain the child handle so Windows PowerShell can read its exit code after exit.
$deadline=[DateTime]::UtcNow.AddSeconds(60)
while(!(Test-Path "$path/output/fields/solution.pvd") -and !$process.HasExited -and [DateTime]::UtcNow -lt $deadline){Start-Sleep -Milliseconds 10;$process.Refresh()}
Assert (!(Test-Path "$path/output/fields/volume/volume_000003.vtu") -and !$process.HasExited) 'IO failure fixture missed injection window'
New-Item -ItemType Directory -Path "$path/output/fields/volume/volume_000003.vtu"|Out-Null
Assert ($process.WaitForExit(60000)) 'IO fixture timed out'
$process.Refresh();Write-Text "$path/exit.txt" "$($process.ExitCode)`n"
Assert ($null -ne $process.ExitCode -and $process.ExitCode -ne 0) 'IO failure swallowed or exit code unavailable'
$failure=Get-Content -Raw -Encoding UTF8 "$path/output/meta/failure.json"|ConvertFrom-Json
$manifest=Get-Content -Raw -Encoding UTF8 "$path/output/meta/manifest.json"|ConvertFrom-Json
$frames=@(Import-Csv "$path/output/meta/frames.csv")
Assert ($failure.stage -eq 'accepted output' -and $failure.accepted_step -eq 3 -and $manifest.accepted_steps -eq 3 -and $manifest.exit_status -eq 'failed' -and
 $frames.Count -eq 3 -and $frames[-1].step -eq '2' -and $failure.reason -match 'volume_000003.vtu') 'IO failure corrupted numerical/frame clocks'
[xml]$pvd=Get-Content -Raw -Encoding UTF8 "$path/output/fields/solution.pvd"
Assert (@($pvd.VTKFile.Collection.DataSet).Count -eq 6) 'Failed payload registered in PVD'
# Exercise the Windows sharing-conflict retry itself. Hold a real read handle
# without delete sharing until the writer reports a retry, then release it.
$body=Fixture 'staticCircleXY';$body['mesh.dict']=$body['mesh.dict'].Replace('resolution 8;','resolution 16;')
$body['run.dict']=$body['run.dict'].Replace('enabled false;','enabled true;')
$path=Case 'atomic_retry' $body
$process=Start-Process -FilePath $Program -ArgumentList @('--case',('"'+$path+'"')) -PassThru -WindowStyle Hidden -RedirectStandardOutput "$path/run.log" -RedirectStandardError "$path/error.log"
$null=$process.Handle
$deadline=[DateTime]::UtcNow.AddSeconds(60)
while(!(Test-Path "$path/output/fields/solution.pvd") -and !$process.HasExited -and [DateTime]::UtcNow -lt $deadline){Start-Sleep -Milliseconds 10;$process.Refresh()}
Assert (!$process.HasExited) 'Atomic retry fixture missed live run'
$lock=[IO.File]::Open("$path/output/meta/manifest.json",[IO.FileMode]::Open,[IO.FileAccess]::Read,[IO.FileShare]::Read)
$sawRetry=$false
try{
 while(!$process.HasExited -and [DateTime]::UtcNow -lt $deadline){
  if((Get-Content -Raw "$path/error.log") -match 'atomic replace retry='){$sawRetry=$true;break}
  Start-Sleep -Milliseconds 10;$process.Refresh()
 }
}finally{$lock.Dispose()}
Assert ($sawRetry) 'Sharing conflict did not exercise bounded retry'
Assert ($process.WaitForExit(60000)) 'Atomic retry fixture timed out'
$process.Refresh();Write-Text "$path/exit.txt" "$($process.ExitCode)`n"
$manifest=Get-Content -Raw -Encoding UTF8 "$path/output/meta/manifest.json"|ConvertFrom-Json
Assert ($process.ExitCode -eq 0 -and $manifest.exit_status -eq 'success' -and $manifest.accepted_steps -eq 4 -and $manifest.accepted_frames -eq 5 -and
 !(Test-Path "$path/output/meta/failure.json") -and !(Test-Path "$path/output/meta/manifest.json.tmp")) 'Atomic retry failed to preserve/replace complete manifest'
foreach($value in @('0','-1','1.5','1000001')){
 $body=Fixture 'staticCircleXY';$body['run.dict']=$body['run.dict'].Replace('enabled false;',"enabled true; writeEverySteps $value;")
 $path=Case "bad_frequency_$value" $body
 $ErrorActionPreference='Continue'; & $Program --case $path --check-case *> "$path/check.log"; $rc=$LASTEXITCODE;$ErrorActionPreference='Stop'
 Write-Text "$path/exit.txt" "$rc`n";Assert ($rc -ne 0 -and !(Test-Path "$path/output")) 'Invalid output frequency accepted'
}
@{status='PASS';frequency_cases=$frequencyCases;frames_verified=$framesChecked;vtk_files_verified=$vtkFiles;interface_samples_verified=$pointSamples;
 generic_ellipse_cases=3;effective_config_replay=$true;numeric_schedule_independence=$true;independent_input_sha256=$true;numeric_rejection=$true;existing_output_preserved=$true;accepted_io_failure=$true;
 atomic_replacement_retry=$true;invalid_frequencies=4}|
 ConvertTo-Json|Set-Content -Encoding UTF8 "$RunRoot/verification.json"
Get-Content -Raw -Encoding UTF8 "$RunRoot/verification.json"
