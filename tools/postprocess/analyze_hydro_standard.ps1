param(
    [string]$Run='',
    [ValidateSet('staticXY','staticRZ','oscillation','hysingTC1','hysingTC2')][string]$Benchmark='staticXY',
    [switch]$AcceptedPrefix,
    [switch]$CompileOnly
)
$ErrorActionPreference='Stop'
Set-StrictMode -Version Latest
$helper=Join-Path $PSScriptRoot 'HydroStandardAnalysis.cs'
Add-Type -Path $helper
if($CompileOnly){Write-Output 'Hydro standard analysis helper compiled; no data were analyzed.';return}
if(-not $Run){throw 'Provide the timestamped case run directory with -Run.'}
$Run=(Resolve-Path -LiteralPath $Run).Path
$project=Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$cases=[IO.Path]::GetFullPath((Join-Path $project 'runs'))+[IO.Path]::DirectorySeparatorChar
if(-not $Run.StartsWith($cases,[StringComparison]::OrdinalIgnoreCase)){throw 'Select a saved run under runs/.'}
$culture=[Globalization.CultureInfo]::InvariantCulture
function Num($v){$n=[double]::Parse([string]$v,$culture);if([double]::IsNaN($n) -or [double]::IsInfinity($n)){throw 'Nonfinite observation input.'};return $n}
function Near([double]$a,[double]$b){[Math]::Abs($a-$b) -le 1e-11*[Math]::Max(1,[Math]::Abs($b))}
function Optional([double]$v){if(-not [double]::IsNaN($v) -and -not [double]::IsInfinity($v)){$v}else{$null}}
$sources=New-Object 'System.Collections.Generic.List[string]'
function Source([string]$relative){$path=Join-Path $Run $relative;if(-not (Test-Path -LiteralPath $path -PathType Leaf)){throw "Missing accepted output: $relative"};$sources.Add($path);return $path}
$record=Get-Content -LiteralPath (Source 'run.json') -Raw -Encoding UTF8 | ConvertFrom-Json
if($record.check_only -or $record.status -notin @('completed','failed')){throw 'Require a terminal physics run, not a preflight, live or launcher-failed run.'}
if($record.status -ne 'completed' -and -not $AcceptedPrefix){throw 'Failed run: use -AcceptedPrefix to analyze only its accepted history.'}
$exe=Source 'runtime/QUANTUM.exe'
if((Get-FileHash -LiteralPath $exe).Hash -ne $record.program_sha256){throw 'Executable snapshot differs from run identity.'}
foreach($inputRecord in $record.inputs){$path=Source $inputRecord.path;if((Get-FileHash -LiteralPath $path).Hash -ne $inputRecord.sha256){throw "Modified input snapshot: $($inputRecord.path)"}}
$meta=Get-Content -LiteralPath (Source 'output/meta/hydro_moment_case.json') -Raw -Encoding UTF8 | ConvertFrom-Json
if($Benchmark -in @('hysingTC1','hysingTC2')){
    $tc=if($Benchmark -eq 'hysingTC1'){1}else{2}
    $rho=if($tc -eq 1){100.}else{1.};$mu=if($tc -eq 1){1.}else{0.1};$sigma=if($tc -eq 1){24.5}else{1.96}
    if($meta.format_version -ne 1 -or $meta.coordinate -ne 'XY' -or $meta.initial_velocity -ne 'zero' -or $meta.convection -ne 'central' -or
        -not (Near $meta.density_negative $rho) -or -not (Near $meta.density_positive 1000) -or
        -not (Near $meta.viscosity_negative $mu) -or -not (Near $meta.viscosity_positive 10) -or
        -not (Near $meta.surface_tension $sigma) -or $meta.acceleration[0] -ne 0 -or -not (Near $meta.acceleration[1] -0.98)){
        throw 'Resolved physics does not match the selected Hysing test case.'
    }
    $box=@(0.,0.,1.,2.);for($i=0;$i -lt 4;$i++){if(-not (Near $meta.domain[$i] $box[$i])){throw 'Hysing requires the original 1 by 2 domain.'}}
    if(@($meta.boundary_patches).Count -ne 2){throw 'Hysing requires two free-slip sidewalls; omitted top/bottom walls are no-slip.'}
    foreach($side in @('west','east')){
        $patch=@($meta.boundary_patches | Where-Object {$_.side -eq $side})
        if($patch.Count -ne 1 -or $patch[0].type -ne 'freeSlip' -or $patch[0].value -ne 0 -or
            -not (Near $patch[0].start 0) -or -not (Near $patch[0].end 2)){throw "Incorrect Hysing wall: $side"}
    }
    $timeBlock=[regex]::Match((Get-Content -LiteralPath (Source 'system/run.dict') -Raw -Encoding UTF8),'Time\s*\{([^{}]*)\}')
    $endMatch=[regex]::Match($timeBlock.Groups[1].Value,'\bendTime\s+([^;\s]+)\s*;')
    if(-not $timeBlock.Success -or -not $endMatch.Success -or $timeBlock.Groups[1].Value -notmatch '\bunits\s+code\s*;'){throw 'Require code time units and endTime in saved run.dict.'}
    $requestedEnd=Num $endMatch.Groups[1].Value
    $energy=@(Import-Csv -LiteralPath (Source 'output/series/hydro_moment_energy.csv'))
    $phase=@(Import-Csv -LiteralPath (Source 'output/series/hydro_moment_phase_integrals.csv') | Where-Object {$_.phase -eq 'negative'})
    if($energy.Count -lt 2 -or $energy.Count -ne $phase.Count){throw 'Incomplete or mismatched accepted Hysing histories.'}
    $time=New-Object 'double[]' $phase.Count
    $values=New-Object 'double[][]' 4;for($k=0;$k -lt 4;$k++){$values[$k]=New-Object 'double[]' $phase.Count}
    $history=New-Object 'System.Collections.Generic.List[object]';$lastTime=-1.;$initialFV=Num $phase[0].volume_fv
    $maxFV=0.;$maxGeometryFV=0.;$maxFraction=0.;$initialArea=[Math]::PI*0.25*0.25
    if($initialFV -le 0 -or [Math]::Abs((Num $phase[0].volume_geometry)/$initialArea-1) -gt 1e-4 -or
        [Math]::Abs((Num $phase[0].circularity)-1) -gt 1e-4 -or [Math]::Abs((Num $phase[0].centroid_x)-0.5) -gt 1e-8 -or
        [Math]::Abs((Num $phase[0].centroid_y)-0.5) -gt 1e-8 -or [Math]::Abs((Num $phase[0].mean_u_y)) -gt 1e-12){throw 'Initial bubble observations differ from the standard resting circle.'}
    for($i=0;$i -lt $phase.Count;$i++){
        $p=$phase[$i];$e=$energy[$i];$t=Num $p.time
        if([int]$p.step -ne $i -or [int]$e.step -ne $i -or -not (Near (Num $e.time) $t) -or $t -le $lastTime -or ($i -eq 0 -and $t -ne 0)){throw 'Hysing accepted clocks have gaps, duplicates or a mismatch.'}
        if($i -gt 0 -and (-not (Near ($t-$lastTime) (Num $e.dt)) -or [int]$e.pressure_available -ne 1 -or
            -not (Near (Num $e.pressure_time) (0.5*($lastTime+$t))))){throw 'Hysing endpoint and midpoint clocks disagree.'}
        if([int]$p.present -ne 1 -or [int]$p.surface_clipped -ne 0){throw 'The benchmark bubble is absent or clipped by the box.'}
        if($tc -eq 1 -and [int]$e.event_count -ne 0){throw 'TC1 does not admit topology events.'}
        $time[$i]=$t;$values[0][$i]=Num $p.volume_geometry;$values[1][$i]=Num $p.circularity;$values[2][$i]=Num $p.centroid_y;$values[3][$i]=Num $p.mean_u_y
        $fv=Num $p.volume_fv;$drift=($fv-$initialFV)/$initialFV;$geometryFV=($values[0][$i]-$fv)/$initialFV
        if(-not (Near (Num $p.geometry_minus_fv) ($values[0][$i]-$fv)) -or -not (Near (Num $p.volume_fv_drift) ($fv-$initialFV))){throw 'Reported phase budgets differ from the saved inventories.'}
        $maxFV=[Math]::Max($maxFV,[Math]::Abs($drift));$maxGeometryFV=[Math]::Max($maxGeometryFV,[Math]::Abs($geometryFV));$maxFraction=[Math]::Max($maxFraction,(Num $e.fraction_error))
        $history.Add([pscustomobject]@{step=$i;time=$t;area_geometry=$values[0][$i];area_fv=$fv;fv_relative_drift=$drift;geometry_minus_fv_relative=$geometryFV;
            circularity=$values[1][$i];centroid_y=$values[2][$i];mean_u_y=$values[3][$i];kinetic_energy=(Num $e.kinetic_energy);surface_energy=(Num $e.surface_energy);
            potential_energy=(Num $e.potential_energy);mechanical_energy=(Num $e.mechanical_energy);mechanical_energy_change=(Num $e.mechanical_energy_change);
            cumulative_viscous_dissipation=(Num $e.cumulative_viscous_dissipation);total_energy_balance=(Num $e.total_energy_balance);
            cumulative_event_dissipation=(Num $e.cumulative_event_dissipation);event_count=[int]$e.event_count;leaf_cells=[int]$e.leaf_cells;fraction_error=(Num $e.fraction_error)})
        $lastTime=$t
    }
    $referenceRoot=Join-Path $project 'cases/hydro_hysing_bubble/reference';$manifestPath=Join-Path $referenceRoot 'sources.json';$sources.Add($manifestPath)
    $manifest=Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
    $references=@($manifest.files | Where-Object {$_.file -like ('c'+$tc+'g*.txt')})
    if($references.Count -ne 3){throw 'Require the three original Hysing reference groups.'}
    $commonBegin=0.;$commonEnd=[Math]::Min(3.,$lastTime)
    foreach($reference in $references){
        $path=Join-Path $referenceRoot $reference.file
        if((Get-FileHash -LiteralPath $path).Hash -ne $reference.sha256){throw "Modified Hysing reference: $($reference.file)"};$sources.Add($path)
        $commonBegin=[Math]::Max($commonBegin,(Num $reference.first_row.value[0]));$commonEnd=[Math]::Min($commonEnd,(Num $reference.last_row.value[0]))
    }
    $comparisons=New-Object 'System.Collections.Generic.List[object]'
    foreach($reference in $references){$result=[HydroStandardAnalysis]::CompareHysing((Join-Path $referenceRoot $reference.file),$time,$values,$commonBegin,$commonEnd)
        $comparisons.Add([pscustomobject]@{file=$reference.file;source=$reference.source;result=$result})}
    $output=Join-Path $Run ('output/analysis/'+(Get-Date -Format 'yyyyMMdd_HHmmss_fff')+'_'+$Benchmark)
    New-Item -ItemType Directory -Path $output -ErrorAction Stop | Out-Null
    $history | Export-Csv -LiteralPath (Join-Path $output 'history.csv') -NoTypeInformation -Encoding UTF8
    $metrics=New-Object 'System.Collections.Generic.List[object]'
    foreach($comparison in $comparisons){
        $r=$comparison.result;$areaAvailable=$r.Metrics[0].Available
        @($r.Aligned | ForEach-Object {[pscustomobject]@{time=$_[0];area=$_[1];circularity=$_[2];centroid_y=$_[3];mean_u_y=$_[4];
            reference_area=$(if($areaAvailable){$_[5]}else{$null});reference_circularity=$_[6];reference_centroid_y=$_[7];reference_mean_u_y=$_[8]}}) |
            Export-Csv -LiteralPath (Join-Path $output ($comparison.file+'.csv')) -NoTypeInformation -Encoding UTF8
        foreach($m in $r.Metrics){$metrics.Add([pscustomobject]@{reference=$comparison.file;quantity=$m.Name;available=$m.Available;begin=$r.Begin;end=$r.End;samples=$r.Samples;
            absolute_time_RMS=$(if($m.Available){$m.AbsoluteRms}else{$null});relative_time_L2=$(if($m.Available){Optional $m.RelativeL2}else{$null});
            maximum_sampled_error=$(if($m.Available){$m.MaximumSampledError}else{$null});
            computed_minimum=$(if($m.Available){$m.ComputedMinimum}else{$null});computed_minimum_time=$(if($m.Available){$m.ComputedMinimumTime}else{$null});
            computed_maximum=$(if($m.Available){$m.ComputedMaximum}else{$null});computed_maximum_time=$(if($m.Available){$m.ComputedMaximumTime}else{$null});
            reference_minimum=$(if($m.Available){$m.ReferenceMinimum}else{$null});reference_minimum_time=$(if($m.Available){$m.ReferenceMinimumTime}else{$null});
            reference_maximum=$(if($m.Available){$m.ReferenceMaximum}else{$null});reference_maximum_time=$(if($m.Available){$m.ReferenceMaximumTime}else{$null});
            computed_last=$(if($m.Available){$m.ComputedLast}else{$null});reference_last=$(if($m.Available){$m.ReferenceLast}else{$null})})}
    }
    $metrics | Export-Csv -LiteralPath (Join-Path $output 'reference_errors.csv') -NoTypeInformation -Encoding UTF8
    $sources.Add($helper);$sources.Add($PSCommandPath)
    $summary=[ordered]@{benchmark=$Benchmark;generated=(Get-Date).ToString('o');run=$Run;program_sha256=$record.program_sha256;run_status=$record.status;
        requested_end=$requestedEnd;standard_end=3.;accepted_end=$lastTime;accepted_steps=$phase.Count-1;full_interval=($record.status -eq 'completed' -and (Near $requestedEnd 3) -and (Near $lastTime 3));
        numerical_acceptance='not_decided_by_single_run';physical_definition_checked=$true;common_comparison_window=@($comparisons[0].result.Begin,$comparisons[0].result.End);
        max_relative_fv_drift=$maxFV;max_relative_geometry_fv_difference=$maxGeometryFV;max_fraction_error=$maxFraction;final_observation=$history[$history.Count-1];
        event_count=[int]$energy[-1].event_count;reference_errors=@($metrics);
        limitations=@('References remain separate; no average is treated as an exact solution.',
            'Comparison uses only common accepted endpoints, linear reference interpolation and trapezoidal time norms; extrema are sampled.',
            'MooNMD zero area is unavailable, not a vanishing bubble.',
            'TC2 late-time topology depends on the resolved interface and declared event model.',
            'Mechanical energy change alone is not a conservation residual; dissipation, event and discrete work histories remain separate.',
            'Mesh/time sensitivity, shapes and predefined numerical acceptance remain required.');
        sources=@($sources | Sort-Object -Unique | ForEach-Object {[ordered]@{path=$_;sha256=(Get-FileHash -LiteralPath $_).Hash}})}
    [IO.File]::WriteAllText((Join-Path $output 'summary.json'),($summary | ConvertTo-Json -Depth 9),[Text.UTF8Encoding]::new($false))
    Write-Output ('Analysis: '+$output);Write-Output ('Full predefined interval: '+$summary.full_interval+'; references remain separate.');return
}
$rz=$Benchmark -eq 'staticRZ';$osc=$Benchmark -eq 'oscillation'
$radius=if($osc){0.1}else{0.4}
$expectedDomain=if($rz){@(0.,-1.,1.,1.)}elseif($osc){@(-0.5,-0.5,0.5,0.5)}else{@(-1.,-1.,1.,1.)}
$expectedCoordinate=if($rz){'RZ'}else{'XY'}
$expectedMu=if($osc){0.}else{[Math]::Sqrt(0.8/120)}
$expectedRhoPositive=if($osc){0.001}else{1.}
if($meta.format_version -ne 1 -or $meta.coordinate -ne $expectedCoordinate -or $meta.arc_basis -ne 'power_in_fraction_0_1' -or
   $meta.initial_velocity -ne 'zero' -or $meta.convection -ne 'central' -or
   -not (Near $meta.density_negative 1) -or -not (Near $meta.density_positive $expectedRhoPositive) -or
   -not (Near $meta.viscosity_negative $expectedMu) -or -not (Near $meta.viscosity_positive $expectedMu) -or
   -not (Near $meta.surface_tension 1) -or $meta.acceleration[0] -ne 0 -or $meta.acceleration[1] -ne 0){throw 'Resolved physics does not match the selected fixed B2/B3 definition.'}
for($i=0;$i -lt 4;$i++){if(-not (Near $meta.domain[$i] $expectedDomain[$i])){throw 'Reference domain differs from the actual domain.'}}
$sides=if($rz){@('south','east','north')}else{@('west','south','east','north')}
if(@($meta.boundary_patches).Count -ne $sides.Count){throw 'Reference requires full free-slip outer walls and the native RZ axis.'}
foreach($side in $sides){
    $patch=@($meta.boundary_patches | Where-Object {$_.side -eq $side})
    $horizontal=$side -in @('south','north');$lo=if($horizontal){$expectedDomain[0]}else{$expectedDomain[1]};$hi=if($horizontal){$expectedDomain[2]}else{$expectedDomain[3]}
    if($patch.Count -ne 1 -or $patch[0].type -ne 'freeSlip' -or $patch[0].value -ne 0 -or -not (Near $patch[0].start $lo) -or -not (Near $patch[0].end $hi)){throw "Incorrect standard wall: $side"}
}
$timeBlock=[regex]::Match((Get-Content -LiteralPath (Join-Path $Run 'system/run.dict') -Raw -Encoding UTF8),'Time\s*\{([^{}]*)\}')
$endMatch=[regex]::Match($timeBlock.Groups[1].Value,'\bendTime\s+([^;\s]+)\s*;')
if(-not $timeBlock.Success -or -not $endMatch.Success -or $timeBlock.Groups[1].Value -notmatch '\bunits\s+code\s*;'){throw 'Require code time units and endTime in saved run.dict.'}
$requestedEnd=Num $endMatch.Groups[1].Value
$standardEnd=if($osc){1.}else{0.64/$expectedMu}
$energy=@(Import-Csv -LiteralPath (Source 'output/series/hydro_moment_energy.csv'))
$phase=@(Import-Csv -LiteralPath (Source 'output/series/hydro_moment_phase_integrals.csv') | Where-Object {$_.phase -eq 'negative'})
if($energy.Count -lt 2 -or $energy.Count -ne $phase.Count){throw 'Incomplete or mismatched accepted energy/phase histories.'}
$history=New-Object 'System.Collections.Generic.List[object]'
$initialVolume=Num $energy[0].liquid_volume;$initialEnergy=(Num $energy[0].kinetic_energy)+(Num $energy[0].surface_energy)
$equilibriumSurface=2*[Math]::Sqrt([Math]::PI*(Num $phase[0].volume_geometry))
$excess=if($osc){$initialEnergy-$equilibriumSurface}else{$initialEnergy}
if($initialVolume -le 0 -or $excess -le 0){throw 'Invalid initial liquid volume or oscillation excess energy.'}
$lastTime=-1.;$maxRaw=0.;$maxInventory=0.;$maxFraction=0.
for($i=0;$i -lt $energy.Count;$i++){
    $e=$energy[$i];$p=$phase[$i];$t=Num $e.time
    if([int]$e.step -ne $i -or [int]$p.step -ne $i -or -not (Near (Num $p.time) $t) -or $t -le $lastTime -or ($i -eq 0 -and $t -ne 0)){throw 'Accepted clocks have a gap, duplicate or mismatch.'}
    if($i -gt 0 -and (-not (Near ($t-$lastTime) (Num $e.dt)) -or [int]$e.pressure_available -ne 1 -or -not (Near (Num $e.pressure_time) (0.5*($lastTime+$t))))){throw 'Endpoint and midpoint clocks disagree.'}
    if([int]$e.event_count -ne 0 -or [int]$p.surface_clipped -ne 0){throw 'B2/B3 does not admit topology events or a clipped free surface.'}
    $raw=(Num $e.kinetic_energy)+(Num $e.surface_energy)-$initialEnergy
    $drift=[Math]::Abs((Num $e.liquid_volume)-$initialVolume)/$initialVolume
    $maxRaw=[Math]::Max($maxRaw,[Math]::Abs($raw)/$excess);$maxInventory=[Math]::Max($maxInventory,$drift);$maxFraction=[Math]::Max($maxFraction,(Num $e.fraction_error))
    $history.Add([pscustomobject]@{step=$i;time=$t;pressure_time=$(if($i){Num $e.pressure_time}else{$null});kinetic_energy=(Num $e.kinetic_energy);surface_energy=(Num $e.surface_energy);
        reference_quadrant_kinetic_energy=$(if($osc){0.5*(Num $e.kinetic_energy)}else{$null});raw_energy_change=$raw;raw_energy_change_over_initial_energy=$raw/$initialEnergy;
        raw_energy_change_over_initial_excess=$(if($osc){$raw/$excess}else{$null});discrete_balance=(Num $e.total_energy_balance);liquid_inventory=(Num $e.liquid_volume);
        relative_inventory_drift=$drift;fraction_error=(Num $e.fraction_error);mean_curvature=(Num $p.mean_curvature);surface_measure=(Num $p.interface_measure);
        curvature_squared_integral=(Num $p.curvature_squared_integral);surface_quadrature_estimate=(Num $p.surface_quadrature_relative_estimate)})
    $lastTime=$t
}
$complete=$record.status -eq 'completed' -and $record.exit_code -eq 0 -and (Near $requestedEnd $standardEnd) -and (Near $lastTime $standardEnd)
if(-not $complete -and -not $AcceptedPrefix){throw 'Full predefined standard interval was not completed; use -AcceptedPrefix to keep that limitation explicit.'}
$polynomials=@(Get-ChildItem -LiteralPath (Join-Path $Run 'output/fields') -Filter '*_arcs_step*.csv' -File | Sort-Object {[long]([regex]::Match($_.BaseName,'step(\d+)$').Groups[1].Value)})
if(-not $polynomials.Count){throw 'Missing raw Hermite frames; sampled contours cannot replace them.'}
$frames=New-Object 'System.Collections.Generic.List[object]';$geometry=New-Object 'System.Collections.Generic.List[object]';$times=New-Object 'System.Collections.Generic.HashSet[int]'
$initialShapeError=$null
foreach($file in $polynomials){
    $sources.Add($file.FullName);$g=[HydroStandardAnalysis]::ReadGeometry($file.FullName,$rz,$radius,$osc)
    if($g.Step -ge $history.Count -or -not (Near $g.Time $history[$g.Step].time) -or -not $times.Add($g.Step)){throw 'A polynomial frame is not a unique accepted time.'}
    $fieldPath=$file.FullName -replace '_arcs_step','_step' -replace '\.csv$','.dat'
    if(-not (Test-Path -LiteralPath $fieldPath)){throw 'An interface frame has no matching fluid frame.'}
    $sources.Add($fieldPath);$f=[HydroStandardAnalysis]::ReadField($fieldPath)
    if(-not (Near $f.Time $g.Time) -or $f.PressureAvailable -ne ($g.Step -gt 0) -or ($f.PressureAvailable -and -not (Near $f.PressureTime $history[$g.Step].pressure_time))){throw 'Field time layers differ from the accepted history.'}
    $s=[HydroStandardAnalysis]::Static($f,$rz,$radius);$capillarySpeed=[Math]::Sqrt(1/$radius);$jump=if($rz){2/$radius}else{1/$radius}
    $frames.Add([pscustomobject]@{step=$g.Step;time=$g.Time;pressure_time=$(if($s.PressureAvailable){$s.PressureTime}else{$null});pressure_available=$s.PressureAvailable;
        mean_pressure_inside=$(if(-not $osc){Optional $s.MeanInside}else{$null});mean_pressure_outside=$(if(-not $osc){Optional $s.MeanOutside}else{$null});
        pressure_jump=$(if(-not $osc){Optional $s.PressureJump}else{$null});pressure_jump_relative_error=$(if(-not $osc){Optional ([Math]::Abs($s.PressureJump/$jump-1))}else{$null});
        inner_cells=$s.InsideCells;outer_cells=$s.OutsideCells;inner_volume=$s.InnerVolume;outer_volume=$s.OuterVolume;maximum_cell_width=$s.MaximumCellWidth;
        speed_center_max=$s.SpeedCenter;speed_sampled_5_max=$s.Speed5;speed_sampled_9_max=$s.Speed9;speed_sampled_9_over_capillary=$s.Speed9/$capillarySpeed;
        speed_sampling_increment_over_capillary=($s.Speed9-$s.Speed5)/$capillarySpeed})
    $geometry.Add([pscustomobject]@{step=$g.Step;time=$g.Time;arcs=$g.Arcs;surface_measure=$g.Measure;geometric_volume=$g.Volume;
        curvature_relative_L2_to_circle_or_sphere=$g.CurvatureRelativeL2;radial_distance_relative_L2_to_circle_or_sphere=$g.RadialDistanceRelativeL2;
        quadrature_absolute_estimate=$g.QuadratureEstimate;polar_available=$g.PolarAvailable;polar_reason=$g.PolarReason;mean_polar_radius=(Optional $g.MeanRadius);
        mode2_cos=(Optional $g.Mode2Cos);mode2_sin=(Optional $g.Mode2Sin);mirror_relative_L2=(Optional $g.MirrorRelativeL2);polar_512_1024_difference=(Optional $g.PolarSamplingDifference)})
    if($g.Step -eq 0){
        if($osc){if(-not $g.PolarAvailable){throw 'Initial oscillation curve has no regular polar representation.'};$square=0.
            for($j=0;$j -lt $g.PolarRadii.Length;$j++){$d=$g.PolarRadii[$j]/$radius-1-0.05*[Math]::Cos(4*[Math]::PI*$j/$g.PolarRadii.Length);$square+=$d*$d}
            $initialShapeError=[Math]::Sqrt($square/$g.PolarRadii.Length)
        }else{$initialShapeError=$g.RadialDistanceRelativeL2}
    }
}
if($null -eq $initialShapeError -or $initialShapeError -gt 1e-3){throw 'Initial geometry is missing or inconsistent with the declared standard shape.'}
$output=Join-Path $Run ('output/analysis/'+(Get-Date -Format 'yyyyMMdd_HHmmss_fff')+'_'+$Benchmark)
New-Item -ItemType Directory -Path $output -ErrorAction Stop | Out-Null
$history | Export-Csv -LiteralPath (Join-Path $output 'history.csv') -NoTypeInformation -Encoding UTF8
$frames | Export-Csv -LiteralPath (Join-Path $output 'frames.csv') -NoTypeInformation -Encoding UTF8
$geometry | Export-Csv -LiteralPath (Join-Path $output 'geometry.csv') -NoTypeInformation -Encoding UTF8
$fit=$null;$omega=$null;$crossingPeriods=@()
if($osc){
    $omega=[Math]::Sqrt(6/((1+0.001)*[Math]::Pow($radius,3)))
    $f=[HydroStandardAnalysis]::FitEnergy([double[]]@($history | ForEach-Object {$_.time}),[double[]]@($history | ForEach-Object {$_.kinetic_energy}),$omega)
    $fit=[ordered]@{available=$f.Available;reason=$f.Reason;window=@($f.Begin,$f.End);search_frequency=@($f.FrequencyLower,$f.FrequencyUpper);search_damping=@($f.DampingLower,$f.DampingUpper);boundary=$f.AtSearchBoundary}
    if($f.Available){
        $fit['omega']=$f.Omega;$fit['frequency_relative_error']=$f.Omega/$omega-1;$fit['damping']=$f.Damping;$fit['amplitude']=$f.Amplitude;$fit['relative_RMS_residual']=$f.RelativeRms
        $fit['shape_amplitude_loss_per_linear_period']=1-[Math]::Exp(-0.5*$f.Damping*2*[Math]::PI/$omega)
        @($history | ForEach-Object {[pscustomobject]@{time=$_.time;kinetic_energy=$_.kinetic_energy;fit=$f.Amplitude*[Math]::Exp(-$f.Damping*($_.time-$f.Begin))*(1-[Math]::Cos(2*$f.Omega*$_.time))}}) |
            Export-Csv -LiteralPath (Join-Path $output 'kinetic_fit.csv') -NoTypeInformation -Encoding UTF8
    }
    # Independent signed-shape timing: pair crossings of the same direction.
    $crossings=@{up=New-Object 'System.Collections.Generic.List[double]';down=New-Object 'System.Collections.Generic.List[double]'}
    $crossingRows=New-Object 'System.Collections.Generic.List[object]'
    for($j=1;$j -lt $geometry.Count;$j++){
        $a=$geometry[$j-1];$b=$geometry[$j];if(-not $a.polar_available -or -not $b.polar_available){continue}
        if(($a.mode2_cos -lt 0 -and $b.mode2_cos -ge 0) -or ($a.mode2_cos -gt 0 -and $b.mode2_cos -le 0)){
            $time=$a.time+($b.time-$a.time)*$a.mode2_cos/($a.mode2_cos-$b.mode2_cos);$direction=if($b.mode2_cos -gt $a.mode2_cos){'up'}else{'down'}
            $ordinal=$crossings[$direction].Count
            $linearTime=($(if($direction -eq 'down'){0.5}else{1.5})+2*$ordinal)*[Math]::PI/$omega
            $crossingRows.Add([pscustomobject]@{direction=$direction;ordinal=$ordinal;time=$time;linear_reference_time=$linearTime;
                timing_difference=$time-$linearTime;phase_lag_radians=$omega*($time-$linearTime);interpolation_interval=$b.time-$a.time})
            $crossings[$direction].Add($time)
        }
    }
    foreach($direction in @('up','down')){for($j=1;$j -lt $crossings[$direction].Count;$j++){$crossingPeriods+= $crossings[$direction][$j]-$crossings[$direction][$j-1]}}
    if($crossingRows.Count){$crossingRows | Export-Csv -LiteralPath (Join-Path $output 'shape_crossings.csv') -NoTypeInformation -Encoding UTF8}
}
$sources.Add($helper);$sources.Add($PSCommandPath)
$referenceVolume=if($rz){4*[Math]::PI*[Math]::Pow($radius,3)/3}elseif($osc){[Math]::PI*$radius*$radius*(1+0.05*0.05/2)}else{[Math]::PI*$radius*$radius}
$summary=[ordered]@{
    benchmark=$Benchmark;generated=(Get-Date).ToString('o');run=$Run;program_sha256=$record.program_sha256;run_status=$record.status;requested_end=$requestedEnd;standard_end=$standardEnd;
    full_interval=$complete;accepted_steps=$energy.Count-1;accepted_end=$lastTime;last_field_step=$geometry[-1].step;last_field_is_last_accepted=($geometry[-1].step -eq $energy.Count-1);
    numerical_acceptance='not_decided_by_single_run';physical_definition_checked=$true;initial_shape_relative_L2=$initialShapeError;initial_geometry_volume_relative_error=$geometry[0].geometric_volume/$referenceVolume-1;
    max_relative_inventory_drift=$maxInventory;max_fraction_error=$maxFraction;max_absolute_raw_energy_change_scaled=$maxRaw;raw_energy_scale=$(if($osc){'initial excess over equal-area circle'}else{'initial total energy'});
    initial_energy=$initialEnergy;initial_excess_energy=$(if($osc){$excess}else{$null});reference_shape_omega=$omega;kinetic_fit=$fit;shape_same_direction_crossing_periods=$crossingPeriods;
    final_saved_frame=$frames[-1];final_saved_geometry=$geometry[-1];max_saved_frame_speed_over_capillary=($frames | Measure-Object -Property speed_sampled_9_over_capillary -Maximum).Maximum;
    max_saved_speed_sampling_increment=($frames | Measure-Object -Property speed_sampling_increment_over_capillary -Maximum).Maximum;
    all_saved_polar_modes_available=$(if($osc){@($geometry | Where-Object {-not $_.polar_available}).Count -eq 0}else{$null});
    max_saved_mirror_relative_L2=$(if($osc){($geometry | Measure-Object -Property mirror_relative_L2 -Maximum).Maximum}else{$null});
    max_saved_polar_sampling_difference=$(if($osc){($geometry | Measure-Object -Property polar_512_1024_difference -Maximum).Maximum}else{$null});
    targets=$(if($osc){[ordered]@{finest_frequency_relative_error=0.05;kinetic_fit_relative_RMS=0.05;amplitude_loss_per_period=0.05;max_raw_energy_over_initial_excess=0.1;mirror_relative_L2=0.001}}else{[ordered]@{pressure_jump_relative_error=0.01;sampled_speed_over_capillary=0.001;speed_sampling_increment_over_capillary=0.0001;finest_curvature_relative_L2=0.01;minimum_curvature_order=2}});
    limitations=@('Targets are project criteria, not literature-wide tolerances; mesh/time comparisons remain required.',
        'Speed maxima cover saved frames and nested leaf samples, not every physical time or a continuous supremum.',
        'Pressure is a midpoint P1 physical-volume mean; fixed reference sampling masks use endpoint pure-cell checks.',
        'Polynomial quadrature/root and angular sampling differences are numerical estimates, not rigorous enclosures.',
        'B3 circle-curvature error includes physical deformation and is not an oscillation accuracy criterion.',
        'Energy fit uses the declared bounded search and fixed zero phase; inspect residual and independent shape crossings.',
        'Linear infinite-domain frequency has finite-amplitude and finite-wall modelling differences.');
    sources=@($sources | Sort-Object -Unique | ForEach-Object {[ordered]@{path=$_;sha256=(Get-FileHash -LiteralPath $_).Hash}})
}
[IO.File]::WriteAllText((Join-Path $output 'summary.json'),($summary | ConvertTo-Json -Depth 9),[Text.UTF8Encoding]::new($false))
Write-Output ('Analysis: '+$output)
Write-Output ('Full predefined interval: '+$complete+'; accepted steps: '+($energy.Count-1)+'. A single-run analysis does not close module acceptance.')
