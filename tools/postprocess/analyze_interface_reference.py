"""Read-only verification of saved FirstJet geometry against analytic references.

Input is a terminal run with a frozen reference/verification.json. This tool
never imports the solver or changes its state. VTP metrics describe exported
samples/chords; they are not full authoritative-geometry convergence norms.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def near(a, b):
    return abs(float(a) - float(b)) <= 2e-10 * max(1.0, abs(float(b)))


class Jet:
    """Second spatial jet; elementary chain rules independent of Core AST."""
    def __init__(self, v, x=0, y=0, xx=0, xy=0, yy=0):
        self.v, self.x, self.y, self.xx, self.xy, self.yy = v, x, y, xx, xy, yy

    @staticmethod
    def cast(value):
        return value if isinstance(value, Jet) else Jet(value)

    def __add__(self, other):
        b = Jet.cast(other)
        return Jet(*(getattr(self, k) + getattr(b, k) for k in ('v', 'x', 'y', 'xx', 'xy', 'yy')))

    __radd__ = __add__

    def __neg__(self):
        return Jet(*(-getattr(self, k) for k in ('v', 'x', 'y', 'xx', 'xy', 'yy')))

    def __sub__(self, other):
        return self + -Jet.cast(other)

    def __rsub__(self, other):
        return Jet.cast(other) + -self

    def __mul__(self, other):
        b = Jet.cast(other)
        a = self
        return Jet(a.v*b.v, a.x*b.v+a.v*b.x, a.y*b.v+a.v*b.y,
                   a.xx*b.v+2*a.x*b.x+a.v*b.xx,
                   a.xy*b.v+a.x*b.y+a.y*b.x+a.v*b.xy,
                   a.yy*b.v+2*a.y*b.y+a.v*b.yy)

    __rmul__ = __mul__

    def unary(self, value, first, second):
        return Jet(value, first*self.x, first*self.y,
                   second*self.x*self.x+first*self.xx,
                   second*self.x*self.y+first*self.xy,
                   second*self.y*self.y+first*self.yy)

    def sin(self):
        return self.unary(np.sin(self.v), np.cos(self.v), -np.sin(self.v))

    def cos(self):
        return self.unary(np.cos(self.v), -np.sin(self.v), -np.cos(self.v))


class RadialTwistCircle:
    def __init__(self, reference):
        self.center = np.asarray(reference['flow_center'], dtype=float)
        self.circle = np.asarray(reference['circle_center'], dtype=float)
        self.radius = float(reference['radius'])
        self.speed = float(reference['speed'])
        self.period = float(reference['period'])
        require(self.radius > 0 and self.period > 0, 'Invalid analytic reference scales')

    def beta(self, t):
        return self.speed*self.period/math.pi*math.sin(math.pi*t/self.period)

    def inverse(self, points, t):
        a = Jet(points[:, 0]-self.center[0], 1)
        b = Jet(points[:, 1]-self.center[1], 0, 1)
        angle = -self.beta(t)*(a*a+b*b)
        c, s = angle.cos(), angle.sin()
        return self.center[0]+c*a-s*b, self.center[1]+s*a+c*b

    def field(self, points, t):
        x, y = self.inverse(points, t)
        a, b = x-self.circle[0], y-self.circle[1]
        return a*a+b*b-self.radius*self.radius

    def geometry(self, points, t):
        f = self.field(points, t)
        g = np.hypot(f.x, f.y)
        require(np.all(g > 1e-14), 'Reference gradient vanishes at interface sample')
        normal = np.column_stack((f.x/g, f.y/g))
        curvature = (f.y*f.y*f.xx-2*f.x*f.y*f.xy+f.x*f.x*f.yy)/(g*g*g)
        return f, g, normal, curvature

    def curve(self, theta, t):
        """Forward circle map and analytic theta derivatives (second reference)."""
        c, s = np.cos(theta), np.sin(theta)
        a = self.circle[0]-self.center[0]+self.radius*c
        b = self.circle[1]-self.center[1]+self.radius*s
        ap, bp, app, bpp = -self.radius*s, self.radius*c, -self.radius*c, -self.radius*s
        beta = self.beta(t)
        angle = beta*(a*a+b*b)
        anglep = 2*beta*(a*ap+b*bp)
        anglepp = 2*beta*(ap*ap+a*app+bp*bp+b*bpp)
        ca, sa = np.cos(angle), np.sin(angle)
        rotate = lambda x, y: np.column_stack((ca*x-sa*y, sa*x+ca*y))
        points = rotate(a, b)+self.center
        first = rotate(ap-anglep*b, bp+anglep*a)
        second = rotate(app-anglepp*b-2*anglep*bp-anglep*anglep*a,
                        bpp+anglepp*a+2*anglep*ap-anglep*anglep*b)
        length = np.hypot(first[:, 0], first[:, 1])
        curvature = (first[:, 0]*second[:, 1]-first[:, 1]*second[:, 0])/length**3
        normal = np.column_stack((first[:, 1]/length, -first[:, 0]/length))
        return points, normal, curvature


def reference_self_test(ref):
    errors = {'phi': 0.0, 'normal': 0.0, 'curvature': 0.0, 'jacobian': 0.0}
    theta = np.linspace(0, 2*math.pi, 128, endpoint=False)
    for fraction in (0, .125, .25, .5, 1, 1.5, 2):
        t = fraction*ref.period
        p, normal, curvature = ref.curve(theta, t)
        f, _, n, k = ref.geometry(p, t)
        x, y = ref.inverse(p, t)
        errors['phi'] = max(errors['phi'], float(np.max(np.abs(f.v))))
        errors['normal'] = max(errors['normal'], float(np.max(np.linalg.norm(n-normal, axis=1))))
        errors['curvature'] = max(errors['curvature'], float(np.max(np.abs(k-curvature))))
        errors['jacobian'] = max(errors['jacobian'], float(np.max(np.abs(x.x*y.y-x.y*y.x-1))))
    require(errors['phi'] < 1e-12 and errors['normal'] < 1e-11 and
            errors['curvature'] < 1e-9 and errors['jacobian'] < 1e-11,
            'Inverse/parametric reference cross-check failed: '+str(errors))
    # Finite differences also check the jet implementation away from the zero set.
    p = np.array([[.58, .51], [.73, .62], [.38, .44]])
    h, t = 2e-5, .31*ref.period
    f = ref.field(p, t)
    finite = []
    for axis, name in ((0, 'x'), (1, 'y')):
        shift = np.zeros_like(p)
        shift[:, axis] = h
        difference = (ref.field(p+shift, t).v-ref.field(p-shift, t).v)/(2*h)
        finite.append(float(np.max(np.abs(difference-getattr(f, name)))))
    require(max(finite) < 2e-7, 'Reference derivative finite-difference check failed')
    errors['finite_difference_gradient'] = max(finite)
    return errors


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def xml_array(parent, name, dtype=float):
    item = parent.find("DataArray[@Name='"+name+"']")
    require(item is not None and item.get('format') == 'ascii', 'Missing ASCII VTK array: '+name)
    data = np.fromstring(item.text or '', sep=' ', dtype=dtype)
    require(np.all(np.isfinite(data)), 'Nonfinite VTK array: '+name)
    return data


def vtk(path, frame):
    root = ET.parse(path).getroot()
    fields = root.find('.//FieldData')
    require(fields is not None, 'Missing VTK identity: '+str(path))
    for key, expected in (('accepted_step', frame['step']), ('geometry_identity', frame['geometry_identity']),
                          ('mesh_revision', frame['mesh_revision']), ('TimeValue', frame['time'])):
        value = xml_array(fields, key)
        require(len(value) == 1 and near(value[0], expected), 'VTK/frame identity mismatch: '+key)
    piece = root.find('.//Piece')
    require(piece is not None, 'Missing VTK Piece')
    points = xml_array(piece.find('Points'), 'Points').reshape(-1, 3)
    require(len(points) == int(piece.get('NumberOfPoints')) and np.all(points[:, 2] == 0), 'Bad XY points')
    return piece, points[:, :2]


def surface(path, frame):
    piece, points = vtk(path, frame)
    ids = xml_array(piece.find('Lines'), 'connectivity', int)
    offsets = xml_array(piece.find('Lines'), 'offsets', int)
    require(len(offsets) == int(piece.get('NumberOfLines')) and len(offsets) > 0 and
            offsets[-1] == len(ids) and np.all(ids >= 0) and np.all(ids < len(points)), 'Bad VTP lines')
    weights = np.zeros(len(points))
    segments = []
    start = 0
    for end in offsets:
        line = ids[start:end]
        require(len(line) >= 2, 'Degenerate VTP line')
        length = np.linalg.norm(points[line[1:]]-points[line[:-1]], axis=1)
        np.add.at(weights, line[:-1], .5*length)
        np.add.at(weights, line[1:], .5*length)
        segments.extend(zip(line[:-1], line[1:]))
        start = end
    data = piece.find('PointData')
    normal = xml_array(data, 'normal').reshape(-1, 3)[:, :2]
    curvature = xml_array(data, 'curvature')
    require(len(normal) == len(points) == len(curvature) and weights.sum() > 0, 'Bad interface samples')
    require(np.max(np.abs(np.linalg.norm(normal, axis=1)-1)) < 1e-10, 'Nonunit numerical normal')
    segments = np.asarray(segments)
    return points, normal, curvature, weights, points[segments[:, 0]], points[segments[:, 1]]


def distance_to_segments(query, a, b):
    delta = b-a
    norm = np.sum(delta*delta, axis=1)
    require(np.all(norm > 0), 'Zero-length contour segment')
    distances = []
    for start in range(0, len(query), 128):
        difference = query[start:start+128, None, :]-a[None, :, :]
        projection = np.clip(np.sum(difference*delta[None, :, :], axis=2)/norm[None, :], 0, 1)
        residual = difference-projection[:, :, None]*delta[None, :, :]
        distances.append(np.sqrt(np.min(np.sum(residual*residual, axis=2), axis=1)))
    return np.concatenate(distances)


def volume_error(path, frame, ref):
    piece, points = vtk(path, frame)
    cells = xml_array(piece.find('Cells'), 'connectivity', int).reshape(-1, 4)
    require(len(cells) == int(piece.get('NumberOfCells')) and np.all(cells >= 0) and
            np.all(cells < len(points)), 'Bad native quadrilateral connectivity')
    centers = np.mean(points[cells], axis=1)
    data = piece.find('CellData')
    weights = xml_array(data, 'cell_capacity')
    phi = xml_array(data, 'phi_center')
    require(len(weights) == len(phi) == len(cells) and np.all(weights > 0) and near(weights.sum(), 1),
            'Bad native cell measures or domain')
    f = ref.field(centers, float(frame['time']))
    errors = phi-f.v
    h = np.max(np.ptp(points[cells], axis=1), axis=1)
    band = np.abs(f.v) <= 2*h*np.hypot(f.x, f.y)
    require(np.any(band), 'Empty reference narrow band')
    return math.sqrt(float(np.sum(weights*errors**2)/weights.sum())), math.sqrt(
        float(np.sum(weights[band]*errors[band]**2)/weights[band].sum()))


def plots(directory, history, conservation, performance, snapshots, ref):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    theta = np.linspace(0, 2*math.pi, 4097)
    for axis, (frame, segments) in zip(axes.flat, snapshots):
        exact = ref.curve(theta, float(frame['time']))[0]
        axis.plot(exact[:, 0], exact[:, 1], color='#db5b33', lw=1.8, label='Exact')
        for j, (a, b) in enumerate(zip(*segments)):
            axis.plot([a[0], b[0]], [a[1], b[1]], color='#155b88', lw=.9,
                      label='Saved Hermite geometry' if j == 0 else None)
        axis.set(xlim=(.18, .82), ylim=(.18, .82), aspect='equal', xlabel='x', ylabel='y',
                 title=f"t={float(frame['time']):.2f}; step {frame['step']}")
        axis.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7, loc='upper left')
    fig.suptitle('Reversible differential rotation: analytic and exported interface')
    for extension in ('png', 'svg'):
        fig.savefig(directory/('interface_evolution.'+extension), dpi=220)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    time = [r['time'] for r in history]
    axes[0, 0].plot([float(r['time']) for r in conservation],
                    [abs(float(r['phase_volume'])/(math.pi*ref.radius**2)-1) for r in conservation])
    axes[0, 0].set(title='Relative area error (all accepted states)', ylabel='|A/A_exact - 1|')
    axes[0, 1].plot(time, [r['position_defect_max_over_radius'] for r in history], label='Normal position defect / R')
    axes[0, 1].plot(time, [r['sampled_chord_hausdorff_over_radius'] for r in history], label='Sampled chord Hausdorff / R')
    axes[0, 1].set(title='Exported interface position'); axes[0, 1].legend(fontsize=8)
    axes[1, 0].plot(time, [r['normal_l2'] for r in history], label='Normal L2')
    axes[1, 0].plot(time, [r['curvature_relative_l2'] for r in history], label='Relative curvature L2')
    axes[1, 0].set(title='Arc-length-weighted sample errors'); axes[1, 0].legend(fontsize=8)
    total = [r for r in performance if r['stage'] == 'step_total']
    axes[1, 1].plot([float(r['time']) for r in total], [float(r['wall_seconds']) for r in total], lw=.7)
    axes[1, 1].set(title='Accepted-step wall time', ylabel='seconds')
    for axis in axes.flat:
        axis.set_xlabel('time'); axis.grid(alpha=.25)
    for extension in ('png', 'svg'):
        fig.savefig(directory/('error_history.'+extension), dpi=220)
    plt.close(fig)


def analyze(run, accepted_prefix=False):
    project = Path(__file__).resolve().parents[2]
    run = run.resolve()
    require(run.is_relative_to((project/'runs').resolve()), 'Select a timestamped run under runs/')
    sources = {}

    def source(relative):
        path = (run/relative).resolve()
        require(path.is_relative_to(run) and path.is_file(), 'Missing or unsafe source: '+relative)
        sources[relative] = sha(path)
        return path

    record = json.loads(source('run.json').read_text(encoding='utf-8-sig'))
    require(not record.get('check_only') and record['status'] in ('completed', 'failed'), 'Require a terminal physics run')
    require(record['status'] == 'completed' or accepted_prefix, 'Failed run; explicitly request --accepted-prefix')
    executable = 'runtime/'+Path(record['program']).name
    require(sha(source(executable)).lower() == record['program_sha256'].lower(), 'Modified executable snapshot')
    for item in record['inputs']:
        require(sha(source(item['path'])).lower() == item['sha256'].lower(), 'Modified input snapshot: '+item['path'])
    reference = json.loads(source('reference/verification.json').read_text(encoding='utf-8-sig'))
    require(reference['schema_version'] == 1 and reference['reference'] == 'radialTwistCircle' and
            reference['coordinates'] == 'xy', 'Unsupported reference model/schema')
    require(len(reference['input_sha256']) == 5, 'Require five frozen reference input hashes')
    for path, expected in reference['input_sha256'].items():
        require(sha(source(path)).lower() == expected.lower(), 'Input differs from frozen study: '+path)
    ref = RadialTwistCircle(reference)
    self_test = reference_self_test(ref)
    manifest = json.loads(source('output/meta/manifest.json').read_text(encoding='utf-8-sig'))
    require(manifest['selected_methods'] == {'Interface': 'hermiteGals', 'transport': 'FirstJet'} and
            manifest['volume_measure'] == 'dx*dy', 'Wrong production method/measure')
    for path, value in manifest['input_file_hashes'].items():
        require(value.lower() == 'sha256:'+sha(source(path)), 'Writer/launcher provenance mismatch')
    conservation = read_csv(source('output/series/conservation.csv'))
    performance = read_csv(source('output/series/performance.csv'))
    events = read_csv(source('output/series/events.csv'))
    solver = read_csv(source('output/series/solver.csv'))
    frames = read_csv(source('output/meta/frames.csv'))
    require(len(conservation) > 1 and len(frames) > 1, 'No accepted physics history')
    require(len(conservation)-1 == manifest['accepted_steps'], 'Accepted series/manifest mismatch')
    for step, row in enumerate(conservation):
        require(int(row['step']) == step and near(row['time'], step*reference['expected_dt']), 'Accepted clock gap or mismatch')
        require(row['inventory'] == row['flux_ledger'] == 'NA' and int(row['curvature_unavailable_samples']) == 0,
                'Invalid pure-route or curvature availability diagnostics')
        require(math.isfinite(float(row['phase_volume'])) and float(row['phase_volume']) > 0, 'Invalid measured area')
        if step:
            require(near(row['dt'], reference['expected_dt']), 'Changed step size')
    require(near(conservation[-1]['time'], manifest['final_time']), 'Final clock mismatch')
    history, snapshots = [], []
    snapshot_steps = {min(frames, key=lambda r: abs(int(r['step'])-fraction*reference['expected_steps']))['step']
                      for fraction in (0, .125, .25, .5, .75, 1)}
    theta = np.linspace(0, 2*math.pi, 4097)
    for index, frame in enumerate(frames):
        step, t = int(frame['step']), float(frame['time'])
        require(int(frame['frame']) == index and 0 <= step < len(conservation) and near(t, conservation[step]['time']) and
                frame['geometry_identity'] == conservation[step]['geometry_identity'], 'Frame/accepted series mismatch')
        if index:
            require(step > int(frames[index-1]['step']), 'Non-increasing accepted frame steps')
        points, normal, curvature, weights, a, b = surface(source('output/'+frame['interface']), frame)
        f, g, exact_normal, exact_curvature = ref.geometry(points, t)
        position = np.abs(f.v)/g
        exact = ref.curve(theta, t)[0]
        sampled_distance = max(float(distance_to_segments(points, exact[:-1], exact[1:]).max()),
                               float(distance_to_segments(exact[:-1], a, b).max()))
        phi_l2, phi_band_l2 = volume_error(source('output/'+frame['volume']), frame, ref)
        w = weights.sum()
        metric = {'step': step, 'time': t, 'mesh_revision': int(frame['mesh_revision']), 'points': len(points),
                  'area_relative_error': abs(float(conservation[step]['phase_volume'])/(math.pi*ref.radius**2)-1),
                  'position_defect_max_over_radius': float(position.max()/ref.radius),
                  'position_defect_l2': math.sqrt(float(np.sum(weights*position**2)/w)),
                  'sampled_chord_hausdorff_over_radius': sampled_distance/ref.radius,
                  'normal_l2': math.sqrt(float(np.sum(weights*np.sum((normal-exact_normal)**2, axis=1))/w)),
                  'normal_max': float(np.max(np.linalg.norm(normal-exact_normal, axis=1))),
                  'curvature_relative_l2': math.sqrt(float(np.sum(weights*(curvature-exact_curvature)**2)/np.sum(weights*exact_curvature**2))),
                  'curvature_absolute_max': float(np.max(np.abs(curvature-exact_curvature))),
                  'phi_center_l2': phi_l2, 'phi_center_narrow_band_l2': phi_band_l2}
        require(all(math.isfinite(v) for v in metric.values()), 'Nonfinite accuracy metric')
        history.append(metric)
        if frame['step'] in snapshot_steps:
            snapshots.append((frame, (a, b)))
    area_peak = max(abs(float(r['phase_volume'])/(math.pi*ref.radius**2)-1) for r in conservation)
    peak = {name: max(r[name] for r in history) for name in
            ('position_defect_max_over_radius', 'sampled_chord_hausdorff_over_radius', 'normal_l2', 'curvature_relative_l2')}
    gates = reference['gates']
    complete = (record['status'] == 'completed' and manifest['exit_status'] == 'success' and
                manifest['accepted_steps'] == reference['expected_steps'] and near(manifest['final_time'], reference['expected_end_time']))
    actual_regrids = len([r for r in events if r['event'] == 'regrid'])
    checks = {'completed_window': complete, 'frame_count': len(frames) == reference['expected_frames'],
              'regrid_count': actual_regrids == reference['expected_regrids'],
              'no_numerical_rejection': all(r['status'] == 'accepted' for r in solver),
              'peak_area': area_peak <= gates['peak_relative_area_error'],
              'position': peak['position_defect_max_over_radius'] <= gates['maximum_position_defect_over_radius'],
              'normal': peak['normal_l2'] <= gates['maximum_normal_l2'],
              'curvature': peak['curvature_relative_l2'] <= gates['maximum_curvature_relative_l2'],
              'final_position': history[-1]['position_defect_max_over_radius'] <= gates['final_position_defect_over_radius']}
    returns = [r for r in history if any(near(r['time'], t) for t in reference['expected_return_times'])]
    checks['return_frames'] = len(returns) == len(reference['expected_return_times'])
    totals = [float(r['wall_seconds']) for r in performance if r['stage'] == 'step_total' and int(r['step']) > 0]
    require(len(totals) == len(conservation)-1 and all(math.isfinite(v) and v >= 0 for v in totals), 'Bad step timing history')
    summary = {'schema_version': 1, 'status': 'PASS' if all(checks.values()) else 'FAIL', 'run': str(run),
               'scope': reference['scope'], 'reference': reference['reference'], 'reference_self_test': self_test,
               'accepted_steps': manifest['accepted_steps'], 'final_time': manifest['final_time'], 'frames_verified': len(frames),
               'regrids': actual_regrids, 'exact_area': math.pi*ref.radius**2,
               'peak_relative_area_error': area_peak, 'peak_sample_metrics': peak, 'return_frames': returns,
               'gates': gates, 'checks': checks, 'launcher_seconds': record['seconds'],
               'step_seconds': {'median': float(np.median(totals)), 'p95': float(np.quantile(totals, .95)), 'total': sum(totals)},
               'sources_sha256': sources, 'analyzer_sha256': sha(Path(__file__).resolve()),
               'reference_curve_segments': 4096,
               'metric_note': 'Position defect is |phi_exact|/|grad(phi_exact)|, not exact distance. Hausdorff is bidirectional sampled point-to-chord distance. L2 uses saved VTP chord arc-length weights. Reference normals and level-set curvature are evaluated at the numerical sample positions. Full geometry and convergence norms are not asserted.'}
    directory = run/'output'/'analysis'/('reference_'+datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f'))
    directory.mkdir(parents=True)
    with (directory/'accuracy_history.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(history[0]))
        writer.writeheader(); writer.writerows(history)
    plots(directory, history, conservation, performance, snapshots, ref)
    (directory/'verification.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    text = ['# Independent Interface reference verification', '', f"Status: **{summary['status']}**",
            f"Accepted steps: {manifest['accepted_steps']}; final time: {manifest['final_time']:.12g}; frames: {len(frames)}.",
            '', '| Observation | Value |', '| --- | ---: |', f'| Peak relative area error | {area_peak:.8e} |']
    text.extend(f'| Peak {name} | {value:.8e} |' for name, value in peak.items())
    text.extend(['', '![Interface evolution](interface_evolution.png)', '', '![Error histories](error_history.png)', '',
                 summary['metric_note'], '', 'This checks a prescribed-velocity Interface run. It does not establish mesh/time convergence or Hydro/EHD/jet-breakup acceptance.'])
    (directory/'REPORT.md').write_text('\n'.join(text)+'\n', encoding='utf-8')
    print(json.dumps({'status': summary['status'], 'analysis': str(directory), 'area_error': area_peak,
                      'peak_metrics': peak, 'failed_checks': [k for k, v in checks.items() if not v]}, indent=2))
    return 0 if summary['status'] == 'PASS' else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path)
    parser.add_argument('--accepted-prefix', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(reference_self_test(RadialTwistCircle({'flow_center': [.5, .5], 'circle_center': [.65, .53],
                                                              'radius': .13, 'speed': 40., 'period': 1.})), indent=2))
        return 0
    require(args.run is not None, 'Provide --run or --self-test')
    return analyze(args.run, args.accepted_prefix)


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError, ET.ParseError) as error:
        print('REFERENCE_ANALYSIS_ERROR: '+str(error), flush=True)
        raise SystemExit(2)
