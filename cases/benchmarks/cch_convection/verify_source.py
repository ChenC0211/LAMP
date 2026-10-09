"""Read-only physical reference for the existing advector/closure/boundary modes.

Usage: python verify_source.py MODE RUN_DIRECTORY [--output result.json]
Uses only the Python standard library; failures exit 1. Does not reconstruct a
velocity field, read patch weights, repair output, or certify the full PDE run.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

NODES = (-math.sqrt(.6), 0., math.sqrt(.6))
SAMPLE_NODES = NODES + (-.5, .5)
WEIGHTS = (5/9, 8/9, 5/9)
TOL = 1e-12


def require(ok, why):
    if not ok:
        raise ValueError(why)


def rows(directory, name):
    with (directory / (name + '.csv')).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)), name + ' header')
        for row in reader:
            require(None not in row and None not in row.values(), name + ' width')
            result = {k: v if k == 'kind' else float(v) for k, v in row.items()}
            require(all(math.isfinite(v) for v in result.values() if isinstance(v, float)), name + ' nonfinite')
            yield result


def key(row, fields):
    return tuple(row[x] for x in fields.split(','))


def unique(table, fields):
    result = {}
    for r in table:
        k = key(r, fields)
        require(k not in result, 'duplicate ' + fields)
        result[k] = r
    return result


def coefficients(row, prefix):
    names = [x for x in row if x.startswith(prefix) and x[len(prefix):].isdigit()]
    require(names == [prefix + str(k) for k in range(5)], 'complete quartic schema ' + prefix)
    return [row[x] for x in names]


def value(c, t):
    return math.fsum(v * t**k for k, v in enumerate(c))


def norm(c, length):
    return math.fsum(abs(v) * length**k for k, v in enumerate(c))


def integral(c, length, radius, rz, axis, reverse=False):
    # GL3 integrates degree <=5 exactly (quartic velocity times linear radius).
    return math.fsum(w * value(c, length*(1+z)/2) * length/2 *
                     (2*math.pi*(radius + (-1 if reverse else 1)*length*(1+z)/2
                                 if axis == 1 else radius) if rz else 1)
                     for z, w in zip(NODES, WEIGHTS))


class Audit:
    def __init__(self):
        self.maxima = defaultdict(float)
        self.counts = defaultdict(int)

    def error(self, a, b, work, label):
        require(work >= 0 and all(math.isfinite(x) for x in (a, b, work)), label + ' finite budget')
        error = abs(a-b)/work if work else (0. if a == b else math.inf)
        require(error <= TOL, f'{label}: {error:.17g} > {TOL}, values={a},{b}, work={work}')
        self.maxima[label] = max(self.maxima[label], error)
        self.counts[label] += 1

    def atom(self, r, both=True):
        length = r['by']-r['ay'] if r['axis'] == 0 else r['bx']-r['ax']
        require(length > 0, 'atom width')
        area = length*(math.pi*(r['ax']+r['bx']) if r['rz'] else 1)
        self.error(r['area'], area, area, 'atom physical area')
        if 'length' in r:
            self.error(r['length'], length, length, 'atom length')
        c = coefficients(r, 'c')
        work = norm(c, length)
        q = integral(c, length, r['ax'], r['rz'], r['axis'])
        if both:
            e = coefficients(r, 'e')
            work = max(work, norm(e, length))
            for z in (0., .173, .5, .827, 1.):
                self.error(value(c, length*z), value(e, length*(1-z)), work, 'double chart')
            for actual in (r['integral'], r['end_integral'], integral(e, length, r['bx'], r['rz'], r['axis'], True)):
                self.error(actual, q, area*work, 'atom GL3 integral')
        return q, area*work, work


def analytic(r, x, y):
    revision = r['revision']
    if r['rz']:
        return ((.13+.05*revision)*x, .2+.11*revision+.155*x*x+.19*y)
    return (.3+.11*revision+.17*x-.23*y, -.2-.07*revision+.31*x+.19*y)


def applicable(r):
    # Original 56/64 analytic fixtures. Preparation diagnostics never select gates.
    return not (r['n'] == 4 and r['selection'] in (2, 3))


def advector(directory, audit):
    fields = 'rz,amr,n,selection,revision'
    summaries = unique(rows(directory, 'advector_summary'), 'rz,amr,n,selection')
    require(len(summaries) == 64, '64 original layouts')
    expected = {(rz, amr, n, sel) for rz in (0, 1) for amr in (0, 1) for n in (4, 8, 16, 32) for sel in range(4)}
    require(set(summaries) == expected, 'original layout matrix')
    for r in summaries.values():
        require(all(r[x] == 1 for x in ('geometry_shared', 'updates_shared', 'fresh_identity')), 'shared owner')
    faces = unique(rows(directory, 'advector_faces'), fields + ',u')
    totals = defaultdict(list)
    area_sums = defaultdict(list)
    atom_counts = defaultdict(int)
    field_work = defaultdict(float)
    for r in rows(directory, 'advector_atoms'):
        g = key(r, fields)
        q, work, amplitude = audit.atom(r)
        field_work[g] = max(field_work[g], amplitude)
        require(r['atom'] == atom_counts[g], 'atom order/coverage')
        atom_counts[g] += 1
        if r['primal'] >= 0:
            f = g+(r['primal'],)
            require(f in faces, 'atom primal ID')
            totals[f].append(q)
            area_sums[f].append(r['area'])
        if applicable(r):
            length = r['length']
            for z in (0., .173, .5, .827, 1.):
                x = r['ax']+(r['bx']-r['ax'])*z
                y = r['ay']+(r['by']-r['ay'])*z
                audit.error(value(coefficients(r, 'c'), length*z), analytic(r, x, y)[int(r['axis'])], 1., 'original affine atom')
    identities = {}
    for f, r in faces.items():
        g = f[:-1]
        require(g[:-1] in summaries and g[-1] in (0, 1), 'face layout')
        require(key(r, 'step,clock_begin,clock_end,position') == (9., .2, .3, 1.), 'original full clock')
        identity = key(r, 'geometry,advector')
        require(identities.setdefault(g, identity) == identity and min(identity) > 0, 'source identity')
        audit.error(math.fsum(area_sums[f]), r['area'], r['area'], 'native atom coverage')
        budget = r['area']*field_work[g]+abs(r['area']*r['input'])
        audit.error(math.fsum(totals[f]), r['area']*r['input'], budget, 'native physical mean')
        audit.error(r['integral'], math.fsum(totals[f]), budget, 'native observed sum')
        if applicable(r):
            for actual, exact in zip((r['point_x'], r['point_y']), analytic(r, r['x'], r['y'])):
                audit.error(actual, exact, 1., 'original affine point')
    require(len(identities) == 128 and len({v[1] for v in identities.values()}) == 128, '128 distinct revisions')
    geometry_ids = defaultdict(set)
    for g, identity in identities.items():
        geometry_ids[g[:3]].add(identity[0])
    require(len(geometry_ids) == 16 and all(len(v) == 1 for v in geometry_ids.values()), '16 geometry identities')
    cell_counts = defaultdict(int)
    for r in rows(directory, 'advector_cells'):
        g = key(r, fields)
        require(r['p'] == cell_counts[g], 'cell ID coverage')
        cell_counts[g] += 1
        exact = 2*(.13+.05*r['revision'])+.19 if r['rz'] else .36
        for name in ('divergence', 'original_divergence'):
            audit.error(r[name], exact, 1., 'original mean divergence')
        if applicable(r):
            for actual, target in zip((r['velocity_x'], r['velocity_y']), analytic(r, r['x'], r['y'])):
                audit.error(actual, target, 1., 'original affine cell')
    for g, count in atom_counts.items():
        require(count == summaries[g[:-1]]['atoms'], 'all actual atoms')
        require(sum(1 for f in faces if f[:-1] == g) == summaries[g[:-1]]['faces'], 'all native faces')
    require(set(cell_counts) == set(identities), 'cell revisions')
    return dict(layouts=64, revisions=128, analytic_layouts=56, faces=len(faces), atoms=sum(atom_counts.values()), cells=sum(cell_counts.values()))


def closure(directory, audit):
    fields = 'rz,amr,basis'
    sources = unique(rows(directory, 'closure_sources'), fields)
    status = unique(rows(directory, 'closure_status'), fields)
    require(len(sources) == 208 and set(sources) == set(status), 'original 208 basis sources')
    require(all(r['accepted'] == 1 for r in status.values()), 'actual source rejection')
    require(len({r['advector'] for r in sources.values()}) == 208, 'unique basis source identities')
    for r in sources.values():
        require(key(r, 'step,begin,end,position') == (0., 0., 0., 0.) and r['geometry'] > 0, 'closure full clock')
    atoms = unique(rows(directory, 'closure_atoms'), fields+',atom')
    faces = unique(rows(directory, 'closure_faces'), fields+',u')
    cells = unique(rows(directory, 'closure_cells'), fields+',p')
    values, work = {}, defaultdict(float)
    for k, r in atoms.items():
        q, w, amplitude = audit.atom(r)
        values[k] = q
        work[k[:-1]] = max(work[k[:-1]], amplitude)
    sample_ids = set()
    for r in rows(directory, 'closure_samples'):
        k = key(r, fields+',atom')
        require(k in atoms and r['node'] in range(5), 'physical sample ID')
        sample_key = k+(r['node'],)
        require(sample_key not in sample_ids, 'duplicate sample')
        sample_ids.add(sample_key)
        a = atoms[k]
        z = (1+SAMPLE_NODES[int(r['node'])])/2
        for name, lo, hi in (('x','ax','bx'), ('y','ay','by')):
            audit.error(r[name], a[lo]+(a[hi]-a[lo])*z, a['length'], 'physical sample coordinate')
        require(0 <= r['pressure'] < sources[k[:-1]]['cells'], 'sample pressure ID')
        audit.error(r['normal'], value(coefficients(a, 'c'), a['length']*z), work[k[:-1]], 'physical source sample')
    require(len(sample_ids) == 5*len(atoms), 'all five physical samples')
    views = defaultdict(list)
    seen = set()
    for r in rows(directory, 'closure_views'):
        g = key(r, fields)
        k = g+(r['kind'], r['id'])
        pair = k+(r['atom'],)
        require(pair not in seen and g+(r['atom'],) in atoms, 'shared view ID')
        seen.add(pair)
        require(r['kind'] in ('primal', 'dual'), 'view kind')
        require(0 <= r['id'] < sources[g]['faces' if r['kind'] == 'primal' else 'segments'], 'view ID range')
        views[k].append(r)
    cell_sums, native_sums, cell_gross, cell_work = (defaultdict(list) for _ in range(4))
    internal = amr_split = zero_gross = 0
    for k, entries in views.items():
        g, kind, identifier = k[:3], k[3], k[4]
        first = entries[0]
        axis = int(first['axis'])
        a, b = ('ay','by') if axis == 0 else ('ax','bx')
        ordered = sorted((atoms[g+(r['atom'],)] for r in entries), key=lambda r:r[a])
        cursor = first[a]
        for row in ordered:
            require(row['axis'] == axis and row[a] == cursor, 'real subface gap/overlap')
            require(row['ax' if axis == 0 else 'ay'] == first['ax' if axis == 0 else 'ay'], 'atom plane')
            cursor = row[b]
        require(cursor == first[b], 'complete view endpoint')
        require(all(all(r[x] == first[x] for x in ('axis','negative','positive','ax','ay','bx','by','area')) for r in entries), 'view metadata consistency')
        audit.error(math.fsum(r['area'] for r in ordered), first['area'], first['area'], 'view physical area')
        flux = math.fsum(values[g+(r['atom'],)] for r in entries)
        if first['negative'] >= 0 and first['positive'] >= 0:
            require(first['negative'] != first['positive'], 'internal distinct neighbors')
            internal += 1
        if kind == 'dual':
            require(all(-1 <= first[side] < sources[g]['faces'] for side in ('negative','positive')), 'dual neighbor ID')
        if g[1] and len(ordered) > 1:
            amr_split += 1
        if kind != 'primal':
            continue
        f = faces[g+(identifier,)]
        require(all(first[x] == f[x] for x in ('axis','negative','positive','ax','ay','bx','by','area')), 'native source view')
        require(f['input'] == (1/f['area'] if identifier == g[2] else 0), 'original basis input')
        expected = f['area']*f['input']
        budget = f['area']*work[g]+abs(expected)
        audit.error(flux, expected, budget, 'basis native physical mean')
        for side, sign in (('negative',1), ('positive',-1)):
            if f[side] >= 0:
                c = g+(f[side],)
                require(c in cells, 'native cell incidence')
                cell = cells[c]
                # Incidence orientation and actual coarse/fine subface geometry
                # are checked without the production graph or divergence code.
                plane = ('bx' if sign == 1 else 'ax') if axis == 0 else ('by' if sign == 1 else 'ay')
                require(f['ax' if axis == 0 else 'ay'] == cell[plane], 'physical outward incidence')
                require(cell[a] <= f[a] < f[b] <= cell[b], 'physical native subface')
                cell_sums[c].append(sign*flux)
                native_sums[c].append(sign*expected)
                cell_gross[c].append(abs(expected))
                cell_work[c].append(budget)
    for k, r in cells.items():
        expected = r['original_divergence']*r['volume']
        gross = math.fsum(cell_gross[k])
        volume = (r['bx']-r['ax'])*(r['by']-r['ay'])*(math.pi*(r['ax']+r['bx']) if r['rz'] else 1)
        audit.error(r['volume'], volume, volume, 'cell physical volume')
        audit.error(expected, math.fsum(native_sums[k]), gross, 'independent native D_I')
        if gross == 0:
            require(r['original_divergence'] == 0, 'zero gross original D_I')
            zero_gross += 1
        audit.error(math.fsum(cell_sums[k]), expected, math.fsum(cell_work[k]), 'basis cell boundary integral')
    for g, r in sources.items():
        for data, count in ((atoms,'atoms'), (faces,'faces'), (cells,'cells')):
            require(sum(k[:-1] == g for k in data) == r[count], 'basis '+count+' coverage')
        require(sum(k[:3] == g for k in views) == r['faces']+r['segments'], 'all primal/dual views')
        require(sum(k[:3] == g for k in seen) == r['view_entries'], 'actual view entries')
        active = {f['u'] for k,f in faces.items() if k[:-1] == g and f['area'] > 0}
        require({k[2] for k in sources if k[:2] == g[:2]} == active, 'every nonaxis basis')
    require(zero_gross > 0 and internal > 0 and amr_split > 0, 'original closure witnesses')
    return dict(bases=208, atoms=len(atoms), samples=len(sample_ids), faces=len(faces), cells=len(cells), view_entries=len(seen), internal_shared_views=internal, amr_split_views=amr_split, zero_gross_cells=zero_gross)


def boundary(directory, audit):
    fields = 'rz,amr,n,selection,bound'
    summaries = unique(rows(directory, 'boundary_summary'), fields)
    require(set(summaries) == {(rz, amr, n, sel, bound) for rz in (0,1) for amr in (0,1) for n in (4,8,16,32) for sel in (0,1) for bound in (0,1)}, 'original 64 boundary groups')
    faces = unique(rows(directory, 'boundary_faces'), fields+',u')
    lines = defaultdict(list)
    totals, field_work, atom_counts = defaultdict(list), defaultdict(float), defaultdict(int)
    for r in rows(directory, 'boundary_atoms'):
        g = key(r, fields)
        q, work, amplitude = audit.atom(r, False)
        field_work[g] = max(field_work[g], amplitude)
        require(r['id'] == atom_counts[g], 'boundary atom IDs')
        atom_counts[g] += 1
        audit.error(r['ql']+r['qg'], q, work, 'phase total')
        audit.error(r['F'], q, work, 'equal density unique F')
        if r['primal'] >= 0:
            totals[g+(r['primal'],)].append(q)
            if r['bound'] and faces[g+(r['primal'],)]['outward']:
                require(all(r[x] == 0 for x in ('ql','qg','F')), 'bound normal atom exactly zero')
        plane = r['ax'] if r['axis'] == 0 else r['ay']
        lines[g+(r['axis'],plane)].append((r,q))
    for group in lines.values():
        group.sort(key=lambda item:item[0]['ay'] if item[0]['axis'] == 0 else item[0]['ax'])
    for k, r in faces.items():
        budget = r['area']*field_work[k[:-1]]+abs(r['area']*r['mean'])
        audit.error(math.fsum(totals[k]), r['area']*r['mean'], budget, 'boundary native mean')
        audit.error(r['F'], math.fsum(totals[k]), budget, 'boundary native F')
    internal = segment_count = 0
    counts = defaultdict(int)
    for r in rows(directory, 'boundary_segments'):
        g = key(r, fields)
        require(r['id'] == counts[g], 'boundary segment IDs')
        counts[g] += 1
        lo, hi, plane = ('ay','by','ax') if r['axis'] == 0 else ('ax','bx','ay')
        chosen = [(a,q) for a,q in lines[g+(r['axis'],r[plane])] if a[lo] >= r[lo] and a[hi] <= r[hi]]
        cursor = r[lo]
        for a,_ in chosen:
            require(a[lo] == cursor, 'dual real subface partition')
            cursor = a[hi]
        require(cursor == r[hi], 'dual segment full coverage')
        audit.error(r['F'], math.fsum(q for _,q in chosen), r['area']*field_work[g], 'dual shared F')
        if r['negative'] >= 0 and r['positive'] >= 0:
            internal += 1
        segment_count += 1
    cell_counts = defaultdict(int)
    for r in rows(directory, 'boundary_cells'):
        g = key(r, fields)
        cell_counts[g] += 1
        # Average divergence is a volume integral, not the old secant point diagnostic.
        audit.error(r['B_div'], r['DI_div'], field_work[g]/min(r['bx']-r['ax'],r['by']-r['ay']), 'boundary mean divergence')
    for g,r in summaries.items():
        require(r['retained_law'] == 1 and key(r,'step,begin,end,position') == (1.,0.,.1,1.), 'boundary provenance')
        require(atom_counts[g] == r['atoms'] and counts[g] == r['segments'] and cell_counts[g] == r['cells'], 'boundary row coverage')
    require(internal > 0, 'internal shared F coverage')
    return dict(revisions=64, atoms=sum(atom_counts.values()), faces=len(faces), segments=segment_count, cells=sum(cell_counts.values()), internal_segments=internal)


def degree_controls():
    # Exact rational monomial antiderivatives independently check GL3. Include
    # genuine cubic/quartic content so a quadratic-only reader must fail.
    from fractions import Fraction as F
    worst = truncated = 0.
    for count in (3, 5):
        c = [F(k+1,8) for k in range(count)]
        for rz in (False, True):
            for axis in (0, 1):
                length, radius = F(3,8), F(1,4)
                primitive = sum(v*length**(k+1)/(k+1) for k,v in enumerate(c))
                radial = sum(v*length**(k+2)/(k+2) for k,v in enumerate(c)) if axis == 1 else 0
                exact = float(radius*primitive+radial)*2*math.pi if rz else float(primitive)
                actual = integral([float(v) for v in c],float(length),float(radius),rz,axis)
                worst = max(worst, abs(actual-exact)/abs(exact))
                if count == 5:
                    bad = integral([float(v) for v in c[:3]],float(length),float(radius),rz,axis)
                    truncated = max(truncated,abs(bad-exact)/abs(exact))
    require(worst < TOL and truncated > 1e-3, 'complete-degree reference control')
    # Shifted P3 and (t+1/4)P3 have zero values at every GL3 node and zero
    # integral with a linear metric. The two extra queries must detect them.
    blind, detected = 0., math.inf
    for c in ((-.05,.6,-1.5,1.),(-.0125,.1,.225,-1.25,1.)):
        blind = max(blind, max(abs(value(c,(1+z)/2)) for z in NODES))
        detected = min(detected,max(abs(value(c,(1+z)/2)) for z in SAMPLE_NODES[3:]))
    require(blind < 1e-15 and detected > 1e-3, 'five-point complete-degree control')
    return dict(positive_relative_error=worst, truncation_relative_error=truncated,
                gl3_invisible_error=blind, five_point_detection=detected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('advector','closure','boundary'))
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    audit = Audit()
    result = dict(status='PASS', mode=args.mode, degree_controls=degree_controls())
    result['coverage'] = globals()[args.mode](args.run_directory/'output/series', audit)
    result['maximum_errors'] = dict(audit.maxima)
    result['checks'] = dict(audit.counts)
    output = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(output, encoding='utf-8')
    print(output)


if __name__ == '__main__':
    main()
