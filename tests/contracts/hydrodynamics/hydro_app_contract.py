#!/usr/bin/env python3
"""VOF App input/clock/output contract; archived FirstJet gates are not reclassified."""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import runpy
import shutil
import subprocess
import sys
import uuid


FILES = ("run.dict", "geometry.dict", "mesh.dict",
         "interface.dict", "materials.dict", "hydro.dict")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def edit_once(path, pattern, replacement):
    text, count = re.subn(pattern, replacement, path.read_text(encoding="utf-8"))
    require(count == 1, "one copied input edit required: " + pattern)
    path.write_text(text, encoding="utf-8")


def input_number(path, key):
    matches = re.findall(r"\b" + re.escape(key) + r"\s+([^;]+);", path.read_text(encoding="utf-8"))
    require(len(matches) == 1, "one actual copied input value: " + key)
    value = float(matches[0])
    require(math.isfinite(value) and value > 0, "finite positive input " + key)
    return value


def field_contract(path, frame, helpers, previous_time, base_cells):
    """Check published data dimensions and clocks, with no geometry/solution oracle."""
    text = path.read_text(encoding="utf-8")
    matches = list(re.finditer(r"(?m)^\s*AUXDATA\s+\w+\s*=.*$", text))
    require(matches, "field provenance missing: " + str(path))
    header, tokens = text[:matches[-1].end()], text[matches[-1].end():].split()
    attribute, auxiliary = helpers["attribute"], helpers["auxiliary"]
    equal = helpers["equal_float"]
    metadata = auxiliary(header)
    require(metadata.get("Primary") == "alpha" and metadata.get("Representation") == "localPlic",
            "field retains the single VOF inventory/PLIC owner")
    require(int(metadata["AcceptedStep"]) == int(frame["step"]) and
            int(metadata["GeometryIdentity"]) == int(frame["geometry_identity"]),
            "field accepted step and geometry agree with authoritative frame")
    time = float(frame["time"])
    equal(float(attribute(header, "SOLUTIONTIME")), time, "field solution time")
    require(attribute(header, "DATAPACKING") == "BLOCK", "field BLOCK format")
    nodes, cells = int(attribute(header, "NODES")), int(attribute(header, "ELEMENTS"))
    require(nodes > 0 and cells > 0, "field real nodes/elements")
    names = helpers["variables"](header)
    kind = attribute(header, "ZONETYPE")
    if kind == "FEQUADRILATERAL":
        require(names == ["X", "Y"] + helpers["VOLUME_FIELDS"], "native volume variables")
        require(cells > base_cells, "initial interface band produces real refined leaf cells")
        count, arity = 2 * nodes + len(helpers["VOLUME_FIELDS"]) * cells, 4
        equal(float(metadata["VTK_TimeValue"]), time, "native accepted time")
        expected = 0.0 if int(frame["step"]) == 0 else (previous_time + time) * .5
        equal(float(metadata["VTK_pressure_sample_time"]), expected, "actual midpoint pressure clock")
    else:
        require(kind == "FELINESEG", "native interface uses independent PLIC segments")
        require(names[:2] == ["X", "Y"], "interface coordinates")
        count, arity = nodes * len(names), 2
    require(len(tokens) == count + arity * cells, "real field block/connectivity dimensions")
    require(all(helpers["number"](value) is not None for value in tokens[:count]), "finite published payload")
    connection = tokens[count:]
    require(all(re.fullmatch(r"\d+", value) and 1 <= int(value) <= nodes for value in connection),
            "integer actual node connectivity bounds")
    if kind == "FEQUADRILATERAL":
        x, y = list(map(float, tokens[:nodes])), list(map(float, tokens[nodes:2 * nodes]))
        widths = set()
        for index in range(cells):
            vertices = [int(v) - 1 for v in connection[4 * index:4 * index + 4]]
            widths.update((max(x[v] for v in vertices) - min(x[v] for v in vertices),
                           max(y[v] for v in vertices) - min(y[v] for v in vertices)))
        require(min(widths) > 0 and len(widths) >= 2, "actual coarse and fine leaf widths")
        equal(max(widths) / min(widths), 2.0, "actual static L1 2:1 grid")


def output_contract(directory, dt, end_time, base_cells, helpers):
    output = directory / "output"
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest["status"] == "COMPLETED" and manifest["accepted_step"] == 2,
            "real two-step completed status")
    helpers["equal_float"](float(manifest["accepted_time"]), end_time, "actual final clock")
    require(manifest.get("primary") == "alpha" and manifest.get("scalar_representation") == "localPlic",
            "manifest declares inventory alpha as the single primary interface")
    require(manifest.get("volume_format") == manifest.get("series_format") == "tecplotDat" and
            manifest.get("frame_index") == "meta/frames.csv", "native DAT-only protocol")
    require(not list(output.rglob("*.vtu")) and not list(output.rglob("*.pvd")) and
            not (output / "meta/tecplot_frames.csv").exists(), "no duplicate VTK/Tecplot index")
    require(not list((output / "series").rglob("*.csv")), "series DAT-only")
    require(not (output / "series/sevof_energy.dat").exists(), "no second area-state history")
    fields, frames = helpers["csv_rows"](output / "meta/frames.csv")
    require(fields == ["frame", "step", "time", "geometry_identity", "volume", "interface"],
            "one authoritative accepted-frame schema")
    require([int(frame["step"]) for frame in frames] == [0, 1, 2] and
            manifest["complete_frames"] == 3, "initial plus both accepted frames")
    conservation = helpers["table_history"](output / "series/conservation.dat", "conservation")
    require([int(float(row["tokens"]["step"])) for row in conservation] == [0, 1, 2],
            "actual accepted history rows")
    previous_time = 0.0
    for index, (frame, row) in enumerate(zip(frames, conservation)):
        require(all(row["available"][name] == 1 for name in ("step", "time", "dt", "geometry_identity")),
                "accepted clock and geometry are available")
        actual_time = float(row["tokens"]["time"])
        helpers["equal_float"](actual_time, index * dt, "two steps retain configured dt")
        helpers["equal_float"](float(row["tokens"]["dt"]), 0.0 if index == 0 else dt,
                               "series actual step length")
        helpers["equal_float"](float(frame["time"]), actual_time, "frame and series accepted clock")
        require(int(frame["geometry_identity"]) == int(float(row["tokens"]["geometry_identity"])),
                "frame and series share accepted geometry")
        for column in ("volume", "interface"):
            relative = Path(frame[column])
            require(not relative.is_absolute() and ".." not in relative.parts and relative.suffix == ".dat",
                    "authoritative frame contains a native relative DAT path")
            field_contract(output / relative, frame, helpers, previous_time, base_cells)
        previous_time = actual_time
    require(int(frames[-1]["geometry_identity"]) == manifest["geometry_identity"],
            "last frame is the final accepted owner")
    for name in FILES:
        require(manifest["input_hashes"]["system/" + name] == "sha256:" + digest(directory / "system" / name),
                "manifest hashes the actual copied dictionary: " + name)
    _, trace = helpers["csv_rows"](output / "meta/configuration_trace.csv")
    require(not any(row["field"].startswith("Time.Adaptive.") for row in trace), "fixed mode must not invent adaptive defaults")
    require(not any(row["field"].startswith("Fields.Constants.") for row in trace), "retired Fields constant trace owner")
    constants = [row for row in trace if row["field"].startswith("Geometry.Constants.")]
    require(all(Path(row["file"]).name == "geometry.dict" and row["used"] == "true" for row in constants),
            "shared constants identify actual geometry.dict source")
    rows = [row for row in trace if row["field"] == "Interface.Initialization.subcellsPerAxis"]
    require(len(rows) == 1, "one initialization trace owner")
    require(rows[0]["effective_value"] == "8" and rows[0]["defaulted"] == "true" and rows[0]["used"] == "true" and
            rows[0]["consumer"] == "Interface.GeometricVofInitialization.fixedSubcells",
            "default eight subcells is actually consumed by the new initialization")
    methods = [row for row in trace if row["field"] == "Interface.method"]
    require(len(methods) == 1 and methods[0]["effective_value"] == "vof" and
            methods[0]["consumer"] == "Interface.VOF.dispatch", "one VOF input dispatch owner")
    effective = output / "meta/effective"
    require({path.name for path in effective.iterdir()} == set(FILES), "same six effective documents")
    interface = (effective / "interface.dict").read_text(encoding="utf-8")
    require(re.search(r"\bmethod\s+vof\s*;", interface) and
            re.search(r"\bInitialization\s*\{\s*subcellsPerAxis\s+8\s*;\s*\}", interface),
            "effective interface retains the consumed default and VOF method")
    info = (output / "meta/build_info.txt").read_text(encoding="utf-8")
    require("interface_method=vof" in info and "primary=alpha" in info,
            "actual observed VOF inventory primary is recorded")
    return {"directory": str(directory), "accepted_steps": 2, "accepted_time": end_time,
            "dt": dt, "frames": len(frames), "subcells_default": 8}


def adaptive_input_contract(task, invoke, helpers):
    """Input/schema/provenance only; no claim about proactive rate selection."""
    template = task / "hydro_rayleigh_taylor/system"
    seed = input_number(template / "run.dict", "dt")
    end_time = 1.5 * seed
    values = {"minDt": seed / 16, "maxDt": 2 * seed, "targetCourant": .25}

    def prepare(label, adaptive, retry=None):
        directory = task / ("adaptive_input_" + label)
        shutil.copytree(template, directory / "system")
        run_path = directory / "system/run.dict"
        source = run_path.read_text(encoding="utf-8")
        source, count = re.subn(r"\bStepRetry\s*\{[^{}]*\}", "", source)
        require(count <= 1 and "Adaptive" not in source, "copied fixed fixture must have one time policy owner")
        source, count = re.subn(r"\bendTime\s+[^;]+;", f"endTime {end_time:.17g};", source)
        require(count == 1, "one actual endTime in copied input")
        blocks = []
        if adaptive is not None:
            blocks.append("Adaptive { " + " ".join(f"{key} {value};" for key, value in adaptive.items()) + " }")
        if retry is not None:
            blocks.append("StepRetry { " + " ".join(f"{key} {value};" for key, value in retry.items()) + " }")
        source, count = re.subn(r"\bTime\s*\{", lambda match: match.group() + "\n" + "\n".join(blocks), source)
        require(count == 1, "one actual Time block")
        run_path.write_text(source, encoding="utf-8")
        return directory

    defaults = prepare("defaults", values)
    checked = invoke(defaults, True)
    require(checked.returncode == 0 and "HYDRO_INPUT_OK" in checked.stdout and "[INIT]" not in checked.stdout and
            not (defaults / "output").exists(), "Adaptive accepts noninteger endTime/seed before compute")
    explicit = dict(values, safety=1, maxGrowth=1)
    combined = prepare("explicit_matching_retry", explicit,
                       {"minDt": values["minDt"], "maxRetries": 4, "recoverAfter": 0})
    checked = invoke(combined, True)
    require(checked.returncode == 0 and "HYDRO_INPUT_OK" in checked.stdout and "[INIT]" not in checked.stdout and
            not (combined / "output").exists(), "explicit adaptive bounds and matching zero-recovery retry are accepted")
    invalid = [
        ("fixed_noninteger", None, None),
        ("min_above_seed", dict(values, minDt=2 * seed), None),
        ("max_below_seed", dict(values, maxDt=seed / 2), None),
        ("target_above_kernel_limit", dict(values, targetCourant=.6), None),
        ("safety_above_one", dict(values, safety=1.01), None),
        ("safety_zero", dict(values, safety=0), None),
        ("growth_below_one", dict(values, maxGrowth=.99), None),
        ("nonfinite_min", dict(values, minDt="nan"), None),
        ("nonfinite_max", dict(values, maxDt="inf"), None),
        ("retry_min_conflict", values, {"minDt": values["minDt"] / 2, "maxRetries": 4, "recoverAfter": 0}),
        ("retry_growth_conflict", values, {"minDt": values["minDt"], "maxRetries": 4, "recoverAfter": 1}),
        ("missing_target", {key: value for key, value in values.items() if key != "targetCourant"}, None),
        ("unknown_adaptive_key", dict(values, unconsumedTimeKnob=1), None),
    ]
    failures = []
    for label, policy, retry in invalid:
        directory = prepare("reject_" + label, policy, retry)
        result = invoke(directory, True)
        require(result.returncode != 0 and re.search(r"run\.dict:\d+:", result.stdout) and
                "HYDRO_INPUT_OK" not in result.stdout and "[INIT]" not in result.stdout and
                not (directory / "output").exists(), "invalid adaptive input must reject before compute: " + label)
        failures.append({"case": label, "directory": str(directory), "rejected_before_compute": True})
    # One tiny copied run observes only the real consumed effective values and
    # default provenance. It neither sets a fixed accepted-step count nor
    # asserts any particular proactive numerical timestep schedule.
    result = invoke(defaults, False)
    require(result.returncode == 0 and "[COMPLETED]" in result.stdout, "adaptive input provenance run failed")
    output = defaults / "output"
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest["status"] == "COMPLETED" and manifest["accepted_step"] > 0,
            "provenance run must preserve actual completion")
    helpers["equal_float"](float(manifest["accepted_time"]), end_time, "noninteger adaptive input final physical time")
    effective = (output / "meta/effective/run.dict").read_text(encoding="utf-8")
    match = re.search(r"\bAdaptive\s*\{([^{}]*)\}", effective)
    require(match is not None, "effective Time.Adaptive block")
    encoded = dict(re.findall(r"\b(\w+)\s+([^;]+);", match.group(1)))
    expected = dict(values, safety=.9, maxGrowth=1.25)
    require(set(encoded) == set(expected), "effective policy explicitly carries all five consumed values")
    for key, value in expected.items():
        require(math.isclose(float(encoded[key]), value, rel_tol=1e-15, abs_tol=0), "effective adaptive value " + key)
    _, trace = helpers["csv_rows"](output / "meta/configuration_trace.csv")
    for key, value in expected.items():
        rows = [row for row in trace if row["field"] == "Time.Adaptive." + key]
        require(len(rows) == 1 and rows[0]["used"] == "true" and rows[0]["consumer"] == "Runtime.VOF.Adaptive",
                "one typed adaptive trace owner: " + key)
        require(math.isclose(float(rows[0]["effective_value"]), value, rel_tol=1e-15, abs_tol=0) and
                rows[0]["defaulted"] == ("true" if key in ("safety", "maxGrowth") else "false") and
                rows[0]["file"].endswith("run.dict"), "adaptive value/default/source provenance: " + key)
    return {"scope": "input and actual provenance; no timestep-policy accuracy claim",
            "defaulted": str(defaults), "explicit_matching_retry": str(combined), "negative_input": failures,
            "accepted_time": manifest["accepted_time"], "accepted_steps": manifest["accepted_step"],
            "effective_values": expected}


def dynamic_regrid_contract(task, invoke, helpers):
    """One tiny real App transaction; run only after D5 receipt output is wired."""
    directory = task / "dynamic_regrid_two_steps"
    shutil.copytree(task / "hydro_rayleigh_taylor/system", directory / "system")
    mesh = directory / "system/mesh.dict"
    edit_once(mesh, r"\bresolution\s*\(\s*\d+\s+\d+\s*\)", "resolution (4 16)")
    edit_once(mesh, r"\binitialMesh\s+interfaceBand\s*;", "initialMesh uniform;")
    edit_once(mesh, r"\bRefinement\s*\{[^{}]*\}", "")
    edit_once(mesh, r"(\bMesh\s*\{)",
              r"\1 Regrid { minLevel 1; maxLevel 2; interfaceBandLayers 1; vorticityTolerance 0.5; coarsenRatio 0.25; }")
    dt = input_number(directory / "system/run.dict", "dt")
    checked = invoke(directory, True)
    require(checked.returncode == 0 and "HYDRO_INPUT_OK" in checked.stdout and
            "[INIT]" not in checked.stdout and not (directory / "output").exists(),
            "dynamic policy is read before mesh/output")
    result = invoke(directory, False)
    require(result.returncode == 0 and re.search(r"\[COMPLETED\]\s+steps=2\b", result.stdout),
            "real two-step dynamic VOF transaction failed: " + str(directory / "run.log"))
    output = directory / "output"
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest["status"] == "COMPLETED" and manifest["accepted_step"] == 2,
            "dynamic transaction has two actual accepts")
    helpers["equal_float"](float(manifest["accepted_time"]), 2 * dt, "dynamic final clock advances once per physical step")
    _, frames = helpers["csv_rows"](output / "meta/frames.csv")
    require([int(row["step"]) for row in frames] == [0, 1, 2], "dynamic frames publish each physical accept once")
    leaves, minimum_widths = [], []
    for index, frame in enumerate(frames):
        helpers["equal_float"](float(frame["time"]), index * dt, "same-End regrid does not add physical time")
        text = (output / frame["volume"]).read_text(encoding="utf-8")
        matches = list(re.finditer(r"(?m)^\s*AUXDATA\s+\w+\s*=.*$", text))
        require(matches, "dynamic actual field provenance")
        header, tokens = text[:matches[-1].end()], text[matches[-1].end():].split()
        nodes, cells = int(helpers["attribute"](header, "NODES")), int(helpers["attribute"](header, "ELEMENTS"))
        require(helpers["variables"](header) == ["X", "Y"] + helpers["VOLUME_FIELDS"], "dynamic native field variables")
        count = 2 * nodes + len(helpers["VOLUME_FIELDS"]) * cells
        require(len(tokens) == count + 4 * cells, "dynamic actual field/connectivity dimensions")
        connection = tokens[count:]
        require(all(re.fullmatch(r"\d+", value) and 1 <= int(value) <= nodes for value in connection),
                "dynamic new native connectivity bounds")
        x, y = list(map(float, tokens[:nodes])), list(map(float, tokens[nodes:2 * nodes]))
        widths = []
        for cell in range(cells):
            vertices = [int(value) - 1 for value in connection[4 * cell:4 * cell + 4]]
            widths.extend((max(x[v] for v in vertices) - min(x[v] for v in vertices),
                           max(y[v] for v in vertices) - min(y[v] for v in vertices)))
        require(min(widths) > 0 and max(widths) / min(widths) <= 2,
                "actual minLevel1/maxLevel2 leaves remain 2:1")
        leaves.append(cells)
        minimum_widths.append(min(widths))
        metadata = helpers["auxiliary"](header)
        require(int(metadata["AcceptedStep"]) == index and int(metadata["GeometryIdentity"]) == int(frame["geometry_identity"]),
                "dynamic volume binds actual accepted target owner")
    require(leaves[0] == 4 * 16 * 4, "positive minLevel creates the actual uniform initial floor")
    require(leaves[1] > leaves[0] and minimum_widths[1] < minimum_widths[0],
            "actual first physical End refines beyond the initial floor using maxLevel allocation")
    history = helpers["table_history"](output / "series/conservation.dat", "conservation")
    require([int(float(row["tokens"]["step"])) for row in history] == [0, 1, 2],
            "same-End regrid does not duplicate conservation rows")
    for index, row in enumerate(history):
        helpers["equal_float"](float(row["tokens"]["time"]), index * dt, "dynamic actual accepted time")
        helpers["equal_float"](float(row["tokens"]["dt"]), dt if index else 0, "dynamic actual dt")
    require((output / "series/regrid_decision.dat").is_file() and
            (output / "series/regrid.dat").is_file(), "D5 publishes actual decisions and changed-tree remap ledger")
    _, trace = helpers["csv_rows"](output / "meta/configuration_trace.csv")
    regrid = [row for row in trace if row["field"].startswith("Mesh.Regrid.")]
    require(len(regrid) == 5 and all(row["used"] == "true" and
            row["consumer"] == "Hydrodynamics.VOF.interfaceFlowRegrid" for row in regrid),
            "one consumed owner for all five dynamic policy values")
    require(not list(output.rglob("*.vtu")) and not list(output.rglob("*.pvd")) and
            not list((output / "series").rglob("*.csv")), "dynamic run preserves DAT-only output")
    return {"directory": str(directory), "accepted_steps": 2, "accepted_time": 2 * dt,
            "leaves": leaves, "minimum_widths": minimum_widths, "initial_min_level": 1, "dynamic_max_level": 2}


def initial_band_level3_contract(task, invoke, helpers):
    """Independent deep initial-tree check; does not weaken static-L1 helpers."""
    directory = task / "initial_band_level3"
    shutil.copytree(task / "hydro_rayleigh_taylor/system", directory / "system")
    (directory / "system/mesh.dict").write_text(
        "Mesh { resolution (4 16); initialMesh interfaceBand; Refinement { level 3; bandCells 2; } }\n", encoding="utf-8")
    dt = input_number(directory / "system/run.dict", "dt")
    checked = invoke(directory, True)
    require(checked.returncode == 0 and "HYDRO_INPUT_OK" in checked.stdout and
            not (directory / "output").exists(), "deep initial-band schema accepted before compute")
    process = invoke(directory, False)
    require(process.returncode == 0 and re.search(r"\[COMPLETED\]\s+steps=2\b", process.stdout),
            "deep initial-band two physical accepts; see " + str(directory / "run.log"))
    output = directory / "output"
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest["status"] == "COMPLETED" and manifest["accepted_step"] == 2 and
            manifest["complete_frames"] == 3, "deep initial-band actual completed prefix")
    helpers["equal_float"](float(manifest["accepted_time"]), 2*dt, "deep initial-band true final clock")
    _, frames = helpers["csv_rows"](output / "meta/frames.csv")
    require([int(x["step"]) for x in frames] == [0,1,2], "deep initial-band unique actual frames")
    saved_rectangles = None;leaf_counts = [];shared_counts = []
    for index, frame in enumerate(frames):
        path = output / frame["volume"];text = path.read_text(encoding="utf-8")
        tags = list(re.finditer(r"(?m)^\s*AUXDATA\s+\w+\s*=.*$", text));require(tags, "deep initial-field provenance")
        header,tokens = text[:tags[-1].end()],text[tags[-1].end():].split()
        nodes,cells = int(helpers["attribute"](header,"NODES")),int(helpers["attribute"](header,"ELEMENTS"))
        names = helpers["variables"](header)
        require(names == ["X","Y"]+helpers["VOLUME_FIELDS"], "deep initial-field exact native schema")
        count = 2*nodes+len(helpers["VOLUME_FIELDS"])*cells
        require(len(tokens) == count+4*cells and all(helpers["number"](x) is not None for x in tokens[:count]),
                "deep initial complete finite native payload")
        arrays = {};offset = 0
        for name in names:
            size = nodes if name in ("X","Y") else cells
            arrays[name] = [float(x) for x in tokens[offset:offset+size]];offset += size
        connection = [int(x)-1 for x in tokens[count:]]
        require(all(0<=x<nodes for x in connection), "deep initial native connectivity bounds")
        rectangles = [];levels = [];vertical = {};horizontal = {}
        for cell in range(cells):
            vertices = connection[4*cell:4*cell+4]
            xs = sorted({arrays["X"][v] for v in vertices});ys = sorted({arrays["Y"][v] for v in vertices})
            require(len(xs)==len(ys)==2 and len(set(vertices))==4, "deep initial actual rectangular native cell")
            x0,x1 = xs;y0,y1 = ys;dx,dy = x1-x0,y1-y0
            require(0<=x0<x1<=1 and 0<=y0<y1<=4 and dx==dy, "deep initial physical RT domain and square leaf")
            level = round(math.log2(.25/dx))
            require(0<=level<=3 and dx==math.ldexp(.25,-level), "deep initial true dyadic leaf level")
            rectangles.append((x0,x1,y0,y1));levels.append(level)
            # Independent native-face interval overlays include real hanging
            # subfaces; global max/min width is deliberately NOT a 2:1 test.
            vertical.setdefault(x0,[[],[]])[0].append((y0,y1,level,cell))
            vertical.setdefault(x1,[[],[]])[1].append((y0,y1,level,cell))
            horizontal.setdefault(y0,[[],[]])[0].append((x0,x1,level,cell))
            horizontal.setdefault(y1,[[],[]])[1].append((x0,x1,level,cell))
        require(set(levels)=={0,1,2,3}, "deep initial tree really contains graded0..3, not a uniform floor workaround")
        require(max(x[1]-x[0] for x in rectangles)==8*min(x[1]-x[0] for x in rectangles), "deep initial far coarse8h")
        shared = 0;hanging = 0
        for edges in (vertical,horizontal):
            for sides in edges.values():
                left,right = sorted(sides[0]),sorted(sides[1]);i = j = 0
                while i<len(left) and j<len(right):
                    a,b = left[i],right[j]
                    if min(a[1],b[1])>max(a[0],b[0]):
                        require(a[3]!=b[3] and abs(a[2]-b[2])<=1, "deep initial actual shared-subface local2:1")
                        shared += 1;hanging += a[2]!=b[2]
                    if a[1]<=b[1]:i += 1
                    if b[1]<=a[1]:j += 1
        require(shared>0 and hanging>0, "deep initial actual hanging-face evidence")
        helpers["equal_float"](math.fsum((x1-x0)*(y1-y0) for x0,x1,y0,y1 in rectangles),4.,"deep initial complete domain measure")
        require(saved_rectangles is None or rectangles==saved_rectangles,"static deep initial mesh remains fixed through physical steps")
        saved_rectangles = rectangles;leaf_counts.append(cells);shared_counts.append(shared)
        if index==0:
            actual_mixed = [i for i,(l,g) in enumerate(zip(arrays["liquid_volume"],arrays["gas_volume"])) if l>0 and g>0]
            require(actual_mixed and all(levels[i]==3 for i in actual_mixed), "deep initial real mixed interface lies on finest cells")
        meta = helpers["auxiliary"](header)
        require(int(meta["AcceptedStep"])==index and int(meta["GeometryIdentity"])==int(frame["geometry_identity"]),"deep initial actual frame owner")
        helpers["equal_float"](float(frame["time"]),index*dt,"deep initial physical clock")
    history = helpers["table_history"](output/"series/conservation.dat","conservation")
    require([int(float(x["tokens"]["step"])) for x in history]==[0,1,2],"deep initial one conservation row per real accept")
    require(not (output/"series/regrid_decision.dat").exists() and not (output/"series/regrid.dat").exists(),"initial band is not a fake dynamic remap event")
    require(not list(output.rglob("*.vtu")) and not list(output.rglob("*.pvd")),"deep initial remains DAT-only")
    return {"directory":str(directory),"accepted_steps":2,"accepted_time":2*dt,"leaves":leaf_counts,
            "levels":[0,1,2,3],"shared_subfaces":shared_counts,"initial_mixed_finest":True}


def main():
    require(len(sys.argv) == 4, "usage: hydro_app_contract.py LAMP.exe projectRoot workRoot")
    executable, project, work = (Path(arg).resolve() for arg in sys.argv[1:])
    require(executable.is_file() and work.drive.upper() == "D:", "formal executable and D: workRoot required")
    sys.dont_write_bytecode = True
    # Reuse only pure parsers. A distinct run_name never invokes that contract's main.
    helpers = runpy.run_path(str(project / "tests/contracts/hydrodynamics/hydro_tecplot_output_contract.py"),
                             run_name="_hydro_app_dat_parsers")
    task = work / ("hydro_app_vof_" + uuid.uuid4().hex[:12])
    task.mkdir(parents=True, exist_ok=False)
    temporary = task / "_temp"
    temporary.mkdir()
    environment = os.environ.copy()
    environment.update(TEMP=str(temporary), TMP=str(temporary), OMP_NUM_THREADS="2", MKL_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    originals = {}
    positives = []

    def invoke(directory, check):
        args = [str(executable), "--threads", "2"] + (["--check-case"] if check else [])
        result = subprocess.run(args, cwd=directory, env=environment, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", timeout=120)
        (directory / ("check.log" if check else "run.log")).write_text(result.stdout, encoding="utf-8")
        return result

    def unchanged():
        for directory, hashes in originals.items():
            require({name: digest(directory / name) for name in FILES} == hashes,
                    "production dictionaries were modified: " + str(directory))

    try:
        for case, shape in (("hydro_static_drop", (8, 8)), ("hydro_rayleigh_taylor", (8, 32))):
            original = project / "cases" / case / "system"
            require({path.name for path in original.iterdir()} == set(FILES), "actual case has six dictionaries")
            originals[original] = {name: digest(original / name) for name in FILES}
            directory = task / case
            shutil.copytree(original, directory / "system")
            # This fixture explicitly verifies fixed mode even when production
            # templates select Adaptive; copied-input edits never touch sources.
            run_path = directory / "system/run.dict"
            fixed_text, adaptive_blocks = re.subn(r"\bAdaptive\s*\{[^{}]*\}", "", run_path.read_text(encoding="utf-8"))
            require(adaptive_blocks <= 1, "at most one copied adaptive policy")
            run_path.write_text(fixed_text, encoding="utf-8")
            require("Initialization" not in (directory / "system/interface.dict").read_text(encoding="utf-8"),
                    "case must exercise the absent initialization block default")
            dt = input_number(directory / "system/run.dict", "dt")
            end_time = 2 * dt
            edit_once(directory / "system/mesh.dict", r"\bresolution\s*\(\s*\d+\s+\d+\s*\)",
                      f"resolution ({shape[0]} {shape[1]})")
            edit_once(directory / "system/run.dict", r"\bendTime\s+[^;]+;", f"endTime {end_time:.17g};")
            edit_once(directory / "system/run.dict", r"\bwriteEverySteps\s+\d+;", "writeEverySteps 1;")
            checked = invoke(directory, True)
            require(checked.returncode == 0 and "HYDRO_INPUT_OK" in checked.stdout,
                    "real six-file VOF CheckOnly failed: " + str(directory / "check.log"))
            require(not (directory / "output").exists() and "[INIT]" not in checked.stdout,
                    "CheckOnly must not create mesh/numerical output")
            result = invoke(directory, False)
            require(result.returncode == 0 and re.search(r"\[COMPLETED\]\s+steps=2\b", result.stdout),
                    "real two-step VOF run failed: " + str(directory / "run.log"))
            positives.append(output_contract(directory, dt, end_time, shape[0] * shape[1], helpers))
            unchanged()
        negative = []
        for kind in ("sevof", "hermiteGals", "unknown_keyword", "obsolete_fields"):
            directory = task / ("reject_" + kind)
            shutil.copytree(task / "hydro_rayleigh_taylor/system", directory / "system")
            if kind == "obsolete_fields":
                (directory / "system/fields.dict").write_text("Fields { Constants {} }\n", encoding="utf-8")
            elif kind == "unknown_keyword":
                edit_once(directory / "system/interface.dict", r"(\bInterface\s*\{)",
                          r"\1 unexpectedContractField 1;")
            else:
                edit_once(directory / "system/interface.dict", r"\bmethod\s+vof\s*;", f"method {kind};")
            result = invoke(directory, True)
            source_pattern = r"fields\.dict:\d+[;:]" if kind == "obsolete_fields" else r"interface\.dict:\d+:"
            require(result.returncode != 0 and re.search(source_pattern, result.stdout),
                    "invalid method/keyword/obsolete document must be rejected at its real input location")
            require("HYDRO_INPUT_OK" not in result.stdout and "[INIT]" not in result.stdout and
                    "[COMPLETED]" not in result.stdout and not (directory / "output").exists(),
                    "rejected input cannot start computation")
            if kind == "obsolete_fields":
                require("unused" in result.stdout and "fields.dict" in result.stdout, "obsolete document reported explicitly")
            elif kind == "unknown_keyword":
                require("unexpectedContractField" in result.stdout, "unknown key reported explicitly")
            else:
                require(kind.lower() in result.stdout.lower() and "vof" in result.stdout.lower(),
                        "removed production method receives explicit VOF-only input error")
            negative.append({"case": kind, "directory": str(directory), "rejected_before_compute": True})
        adaptive_inputs = adaptive_input_contract(task, invoke, helpers)
        dynamic_regrid = dynamic_regrid_contract(task, invoke, helpers)
        initial_band_level3 = initial_band_level3_contract(task, invoke, helpers)
        unchanged()
        record = {"status": "PASS", "scope": "VOF App wiring, not archived FirstJet gates or physical accuracy",
                  "positives": positives, "negative_input": negative, "adaptive_input": adaptive_inputs, "dynamic_regrid": dynamic_regrid, "initial_band_level3": initial_band_level3, "production_case_sha_unchanged": True,
                  "original_case_sha256": {str(path): hashes for path, hashes in originals.items()}}
        (task / "result.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        print("HYDRO_APP_CONTRACT_PASS: " + str(task))
    finally:
        unchanged()


if __name__ == "__main__":
    main()
