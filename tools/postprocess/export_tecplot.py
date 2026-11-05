#!/usr/bin/env python3
"""Export saved Hydro ASCII VTU fields and CSV curves without changing sources.

Only saved values are serialized: no solver, mesh interpolation or reconstruction.
A run is selected directly; --retrofit explicitly enables recursive run discovery.
"""
import argparse
import codecs
import csv
import io
import json
from decimal import Decimal, InvalidOperation
import os
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?\Z")
INTEGER = re.compile(r"[+-]?\d+\Z")
VTK_TYPES = {"Float32", "Float64", "Int8", "UInt8", "Int16", "UInt16",
             "Int32", "UInt32", "Int64", "UInt64"}
MISSING = {"", "na", "nan", "+nan", "-nan",
           "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"}
CATEGORIES = {"stage", "status", "module", "event", "reason", "rejected_reason", "kind",
              "phase", "species", "side", "boundary", "boundary_name", "boundary_id",
              "species_id", "phase_id", "target", "type", "pol_solver_mode", "top1_name",
              "top2_name", "top3_name", "dominant_source", "dominant_boundary",
              "dominant_reaction_channel", "curvature_measure", "force_discretization", "area_update"}


def missing(token):
    key = token.strip().lower()
    return key in MISSING or bool(re.fullmatch(r"[+-]?nan\([^)]*\)", key))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def quote(value):
    """One Tecplot quoted string, including escaped source metadata."""
    value = str(value)
    value = "".join({"\n": r"\n", "\r": r"\r", "\t": r"\t"}.get(c,
                    c if ord(c) >= 32 and ord(c) != 127 else "\\x" + format(ord(c), "02x")) for c in value)
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def scalar(token, boolean=False):
    """Return a finite Tecplot token without converting through binary64."""
    token = token.strip()
    lower = token.lower()
    if boolean and lower in ("true", "false"):
        return "1" if lower == "true" else "0"
    if missing(token):
        return None
    if not NUMBER.fullmatch(token):
        return None
    try:
        value = Decimal(token.replace("d", "e").replace("D", "E"))
    except InvalidOperation:
        return None
    return token.replace("d", "e").replace("D", "E") if value.is_finite() else None


def integer(token, description):
    require(bool(INTEGER.fullmatch(token)), description + " must be an integer")
    return int(token)


def write_atomic(path, content, overwrite=False):
    if path.exists() and not overwrite:
        return "existing"
    if path.exists() and overwrite:
        require(not path.is_symlink(), "cannot overwrite a symlink")
        expected = [line for line in content.splitlines() if line.startswith("DATASETAUXDATA Generator=")
                    or line.startswith("DATASETAUXDATA SourceCsv=") or line.startswith("DATASETAUXDATA SourceVTU=")]
        old_header = path.read_text(encoding="utf-8").split("ZONE", 1)[0].splitlines()
        require(len(expected) == 2 and all(line in old_header for line in expected),
                "existing DAT is not a matching generated companion: " + str(path))
    temporary = path.with_name(path.name + ".tmp." + str(os.getpid()))
    require(not temporary.exists(), "temporary export already exists: " + str(temporary))
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as output:
            output.write(content)
        if path.exists() and not overwrite:
            temporary.unlink()
            return "existing"
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return "written"


def array_tokens(node, expected=None):
    require(node.tag == "DataArray" and node.get("format", "ascii") == "ascii",
            "only inline ASCII DataArray is supported")
    kind = node.get("type")
    require(kind in VTK_TYPES, "unsupported VTK array type: " + str(kind))
    require(not list(node) and not any(k in node.attrib for k in ("offset", "RangeMinOffset")),
            "appended/nested VTK arrays are unsupported")
    components = integer(node.get("NumberOfComponents", "1"), "component count")
    require(components in (1, 3), "only scalar and three-component arrays are supported")
    tokens = (node.text or "").split()
    if expected is not None:
        require(len(tokens) == expected * components,
                "array count mismatch: " + str(node.get("Name")))
        if node.get("NumberOfTuples") is not None:
            require(integer(node.get("NumberOfTuples"), "tuple count") == expected,
                    "declared tuple count mismatch")
    else:
        require(len(tokens) % components == 0, "incomplete VTK tuple")
        if node.get("NumberOfTuples") is not None:
            require(integer(node.get("NumberOfTuples"), "tuple count") * components == len(tokens),
                    "declared tuple count mismatch")
    for token in tokens:
        require(scalar(token) is not None or (missing(token) and token.lower() not in {"", "na"}),
                "nonnumeric VTK array token: " + token)
        if kind.startswith(("Int", "UInt")):
            number = integer(token, "VTK integer array token")
            bits = int(re.search(r"\d+", kind).group())
            low, high = (0, 2**bits - 1) if kind.startswith("UInt") else (-2**(bits-1), 2**(bits-1)-1)
            require(low <= number <= high, "VTK integer array token outside declared range")
    return components, tokens


def blocks(values):
    return "".join(" ".join(values[i:i+6]) + "\n" for i in range(0, len(values), 6))


def coordinate_names(run):
    for path in (run / "output/meta/effective/geometry.dict", run / "system/geometry.dict"):
        if path.is_file():
            body = re.sub(r"//[^\n]*|/\*.*?\*/", "", path.read_text(encoding="utf-8-sig"), flags=re.S)
            match = re.search(r"\bcoordinates\s+(xy|rz)\s*;", body)
            if match:
                return ("R", "Z") if match.group(1) == "rz" else ("X", "Y")
    return ("X", "Y")


def vtu_text(source, coordinates=("X", "Y")):
    raw = source.read_bytes()
    require(b"<!DOCTYPE" not in raw.upper() and b"<!ENTITY" not in raw.upper(),
            "VTK declarations/entities are unsupported")
    document = ET.fromstring(raw)
    require(document.tag == "VTKFile" and document.get("type") == "UnstructuredGrid",
            "requires VTK UnstructuredGrid")
    require(not document.findall(".//AppendedData"), "appended VTK data is unsupported")
    grids = document.findall("UnstructuredGrid")
    require(len(grids) == 1, "requires exactly one UnstructuredGrid")
    grid = grids[0]
    pieces = grid.findall("Piece")
    require(len(pieces) == 1, "requires one complete unstructured Piece")
    piece = pieces[0]
    require({n.tag for n in piece} <= {"Points", "Cells", "PointData", "CellData", "FieldData"},
            "unsupported VTK Piece section")
    nodes = integer(piece.get("NumberOfPoints", ""), "point count")
    cells = integer(piece.get("NumberOfCells", ""), "cell count")
    require(nodes > 0 and cells > 0, "empty VTK mesh")
    require(len(piece.findall("Points")) == 1, "requires one Points section")
    points = piece.findall("Points/DataArray")
    require(len(points) == 1, "requires one Points array")
    components, values = array_tokens(points[0], nodes)
    require(components == 3 and all(scalar(x) is not None for x in values),
            "requires finite three-component VTK points")
    require(all(Decimal(scalar(values[i+2])) == 0 for i in range(0, len(values), 3)),
            "only planar two-coordinate fields are supported; nonzero third coordinate")
    xy = [[scalar(values[3*i+j]) for i in range(nodes)] for j in range(2)]
    sections = piece.findall("Cells")
    require(len(sections) == 1, "requires one Cells section")
    arrays = sections[0].findall("DataArray")
    connectivity_arrays = {n.get("Name"): n for n in arrays}
    require(len(arrays) == 3 and set(connectivity_arrays) == {"connectivity", "offsets", "types"},
            "requires unique connectivity/offsets/types arrays")
    _, connections = array_tokens(connectivity_arrays["connectivity"])
    _, offsets = array_tokens(connectivity_arrays["offsets"], cells)
    _, types = array_tokens(connectivity_arrays["types"], cells)
    require(all(n.get("NumberOfComponents", "1") == "1" for n in arrays),
            "cell topology arrays must be scalar")
    connectivity = [integer(x, "connectivity") for x in connections]
    groups = {3: [], 9: []}
    start = 0
    for index, (offset, cell_type) in enumerate(zip(offsets, types)):
        end, kind = integer(offset, "offset"), integer(cell_type, "cell type")
        require(kind in groups, "only VTK_LINE=3 and VTK_QUAD=9 are supported")
        require(end > start and end <= len(connectivity), "invalid connectivity offset")
        item = connectivity[start:end]
        require(len(item) == (2 if kind == 3 else 4) and len(set(item)) == len(item),
                "invalid line/quad connectivity arity or repeated node")
        require(all(0 <= x < nodes for x in item), "connectivity index outside source points")
        groups[kind].append((index, item))
        start = end
    require(start == len(connectivity), "unconsumed connectivity values")
    fields = []
    names = set(coordinates)
    for section in ("Points", "Cells", "PointData", "CellData", "FieldData"):
        for parent in piece.findall(section):
            require(all(child.tag == "DataArray" for child in parent), "non-DataArray child in " + section)
    for section, expected, centered in (("PointData", nodes, False), ("CellData", cells, True)):
        require(len(piece.findall(section)) <= 1, "duplicate " + section)
        for node in piece.findall(section + "/DataArray"):
            base = node.get("Name", "")
            require(base, "unnamed source field")
            count, tokens = array_tokens(node, expected)
            for component in range(count):
                name = base if count == 1 else base + "_" + "XYZ"[component]
                require(name not in names, "duplicate exported variable: " + name)
                names.add(name)
                original = tokens[component::count]
                finite = [scalar(token) for token in original]
                fields.append((name, centered, [token if token is not None else "0" for token in finite]))
                if any(token is None for token in finite):
                    flag = name + "Available"
                    require(flag not in names, "availability variable collides with source: " + flag)
                    names.add(flag)
                    fields.append((flag, centered, ["1" if token is not None else "0" for token in finite]))
    metadata = {}
    for node in grid.findall("FieldData/DataArray") + piece.findall("FieldData/DataArray"):
        name = node.get("Name", "")
        require(name and name not in metadata, "empty/duplicate FieldData name")
        _, tokens = array_tokens(node)
        metadata[name] = tokens
    require(len(grid.findall("FieldData")) <= 1 and len(piece.findall("FieldData")) <= 1,
            "duplicate FieldData section")
    variables = list(coordinates) + [field[0] for field in fields]
    output = ["TITLE=" + quote("Saved LAMP field; no interpolation") + "\nVARIABLES=\n"]
    output.extend(quote(name) + "\n" for name in variables)
    output.extend(('DATASETAUXDATA Generator="LAMP_AsciiVtu_v1"\n',
                   "DATASETAUXDATA SourceVTU=" + quote(source.resolve().as_posix()) + "\n"))
    aux = {"SourceVTU": str(source), "CoordinateNames": "/".join(coordinates),
           "MissingValues": "0 is a serialization placeholder only when the matching Available variable is 0"}
    for name, tokens in metadata.items():
        key = "VTK_" + re.sub(r"[^A-Za-z0-9_]", "_", name)
        require(key not in aux, "FieldData AUX name collision")
        aux[key] = " ".join(tokens)
    aliases = {"accepted_step": "AcceptedStep", "geometry_identity": "GeometryIdentity"}
    for source_name, alias in aliases.items():
        if source_name in metadata:
            require(len(metadata[source_name]) == 1, source_name + " must be scalar")
            aux[alias] = metadata[source_name][0]
    clock = metadata.get("TimeValue")
    if clock is not None:
        require(len(clock) == 1 and scalar(clock[0]) is not None, "TimeValue must be finite and scalar")
    for kind, group in groups.items():
        if not group:
            continue
        output.append("ZONE T=" + quote(source.stem + (" lines" if kind == 3 else " quads")) + "\n")
        output.append("NODES=" + str(nodes) + ", ELEMENTS=" + str(len(group))
                      + ", ZONETYPE=" + ("FELINESEG" if kind == 3 else "FEQUADRILATERAL")
                      + ", DATAPACKING=BLOCK\n")
        if clock is not None:
            output.append("SOLUTIONTIME=" + scalar(clock[0]) + ", STRANDID=1\n")
        output.append("DT=(" + " ".join(["DOUBLE"] * len(variables)) + ")\n")
        locations = [str(i+3) for i, field in enumerate(fields) if field[1]]
        if locations:
            output.append("VARLOCATION=([" + ",".join(locations) + "]=CELLCENTERED)\n")
        output.extend("AUXDATA " + key + "=" + quote(value) + "\n" for key, value in aux.items())
        output.extend(blocks(values) for values in xy)
        for _, centered, values in fields:
            output.append(blocks([values[index] for index, _ in group] if centered else values))
        output.extend(" ".join(str(index+1) for index in item) + "\n" for _, item in group)
    return "".join(output)


class IncompleteCsv(ValueError):
    """The single live CSV snapshot ends inside an uncommitted record."""


def csv_encoding(value):
    """Validate an explicit text codec without guessing a file's encoding."""
    try:
        codec = codecs.lookup(value)
        if not getattr(codec, "_is_text_encoding", False):
            raise LookupError("not a text codec")
        name = codec.name
        b"".decode(name, errors="strict")
    except (LookupError, UnicodeError) as error:
        raise argparse.ArgumentTypeError("invalid CSV text encoding: " + str(value)) from error
    return name


def csv_text(source, encoding="utf-8-sig"):
    encoding = csv_encoding(encoding)
    snapshot = source.read_text(encoding=encoding, errors="strict")
    if not snapshot:
        return None
    if not snapshot.endswith("\n"):
        raise IncompleteCsv("snapshot has no complete final newline: " + str(source))
    reader = csv.reader(io.StringIO(snapshot, newline=""), strict=True)
    try:
        columns = next(reader)
    except StopIteration:
        return None
    columns = [name.strip() for name in columns]
    require(columns and all(columns) and len(set(columns)) == len(columns), "empty/duplicate CSV header")
    rows = list(reader)
    require(all(len(row) == len(columns) for row in rows), "CSV record width differs from header")
    if not rows:
        return None
    numeric = []
    categories = []
    for index in range(len(columns)):
        values = [row[index].strip() for row in rows]
        has_text = any(not missing(value) and scalar(value, True) is None for value in values)
        has_numeric = any(scalar(value, True) is not None for value in values)
        require(not (has_text and has_numeric) or columns[index].lower() in CATEGORIES,
                "CSV column mixes numeric and categorical values: " + columns[index])
        if has_text or columns[index].lower() in CATEGORIES:
            categories.append(index)
        elif any(scalar(value, True) is not None for value in values):
            numeric.append(index)
    lower = {name.lower(): i for i, name in enumerate(columns)}
    time = next((lower[name] for name in ("time", "accepted_time", "time_s", "time_code", "t_si_s",
                                          "time_si_s", "physical_time", "t") if name in lower), None)
    step = next((lower[name] for name in ("step", "accepted_step", "timestep", "attempted_step")
                 if name in lower), None)
    axis = time if time is not None else step
    require(all(index is None or index not in categories for index in (axis, step)), "CSV clock column contains text")
    axis_name = columns[axis] if axis is not None else "Record"
    third = "Step" if step is not None else "Record"
    if third.lower() == axis_name.lower():
        third = "Source" + third
    metrics = [index for index in numeric if index not in (axis, step)]
    zones = []
    for metric in metrics:
        segments = {}
        previous = {}
        for record, row in enumerate(rows, 1):
            classification = tuple(row[i].strip() for i in categories)
            value = scalar(row[metric], True)
            x = scalar(row[axis], True) if axis is not None else str(record)
            z = scalar(row[step], True) if step is not None else str(record)
            active = segments.setdefault(classification, [])
            if value is None or x is None or z is None:
                previous.pop(classification, None)
                continue
            current = previous.get(classification)
            if current is None or Decimal(x) <= Decimal(current[-1][0]):
                current = []
                active.append(current)
                previous[classification] = current
            current.append((x, value, z))
        for classification, parts in segments.items():
            for part_number, values in enumerate(parts, 1):
                if not values:
                    continue
                labels = [columns[index] + "=" + value for index, value in zip(categories, classification)]
                title = columns[metric] + (" | " + ", ".join(labels) if labels else "") + " | segment " + str(part_number)
                zones.append((title[:128], columns[metric], classification, values))
    if not zones:
        return None
    source_name = source.resolve().as_posix()
    meaning = ("source logical data record (1-based); no physical clock column" if axis is None else
               "actual source step column; no physical time column" if time is None else
               "actual source time column; no resampling")
    output = ["VARIABLES=" + " ".join(quote(name) for name in (axis_name, "Value", third)) + "\n",
              'DATASETAUXDATA Generator="LAMP_CsvHistory_v1"\n',
              "DATASETAUXDATA SourceCsv=" + quote(source_name) + "\n"]
    if encoding != "utf-8-sig":
        output.append("DATASETAUXDATA CsvEncoding=" + quote(encoding) + "\n")
    for title, metric, classification, values in zones:
        output.append("ZONE T=" + quote(title) + ", I=" + str(len(values)) + ", ZONETYPE=ORDERED, DATAPACKING=POINT\nDT=(DOUBLE DOUBLE DOUBLE)\n")
        aux = {"SourceCsv": source_name, "Column": metric, "AxisColumn": axis_name,
               "XMeaning": meaning, "StepColumn": columns[step] if step is not None else "Record",
               "MissingPolicy": "NA, empty and nonfinite omitted; gaps and nonincreasing X start new zones",
               "BooleanValues": "source true/false flags are serialized as 1/0"}
        if encoding != "utf-8-sig":
            aux["CsvEncoding"] = encoding
        if step is None:
            aux["RecordMeaning"] = "source logical data record (1-based); not a physical step"
        for number, (index, category) in enumerate(zip(categories, classification)):
            aux["Category" + str(number) + "Name"] = columns[index]
            aux["Category" + str(number) + "Value"] = category
        output.extend("AUXDATA " + key + "=" + quote(value) + "\n" for key, value in aux.items())
        output.extend(" ".join(value) + "\n" for value in values)
    return "".join(output)


def selected_run(run):
    if any(part.lower() == "benchmarks" for part in run.parts):
        return False
    config = run / "system/run.dict"
    if not config.is_file() or not (run / "output").is_dir():
        return False
    body = re.sub(r"//[^\n]*|/\*.*?\*/", "", config.read_text(encoding="utf-8-sig"), flags=re.S)
    return bool(re.search(r"\bapp\s+Hydrodynamics\s*;", body))


def series_run_identity(run):
    """Require a launcher run identity and matching standalone App metadata."""
    if any(part.lower() == "benchmarks" for part in run.parts):
        return False, "benchmark path is excluded"
    if not (run / "output/series").is_dir():
        return False, "no saved output/series directory"
    try:
        configs = [path for path in (run / "system/run.dict", run / "output/meta/effective/run.dict") if path.is_file()]
        if not configs:
            return False, "no original or effective Run.app input"
        declared = set()
        for config in configs:
            body = re.sub(r"//[^\n]*|/\*.*?\*/", "", config.read_text(encoding="utf-8-sig"), flags=re.S)
            apps = re.findall(r"\bapp\s+([A-Za-z][A-Za-z0-9_]*)\s*;", body)
            if len(apps) != 1:
                return False, "missing or ambiguous Run.app"
            declared.add(apps[0])
        if len(declared) != 1:
            return False, "original/effective Run.app identities conflict"
        app = next(iter(declared))
        if app not in {"Hydrodynamics", "Interface"}:
            return False, "App " + app + " is excluded (only non-electric Hydro/Interface series)"
        if any((run / relative).is_file() for relative in ("system/electric.dict", "output/meta/effective/electric.dict")):
            return False, "electric coupling dictionary is present"
        record_path = run / "run.json"
        manifest_path = run / "output/meta/manifest.json"
        if not record_path.is_file() or not manifest_path.is_file():
            return False, "missing launcher run.json or output/meta/manifest.json identity"
        record = json.loads(record_path.read_text(encoding="utf-8-sig"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        if not isinstance(record, dict) or not isinstance(manifest, dict):
            return False, "run/manifest identity must be an object"
        if not isinstance(record.get("run"), str) or Path(record["run"]).resolve() != run.resolve():
            return False, "launcher run path is absent or differs from this directory"
        if manifest.get("app") != app or record.get("app", app) != app:
            return False, "Run.app and run/manifest App identities conflict"
        if str(record.get("kind", "")).lower() == "benchmark":
            return False, "launcher identifies a benchmark run"
        if isinstance(record.get("case"), str) and any(part.lower() == "benchmarks" for part in Path(record["case"]).parts):
            return False, "source case is under benchmarks"
    except (OSError, ValueError) as error:
        return False, "unreadable or malformed run identity: " + str(error)
    return True, app


def series_run_candidates(root):
    """Discover run-level output directories without scanning field payloads."""
    found = set()
    def inaccessible(error):
        print("SKIP unreadable subtree: " + str(error))
    for directory, children, _ in os.walk(root, followlinks=False, onerror=inaccessible):
        current = Path(directory)
        children[:] = [name for name in children if name.lower() not in
                       {"benchmarks", ".git", "runtime", "fields", "build", "cache", "_temp"}]
        if current.name == "output":
            if (current / "series").is_dir():
                found.add(current.parent)
            children[:] = []
    return sorted(found)


NATIVE_VOLUME_FIELDS = (
    "cell_id", "cell_capacity", "phase_volume_fraction", "liquid_volume", "gas_volume",
    "plic_normal_X", "plic_normal_Y", "plic_normal_Z", "plic_offset", "plic_resolved",
    "velocity_X", "velocity_Y", "velocity_Z", "pressure", "kappa_HF_force", "kappa_HF_forceAvailable")


def validate_native_volume(source, frame):
    """Authenticate one published native DAT without converting or replacing it."""
    text = source.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    auxiliaries = {}
    last_auxiliary = -1
    pattern = re.compile(r'^\s*AUXDATA\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*"((?:[^"\\]|\\.)*)"\s*$')
    for line_number, line in enumerate(lines):
        match = pattern.fullmatch(line)
        if match:
            key, value = match.groups()
            require(key not in auxiliaries, "duplicate native volume AUXDATA: " + key)
            auxiliaries[key] = value
            last_auxiliary = line_number
    require(last_auxiliary >= 0, "native volume requires accepted-state AUXDATA")
    header = "\n".join(lines[:last_auxiliary+1])
    control = re.sub(r'"(?:[^"\\]|\\.)*"', '""', header)
    require(len(re.findall(r"\bZONE\b", control)) == 1, "native volume requires one complete zone")
    variables = re.search(r"\bVARIABLES\s*=\s*(.*?)\bZONE\b", header, re.S)
    require(variables is not None, "native volume lacks VARIABLES header")
    names = re.findall(r'"([^"\\]+)"', variables.group(1))
    require(tuple(names[:2]) in (("X", "Y"), ("R", "Z"))
            and tuple(names[2:]) == NATIVE_VOLUME_FIELDS,
            "native volume must retain XY/RZ and the canonical 16 cell-centered fields")

    def attribute(name):
        values = re.findall(r"\b" + name + r"\s*=\s*([^\s,()]+)", control)
        require(len(values) == 1, "native volume requires one " + name)
        return values[0]

    nodes = integer(attribute("NODES"), "native node count")
    elements = integer(attribute("ELEMENTS"), "native element count")
    require(nodes > 0 and elements > 0, "native volume has empty dimensions")
    require(attribute("ZONETYPE") == "FEQUADRILATERAL" and attribute("DATAPACKING") == "BLOCK",
            "native volume must be FEQUADRILATERAL/BLOCK")
    types = re.findall(r"\bDT\s*=\s*\(([^)]*)\)", control, re.S)
    require(len(types) == 1 and types[0].split() == ["DOUBLE"] * len(names),
            "native volume requires explicit DOUBLE for every variable")
    locations = re.findall(r"\bVARLOCATION\s*=\s*\(\s*\[([^]]+)\]\s*=\s*CELLCENTERED\s*\)", control, re.S)
    require(len(locations) == 1, "native volume requires one cell-centered VARLOCATION")
    centered = []
    for item in locations[0].split(","):
        match = re.fullmatch(r"\s*(\d+)\s*(?:-\s*(\d+)\s*)?", item)
        require(match is not None, "invalid native cell-centered variable index")
        first, last = int(match.group(1)), int(match.group(2) or match.group(1))
        require(1 <= first <= last <= len(names), "native VARLOCATION index outside variables")
        centered.extend(range(first, last+1))
    require(len(centered) == 16 and set(centered) == set(range(3, 19)),
            "native volume fields must all remain cell-centered")
    payload = " ".join(lines[last_auxiliary+1:]).split()
    numeric_count = 2*nodes + 16*elements
    require(len(payload) == numeric_count + 4*elements, "native volume BLOCK/connectivity dimensions differ from header")
    require(all(scalar(token) is not None for token in payload[:numeric_count]),
            "native volume contains nonfinite or nonnumeric field data")
    connectivity = [integer(token, "native quadrilateral connectivity") for token in payload[numeric_count:]]
    require(all(1 <= node <= nodes for node in connectivity), "native connectivity index outside nodes")
    require(all(len(set(connectivity[offset:offset+4])) == 4 for offset in range(0, len(connectivity), 4)),
            "native quadrilateral repeats a node")
    step = auxiliaries.get("AcceptedStep", "")
    geometry = auxiliaries.get("GeometryIdentity", "")
    require(integer(step, "native accepted step") >= 0 and integer(geometry, "native geometry identity") > 0,
            "native volume has invalid accepted provenance")
    time = scalar(attribute("SOLUTIONTIME"))
    require(time is not None, "native volume requires finite SOLUTIONTIME")
    require(all(key in frame for key in ("step", "time", "geometry_identity")),
            "native DAT frame index must include step/time/geometry_identity")
    frame_step = integer(frame["step"], "frame accepted step")
    frame_geometry = integer(frame["geometry_identity"], "frame geometry identity")
    frame_time = scalar(frame["time"])
    require(frame_time is not None and int(step) == frame_step and int(geometry) == frame_geometry
            and Decimal(time) == Decimal(frame_time), "native volume accepted step/time/geometry differ from published frame")
    for key, value in (("VTK_accepted_step", step), ("VTK_geometry_identity", geometry), ("VTK_TimeValue", time)):
        if key in auxiliaries:
            actual = scalar(auxiliaries[key])
            require(actual is not None and Decimal(actual) == Decimal(value),
                    "native volume " + key + " differs from accepted provenance")


def published_volumes(run):
    """An atomically published frame index is authoritative, including live runs."""
    index = run / "output/meta/frames.csv"
    if not index.is_file():
        return sorted((run / "output/fields").rglob("*.vtu"))
    snapshot = index.read_text(encoding="utf-8-sig")
    require(snapshot.endswith("\n"), "published frame index has an incomplete final record")
    reader = csv.reader(io.StringIO(snapshot, newline=""), strict=True)
    try:
        columns = next(reader)
    except StopIteration:
        raise ValueError("published frame index is empty")
    require(len(columns) == len(set(columns)) and "volume" in columns, "frame index requires a unique volume column")
    column = columns.index("volume")
    sources = []
    output = (run / "output").resolve()
    seen = set()
    for row in reader:
        require(len(row) == len(columns), "published frame index record width differs from header")
        if not row[column]:
            continue
        source = (output / row[column]).resolve()
        require(source.is_relative_to(output) and source.suffix.lower() in {".vtu", ".dat"} and source.is_file(),
                "published volume is absent, unsupported or outside run output: " + row[column])
        require(source not in seen, "duplicate published volume in frame index")
        if source.suffix.lower() == ".dat":
            validate_native_volume(source, dict(zip(columns, row)))
        seen.add(source)
        sources.append(source)
    return sources


def export_run(run, overwrite=False, series_only=False, csv_encoding_name="utf-8-sig"):
    if series_only:
        selected, reason = series_run_identity(run)
        require(selected, "series-only run excluded: " + reason)
    else:
        require(selected_run(run), "requires a saved standalone Hydrodynamics run; electric/benchmark runs are excluded")
    counts = {"written": 0, "existing": 0, "no_numeric_data": 0}
    sources = ([] if series_only else published_volumes(run)) + sorted((run / "output/series").glob("*.csv"))
    for source in sources:
        if source.suffix.lower() == ".dat":
            # Published native fields were validated against their index above.
            # Never feed them to a converter or replace their authoritative bytes.
            counts["existing"] += 1
            counts["native_existing"] = counts.get("native_existing", 0) + 1
            continue
        target = source.with_suffix(".dat")
        if target.exists() and not overwrite:
            counts["existing"] += 1
            continue
        try:
            content = vtu_text(source, coordinate_names(run)) if source.suffix.lower() == ".vtu" else csv_text(source, csv_encoding_name)
        except IncompleteCsv as error:
            print("SKIP incomplete CSV snapshot: " + str(error))
            counts["incomplete_csv"] = counts.get("incomplete_csv", 0) + 1
            continue
        if content is None:
            counts["no_numeric_data"] += 1
            continue
        counts[write_atomic(target, content, overwrite)] += 1
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="one saved Hydro run, or a root with --retrofit")
    parser.add_argument("--retrofit", action="store_true", help="recursively discover saved Hydro runs under path")
    parser.add_argument("--overwrite", action="store_true", help="replace only existing derived .dat companions")
    parser.add_argument("--series-only", action="store_true", help="export only curves from identified non-electric Hydro/Interface runs; no field traversal")
    parser.add_argument("--csv-encoding", type=csv_encoding, default="utf-8-sig",
                        help="explicit strict text codec for source CSV only (default: utf-8-sig); no encoding guessing")
    args = parser.parse_args(argv)
    root = args.path.resolve()
    require(root.is_dir() and root.drive.upper() == "D:", "input/output root must be an existing directory on D:")
    if args.series_only:
        candidates = series_run_candidates(root) if args.retrofit else [root]
        runs = []
        for candidate in candidates:
            selected, reason = series_run_identity(candidate)
            if selected:
                runs.append(candidate)
            else:
                print("SKIP run " + str(candidate) + ": " + reason)
        require(runs, "no identified non-electric Hydro/Interface runs with saved series found")
    elif args.retrofit:
        runs = sorted({config.parent.parent for config in root.rglob("system/run.dict")
                       if selected_run(config.parent.parent)})
        if selected_run(root):
            runs = sorted(set(runs) | {root})
    else:
        require(selected_run(root), "path is not a selected Hydro run; use --retrofit for recursive discovery")
        runs = [root]
    require(runs, "no standalone Hydrodynamics runs found; electric/benchmark runs are excluded")
    for run in runs:
        print(str(run) + ": " + ", ".join(key + "=" + str(value) for key, value in export_run(run, args.overwrite, args.series_only, args.csv_encoding).items()))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, ET.ParseError, csv.Error) as error:
        print("Tecplot export failed: " + str(error), file=sys.stderr)
        sys.exit(2)
