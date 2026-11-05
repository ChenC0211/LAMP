#!/usr/bin/env python3
"""Formal fixed-geometry Electric App input/clock/output regression.

Usage: electric_app_contract.py LAMP.exe projectRoot workRoot
Creates five-file test data only under workRoot on D:, invokes read-only
preflight/rejections and two three-step formal numerical windows, then reads
the production receipts. This is a wiring/serialization contract, with no
analytic oracle, reference solver, tolerance change or kernel implementation.
Exit 0 means this App scope passed; any exception exits 1 and preserves evidence.
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import queue
import re
import subprocess
import sys
import threading
import time
import uuid


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def marker(path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8")
    return path.stat().st_birthtime_ns


INPUTS = {
    "run.dict": """Run { app Electric; }
Time { units code; dt 0.06; endTime 0.18; }
Parallel { threads 8; }
Output {
    prefix electric;
    Series { profile diagnostic; }
    Fields {
        enabled true; format tecplot; profile physical;
        writeControl stepInterval; writeInterval 2;
        writeInitial true; writeFirst false; writeFinal true;
    }
}
Log { everySteps 1; level summary; }
Restart { model none; }
""",
    "geometry.dict": """Geometry {
    coordinates rz;
    Domain { min (0 0); max (1 1); }
    Phases { negative 1; positive 2; }
    Initial { phi "r^2+(z-0.5)^2-0.213^2"; }
}
""",
    "mesh.dict": """Mesh { resolution (8 8); initialMesh uniform; }
""",
    "materials.dict": """Materials {
    Phases {
        liquid { phaseId 1; permittivity 5; conductivity 0.7; }
        gas { phaseId 2; permittivity 1; conductivity 0.03; }
    }
}
""",
    "electric.dict": """Electric {
    method HHO;
    interfaceMotion fixed;
    Initial { surfaceChargeDensity "0.05+0.01*z"; }
    Boundary {
        W { type Neumann; value "0"; }
        S { type Dirichlet; value "0"; }
        E { type Dirichlet; value "0.05*z"; }
        N { type Dirichlet; value "0.05"; }
    }
    TimeIntegration { method variableBdf2; }
}
""",
}


def replace(directory, filename, old, new):
    path = directory / "system" / filename
    content = path.read_text(encoding="utf-8")
    require(content.count(old) == 1, "test mutation must identify exactly one field")
    path.write_text(content.replace(old, new), encoding="utf-8")


def source_snapshot(project):
    paths = [
        "src/Apps/Electric/ElectricApp.cpp",
        "include/Apps/Electric/IO/ElectricCaseInput.h",
        "src/Apps/Electric/IO/ElectricCaseInput.cpp",
        "include/Apps/Electric/IO/ElectricCaseOutput.h",
        "src/Apps/Electric/IO/ElectricCaseOutput.cpp",
        "include/Modules/Electric/ElectricPotential.h",
        "src/Modules/Electric/ElectricPotential.cpp",
        "include/Modules/Electric/SurfaceChargeField.h",
        "src/Modules/Electric/SurfaceChargeField.cpp",
        "tests/CMakeLists.txt",
        "tests/contracts/electric/electric_app_contract.py",
    ]
    return {name: sha(project / name) for name in paths}


def invoke(executable, directory, environment, check):
    arguments = [str(executable), "--case", str(directory), "--threads", "2"]
    if check:
        arguments.append("--check-case")
    start = marker(directory / "process_started.json", {"arguments": arguments})
    observed = []
    lines = []
    first_accepted = completed = None
    events = queue.Queue()
    deadline = time.monotonic() + 180
    process = subprocess.Popen(arguments, cwd=directory, env=environment,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, encoding="utf-8", errors="replace",
                               bufsize=1)

    def read_stream():
        try:
            for line in process.stdout:
                events.put((time.perf_counter_ns(), line))
        finally:
            events.put(None)

    reader = threading.Thread(target=read_stream, daemon=True)
    reader.start()
    ended = False
    try:
        with (directory / "process.log").open("w", encoding="utf-8") as log:
            while not ended:
                require(time.monotonic() < deadline, "formal App invocation exceeded 180 seconds")
                try:
                    event = events.get(timeout=0.1)
                except queue.Empty:
                    continue
                if event is None:
                    ended = True
                    continue
                received, line = event
                lines.append(line)
                observed.append({"received_perf_counter_ns": received, "text": line.rstrip()})
                log.write(line)
                log.flush()
                if first_accepted is None and re.match(r"\[Step 1\]", line):
                    first_accepted = marker(directory / "first_accepted_observed.json",
                                            {"received_perf_counter_ns": received, "text": line.rstrip()})
                if line.startswith("[COMPLETED]"):
                    completed = marker(directory / "completed_observed.json",
                                       {"received_perf_counter_ns": received, "text": line.rstrip()})
        code = process.wait(timeout=max(1.0, deadline-time.monotonic()))
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=10)
        reader.join(timeout=2)
        process.stdout.close()
    finish = marker(directory / "process_finished.json", {"returncode": code})
    observation = {
        "arguments": arguments, "returncode": code,
        "start_file_birthtime_ns": start, "finish_file_birthtime_ns": finish,
        "file_creation_delta_seconds": (finish-start)/1e9,
        "first_accepted_file_birthtime_ns": first_accepted,
        "completed_file_birthtime_ns": completed,
        "stream_observations": observed,
    }
    (directory / "process_receipt.json").write_text(json.dumps(observation, indent=2)+"\n",
                                                   encoding="utf-8")
    return code, "".join(lines), observation


def finite_number(value, description):
    result = float(value)
    require(math.isfinite(result), description + " is nonfinite")
    return result


def tecplot(path, expected_variables):
    lines = path.read_text(encoding="utf-8").splitlines()
    start = next((i for i, line in enumerate(lines)
                  if re.match(r"^\s*[+-]?(?:\d|\.\d)", line)), None)
    require(start is not None, str(path) + " has no numerical payload")
    header = "\n".join(lines[:start])
    variables_header = re.search(r"\bVARIABLES\s*=\s*(.*?)\bZONE\b", header, flags=re.DOTALL)
    require(variables_header is not None, "missing Tecplot VARIABLES")
    variables = re.findall(r'"([^"]+)"', variables_header.group(1))
    require(variables == expected_variables, path.name + " has the wrong public variable list")
    require(len(re.findall(r"^\s*ZONE\b", header, flags=re.MULTILINE)) == 1,
            "field must publish exactly one zone")
    require(re.search(r"\bDATAPACKING\s*=\s*BLOCK", header), "field must use BLOCK packing")
    types = re.search(r"\bDT\s*=\s*\((.*?)\)", header, flags=re.DOTALL)
    require(types is not None and types.group(1).split() == ["DOUBLE"]*len(variables),
            "field variables must be explicit DOUBLE")
    require(re.search(r"\bZONETYPE\s*=\s*ORDERED", header), "electric field must use actual point samples")
    count = re.search(r"\bI\s*=\s*(\d+)", header)
    require(count is not None and int(count.group(1)) > 0, "field has no source points")
    count = int(count.group(1))
    tokens = " ".join(lines[start:]).split()
    require(len(tokens) == count*len(variables), "field payload length differs from declared points")
    numbers = [finite_number(token, path.name) for token in tokens]
    auxiliary = dict(re.findall(r'^\s*AUXDATA\s+(\w+)\s*=\s*"([^"]*)"\s*$',
                                header, flags=re.MULTILINE))
    solution = re.search(r"\bSOLUTIONTIME\s*=\s*([^\s,]+)", header)
    require(solution is not None, "field omitted physical solution time")
    return {"variables": variables, "points": count, "auxiliary": auxiliary,
            "time": finite_number(solution.group(1), path.name),
            "blocks": {name: numbers[index*count:(index+1)*count]
                       for index, name in enumerate(variables)}}


def tecplot_contour(path, expected_variables):
    """Parse the public FE serialization independently of the C++ writer."""
    lines = path.read_text(encoding="utf-8").splitlines()
    zones = [index for index, line in enumerate(lines) if re.match(r"^\s*ZONE\b", line)]
    require(len(zones) == 2, "one contour file must contain both phase zones")
    require(len(re.findall(r"^\s*VARIABLES\s*=", "\n".join(lines), flags=re.MULTILINE)) == 1,
            "contour must publish one common VARIABLES header")
    variables_header = re.search(r"\bVARIABLES\s*=\s*(.*?)\bZONE\b",
                                 "\n".join(lines), flags=re.DOTALL)
    require(variables_header is not None, "contour has no VARIABLES header")
    variables = re.findall(r'"([^"]+)"', variables_header.group(1))
    require(variables == expected_variables, "contour has the wrong physical variable list")
    result = []
    for zone_index, begin in enumerate(zones):
        end = zones[zone_index+1] if zone_index+1 < len(zones) else len(lines)
        segment = lines[begin:end]
        data = next((index for index, line in enumerate(segment)
                     if re.match(r"^\s*[+-]?(?:\d|\.\d)", line)), None)
        require(data is not None, "contour zone has no numerical payload")
        header = "\n".join(segment[:data])
        require(re.search(r"\bZONETYPE\s*=\s*FETRIANGLE\b", header),
                "contour must carry FE triangle connectivity")
        require(re.search(r"\bDATAPACKING\s*=\s*BLOCK\b", header),
                "contour must use BLOCK packing")
        require("VARLOCATION" not in header,
                "contour physical values must be original element-local nodal samples")
        types = re.search(r"\bDT\s*=\s*\((.*?)\)", header, flags=re.DOTALL)
        require(types is not None and types.group(1).split() == ["DOUBLE"]*len(variables),
                "contour variables must be explicit DOUBLE")
        nodes = re.search(r"\bNODES\s*=\s*(\d+)", header)
        elements = re.search(r"\bELEMENTS\s*=\s*(\d+)", header)
        require(nodes is not None and elements is not None,
                "contour omitted FE node/element counts")
        nodes, elements = int(nodes.group(1)), int(elements.group(1))
        require(nodes > 0 and elements > 0, "contour phase has no mesh")
        tokens = " ".join(segment[data:]).split()
        field_size = nodes*len(variables)
        require(len(tokens) == field_size+3*elements,
                "contour payload length differs from FE declaration")
        numbers = [finite_number(token, path.name) for token in tokens[:field_size]]
        raw_connectivity = tokens[field_size:]
        require(all(re.fullmatch(r"[1-9]\d*", token) for token in raw_connectivity),
                "contour connectivity must contain positive integer indices")
        connectivity = [tuple(int(token)-1 for token in raw_connectivity[index:index+3])
                        for index in range(0, len(raw_connectivity), 3)]
        require(all(len(set(triangle)) == 3 and all(0 <= node < nodes for node in triangle)
                    for triangle in connectivity),
                "contour connectivity is degenerate or crosses its phase zone")
        auxiliary = dict(re.findall(r'^\s*AUXDATA\s+(\w+)\s*=\s*"([^"]*)"\s*$',
                                    header, flags=re.MULTILINE))
        solution = re.search(r"\bSOLUTIONTIME\s*=\s*([^\s,]+)", header)
        require(solution is not None, "contour omitted physical solution time")
        result.append({"variables": variables, "points": nodes, "triangles": elements,
                       "auxiliary": auxiliary, "connectivity": connectivity,
                       "time": finite_number(solution.group(1), path.name),
                       "blocks": {name: numbers[index*nodes:(index+1)*nodes]
                                  for index, name in enumerate(variables)}})
    return result


def verify_field_source(zone, step, identity):
    """Keep the original method, immutable geometry and full electric-clock gates."""
    aux = zone["auxiliary"]
    require(aux["Method"] == "HHO" and aux["InterfaceMotion"] == "fixed"
            and aux["CoordinateSystem"] == "axisymmetricRZ",
            "field public method/geometry source differs")
    require(int(aux["AcceptedStep"]) == step and int(aux["GeometryIdentity"]) == identity,
            "phase/interface output does not share the actual accepted source")
    require(aux["ElectricClockPosition"] == ("Initial" if step == 0 else "End"),
            "output physical clock position differs from accepted state")
    near(float(aux["ElectricClockBegin"]), 0 if step == 0 else (step-1)*.06,
         "output physical interval begin")
    near(float(aux["ElectricClockEnd"]), step*.06, "output physical interval end")
    near(zone["time"], step*.06, "output solution time")
    require(aux["GeometrySourceStep"] == "0" and aux["GeometrySourceClockPosition"] == "Initial"
            and float(aux["GeometrySourceClockBegin"]) == 0
            and float(aux["GeometrySourceClockEnd"]) == 0,
            "immutable Initial geometry source was rewritten")


def near(actual, expected, description, relative=5e-13):
    require(abs(actual-expected) <= relative*(1+abs(expected)), description)


def verify_performance(directory, history_rows, frame_steps):
    """Check passive timing ownership and clocks, without a wall-time gate."""
    path = directory / "output/series/performance.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames
        require(columns is not None and len(columns) == len(set(columns)),
                "performance receipt has missing or duplicate column names")
        serialized = list(reader)
    initial_fields = [
        "input_s", "output_startup_s", "mesh_s", "publication_s", "phase_geometry_s",
        "geometry_phase_quadrature_s", "geometry_cartesian_faces_s",
        "geometry_interface_faces_s", "geometry_local_operators_s",
        "geometry_matrix_comparison_s", "geometry_remaining_s", "panel_charge_s",
        "moment_prepare_s", "initial_electric_s", "electric_boundary_s",
        "electric_local_operators_s", "electric_local_loads_s", "electric_assembly_s",
        "electric_factor_s", "electric_solve_s", "electric_recovery_s",
        "electric_gauss_affine_s", "charge_prepare_s", "initial_current_s",
    ]
    step_fields = [
        "charge_advance_s", "charge_schur_s", "charge_schur_inner_electric_s",
        "charge_moment_load_s", "charge_final_electric_s", "charge_recovery_s",
        "charge_current_s", "charge_balance_s",
    ]
    output_fields = [
        "output_cloud_build_s", "output_cloud_write_s", "output_samples_s",
    ]
    fixed_counts = [
        "phase_elements", "free_skeleton_dofs", "panels", "volume_points",
        "surface_points", "quadrature_refinement",
    ]
    required = set(initial_fields + step_fields + output_fields + fixed_counts + [
        "step", "time", "schur_iterations", "electric_solves", "total_s",
        "maxwell_s", "commit_s", "output_s", "output_validation_s",
        "output_series_s", "app_unclassified_s",
    ])
    require(required.issubset(columns), "performance receipt omitted a public stage")
    require(all(set(row) == set(columns) and None not in row.values() for row in serialized),
            "performance receipt row width differs from its header")
    rows = [{name: finite_number(value, "performance " + name)
             for name, value in row.items()} for row in serialized]
    require([row["step"] for row in rows] == [0, 1, 2, 3],
            "performance receipt must follow every accepted physical step")
    console = {}
    for match in re.finditer(r"^\[Step (\d+)\][^\r\n]*?iterations=(\d+) electricSolves=(\d+)",
                             (directory / "process.log").read_text(encoding="utf-8"),
                             flags=re.MULTILINE):
        require(int(match.group(1)) not in console, "duplicated console iteration receipt")
        console[int(match.group(1))] = (int(match.group(2)), int(match.group(3)))
    require(set(console) == {1, 2, 3}, "console omitted actual iteration/solve counts")
    parents = {
        "phase_geometry_s": ["geometry_phase_quadrature_s", "geometry_cartesian_faces_s",
                             "geometry_interface_faces_s", "geometry_local_operators_s",
                             "geometry_matrix_comparison_s", "geometry_remaining_s"],
        "initial_electric_s": ["electric_boundary_s", "electric_local_operators_s",
                              "electric_local_loads_s", "electric_assembly_s",
                              "electric_factor_s", "electric_solve_s",
                              "electric_recovery_s", "electric_gauss_affine_s"],
        "charge_advance_s": ["charge_schur_s", "charge_moment_load_s",
                             "charge_final_electric_s", "charge_recovery_s",
                             "charge_current_s", "charge_balance_s"],
        "charge_schur_s": ["charge_schur_inner_electric_s"],
        "output_s": ["output_validation_s", "output_cloud_build_s",
                     "output_cloud_write_s", "output_samples_s", "output_series_s"],
    }
    app_parts = [
        "input_s", "output_startup_s", "mesh_s", "publication_s", "phase_geometry_s",
        "panel_charge_s", "moment_prepare_s", "initial_electric_s", "charge_prepare_s",
        "initial_current_s", "charge_advance_s", "maxwell_s", "commit_s", "output_s",
        "app_unclassified_s",
    ]
    for step, (timing, history) in enumerate(zip(rows, history_rows)):
        near(timing["time"], history["time"], "performance physical clock differs from history")
        for name in columns:
            if name.endswith("_s"):
                require(timing[name] >= 0, "negative passive stage timing " + name)
        require(all(timing[name] == rows[0][name]
                    and timing[name] >= (0 if name == "quadrature_refinement" else 1)
                    and timing[name].is_integer() for name in fixed_counts),
                "fixed geometry/count receipt changed between accepted states")
        require(all(timing[name].is_integer() and timing[name] >= 0
                    for name in ("schur_iterations", "electric_solves")),
                "performance iteration/solve counters are not actual nonnegative integers")
        near(timing["schur_iterations"], history["Schur_iterations"],
             "performance Schur iterations differ from physical history")
        if step == 0:
            require(timing["schur_iterations"] == 0 and timing["electric_solves"] == 0
                    and all(timing[name] == 0 for name in step_fields),
                    "initial receipt was presented as an accepted charge advance")
        else:
            require(all(timing[name] == 0 for name in initial_fields),
                    "accepted-step receipt inherited initialization timings")
            require(timing["schur_iterations"] > 0
                    and timing["electric_solves"] >= timing["schur_iterations"] + 2,
                    "nonzero evolving charge omitted actual Schur/electric work")
            require(console[step] == (timing["schur_iterations"], timing["electric_solves"]),
                    "performance counters differ from real-time console")
        for name in output_fields:
            require(timing[name] > 0 if step in frame_steps else timing[name] == 0,
                    "field timing does not follow the formal sparse output schedule")
        for parent, children in parents.items():
            # Inclusive child clocks are bounded by their own parent. They are
            # never added a second time to the top-level application duration.
            require(math.fsum(timing[name] for name in children)
                    <= timing[parent] + 1e-8*(1 + timing[parent]),
                    "passive child timing escaped its parent " + parent)
        near(math.fsum(timing[name] for name in app_parts), timing["total_s"],
             "exclusive App timings plus unclassified time do not close")
    return {"path": "output/series/performance.csv", "rows": rows,
            "clock_aligned_with_history": True, "sparse_field_timings": frame_steps,
            "inclusive_parents_checked": True, "wall_time_speed_gate": False}


def verify_output(directory, series_profile="diagnostic"):
    meta = directory / "output/meta"
    history_path = directory / "output/series/electric_electric.dat"
    lines = history_path.read_text(encoding="utf-8").splitlines()
    start = next((i for i, line in enumerate(lines)
                  if re.match(r"^\s*[+-]?(?:\d|\.\d)", line)), None)
    require(start is not None, "accepted physical history has no rows")
    header = "\n".join(lines[:start])
    match = re.search(r"\bVARIABLES\s*=\s*(.*?)\bZONE\b", header, flags=re.DOTALL)
    require(match is not None, "history has no VARIABLES header")
    columns = re.findall(r'"([^"]+)"', match.group(1))
    expected = ["step", "time", "geometry_source_step", "geometry_source_time",
                "total_Q", "panel_current", "exterior_current", "bulk_current",
                "current_balance", "current_balance_relative", "linear_residual",
                "local_Gauss_relative", "global_Gauss", "global_Gauss_relative",
                "time_update_present", "temporal_order", "dt", "previous_dt",
                "a0", "a1", "a2", "Schur_iterations", "Schur_relative",
                "panel_equation_relative", "charge_equation", "charge_equation_relative",
                "Maxwell_charge_jump_relative"]
    require(columns == expected, "physical history public columns differ")
    rows = []
    for line in lines[start:]:
        tokens = line.split()
        require(len(tokens) == len(columns), "history row width differs from its header")
        rows.append(dict(zip(columns, (finite_number(value, "history "+column)
                                       for column, value in zip(columns, tokens)))))
    require([row["step"] for row in rows] == [0, 1, 2, 3], "history did not retain all accepted physical steps")
    require([row["temporal_order"] for row in rows] == [0, 1, 2, 2], "history did not retain actual BE/BDF2 order")
    require([row["time_update_present"] for row in rows] == [0, 1, 1, 1],
            "initial placeholders were presented as a time equation")
    for index, row in enumerate(rows):
        near(row["time"], index*.06, "history physical time")
        require(row["geometry_source_step"] == 0 and row["geometry_source_time"] == 0,
                "fixed geometry source was retagged as a physical electric state")
        for name, limit in (("linear_residual", 1e-11), ("local_Gauss_relative", 1e-10),
                            ("global_Gauss_relative", 1e-10), ("current_balance_relative", 1e-10),
                            ("Maxwell_charge_jump_relative", 1e-10)):
            require(0 <= row[name] <= limit, "original public "+name+" gate failed")
        current_defect = row["panel_current"] + row["exterior_current"] - row["bulk_current"]
        current_scale = abs(row["panel_current"])+abs(row["exterior_current"])+abs(row["bulk_current"])
        require(abs(current_defect) <= 1e-10*current_scale,
                "published original Ohmic current ledger does not close")
        if index:
            for name, limit in (("Schur_relative", 1e-11), ("panel_equation_relative", 1e-10),
                                ("charge_equation_relative", 1e-10)):
                require(0 <= row[name] <= limit, "original public "+name+" gate failed")
            near(row["dt"], .06, "accepted physical dt")
            near(row["previous_dt"], 0 if index == 1 else .06, "BDF2 previous physical dt")
            coefficients = (1., -1., 0.) if index == 1 else (1.5, -2., .5)
            for name, expected_value in zip(("a0", "a1", "a2"), coefficients):
                near(row[name], expected_value, "actual BE/BDF2 coefficient "+name)
            terms = [row["a0"]*row["total_Q"], row["a1"]*rows[index-1]["total_Q"],
                     row["a2"]*rows[index-2]["total_Q"] if index > 1 else 0.,
                     -row["dt"]*row["panel_current"]]
            require(abs(math.fsum(terms)) <= 1e-10*math.fsum(abs(value) for value in terms),
                    "published Q and raw same-field current do not close the time equation")
    require(rows[0]["total_Q"] > 0 and any(row["total_Q"] != rows[0]["total_Q"] for row in rows[1:]),
            "nonzero configured charge did not actually evolve")

    with (meta / "frames.csv").open(encoding="utf-8", newline="") as stream:
        frame_reader = csv.DictReader(stream)
        require(frame_reader.fieldnames == ["step", "time", "geometry_identity", "liquid",
                                           "gas", "interface", "cloud"],
                "frame links must retain original sample columns and append the contour")
        frames = list(frame_reader)
    require([int(frame["step"]) for frame in frames] == [0, 2, 3],
            "field schedule must use accepted steps 0/2/3 rather than frame indices")
    identity = int(frames[0]["geometry_identity"])
    require(identity > 0 and all(int(frame["geometry_identity"]) == identity for frame in frames),
            "fixed output geometry identity changed")
    volume_variables = ["r [code]", "z [code]", "V [code]", "epsilon [code]",
                        "sigma [code]", "Er [code]", "Ez [code]"]
    interface_variables = ["r [code]", "z [code]", "nr", "nz", "V_Gamma [code]",
                           "q_h [code]", "Et [code]", "Dn_liquid [code]", "Dn_gas [code]",
                           "Maxwell_normal [code]", "Maxwell_tangent [code]",
                           "Maxwell_r [code]", "Maxwell_z [code]"]
    contour_variables = ["r [code]", "z [code]", "V [code]", "Er [code]", "Ez [code]",
                         "E_mag [code]", "phase_id"]
    expected_files = set()
    point_counts = []
    contour_meshes = []
    for frame in frames:
        step = int(frame["step"])
        suffix = str(step).zfill(6)+".dat"
        near(float(frame["time"]), step*.06, "frame receipt time")
        for kind in ("liquid", "gas", "interface"):
            filename = ("interface/electric_interface_samples_step"+suffix if kind == "interface"
                        else "volume/samples/electric_phase_"+kind+"_samples_step"+suffix)
            require(frame[kind] == filename, "field filename does not encode the actual accepted step")
            expected_files.add(filename)
            zone = tecplot(directory / "output/fields" / filename,
                           interface_variables if kind == "interface" else volume_variables)
            verify_field_source(zone, step, identity)
            aux = zone["auxiliary"]
            if kind != "interface":
                epsilon, sigma = (5., .7) if kind == "liquid" else (1., .03)
                require(all(value == epsilon for value in zone["blocks"]["epsilon [code]"])
                        and all(value == sigma for value in zone["blocks"]["sigma [code]"]),
                        "one-sided phase output mixed or ignored actual material coefficients")
                require(aux["PhaseId"] == ("1" if kind == "liquid" else "2"),
                        "phase output differs from geometry/material phase IDs")
            else:
                require(any(value != 0 for value in zone["blocks"]["q_h [code]"]),
                        "interface output lost actual nonzero q_h")
            point_counts.append({"step": step, "kind": kind, "points": zone["points"]})

        cloud = "volume/electric_fields_step"+suffix
        require(frame["cloud"] == cloud, "contour filename does not encode the accepted step")
        expected_files.add(cloud)
        zones = tecplot_contour(directory / "output/fields" / cloud, contour_variables)
        require(sorted(zone["auxiliary"].get("PhaseId") for zone in zones) == ["1", "2"],
                "one contour file must retain each physical phase exactly once")
        phase_areas = []
        for zone in zones:
            verify_field_source(zone, step, identity)
            phase = int(zone["auxiliary"]["PhaseId"])
            blocks = zone["blocks"]
            require(all(value == phase for value in blocks["phase_id"]),
                    "contour triangles interpolate or mix the two phase IDs")
            require(all(-5e-14 <= value <= 1+5e-14
                        for name in ("r [code]", "z [code]") for value in blocks[name]),
                    "contour nodes escaped the actual fixed domain")
            for radial, axial, magnitude in zip(blocks["Er [code]"], blocks["Ez [code]"],
                                                 blocks["E_mag [code]"]):
                require(magnitude >= 0, "electric field magnitude must be nonnegative")
                near(magnitude, math.hypot(radial, axial),
                     "contour magnitude differs from the same one-sided electric field")
            r, z = blocks["r [code]"], blocks["z [code]"]
            areas = []
            for i, j, k in zone["connectivity"]:
                signed_area = .5*((r[j]-r[i])*(z[k]-z[i])-(z[j]-z[i])*(r[k]-r[i]))
                require(signed_area > 0, "contour mesh has an inverted or zero-area triangle")
                areas.append(signed_area)
            area = math.fsum(areas)
            require(area > 0, "contour phase has zero meridional area")
            phase_areas.append(area)
            contour_meshes.append({"step": step, "phase_id": phase, "points": zone["points"],
                                   "triangles": zone["triangles"], "meridional_area": area})
        # This is the display mesh's planar coverage, not an axisymmetric volume,
        # an integration rule, or an analytic phase-sign test on curved chords.
        near(math.fsum(phase_areas), 1.,
             "the two contour phase meshes do not cover the original rectangular domain")
    fields = directory / "output/fields"
    require({path.relative_to(fields).as_posix() for path in fields.rglob("*.dat")} == expected_files,
            "unexpected or missing physical fields")
    require((fields / "README.txt").is_file(),
            "field output omitted the entry point and sample/contour explanation")
    hashes = json.loads((meta / "input_hashes.json").read_text(encoding="utf-8"))
    require(set(hashes) == {"system/"+filename for filename in INPUTS},
            "source/effective provenance did not cover exactly the five input owners")
    for filename in INPUTS:
        entry = hashes["system/"+filename]
        require(entry["source_sha256"] == sha(directory / "system" / filename)
                == sha(meta / "input" / filename),
                "input source SHA differs from the saved immutable source")
        require(entry["effective_sha256"] == sha(meta / "effective" / filename),
                "resolved effective input SHA differs from the saved document")
    effective_run = (meta / "effective/run.dict").read_text(encoding="utf-8")
    require(re.search(r"Parallel\s*\{\s*threads\s+2\s*;", effective_run),
            "effective dictionary lost CLI thread override")
    require(re.search(r"Series\s*\{\s*profile\s+"+series_profile+r"\s*;", effective_run),
            "effective dictionary lost the selected series profile")
    with (meta / "configuration_trace.csv").open(encoding="utf-8", newline="") as stream:
        entries = list(csv.DictReader(stream))
    thread = [entry for entry in entries if entry["field"] == "Parallel.threads"]
    require(len(thread) == 1 and thread[0]["effective_value"] == "2"
            and "--threads" in thread[0]["consumer"] and thread[0]["used"] == "true",
            "configuration trace lost the resolved CLI thread ownership")
    require("build_type=Release" in (meta / "build_info.txt").read_text(encoding="utf-8"),
            "formal numerical evidence must use Release")
    performance = None
    if series_profile == "diagnostic":
        performance = verify_performance(directory, rows, [0, 2, 3])
    else:
        require(not (directory / "output/series/performance.csv").exists(),
                "physical series profile published development performance output")
    return {"history_rows": rows, "field_steps": [0, 2, 3], "geometry_identity": identity,
            "series_profile": series_profile,
            "performance": performance,
            "point_counts": point_counts, "contour_meshes": contour_meshes,
            "source_and_effective_SHA_matched": True,
            "thread_override": 2, "clock_scope": "matching Q/V/I/Maxwell; immutable Initial geometry",
            "charge_density_integral_scope": "writer verifies original weighted Q moments; no extra test quadrature"}


def main():
    require(len(sys.argv) == 4,
            "usage: electric_app_contract.py LAMP.exe projectRoot workRoot")
    executable, project, work = (Path(arg).resolve() for arg in sys.argv[1:])
    require(executable.is_file(), "missing formal LAMP executable")
    require(work.drive.upper() == "D:", "Electric App contract workRoot must be on D:")
    task = work / ("electric_app_" + uuid.uuid4().hex[:12])
    task.mkdir(parents=True, exist_ok=False)
    (work / "last_run.txt").write_text(str(task)+"\n", encoding="utf-8")
    temporary = task / "_temp"
    temporary.mkdir()
    environment = os.environ.copy()
    environment.update(TEMP=str(temporary), TMP=str(temporary), OMP_NUM_THREADS="2",
                       OMP_THREAD_LIMIT="2", OMP_DYNAMIC="false",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       PYTHONDONTWRITEBYTECODE="1")
    binary_before = sha(executable)
    source_before = source_snapshot(project)
    result = {"status": "RUNNING", "scope": "formal_fixed_RZ_electric_App_wiring",
              "task": str(task), "executable": str(executable),
              "executable_sha256": binary_before, "source_sha256": source_before,
              "actual_numerical_windows": 0, "checks": []}

    def fresh(name):
        directory = task / name
        (directory / "system").mkdir(parents=True)
        for filename, content in INPUTS.items():
            (directory / "system" / filename).write_text(content, encoding="utf-8")
        return directory

    def reject(name, mutation, filename, field):
        directory = fresh(name)
        mutation(directory)
        code, text, receipt = invoke(executable, directory, environment, True)
        require(code != 0, name + " was silently accepted")
        require(not (directory / "output").exists() and not (directory / "cache").exists(),
                name + " preflight produced computational state")
        require("[INIT]" not in text and "[Step " not in text,
                name + " preflight constructed mesh or advanced a solver")
        require(re.search(re.escape(filename)+r":\d+:", text),
                name + " rejection omitted actual file:line")
        require(field in text, name + " rejection omitted failing field")
        result["checks"].append({"name": name, "status": "PASS", "process": receipt})

    try:
        directory = fresh("valid_preflight")
        code, text, receipt = invoke(executable, directory, environment, True)
        require(code == 0 and "ELECTRIC_INPUT_OK" in text, "valid five-file preflight failed")
        require("[INIT]" not in text and "[Step " not in text
                and not (directory / "output").exists() and not (directory / "cache").exists(),
                "read-only preflight constructed mesh, solved or output")
        require(re.search(r"\[Runtime\].*\bthreads=2\b.*\bsource=command_line\b", text),
                "preflight did not honor CLI threads over dictionary threads=8")
        result["checks"].append({"name": "valid_preflight", "status": "PASS", "process": receipt})
        reject("old_method", lambda d: replace(d, "electric.dict", "method HHO;",
                                               "method sharpHermiteP3;"), "electric.dict", "method")
        reject("old_method_name", lambda d: replace(d, "electric.dict", "method HHO;",
                                                    "method sharpHHO2;"), "electric.dict", "method")
        reject("unknown_key", lambda d: replace(d, "electric.dict", "interfaceMotion fixed;",
                                               "interfaceMotion fixed;\n    unknownKey true;"),
               "electric.dict", "unknownKey")
        reject("misplaced_material", lambda d: replace(d, "electric.dict", "interfaceMotion fixed;",
                                                      "interfaceMotion fixed;\n    permittivity 5;"),
               "electric.dict", "permittivity")
        reject("unanchored", lambda d: replace(d, "electric.dict",
               'S { type Dirichlet; value "0"; }\n        E { type Dirichlet; value "0.05*z"; }\n        N { type Dirichlet; value "0.05"; }',
               'S { type Neumann; value "0"; }\n        E { type Neumann; value "0"; }\n        N { type Neumann; value "0"; }'),
               "electric.dict", "Boundary")
        reject("wrong_axis", lambda d: replace(d, "electric.dict", 'W { type Neumann; value "0"; }',
                                               'W { type Dirichlet; value "0"; }'),
               "electric.dict", "W")
        reject("restart", lambda d: replace(d, "run.dict", "Restart { model none; }",
                                            "Restart { model read; }"), "run.dict", "model")
        reject("unsupported_amr", lambda d: replace(d, "mesh.dict", "initialMesh uniform;",
                                                    "initialMesh amr;"), "mesh.dict", "initialMesh")
        reject("unsupported_motion", lambda d: replace(d, "electric.dict", "interfaceMotion fixed;",
                                                       "interfaceMotion moving;"),
               "electric.dict", "interfaceMotion")

        directory = fresh("formal_window")
        result["actual_numerical_windows"] = 1
        code, text, receipt = invoke(executable, directory, environment, False)
        require(code == 0 and "[COMPLETED] Electric steps=3" in text,
                "one formal three-step electric window failed")
        require(re.search(r"\[Runtime\].*\bthreads=2\b.*\bsource=command_line\b", text),
                "formal window did not honor CLI threads over dictionary threads=8")
        summaries = [line for line in text.splitlines() if re.match(r"\[Step \d+\]", line)]
        require([int(re.search(r"\[Step (\d+)\]", line).group(1)) for line in summaries] == [1, 2, 3],
                "console accepted-step progress is missing, repeated or numbered as frames")
        require([int(re.search(r"\border=(\d+)", line).group(1)) for line in summaries] == [1, 2, 2],
                "formal charge time method did not perform BE startup then true BDF2")
        observations = receipt["stream_observations"]
        first_index = next(i for i, row in enumerate(observations) if re.match(r"\[Step 1\]", row["text"]))
        final_index = next(i for i, row in enumerate(observations) if row["text"].startswith("[COMPLETED]"))
        require(first_index < final_index and receipt["first_accepted_file_birthtime_ns"] is not None
                and receipt["completed_file_birthtime_ns"] > receipt["first_accepted_file_birthtime_ns"],
                "accepted progress was not observed before completion")
        output = verify_output(directory)
        for filename in INPUTS:
            require(sha(directory / "system" / filename) == sha(directory / "output/meta/input" / filename),
                    "formal output did not preserve exact five-file source input")
        effective_directory = fresh("effective_recheck")
        for filename in INPUTS:
            (effective_directory / "system" / filename).write_bytes(
                (directory / "output/meta/effective" / filename).read_bytes())
        effective_code, effective_text, effective_receipt = invoke(
            executable, effective_directory, environment, True)
        require(effective_code == 0 and "ELECTRIC_INPUT_OK" in effective_text
                and "[INIT]" not in effective_text and not (effective_directory / "output").exists(),
                "saved effective five-file configuration could not be read back without computation")
        result["checks"].append({"name": "effective_recheck", "status": "PASS", "process": effective_receipt})
        result["checks"].append({"name": "formal_window", "status": "PASS",
                                 "process": receipt, "output": output})

        physical_directory = fresh("physical_series_window")
        replace(physical_directory, "run.dict", "Series { profile diagnostic; }",
                "Series { profile physical; }")
        result["actual_numerical_windows"] += 1
        physical_code, physical_text, physical_receipt = invoke(
            executable, physical_directory, environment, False)
        require(physical_code == 0 and "[COMPLETED] Electric steps=3" in physical_text,
                "physical series replay of the formal electric window failed")
        physical_output = verify_output(physical_directory, "physical")
        physical_paths = [Path("series/electric_electric.dat")]
        physical_paths.extend(path.relative_to(directory / "output")
                              for path in sorted((directory / "output/fields").rglob("*.dat")))
        for path in physical_paths:
            require((directory / "output" / path).read_bytes()
                    == (physical_directory / "output" / path).read_bytes(),
                    "series profile changed physical output " + path.as_posix())
        result["checks"].append({"name": "physical_series_window", "status": "PASS",
                                 "process": physical_receipt, "output": physical_output,
                                 "identical_physical_files": [path.as_posix() for path in physical_paths]})
        result["source_after_sha256"] = source_snapshot(project)
        result["executable_after_sha256"] = sha(executable)
        require(result["source_after_sha256"] == source_before
                and result["executable_after_sha256"] == binary_before,
                "sources or the built executable changed during this frozen check")
        result["status"] = "PASS_WITHIN_FORMAL_FIXED_RZ_electric_APP_SCOPE"
        result["not_proven"] = ["coarse_P5_Maxwell_accuracy", "moving_charge_convergence",
                               "Maxwell_Hydro_coupling", "AMR_small_cut", "electrospray_completion"]
    except Exception as error:
        result["status"] = "FAIL"
        result["error"] = str(error)
        raise
    finally:
        (task / "receipt.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n",
                                         encoding="utf-8")
        print("Electric electric App contract receipt: " + str(task / "receipt.json"), flush=True)
    print(result["status"], flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("Electric electric App contract FAIL: " + str(error), file=sys.stderr, flush=True)
        sys.exit(1)
