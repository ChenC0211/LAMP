#!/usr/bin/env python3
"""Formal native RZ sharp-interface electrostatics/Hydro AC App wiring, clock and output contract.

Exactly one two-step production window follows read-only eight-file preflight.
No analytic oracle or numerical gate changes. Generated evidence stays on D:.
"""
import csv
from datetime import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
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


def marker(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    return path.stat().st_birthtime_ns


def helper(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


INPUTS = {
    "run.dict": """Run { app Electrohydrodynamics; }
Parallel { threads 8; }
Time { units code; dt .00025; endTime .0005; }
Output { Fields { enabled true; writeEverySteps 2; } }
Log { everySteps 2; level detail; }
""",
    "geometry.dict": """Geometry {
    coordinates rz;
    Domain { min (0 0); max (1 1); }
    Phases { negative 1; positive 2; }
    Initial { phi "r^2+(z-centerZ)^2-radius^2"; }
}
""",
    "mesh.dict": "Mesh { resolution (4 4); initialMesh uniform; }\n",
    "fields.dict": """Fields {
    Constants { centerZ .51; radius .27; electricStrength .05; chargeAmplitude .1; }
}
""",
    "interface.dict": "Interface { method hermiteGals; Boundary { type none; } }\n",
    "materials.dict": """Materials {
    Phases {
        Liquid { phaseId 1; density 1000; dynamicViscosity 0; permittivity 5; conductivity .7; }
        Gas { phaseId 2; density 1; dynamicViscosity 0; permittivity 1; conductivity .03; }
    }
    Interfaces { Shared { phases (1 2); surfaceTension 0; } }
}
""",
    "hydro.dict": """Hydrodynamics {
    method staggered; densityRepresentation sharpGeometry; capillaryModel balancedCsf;
    Initial { velocity (0 0); }
    Boundary { type impermeable; }
    BodyForce { acceleration (0 0); }
    Numerics { characteristicSubsteps 1; materialInitialPanels 1; }
    TimeIntegration { method adaptiveLinear; }
}
""",
    "electric.dict": """Electric {
    method HHO; interfaceMotion coupled;
    Initial { surfaceChargeDensity "chargeAmplitude*(sin(1.7*(z-centerZ))+.05*(z-centerZ))"; }
    Boundary {
        W { type Neumann; value "0"; }
        S { type Dirichlet; value "0"; }
        E { type Dirichlet; value "electricStrength*z"; }
        N { type Dirichlet; value "electricStrength"; }
    }
    TimeIntegration { method variableBdf2; }
}
""",
}


def replace(directory, filename, old, new):
    path = directory / "system" / filename
    text = path.read_text(encoding="utf-8")
    require(text.count(old) == 1, "mutation must identify exactly one record: "+old)
    path.write_text(text.replace(old, new), encoding="utf-8")


def source_snapshot(project):
    paths = set()
    for owner in ("include", "src"):
        paths.update(path for path in (project / owner).rglob("*")
                     if path.is_file() and path.suffix in (".h", ".hpp", ".cpp", ".c"))
    paths.update(project / name for name in (
        "CMakeLists.txt", "tests/CMakeLists.txt",
        "tests/contracts/hydrodynamics/hydro_electric_app_contract.py",
        "tests/contracts/hydrodynamics/hydro_app_contract.py",
        "tests/contracts/electric/electric_app_contract.py"))
    return {path.relative_to(project).as_posix(): sha(path) for path in sorted(paths)}


def invoke(executable, directory, environment, check):
    arguments = [str(executable), "--case", str(directory), "--threads", "2"]
    if check:
        arguments.append("--check-case")
    started = marker(directory / "process_started.json", {"arguments": arguments})
    deadline = time.monotonic()+180
    process = subprocess.Popen(arguments, cwd=directory, env=environment,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, encoding="utf-8", errors="replace", bufsize=1)
    events = queue.Queue()

    def receive():
        try:
            for line in process.stdout:
                events.put((time.perf_counter_ns(), line, process.poll() is None))
        finally:
            events.put(None)

    reader = threading.Thread(target=receive, daemon=True)
    reader.start()
    observations = []
    lines = []
    try:
        with (directory / "process.log").open("w", encoding="utf-8") as log:
            while True:
                require(time.monotonic() < deadline, "formal App invocation exceeded 180 seconds")
                try:
                    event = events.get(timeout=.1)
                except queue.Empty:
                    continue
                if event is None:
                    break
                received, line, live = event
                lines.append(line)
                observations.append({"received_perf_counter_ns": received,
                                     "process_live_when_received": live, "text": line.rstrip()})
                log.write(line)
                log.flush()
                if re.match(r"\[Step 1\]", line):
                    marker(directory / "first_accepted_observed.json",
                           {"received_perf_counter_ns": received, "process_live": live})
            code = process.wait(timeout=max(1., deadline-time.monotonic()))
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=10)
        reader.join(timeout=2)
        process.stdout.close()
    finished = marker(directory / "process_finished.json", {"returncode": code})
    receipt = {"arguments": arguments, "returncode": code,
               "start_file_birthtime_ns": started, "finish_file_birthtime_ns": finished,
               "file_creation_delta_seconds": (finished-started)/1e9,
               "stream_observations": observations}
    marker(directory / "process_receipt.json", receipt)
    return code, "".join(lines), receipt


def near(actual, expected, message):
    require(math.isfinite(actual) and abs(actual-expected) <= 5e-13*(1+abs(expected)), message)


def rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def output_path(directory, reference):
    path = (directory / "output" / reference).resolve()
    require(path.is_relative_to((directory / "output").resolve()) and path.is_file(),
            "output receipt references a missing file or escapes output: "+reference)
    return path


def verify_output(directory, hydro_parser, electric_parser):
    meta = directory / "output/meta"
    manifest = json.loads((meta / "manifest.json").read_text(encoding="utf-8"))
    require(manifest["status"] == "COMPLETED" and manifest["accepted_step"] == 2,
            "joint App must publish the actual two accepted steps")
    near(float(manifest["accepted_time"]), .0005, "joint App accepted End time")
    require("build_type=Release" in (meta / "build_info.txt").read_text(encoding="utf-8"),
            "formal performance evidence must use the Release executable")
    require(manifest["capillary_method"] == "balancedCsf"
            and manifest["density_representation"] == "sharpGeometry",
            "joint output claims a different consumed spatial route")
    require(manifest["complete_frames"] == 2, "sparse output must contain only Initial and final End")
    require(set(manifest["input_hashes"]) == {"system/"+name for name in INPUTS},
            "input provenance must cover exactly all eight owners")
    for relative, claimed in manifest["input_hashes"].items():
        require(claimed == "sha256:"+sha(directory / relative), "source SHA differs from manifest")
    effective = meta / "effective"
    require({path.name for path in effective.iterdir()} == set(INPUTS), "effective configuration is not eight files")
    trace = rows(meta / "configuration_trace.csv")
    require(trace and all(row["used"] == "true" for row in trace), "accepted configuration has unused options")
    require({Path(row["file"]).name for row in trace} == set(INPUTS), "trace lost an input owner")
    thread = [row for row in trace if row["field"] == "Parallel.threads"]
    require(len(thread) == 1 and thread[0]["raw_value"] == "8"
            and thread[0]["effective_value"] == "2", "trace lost actual CLI thread override")
    require(re.search(r"Parallel\s*\{\s*threads\s+2\s*;", (effective / "run.dict").read_text(encoding="utf-8")),
            "effective input lost threads override")
    require("interfaceMotion coupled;" in (effective / "electric.dict").read_text(encoding="utf-8"),
            "effective input lost coupled motion")
    frames = rows(meta / "frames.csv")
    require([int(row["step"]) for row in frames] == [0, 2], "cadence two must output steps 0/2, never step 1")
    geometries = {}
    for frame in frames:
        step = int(frame["step"])
        geometry = int(frame["geometry_identity"])
        require(geometry > 0, "missing real geometry identity")
        geometries[step] = geometry
        near(float(frame["time"]), step*.00025, "frame time differs from accepted End")
        require("native_faces" not in frame, "retired native face output must not remain in frames.csv")
        for group, variables in (("interface", ["R", "Z", "NR", "NZ", "Kappa"]),):
            path = output_path(directory, frame[group])
            require(path.stem.endswith("_"+format(step, "06d")), "field filename must encode accepted step")
            zone = hydro_parser.tecplot_data(path, variables, group)
            require(zone["step"] == step and zone["geometry"] == geometry, "Hydro/interface field source mismatch")
            near(zone["time"], step*.00025, "Hydro/interface physical clock")
    require([int(row["step"]) for row in rows(directory / "output/series/conservation.csv")] == [0, 1, 2],
            "sparse field output erased actual physical history")
    electric = verify_electric_output(directory, frames, geometries, electric_parser)
    require(electric["driven_native_velocity"], "sigma=0 electric never changed native velocity")
    return {"accepted_steps": 2, "frame_steps": [0, 2], "geometry_identity_by_step": geometries,
            "electric": electric, "effective": str(effective)}


def verify_electric_output(directory, frames, geometries, electric_parser):
    manifest = json.loads((directory / "output/meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest["electric"] == "HHO" and manifest["electric_interface_motion"] == "moving",
            "manifest did not identify the actual native moving electric path")
    history = rows(directory / "output/series/electric.csv")
    require([int(row["step"]) for row in history] == [0, 1, 2], "joint End charge history lost accepted steps")
    require([row["time_update_present"] for row in history] == ["false", "true", "true"],
            "Initial placeholders were presented as a time equation")
    require(history[0]["temporal_order"] == "NA" and [int(row["temporal_order"]) for row in history[1:]] == [1, 2],
            "joint charge did not use actual BE startup followed by accepted-history BDF2")
    for index, row in enumerate(history):
        near(float(row["time"]), index*.00025, "electric accepted physical clock")
        require(int(row["geometry_identity"]) > 0, "electric accepted geometry identity missing")
        if index in geometries:
            require(int(row["geometry_identity"]) == geometries[index], "electric/Hydro accepted geometry mismatch")
        for field, limit in (("linear_residual", 1e-11), ("local_Gauss_relative", 1e-10),
                             ("global_Gauss_relative", 1e-10), ("current_balance_relative", 1e-10),
                             ("Maxwell_charge_jump_relative", 1e-10)):
            value = float(row[field])
            require(math.isfinite(value) and 0 <= value <= limit, "original public "+field+" gate failed")
        if index:
            near(float(row["dt"]), .00025, "actual charge physical dt")
            near(float(row["previous_dt"]), 0 if index == 1 else .00025, "actual BDF2 previous accepted dt")
            for field, value in zip(("a0", "a1", "a2"), (1., -1., 0.) if index == 1 else (1.5, -2., .5)):
                near(float(row[field]), value, "charge coefficients did not come from accepted physical history")
            for field, limit in (("Schur_relative", 1e-11), ("panel_equation_relative", 1e-10),
                                 ("charge_equation_relative", 1e-10)):
                value = float(row[field])
                require(math.isfinite(value) and 0 <= value <= limit, "original public "+field+" gate failed")
        else:
            require(all(row[field] == "NA" for field in ("dt", "previous_dt", "a0", "a1", "a2",
                    "Schur_iterations", "Schur_relative", "panel_equation_relative",
                    "charge_equation", "charge_equation_relative")),
                    "Initial electric output fabricated unavailable time equation data")
    phase_variables = ["R", "Z", "V", "Er", "Ez"]
    interface_variables = ["R", "Z", "NR", "NZ", "V_Gamma", "q_h", "Et", "Dn_liquid", "Dn_gas",
                           "Maxwell_normal", "Maxwell_tangent", "Maxwell_r", "Maxwell_z"]
    native_velocity = {int(row["step"]): float(row["native_max_speed"])
                       for row in rows(directory / "output/series/conservation.csv")}
    require(all(math.isfinite(value) and value >= 0 for value in native_velocity.values()),
            "actual native speed diagnostics must remain finite and nonnegative")
    require(not any((directory / "output/fields").rglob("native_faces_*.dat")),
            "retired native face field must not be written")
    signed_initial = False
    electric_fields = []
    for frame in frames:
        step = int(frame["step"])
        for kind in ("liquid", "gas", "interface"):
            column = "electric_"+kind
            path = output_path(directory, frame[column])
            require(path.name == column+"_"+format(step, "06d")+".dat", "electric filename lost actual accepted step")
            zone = electric_parser.tecplot(path, interface_variables if kind == "interface" else phase_variables)
            auxiliary = zone["auxiliary"]
            require(auxiliary["Method"] == "HHO" and auxiliary["InterfaceMotion"] == "moving",
                    "electric field claims a different path")
            require(int(auxiliary["AcceptedStep"]) == step
                    and int(auxiliary["GeometryIdentity"]) == geometries[step],
                    "electric field is not coeval with accepted Hydro/interface")
            require(auxiliary["ClockPosition"] == ("Initial" if step == 0 else "End"),
                    "candidate Midpoint was output as accepted End")
            near(float(auxiliary["ClockBegin"]), 0 if step == 0 else (step-1)*.00025, "electric accepted interval begin")
            near(float(auxiliary["ClockEnd"]), step*.00025, "electric accepted interval end")
            near(zone["time"], step*.00025, "electric field physical solution time")
            if kind != "interface":
                epsilon, conductivity, phase = (5., .7, 1) if kind == "liquid" else (1., .03, 2)
                require(float(auxiliary["Epsilon"]) == epsilon
                        and float(auxiliary["Conductivity"]) == conductivity
                        and int(auxiliary["PhaseId"]) == phase, "electric output mixed or ignored material coefficients")
                require(any(value != 0 for field in ("Er", "Ez") for value in zone["blocks"][field]),
                        "configured nonzero boundary E was ignored")
            elif step == 0:
                charge = zone["blocks"]["q_h"]
                signed_initial = min(charge) < 0 < max(charge) and len(set(charge)) > 1
            electric_fields.append(path.relative_to(directory / "output").as_posix())
    require(signed_initial, "initial charge output lost configured signed nonuniform charge")
    require(native_velocity[0] == 0, "zero initial native velocity was not consumed")
    actual_fields = {path.relative_to(directory / "output").as_posix()
                     for path in (directory / "output/fields").rglob("electric_*.dat")}
    require(actual_fields == set(electric_fields), "electric output wrote a candidate/unscheduled frame or lost an accepted frame")
    energy = rows(directory / "output/series/energy.csv")
    require([int(row["step"]) for row in energy] == [0, 1, 2], "accepted energy history missing")
    require(all(row["interface_power_available"] == "true" and math.isfinite(float(row["interface_power"]))
                and row["capillary_power_available"] == "false" and row["capillary_power"] == "NA"
                for row in energy[1:]),
            "joint output mislabeled combined interface work as pure capillary work")
    return {"temporal_orders": [1, 2], "signed_initial_charge": signed_initial,
            "native_velocity_max_by_step": native_velocity,
            "driven_native_velocity": native_velocity[2] > 0, "electric_fields": electric_fields}


def main():
    require(len(sys.argv) == 4, "usage: hydro_electric_app_contract.py LAMP.exe projectRoot workRoot")
    executable, project, work = (Path(argument).resolve() for argument in sys.argv[1:])
    require(executable.is_file() and work.drive.upper() == "D:", "formal executable and D: workRoot are required")
    task = work / ("hydro_electric_"+datetime.now().strftime("%Y%m%d_%H%M%S_%f")+"_"+uuid.uuid4().hex[:8])
    (task / "_temp").mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    environment.update(TEMP=str(task / "_temp"), TMP=str(task / "_temp"), OMP_NUM_THREADS="2",
                       OMP_THREAD_LIMIT="2", OMP_DYNAMIC="false", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    hydro_parser = helper(project / "tests/contracts/hydrodynamics/hydro_app_contract.py", "hydro_app_parser")
    electric_parser = helper(project / "tests/contracts/electric/electric_app_contract.py", "electric_app_parser")
    before = source_snapshot(project)
    binary = sha(executable)
    result = {"status": "RUNNING", "scope": "formal_native_RZ_joint_AC_App_wiring",
              "task": str(task), "source_sha256": before, "executable_sha256": binary,
              "actual_numerical_windows": 0, "checks": []}

    def fresh(name):
        directory = task / name
        (directory / "system").mkdir(parents=True)
        for filename, content in INPUTS.items():
            (directory / "system" / filename).write_text(content, encoding="utf-8")
        return directory

    def preflight(name, mutate=None, expected_file=None, expected_field=None):
        directory = fresh(name)
        if mutate:
            mutate(directory)
        code, text, process = invoke(executable, directory, environment, True)
        if mutate:
            require(code != 0, name+" silently accepted invalid input")
            require(expected_file in text, name+" omitted failing input owner")
            if expected_field:
                require(re.search(re.escape(expected_file)+r":\d+:", text), name+" omitted true file:line")
                require(expected_field in text, name+" omitted failing field")
        else:
            require(code == 0 and "HYDRO_ELECTRIC_INPUT_OK" in text, "valid eight-file preflight failed")
            require(re.search(r"\[Runtime\].*\bthreads=2\b.*\bsource=command_line\b", text),
                    "read-only preflight ignored CLI thread override")
        require(not (directory / "output").exists() and not (directory / "cache").exists()
                and "[INIT]" not in text and "[Step " not in text,
                name+" read-only preflight constructed mesh, solver or output")
        result["checks"].append({"name": name, "status": "PASS", "process": process})

    try:
        preflight("valid_eight_files")
        preflight("missing_electric_file", lambda d: (d / "system/electric.dict").unlink(), "electric.dict")
        preflight("obsolete_extra_file", lambda d: (d / "system/hydro_fluid.dict").write_text("shape circle;\n", encoding="utf-8"),
                  "hydro_fluid.dict")
        preflight("unknown_electric_field", lambda d: replace(d, "electric.dict", "interfaceMotion coupled;",
                  "interfaceMotion coupled; unknownKey true;"), "electric.dict", "unknownKey")
        preflight("zero_permittivity", lambda d: replace(d, "materials.dict", "permittivity 5;", "permittivity 0;"),
                  "materials.dict", "permittivity")
        preflight("negative_conductivity", lambda d: replace(d, "materials.dict", "conductivity .7;", "conductivity -.7;"),
                  "materials.dict", "conductivity")
        preflight("wrong_axis", lambda d: replace(d, "electric.dict", 'W { type Neumann; value "0"; }',
                  'W { type Dirichlet; value "0"; }'), "electric.dict", "W")
        preflight("unanchored", lambda d: replace(d, "electric.dict",
                  'S { type Dirichlet; value "0"; }\n        E { type Dirichlet; value "electricStrength*z"; }\n        N { type Dirichlet; value "electricStrength"; }',
                  'S { type Neumann; value "0"; }\n        E { type Neumann; value "0"; }\n        N { type Neumann; value "0"; }'),
                  "electric.dict", "Boundary")
        preflight("sharp_external_route", lambda d: replace(d, "hydro.dict", "capillaryModel balancedCsf;",
                  "capillaryModel sharpNormalStress;"), "hydro.dict", "capillaryModel")
        preflight("unsupported_xy_joint", lambda d: replace(d, "geometry.dict", "coordinates rz;", "coordinates xy;"),
                  "geometry.dict", "coordinates")
        preflight("fixed_in_joint", lambda d: replace(d, "electric.dict", "interfaceMotion coupled;",
                  "interfaceMotion fixed;"), "electric.dict", "interfaceMotion")

        def disabled(directory):
            replace(directory, "run.dict", "app Electrohydrodynamics;", "app Hydrodynamics;")
            (directory / "system/electric.dict").unlink()
        preflight("electric_material_without_coupling", disabled, "materials.dict", "permittivity")
        preflight("retired_coupling_switch", lambda d: replace(d, "run.dict", "Parallel { threads 8; }",
                  "Coupling { electric electrostatic; }\nParallel { threads 8; }"), "run.dict", "Coupling")
        preflight("fluid_app_with_electric_file", lambda d: replace(d, "run.dict", "app Electrohydrodynamics;",
                  "app Hydrodynamics;"), "electric.dict")
        directory = fresh("formal_two_step_window")
        result["actual_numerical_windows"] = 1
        code, text, process = invoke(executable, directory, environment, False)
        require(code == 0 and re.search(r"\[COMPLETED\].*\bsteps=2\b.*\bcoupling=sharp-interface electrostatics\b", text),
                "formal joint two-step window failed: "+text[-6000:])
        require("OMP: Info" not in text and "H14P" not in text, "console leaked deprecated runtime/research output")
        progress = process["stream_observations"]
        for pattern in (r"\[INIT\]", r"\[ADVANCE\].*step=1", r"\[Step 1\]"):
            require(any(re.match(pattern, row["text"]) and row["process_live_when_received"] for row in progress),
                    "progress was not flushed while App was running: "+pattern)
        summaries = [row["text"] for row in progress if re.match(r"\[Step \d+\]", row["text"])]
        require([int(re.search(r"\[Step (\d+)\]", line).group(1)) for line in summaries] == [0, 1, 2],
                "logEvery=2 must still show first accepted progress and final step")
        output = verify_output(directory, hydro_parser, electric_parser)
        result["checks"].append({"name": "formal_two_step_window", "status": "PASS", "process": process, "output": output})
        readback = task / "effective_readback"
        shutil.copytree(Path(output["effective"]), readback / "system")
        code, text, process = invoke(executable, readback, environment, True)
        require(code == 0 and "HYDRO_ELECTRIC_INPUT_OK" in text and "[INIT]" not in text
                and not (readback / "output").exists() and not (readback / "cache").exists(),
                "saved effective eight-file input failed read-only roundtrip")
        result["checks"].append({"name": "effective_readback", "status": "PASS", "process": process})
        result["source_after_sha256"] = source_snapshot(project)
        result["executable_after_sha256"] = sha(executable)
        require(result["source_after_sha256"] == before and result["executable_after_sha256"] == binary,
                "sources or executable changed during the frozen check")
        result["status"] = "PASS_WITHIN_FORMAL_NATIVE_RZ_JOINT_AC_APP_SCOPE"
        result["not_proven"] = ["spatial_temporal_convergence", "coarse_P5_Maxwell_accuracy", "curvature_stability",
                               "AMR_small_cut", "electrospray_nozzle", "large_mesh_performance", "electrospray_completion"]
    except Exception as error:
        result["status"] = "FAIL"
        result["error"] = str(error)
        raise
    finally:
        marker(task / "receipt.json", result)
        print("Joint Hydro/Electric App receipt: "+str(task / "receipt.json"), flush=True)
    print(result["status"], flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("Joint Hydro/Electric App FAIL: "+str(error), file=sys.stderr, flush=True)
        sys.exit(1)
