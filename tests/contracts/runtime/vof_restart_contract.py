#!/usr/bin/env python3
"""Actual short VOF restart trajectories; no launcher or plotting reconstruction."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
import struct
import subprocess
import sys
import threading
import time


def require(ok, reason):
    if not ok:
        raise AssertionError(reason)


class Archive:
    """Independent reader of the declared native binary64 checkpoint ABI."""
    def __init__(self, data, tag):
        require(len(data) >= 64, "truncated archive")
        require(hashlib.sha256(data[:-64]).hexdigest().encode() == data[-64:], "archive SHA256")
        self.data, self.end, self.at = data, len(data) - 64, 0
        require(self.text().decode() == tag, "archive kind/version " + tag)
        require(self.scalar("I") == 0x01020304, "little-endian archive required")
        require([self.scalar("I") for _ in range(6)] == [4, 4, 8, 8, 8, 53], "native binary64 ABI")
        require(self.scalar("i") == -1021 and self.scalar("i") == 1024, "native exponent ABI")

    def take(self, size):
        require(size >= 0 and self.at + size <= self.end, "archive bounds")
        out = self.data[self.at:self.at + size]
        self.at += size
        return out

    def scalar(self, code):
        return struct.unpack("<" + code, self.take(struct.calcsize("<" + code)))[0]

    def boolean(self):
        out = self.scalar("B")
        require(out in (0, 1), "archive boolean")
        return bool(out)

    def text(self):
        return self.take(self.scalar("Q"))

    def values(self, code):
        return self.take(self.scalar("Q") * struct.calcsize("<" + code))

    def clock(self, position_code="I"):
        return (self.scalar("q"), self.scalar("d"), self.scalar("d"), self.scalar(position_code))

    def finish(self):
        require(self.at == self.end, "archive trailing fields")


def publication(data):
    a = Archive(data, "LAMP.LocalPlicPublication.1")
    a.values("Q")
    radial = a.boolean()
    a.values("i"); a.values("i"); a.values("d")
    clock = a.clock()
    a.take(16); a.values("d")
    a.take(a.scalar("Q") * 16)
    cells = a.scalar("Q")
    stocks, rectangles = [], []
    for _ in range(cells):
        rectangles.append(struct.unpack("<4d", a.take(32)))
        a.take(24)
        stocks.append(a.take(16))
        a.take(16); a.boolean(); a.values("d")
        a.take(4 + 32); a.values("Q")
    # The entire sealed publication is compared too, including stored faces,
    # actual planes/measures, exact certificates and retained native support.
    return {"raw": data, "clock": clock, "stocks": b"".join(stocks),
            "cells": cells, "radial": radial, "rectangles": rectangles}


def hydro(data):
    a = Archive(data, "LAMP.HydroVofState.1")
    a.values("Q"); a.boolean()
    a.values("i"); a.values("i"); a.values("i")
    flow = a.text()
    current = publication(a.text())
    velocity, pressure = a.values("d"), a.values("d")
    has_history = a.boolean()
    previous = None
    if has_history:
        previous = publication(a.text())
        old_velocity, old_pressure = a.values("d"), a.values("d")
        require(previous["clock"][0] + 1 == current["clock"][0], "genuine immediate predecessor step")
        require(previous["clock"][2] == current["clock"][1], "genuine predecessor interval")
        require(len(old_velocity) == len(velocity) and len(old_pressure) == len(pressure), "history native layout")
    if a.boolean():
        a.text()
    a.finish()
    return {"raw": data, "publication": current, "velocity": velocity,
            "pressure": pressure, "flow": flow, "has_history": has_history, "previous": previous}


def observer(data, app):
    a = Archive(data, "HydrodynamicsVofOutput.v1" if app == "Hydrodynamics" else "VofElectrohydrodynamicsOutput.A.v1")
    clock = a.clock()
    raw_totals = a.take(24)
    totals = struct.unpack("<3d", raw_totals)
    frames = a.scalar("Q")
    a.finish()
    require(all(math.isfinite(x) for x in totals), "finite output reference totals")
    return {"raw": data, "clock": clock, "totals_raw": raw_totals,
            "initialVolume": totals[0], "bodyWork": totals[1], "viscousLoss": totals[2], "frames": frames}


def checkpoint(directory, app):
    data = (directory / "state.bin").read_bytes()
    a = Archive(data, "LAMP/VOF/accepted-checkpoint/1")
    require(a.text().decode() == app, "checkpoint App")
    compatibility = a.text()
    clock = a.clock("B")
    next_dt, without_retry, transaction = a.scalar("d"), a.scalar("i"), a.text()
    a.finish()
    t = Archive(transaction, "LAMP/" + app + "/transaction/1")
    fluid = hydro(t.text())
    charge, flow_means = b"", b""
    has_flow = False
    if app == "VofElectrohydrodynamics":
        charge = t.values("d")
        has_flow = t.boolean()
        if has_flow:
            flow_means = t.values("d")
    output = observer(t.text(), app)
    t.finish()
    require(clock == fluid["publication"]["clock"] == output["clock"], "one complete physical/controller/output clock")
    return {"raw": data, "clock": clock, "next_dt": next_dt, "without_retry": without_retry,
            "compatibility": compatibility, "fluid": fluid, "charge": charge, "has_flow": has_flow,
            "flow_means": flow_means, "observer": output}


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


DT, STEPS, SPLIT = 1 / 1024, 6, 3


def fixture(directory, app, end_steps, restart=None, open_flow=False):
    rz = app == "VofElectrohydrodynamics"
    require(not open_flow or rz, "open restart scene belongs to A")
    run = f"""Run {{ app {app}; }}
Parallel {{ threads 1; }}
Time {{ units code; dt {DT:.17g}; endTime {end_steps * DT:.17g};
 Adaptive {{ minDt {DT:.17g}; maxDt {DT:.17g}; targetCourant .4; safety .9; maxGrowth 1.25; }}
 StepRetry {{ minDt {DT:.17g}; maxRetries 2; recoverAfter 0; }} }}
Output {{ Fields {{ enabled true; writeEverySteps 1; }} }}
Checkpoint {{ everySteps 1; }}
Log {{ everySteps 1; level summary; }}
"""
    if restart is not None:
        run += 'Restart { directory "' + restart.as_posix() + '"; }\n'
    geometry = """Geometry {
 coordinates COORDS; Domain { min (0 0); max (1 2); }
 Phases { negative 1; positive 2; }
 Initial { phi "FORMULA"; }
}
""".replace("COORDS", "rz" if rz else "xy").replace("FORMULA", "z-.8" if open_flow else "r^2+(z-.8)^2-.2^2" if rz else ".8+.02*cos(2*pi*x)-y")
    materials = """Materials { Phases {
 Negative { phaseId 1; density 3; dynamicViscosity .04; ELECTRIC_INNER }
 Positive { phaseId 2; density 1; dynamicViscosity .01; ELECTRIC_OUTER }
} Interfaces { Shared { phases (1 2); surfaceTension 0; } } }
""".replace("ELECTRIC_INNER", "permittivity 2; conductivity .1;" if rz else "").replace("ELECTRIC_OUTER", "permittivity 1; conductivity .05;" if rz else "")
    texts = {"run.dict": run, "geometry.dict": geometry, "materials.dict": materials,
             "mesh.dict": """Mesh { resolution (4 8); initialMesh uniform;
 Regrid { minLevel 0; maxLevel 1; interfaceBandLayers 1; vorticityTolerance .5; coarsenRatio .25; } }
""", "interface.dict": """Interface { method vof; Boundary { type none; }
 Initialization { subcellsPerAxis 4; } Reconstruction { relativeVolumeTolerance 1e-12; maximumVolumeIterations 80; } }
""", "hydro.dict": """Hydrodynamics { method staggered; densityRepresentation cellAverageProjection;
 capillaryModel heightFunction; Initial { velocity (0 0); } Boundary { type impermeable; }
 BodyForce { acceleration (0 -.1); } Numerics { maximumLocalCourant .5; relativeConservationTolerance 1e-12; }
 TimeIntegration { method adaptiveLinear; interfaceMismatchTolerance .01; maximumCorrectionRelative .25;
 correctionReferenceVelocity 1e-8; capillaryCorrectionCourant .1; maximumCapillaryCourant .5; } }
"""}
    if rz:
        texts["electric.dict"] = """Electric { method vofFiniteVolume; Advection { scheme limitedLinearFct; INLETS }
 Initial { chargeDensity ".01*(z-.8)"; } Boundary {
 W { type Neumann; value 0; } S { type Dirichlet; value .2; }
 E { type Neumann; value 0; } N { type Dirichlet; value -.2; } }
 Solver { relativeLinearTolerance 1e-11; relativeConservationTolerance 1e-10; } }
""".replace("INLETS", "Inlets { feed { chargeDensity -.008; } }" if open_flow else "")
    if open_flow:
        texts["hydro.dict"] = texts["hydro.dict"].replace("Boundary { type impermeable; }", """Boundary {
 type prescribedFlow; Backflow { model ambient; velocity (0 0); } Patches {
 feed { side S; interval (0 1); type prescribedVelocity; phase liquid; velocity (0 .02); }
 drain { side N; interval (0 1); type pressureOutlet; phase gas; pressure 0; } } }""")
    for name, text in texts.items():
        write(directory / "system" / name, text)


def tree(path):
    return {p.relative_to(path).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*") if p.is_file()}


def latest(directory):
    root = directory / "output" / "checkpoints"
    name = (root / "latest.txt").read_text(encoding="utf-8").strip()
    require(re.fullmatch(r"checkpoint_\d+", name) is not None, "one committed latest generation")
    return root / name


def metadata(directory):
    return json.loads((directory / "output" / "meta" / "manifest.json").read_text(encoding="utf-8"))


def frame_steps(directory):
    with (directory / "output" / "meta" / "frames.csv").open(newline="", encoding="utf-8") as f:
        return [int(row["step"]) for row in csv.DictReader(f)]


def raw_equal(label, a, b):
    if a != b:
        first = next((i for i, values in enumerate(zip(a, b)) if values[0] != values[1]), min(len(a), len(b)))
        raise AssertionError(f"{label} differs at raw byte {first}, lengths {len(a)}/{len(b)}")


def main():
    require(len(sys.argv) == 3, "usage: vof_restart_contract.py LAMP.exe work_directory")
    executable, work = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    require(executable.is_file(), "actual executable missing")
    if os.name == "nt":
        require(work.drive.upper() == "D:", "contract work must stay on D:")
    root = work / ("runtime_" + str(time.time_ns()))
    root.mkdir(parents=True)
    temp = root / "temp"
    temp.mkdir()
    env = dict(os.environ, TEMP=str(temp), TMP=str(temp), PYTHONDONTWRITEBYTECODE="1",
               OMP_NUM_THREADS="1", OMP_DYNAMIC="false", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    evidence = {"root": str(root), "executable": str(executable), "scenes": [], "guards": [], "interrupted": []}

    def run(directory, success=True, reason=None, check_only=False):
        command = [str(executable), "--case", str(directory)]
        if check_only:
            command.append("--check-case")
        result = subprocess.run(command, cwd=directory, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                encoding="utf-8", errors="replace", timeout=180)
        write(directory / "contract_run.log", result.stdout)
        require((result.returncode == 0) == success, "unexpected runtime result: " + str(directory / "contract_run.log"))
        if reason is not None:
            require(reason.lower() in result.stdout.lower(), "missing rejection reason: " + result.stdout)
        if check_only:
            require(not (directory / "output").exists(), "check-only created physical output")
        if success and not check_only:
            require("[COMPLETED]" in result.stdout, "runtime did not complete")
        elif not success:
            require("[COMPLETED]" not in result.stdout and "[INIT]" not in result.stdout,
                    "invalid restart silently initialized or completed")
            evidence["guards"].append({"directory": str(directory), "returncode": result.returncode,
                                       "reason": reason, "check_only": check_only})
        return result.stdout

    def kill_at_checkpoint(directory):
        process = subprocess.Popen([str(executable), "--case", str(directory)], cwd=directory, env=env,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   encoding="utf-8", errors="replace", bufsize=1)
        lines, pending = [], queue.Queue()

        def capture():
            try:
                for line in process.stdout:
                    lines.append(line)
                    pending.put(line)
            finally:
                pending.put(None)

        reader = threading.Thread(target=capture, daemon=True)
        reader.start()
        stopped = False
        deadline = time.monotonic() + 180
        try:
            while time.monotonic() < deadline:
                line = pending.get(timeout=max(.01, deadline - time.monotonic()))
                require(line is not None, "solver exited before the requested real interruption")
                found = re.search(r"\[CHECKPOINT\] step=(\d+)\b", line)
                if found and int(found.group(1)) >= SPLIT:
                    require(process.poll() is None, "solver already exited before terminate")
                    process.terminate()
                    stopped = True
                    break
            require(stopped, "no accepted checkpoint available for real interruption")
            code = process.wait(timeout=30)
            reader.join(timeout=5)
            require(not reader.is_alive(), "interrupted stdout reader did not drain")
            log = "".join(lines)
            require(code != 0 and "[COMPLETED]" not in log, "requested interruption completed normally")
            write(directory / "contract_run.log", log)
            return code
        finally:
            if process.poll() is None:
                process.kill()
                process.wait(timeout=30)
            reader.join(timeout=5)
            process.stdout.close()

    try:
        for cadence in (-1, 1000001):
            invalid = root / ("reject_cadence_" + str(cadence))
            fixture(invalid, "Hydrodynamics", 1)
            path = invalid / "system" / "run.dict"
            write(path, path.read_text(encoding="utf-8").replace("Checkpoint { everySteps 1; }", "Checkpoint { everySteps " + str(cadence) + "; }"))
            run(invalid, False, "accepted-step interval", True)
        for restart in (False, True):
            invalid = root / ("reject_b1_restart" if restart else "reject_b1_checkpoint")
            fixture(invalid, "VofElectrohydrodynamics", STEPS, root / "unread_B1_checkpoint" if restart else None)
            write(invalid / "system" / "geometry.dict", 'Geometry { coordinates xy; Domain { min (0 0); max (1 2); } Phases { negative 1; positive 2; } Initial { phi "(x-.5)^2+(y-1)^2-.2^2"; } }\n')
            write(invalid / "system" / "mesh.dict", "Mesh { resolution (4 8); initialMesh uniform; }\n")
            electric = invalid / "system" / "electric.dict"
            text = re.sub(r"Advection\s*\{\s*scheme\s+limitedLinearFct;\s*\}", "", electric.read_text(encoding="utf-8"))
            write(electric, text.replace("vofFiniteVolume", "vofSharpSurface").replace("chargeDensity", "surfaceChargeDensity").replace("(z-.8)", "(y-.8)"))
            if restart:
                path = invalid / "system" / "run.dict"
                write(path, path.read_text(encoding="utf-8").replace("Checkpoint { everySteps 1; }", "Checkpoint { everySteps 0; }"))
            run(invalid, False, "B1 restart is not implemented" if restart else "0 for B1", True)
        disabled = root / "checkpoint_disabled"
        fixture(disabled, "Hydrodynamics", 1)
        path = disabled / "system" / "run.dict"
        write(path, path.read_text(encoding="utf-8").replace("Checkpoint { everySteps 1; }", "Checkpoint { everySteps 0; }"))
        run(disabled)
        require(not (disabled / "output" / "checkpoints").exists(), "zero cadence published checkpoint files")
        evidence["zero_cadence"] = {"accepted_steps": 1, "checkpoint_directory_absent": True}
        scenes = (("Hydrodynamics", False), ("VofElectrohydrodynamics", False),
                  ("VofElectrohydrodynamics", True))
        for app, open_flow in scenes:
            scene = root / (app + ("_open" if open_flow else ""))
            full, prefix, resumed = scene / "full", scene / "prefix", scene / "resumed"
            fixture(full, app, STEPS, open_flow=open_flow); fixture(prefix, app, SPLIT, open_flow=open_flow)
            run(full); run(prefix)
            old = latest(prefix)
            old_tree = tree(prefix)
            start = checkpoint(old, app)
            require(start["clock"][0] == SPLIT and start["clock"][2] == SPLIT * DT, "prefix actual accepted End")
            require(start["has_flow"] == open_flow, "prefix physical port capability")
            fixture(resumed, app, STEPS, old, open_flow)
            resumed_log = run(resumed)
            require("[RESTART]" in resumed_log and "[INIT]" not in resumed_log, "restart reran initial state")
            require(tree(prefix) == old_tree, "restart changed its immutable parent run")
            a, b = checkpoint(latest(full), app), checkpoint(latest(resumed), app)
            require(a["clock"] == b["clock"] and a["clock"][0] == STEPS and a["clock"][2] == STEPS * DT, "global final clock")
            require(a["next_dt"] == b["next_dt"] == DT and a["without_retry"] == b["without_retry"], "adaptive/retry continuation")
            raw_equal("phase stocks", a["fluid"]["publication"]["stocks"], b["fluid"]["publication"]["stocks"])
            for key in ("velocity", "pressure", "raw"):
                raw_equal("Hydro " + key, a["fluid"][key], b["fluid"][key])
            raw_equal("signed Q", a["charge"], b["charge"])
            require(a["has_flow"] == b["has_flow"] == open_flow, "restored actual flow capability")
            raw_equal("actual prescribed native normal means", a["flow_means"], b["flow_means"])
            if open_flow:
                for saved in (start, a, b):
                    means = struct.unpack("<" + str(len(saved["flow_means"]) // 8) + "d", saved["flow_means"])
                    require(all(math.isfinite(x) for x in means) and any(x != 0 for x in means),
                            "actual native prescribed port normals were not retained")
            raw_equal("initialVolume/bodyWork/viscousLoss", a["observer"]["totals_raw"], b["observer"]["totals_raw"])
            require(a["observer"]["frames"] == b["observer"]["frames"] == STEPS + 1, "cumulative field frame schedule")
            if not open_flow:
                require(a["observer"]["bodyWork"] != 0 and a["observer"]["viscousLoss"] > 0, "nontrivial cumulative work fixture")
            raw_equal("complete accepted checkpoint", a["raw"], b["raw"])
            require(frame_steps(full) == list(range(STEPS + 1)) and frame_steps(prefix) == list(range(SPLIT + 1))
                    and frame_steps(resumed) == list(range(SPLIT + 1, STEPS + 1)), "new segment duplicated or reset accepted frames")
            m = metadata(resumed)
            require(m["status"] == "COMPLETED" and m["accepted_step"] == STEPS and m["restored_checkpoint_step"] == SPLIT,
                    "resume metadata global step/provenance")
            require(m["frame_offset"] == SPLIT + 1 and m["cumulative_complete_frames"] == STEPS + 1, "resume metadata cumulative frame count")
            require(a["fluid"]["has_history"] and b["fluid"]["has_history"], "real multistep history not retained at stable final mesh")
            p = a["fluid"]["publication"]
            widths = {round(rect[2] - rect[0], 14) for rect in p["rectangles"]}
            require(p["cells"] > 32 and len(widths) >= 2, "actual graded AMR change was not exercised")
            evidence["scenes"].append({"app": app, "open_flow": open_flow, "final_step": STEPS, "prefix_step": SPLIT,
                "time": a["clock"][2], "cells": p["cells"], "widths": sorted(widths),
                "history": True, "raw_final_state_equal": True, "bodyWork": a["observer"]["bodyWork"],
                "viscousLoss": a["observer"]["viscousLoss"], "checkpoint_sha256": hashlib.sha256(a["raw"]).hexdigest()})
            if open_flow:
                continue
            bad = scene / "damaged_checkpoint"
            shutil.copytree(old, bad)
            state = bytearray((bad / "state.bin").read_bytes()); state[len(state) // 2] ^= 1
            (bad / "state.bin").write_bytes(state)
            bad_run = scene / "reject_corrupt"
            fixture(bad_run, app, STEPS, bad); run(bad_run, False, "mismatch")
            mismatch = scene / "reject_physics"
            fixture(mismatch, app, STEPS, old)
            material = mismatch / "system" / "materials.dict"
            material.write_text(material.read_text(encoding="utf-8").replace("density 3;", "density 3.1;"), encoding="utf-8")
            run(mismatch, False, "differ")
            require(tree(prefix) == old_tree, "negative restart changed parent checkpoint or outputs")
            if app == "Hydrodynamics":
                for name in ("mesh", "manifest_time"):
                    damaged = scene / ("damaged_" + name)
                    shutil.copytree(old, damaged)
                    if name == "mesh":
                        path = damaged / "mesh_cache.bin"
                        blob = bytearray(path.read_bytes()); blob[len(blob) // 2] ^= 1
                        path.write_bytes(blob)
                    else:
                        path = damaged / "manifest.txt"
                        original = path.read_bytes().decode("utf-8")
                        path.write_bytes(re.sub(r"(?m)^time=.*$", "time=9", original).encode("utf-8"))
                    rejected = scene / ("reject_" + name)
                    fixture(rejected, app, STEPS, damaged)
                    run(rejected, False, "mesh size/checksum mismatch" if name == "mesh" else "manifest/source time mismatch", True)
                require(tree(prefix) == old_tree, "check-only damage guards modified the original source")
                interrupted, continued = scene / "interrupted", scene / "continued_after_kill"
                fixture(interrupted, app, STEPS)
                exit_code = kill_at_checkpoint(interrupted)
                require(metadata(interrupted)["status"] == "RUNNING", "terminated source must retain RUNNING provenance")
                killed_checkpoint = latest(interrupted)
                saved = checkpoint(killed_checkpoint, app)
                saved_step = saved["clock"][0]
                require(SPLIT <= saved_step < STEPS, "real interruption did not leave an intermediate accepted checkpoint")
                frozen = tree(interrupted)
                fixture(continued, app, STEPS, killed_checkpoint)
                log = run(continued)
                require("[RESTART]" in log and "[INIT]" not in log, "interrupted continuation reran initialization")
                finished = checkpoint(latest(continued), app)
                raw_equal("actual kill/resume complete accepted checkpoint", a["raw"], finished["raw"])
                require(frame_steps(continued) == list(range(saved_step + 1, STEPS + 1)), "interrupted continuation reset field schedule")
                cm = metadata(continued)
                require(cm["restored_checkpoint_step"] == saved_step and cm["accepted_step"] == STEPS
                        and cm["frame_offset"] == saved["observer"]["frames"]
                        and cm["cumulative_complete_frames"] == STEPS + 1, "actual kill continuation provenance")
                require(tree(interrupted) == frozen, "actual kill continuation modified the interrupted source")
                evidence["interrupted"].append({"app": app, "exit_code": exit_code, "source_status": "RUNNING",
                    "restored_step": saved_step, "final_step": STEPS, "same_original_endTime": True,
                    "raw_final_state_equal": True, "parent_immutable": True})
        evidence["status"] = "PASS"
        print("VOF_RESTART_APP_PASS scenes=" + str(len(evidence["scenes"])) + " guards=" + str(len(evidence["guards"]))
              + " killed=" + str(len(evidence["interrupted"]))
              + " actual_Hydro_A_XY_RZ_graded_AMR open_native_normals raw_LG_Q_Up_history controller observers parent_immutable")
    except BaseException as error:
        evidence["status"] = "FAIL"
        evidence["failure"] = str(error)
        raise
    finally:
        write(root / "results.json", json.dumps(evidence, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
