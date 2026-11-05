#!/usr/bin/env python3
"""Exercise the real Windows launcher with immutable restart fixtures on D:.

Only PrepareOnly is used: the EXE is a snapshot marker, not another solver.
Runtime checkpoint byte/schema validation belongs to the C++ restart contract.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


def require(ok, reason):
    if not ok:
        raise AssertionError(reason)


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def tree(path):
    return {p.relative_to(path).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*") if p.is_file()}


def main():
    require(os.name == "nt", "restart launcher contract requires Windows PowerShell")
    source = Path(sys.argv[1]).resolve()
    requested = Path(sys.argv[2]).resolve()
    require(requested.drive.upper() == "D:", "launcher contract work root must be on D:")
    root = requested / ("launcher_" + str(time.time_ns()))
    root.mkdir(parents=True)
    project = root / "project with spaces"
    (project / "tools").mkdir(parents=True)
    for name in ("run_case.ps1", "run_case.cmd"):
        shutil.copy2(source / "tools" / name, project / "tools" / name)
    write(project / "bin" / "LAMP.exe", "current runtime snapshot marker\n")
    write(project / "bin" / "runtime.dll", "current DLL marker\n")
    case = "restart_case"
    write(project / "cases" / case / "system" / "run.dict", "CURRENT CASE MUST NOT REPLACE FROZEN INPUT\n")
    old = project / "runs" / case / "old run 原始"
    frozen_run = '''Run { app Hydrodynamics; }
Time { units code; dt .01; endTime /*not a Core comment*/ .2; }
Output { Fields { enabled false; } }
# Restart { directory "D:/comment-is-preserved"; }
// Time { endTime 999; } must remain a comment.
Restart // the prior continuation block is replaced, not duplicated
{ directory "D:/older-parent/checkpoint_00000001"; }
Log { everySteps 1; level summary; }
'''.replace("/*not a Core comment*/ ", "")
    write(old / "system" / "run.dict", frozen_run)
    write(old / "system" / "geometry.dict", "frozen geometry bytes\n")
    write(old / "input" / "sample.dat", "frozen input bytes\n")
    write(old / "runtime" / "LAMP.exe", "old runtime must not be selected by default\n")
    write(old / "run.json", json.dumps({"case_name": case, "status": "failed"}))
    checkpoint_root = old / "output" / "checkpoints"
    first = checkpoint_root / "checkpoint_00000010"
    latest = checkpoint_root / "checkpoint_00000020"
    for directory in (first, latest):
        write(directory / "manifest.txt", "LAMP_VOF_CHECKPOINT 1\nfixture validation belongs to C++\n")
    write(checkpoint_root / "latest.txt", latest.name + "\n")
    temp = root / "temp"
    temp.mkdir()
    env = dict(os.environ, TEMP=str(temp), TMP=str(temp), PYTHONDONTWRITEBYTECODE="1")
    checks = 0
    results = []

    def launch(arguments, success=True, reason="", user_input=None, cancelled=False):
        nonlocal checks
        before = tree(old)
        before_runs = sorted(p.name for p in (project / "runs" / case).iterdir())
        command = ["powershell.exe", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
                   str(project / "tools" / "run_case.ps1"), case] + list(arguments)
        done = subprocess.run(command, cwd=project, env=env, input=user_input, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, encoding="utf-8", errors="replace", timeout=45)
        results.append({"arguments": list(arguments), "returncode": done.returncode, "output": done.stdout})
        require((done.returncode == 0) == success, "unexpected launcher status: " + done.stdout)
        if reason:
            require(reason in done.stdout, "launcher rejection lost reason: " + done.stdout)
        require(tree(old) == before, "restart launcher changed the old failed run")
        checks += 1
        if cancelled or not success:
            require(sorted(p.name for p in (project / "runs" / case).iterdir()) == before_runs,
                    "cancelled or rejected selection created a run")
        if cancelled:
            require("Restart cancelled" in done.stdout, "interactive cancellation lost its result")
            return None
        if success:
            newest = Path((project / "runs" / case / "latest.txt").read_text(encoding="utf-8"))
            record = json.loads((newest / "run.json").read_text(encoding="utf-8-sig"))
            require(record["status"] == "prepared", "PrepareOnly entered a physical/check run")
            require(not (newest / "output").exists(), "PrepareOnly created physical output")
            require(Path(record["input_source"]) == old, "restart did not use the original frozen input root")
            require((newest / "system" / "geometry.dict").read_bytes() == (old / "system" / "geometry.dict").read_bytes(),
                    "frozen geometry changed")
            require((newest / "input" / "sample.dat").read_bytes() == (old / "input" / "sample.dat").read_bytes(),
                    "frozen ancillary input changed")
            require((newest / "runtime" / "LAMP.exe").read_bytes() == (project / "bin" / "LAMP.exe").read_bytes(),
                    "default restart program did not use current bin")
            return newest, record
        return None

    run, record = launch(["-Restart", os.path.relpath(old, project), "-PrepareOnly"])
    require(Path(record["restart"]["source_run"]) == old and Path(record["restart"]["checkpoint_directory"]) == latest,
            "old-run restart did not select the latest checkpoint")
    text = (run / "system" / "run.dict").read_text(encoding="utf-8-sig")
    require('Restart { directory "' + latest.as_posix() + '"; }' in text
            and "D:/older-parent/checkpoint_00000001" not in text,
            "actual prior Restart block was not replaced")
    require('# Restart { directory "D:/comment-is-preserved"; }' in text
            and "// Time { endTime 999; } must remain a comment." in text and "endTime .2;" in text,
            "restart altered comments or the frozen endTime without an override")
    require("end_time_override" not in record["restart"], "omitted EndTime manufactured an override")

    run, record = launch(["-Restart", str(first), "-EndTime", ".4", "-PrepareOnly"])
    require(record["restart"]["explicit_checkpoint"] and Path(record["restart"]["checkpoint_directory"]) == first,
            "explicit older checkpoint selection was ignored")
    text = (run / "system" / "run.dict").read_text(encoding="utf-8-sig")
    require(re.search(r"Time \{ units code; dt \.01; endTime 0\.4; \}", text) is not None,
            "explicit EndTime did not replace actual Time.endTime")
    require("// Time { endTime 999; } must remain a comment." in text
            and record["restart"]["end_time_override"] == {"old": .2, "new": .4, "units": "code_time"},
            "EndTime override lost old/new provenance or changed comment text")

    run, record = launch(["-Resume", "-PrepareOnly"], user_input="1\n")
    require(Path(record["restart"]["source_run"]) == old and Path(record["restart"]["checkpoint_directory"]) == latest,
            "interactive single choice did not use its indexed checkpoint")
    require("Enter a number or saved run/checkpoint path" in results[-1]["output"], "interactive prompt is invisible")
    second = old.parent / "zz another indexed run"
    shutil.copytree(old, second)
    before_second = tree(second)
    run, record = launch(["-Resume", "-PrepareOnly"], user_input="2\n")
    require(Path(record["restart"]["source_run"]) == old and tree(second) == before_second,
            "interactive multi-choice selected or changed another saved run")
    run, record = launch(["-Resume", "-EndTime", ".4", "-PrepareOnly"], user_input='"' + str(first) + '"\n')
    require(Path(record["restart"]["checkpoint_directory"]) == first and record["restart"]["end_time_override"]["new"] == .4,
            "interactive quoted Unicode checkpoint path or endTime override failed")
    launch(["-Resume", "-PrepareOnly"], user_input="q\n", cancelled=True)
    for answer in ("", "\n"):
        launch(["-Resume", "-PrepareOnly"], False, "No restart source selected", user_input=answer)
    launch(["-Resume", "-PrepareOnly"], False, "outside the displayed range", user_input="0\n")
    launch(["-Resume", "-Restart", str(old), "-PrepareOnly"], False, "Use either -Resume or -Restart")
    help_result = subprocess.run(["powershell.exe", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
                                 str(project / "tools" / "run_case.ps1"), "-Help"], cwd=project, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8", timeout=20)
    require(help_result.returncode == 0 and all(word in help_result.stdout for word in ("-Resume", "-Restart <path>", "Checkpoint", "everySteps")),
            "case-free help did not expose restart and checkpoint input")
    checks += 1

    # The launcher only selects the requested generation; the runtime owns
    # its integrity check. A present corrupt manifest must not select first.
    write(latest / "manifest.txt", "corrupt newest committed fixture\n")
    _, record = launch(["-Restart", str(old), "-PrepareOnly"])
    require(Path(record["restart"]["checkpoint_directory"]) == latest, "launcher silently rolled back a corrupt newest checkpoint")
    missing = checkpoint_root / "checkpoint_00000030"
    missing.mkdir()
    write(checkpoint_root / "latest.txt", missing.name)
    launch(["-Restart", str(old), "-PrepareOnly"], False, "manifest.txt is missing")
    write(checkpoint_root / "latest.txt", "../checkpoint_00000010")
    launch(["-Restart", str(old), "-PrepareOnly"], False, "only a checkpoint_<step>")
    write(checkpoint_root / "latest.txt", latest.name)
    launch(["-Restart", str(root), "-PrepareOnly"], False, "inside this project runs/")
    launch(["-Restart", "", "-PrepareOnly"], False, "nonempty saved run")
    for invalid in ("0", "-.1", "nan", "inf", "1e309"):
        launch(["-Restart", str(old), "-EndTime", invalid, "-PrepareOnly"], False, "finite positive code time")
    launch(["-EndTime", ".4", "-PrepareOnly"], False, "only supported together with -Restart")
    write(old / "run.json", json.dumps({"case_name": "a_different_case", "status": "failed"}))
    launch(["-Restart", str(old), "-PrepareOnly"], False, "case_name differs")
    write(old / "run.json", json.dumps({"case_name": case, "status": "failed"}))

    # Junction creation does not require the symlink privilege on Windows.
    link = project / "runs" / case / "linked-old"
    create = subprocess.run(["powershell.exe", "-NoProfile", "-Command",
                             "New-Item -ItemType Junction -Path '" + str(link).replace("'", "''")
                             + "' -Target '" + str(old).replace("'", "''") + "' | Out-Null"],
                            env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            encoding="utf-8", errors="replace", timeout=20)
    skipped = []
    if create.returncode == 0:
        launch(["-Restart", str(link), "-PrepareOnly"], False, "directory links are not supported")
    else:
        skipped.append({"test": "directory_junction_rejection", "reason": "host denied junction fixture creation",
                        "setup_output": create.stdout})
    write(root / "results.json", json.dumps({"checks": checks, "results": results, "skipped": skipped}, ensure_ascii=False, indent=2))
    print("VOF_RESTART_LAUNCHER_PASS checks=" + str(checks) + " skipped=" + str(len(skipped))
          + " newrun frozen_inputs latest explicit_older immutable_parent endTime path_guards")


if __name__ == "__main__":
    main()
