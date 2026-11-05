#!/usr/bin/env python3
"""Verify the portable entry and Windows menu without running another solver.

The fixture launcher records delegation. Physical continuity belongs to the
real package acceptance, not this menu/asset-integrity contract.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def require(ok, reason):
    if not ok:
        raise AssertionError(reason)


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data, encoding="utf-8-sig" if path.suffix == ".ps1" else "utf-8")


def main():
    require(os.name == "nt", "portable entry contract requires Windows")
    source = Path(sys.argv[1]).resolve()
    work = Path(sys.argv[2]).resolve() / ("bootstrap_" + str(time.time_ns()))
    require(work.drive.upper() == "D:", "contract outputs must be on D:")
    package = work / "portable moved 中文 with spaces"
    package.mkdir(parents=True)
    for old, new in (("start.ps1", "Start.ps1"), ("start.cmd", "Start.cmd")):
        shutil.copy2(source / "tools/portable" / old, package / new)
    write(package / "bin/LAMP.exe", "test runtime marker\n")
    write(package / "bin/runtime.dll", "test runtime dependency marker\n")
    write(package / "tools/run_case.cmd", "@exit /b 0\n")
    write(package / "tools/run_case.ps1", r'''param(
 [string]$Case,[string]$Program,[string]$Restart,[string]$EndTime,
 [int]$Threads,[switch]$Resume,[switch]$CheckOnly)
$root=Split-Path -Parent $PSScriptRoot
if($Case -eq 'failing'){[Console]::Error.WriteLine('fixture solver rejected');exit 7}
$record=[ordered]@{case=$Case;program=$Program;restart=$Restart;end_time=$EndTime;
 threads=$Threads;resume=[bool]$Resume;check=[bool]$CheckOnly;
 caller_cwd=(Get-Location).Path;threads_explicit=$PSBoundParameters.ContainsKey('Threads')}
if($Resume){
 [Console]::Write('Fixture restart source (q cancels): ')
 $answer=[Console]::ReadLine()
 if($null -eq $answer -or $answer -eq 'q'){[Console]::WriteLine('Restart cancelled');exit 0}
 $record.restart=$answer
}
$runs=Join-Path $root 'runs'
[IO.Directory]::CreateDirectory($runs)|Out-Null
[IO.File]::WriteAllText((Join-Path $runs 'delegation.json'),($record|ConvertTo-Json))
exit 0
''')
    for case in ("alpha", "nested/液滴 case", "failing"):
        write(package / "cases" / case / "system/run.dict", "Time { dt .01; endTime .1; }\n")
    write(package / "cases/retired/system/run.dict", "retired fixture\n")
    write(package / "cases/retired/retired.json", "{}\n")
    entries = []
    for item in (package / "bin").iterdir():
        entries.append({"path": item.relative_to(package).as_posix(),
                        "sha256": hashlib.sha256(item.read_bytes()).hexdigest(), "bytes": item.stat().st_size})
    for name in ("run_case.cmd", "run_case.ps1"):
        item = package / "tools" / name
        entries.append({"path": "tools/" + name, "sha256": hashlib.sha256(item.read_bytes()).hexdigest(),
                        "bytes": item.stat().st_size})
    manifest = package / "portable_manifest.json"
    write(manifest, json.dumps({"schema_version": 1, "files": entries}))
    original_manifest = manifest.read_bytes()
    powershell = Path(os.environ["SystemRoot"]) / "System32/WindowsPowerShell/v1.0/powershell.exe"
    cmd = Path(os.environ["SystemRoot"]) / "System32/cmd.exe"
    temp = work / "temp"
    temp.mkdir()
    env = dict(os.environ, PATH=str(work / "empty path"), TEMP=str(temp), TMP=str(temp),
               PYTHONDONTWRITEBYTECODE="1")
    record_path = package / "runs/delegation.json"
    results = []

    def run(args=(), user_input=None, expected=0, reason=None, launched=True, via_cmd=False):
        if record_path.exists():
            record_path.unlink()
        if via_cmd:
            command = [str(cmd), "/d", "/c", str(package / "Start.cmd")] + list(args)
        else:
            command = [str(powershell), "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
                       str(package / "Start.ps1")] + list(args)
        done = subprocess.run(command, cwd=work, env=env, input=user_input if user_input is not None else "",
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              encoding="utf-8", errors="replace", timeout=30)
        results.append({"arguments": list(args), "input": user_input, "returncode": done.returncode,
                        "cmd": via_cmd, "output": done.stdout})
        require(done.returncode == expected, "unexpected portable status: " + done.stdout)
        if reason:
            require(reason in done.stdout, "portable entry lost failure/cancellation reason: " + done.stdout)
        require(record_path.exists() == launched, "portable menu launched unexpectedly or failed to delegate: " + done.stdout)
        return json.loads(record_path.read_text(encoding="utf-8-sig")) if launched else None

    record = run(["-Case", "alpha", "-Action", "new", "-NoPause"])
    require(record["program"] == str(package / "bin/LAMP.exe"), "runtime path is not relative to the moved package")
    require(not record["threads_explicit"], "default entry overwrote dictionary thread selection")
    require(not record["resume"] and not record["check"], "new action delegated incorrect mode")
    record = run(["-Case", "nested/液滴 case", "-Action", "check", "-Threads", "6"])
    require(record["check"] and record["threads"] == 6, "check action/thread override was lost")
    record = run(["-Case", "alpha", "-Action", "resume", "-Restart", "runs/alpha/old run", "-EndTime", ".2"])
    require(record["restart"] == str(package / "runs/alpha/old run") and record["end_time"] == ".2",
            "relative restart path or end time was not preserved")
    record = run(user_input="1\n1\n")
    require(record["case"] == "alpha" and not record["check"], "number menu selected the wrong case/action")
    record = run(user_input="nested/液滴 case\n3\n")
    require(record["case"] == "nested/液滴 case" and record["check"], "Unicode case menu lost selection")
    run(user_input="q\n", launched=False, reason="Cancelled")
    run(user_input="", launched=False, reason="Cancelled")
    run(user_input="1\n", launched=False, reason="Cancelled")
    run(user_input="1\nq\n", launched=False, reason="Cancelled")
    run(user_input="999\n", expected=2, launched=False, reason="outside the displayed range")
    run(user_input="1\n4\n", expected=2, launched=False, reason="Action selection")
    run(["-Case", "missing", "-Action", "new"], expected=2, launched=False, reason="Unknown packaged case")
    run(["-Case", "retired", "-Action", "new"], expected=2, launched=False, reason="Unknown packaged case")
    run(["-Case", "alpha", "-Action", "resume"], user_input="q\n", launched=False, reason="Restart cancelled")
    run(["-Case", "alpha", "-Action", "new", "-Restart", "runs/alpha/old"], expected=2,
        launched=False, reason="require -Action resume")
    run(["-Case", "failing", "-Action", "new"], expected=7, launched=False, reason="fixture solver rejected")
    # .cmd must find PowerShell with an empty PATH, preserve the child's exit
    # status and accept -NoPause without waiting for a console key.
    run(["-Case", "alpha", "-Action", "check", "-NoPause"], via_cmd=True)
    run(["-Case", "failing", "-Action", "new", "-NoPause"], expected=7, launched=False,
        via_cmd=True, reason="fixture solver rejected")
    write(package / "cases/alpha/system/run.dict", "user-edited physical dictionary\n")
    run(["-Case", "alpha", "-Action", "check"])
    dll = package / "bin/runtime.dll"
    original_dll = dll.read_bytes()
    dll.write_bytes(original_dll + b"changed")
    run(["-Case", "alpha", "-Action", "check"], expected=2, launched=False, reason="damaged or changed")
    dll.write_bytes(original_dll)
    dll.unlink()
    run(["-Case", "alpha", "-Action", "check"], expected=2, launched=False, reason="asset is missing")
    dll.write_bytes(original_dll)
    write(package / "bin/unlisted.dll", "unexpected loader input")
    run(["-Case", "alpha", "-Action", "check"], expected=2, launched=False, reason="unlisted executable asset")
    (package / "bin/unlisted.dll").unlink()
    manifest.unlink()
    run(["-Case", "alpha", "-Action", "check"], expected=2, launched=False, reason="manifest.json is missing")
    manifest.write_bytes(original_manifest)
    bad = json.loads(original_manifest)
    bad["files"] = [item for item in bad["files"] if item["path"] != "bin/LAMP.exe"]
    write(manifest, json.dumps(bad))
    run(["-Case", "alpha", "-Action", "check"], expected=2, launched=False, reason="no required executable asset")
    manifest.write_bytes(original_manifest)
    write(work / "verification.json", json.dumps({"checks": len(results), "results": results}, indent=2, ensure_ascii=False))
    print("PASS portable entry: %d checks; %s" % (len(results), work / "verification.json"))


if __name__ == "__main__":
    main()
