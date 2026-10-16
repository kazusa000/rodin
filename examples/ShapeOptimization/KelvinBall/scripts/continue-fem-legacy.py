#!/usr/bin/env python3
"""Resume a pinned FEM geometry on 28h, with a separate serial BEM checkpoint queue."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import time

EXAMPLE = Path(__file__).resolve().parents[1]
ROOT = EXAMPLE.parents[2]
NATIVE_COMMIT = "df602b5afb0504449c0822e0c796fe21b4b08f2a"
BEM_COMMIT = "12f38a9a33dbbb0a1b7c0d0affcfff91620b927a"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["output", "native-source", "build-dir", "build-record", "pilot", "verifier"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--native-commit", default=NATIVE_COMMIT)
    parser.add_argument("--restart-campaign", type=Path)
    parser.add_argument("--restart-step", type=int, default=5, help="Global fully evaluated state to load.")
    parser.add_argument("--updates", type=int, default=100)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("evaluate_fem_bem", EXAMPLE / "scripts/evaluate-fem-bem.py")
    evaluator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(evaluator)
    git, digest = evaluator.git, evaluator.digest
    if platform.node() != "28h" or args.updates < 1:
        raise ValueError("Require 28h and a positive finite update budget.")
    for source, commit in [(ROOT, args.commit), (args.native_source, args.native_commit),
                           (args.verifier, BEM_COMMIT)]:
        if git(source, "rev-parse", "HEAD") != commit or git(source, "status", "--porcelain", "--untracked-files=all"):
            raise ValueError(f"Require exact clean pinned source: {source}")
    changes = git(ROOT, "diff", "--name-only", args.native_commit, args.commit).splitlines()
    allowed = (str(EXAMPLE.relative_to(ROOT) / "scripts") + "/",
               str(EXAMPLE.relative_to(ROOT) / "reproduction") + "/")
    if any(not name.startswith(allowed) and name != str(EXAMPLE.relative_to(ROOT) / "README.md")
           for name in changes):
        raise ValueError("Reused binary requires unchanged C++ and build inputs.")
    pilot = json.loads((args.pilot / "campaign.json").read_text())
    if pilot["commit"] != NATIVE_COMMIT or pilot["real_updates"] != 5:
        raise ValueError("Require the preserved five-update historical FEM pilot.")
    if next(s for s in pilot["stages"] if s["name"] == "optimization")["exit_code"] != 0:
        raise ValueError("Pilot native optimization did not complete.")
    binary = args.build_dir / "examples/ShapeOptimization/KelvinBall/KelvinBall"
    build = json.loads(args.build_record.read_text())
    expected_binary = pilot["binary_sha256"] if args.native_commit == NATIVE_COMMIT else build.get("binary_sha256")
    if args.native_commit != NATIVE_COMMIT and (build.get("commit") != args.native_commit or build.get("exit_code") != 0):
        raise ValueError("A changed native source requires a successful build record for that exact commit.")
    if digest(binary) != expected_binary:
        raise ValueError("Native binary hash differs from its pinned build or pilot.")
    restart = args.restart_campaign or args.pilot
    restart_record = json.loads((restart / "campaign.json").read_text())
    local_step = args.restart_step - restart_record.get("restart_global_step", 0)
    with (restart / "optimization/kelvin-ball.csv").open() as stream:
        complete_steps = {int(row["iteration"]) for row in csv.DictReader(stream)
                          if None not in row and row.get("omega_z") not in (None, "")}
    if local_step < 0 or local_step not in complete_steps:
        raise ValueError("Restart must identify a fully evaluated native CSV state.")
    if args.output.exists() or args.output.resolve().is_relative_to(ROOT):
        raise ValueError("Require a fresh output directory outside source.")
    if shutil.disk_usage(args.output.parent).free < 30 * 1024**3:
        raise ValueError("Require at least 30 GiB of free checkpoint storage.")
    preset = json.loads((EXAMPLE / "reproduction/fem-domain-L15-legacy.json").read_text())
    # The native target_volume is a reporting reference only; its derivative
    # constraint remains D volume[theta]=0 after a geometry restart.
    record = dict(controller_commit=args.commit, commit=args.native_commit, source=str(ROOT),
        native_source=str(args.native_source), binary=str(binary), binary_sha256=digest(binary),
        host=platform.node(), pilot=str(args.pilot), restart_campaign=str(restart),
        restart_local_step=local_step, restart_global_step=args.restart_step,
        requested_real_updates=args.updates, requested_final_global_step=args.restart_step+args.updates,
        preset=preset, original_target_volume=float(pilot["fem_history"][0]["target_volume"]),
        real_updates=0, state="validated", stages=[], bem=[], started_at=time.time(),
        fem_threads=6, fem_blas_threads=1, bem_blas_threads=2,
        optimization_timeout_seconds=259200, bem_checkpoint_timeout_seconds=43320)
    if args.validate_only:
        print(json.dumps(record, indent=2))
        return 0
    args.output.mkdir(parents=True)
    optimization = args.output / "optimization"
    optimization.mkdir()
    checkpoint = args.output / "restart-original-double.mesh"
    record["restart_input"] = evaluator.prepare(
        restart / f"optimization/KelvinBall.Chamber.mesh.{local_step:06d}.h5",
        restart / f"optimization/KelvinBallState-{local_step:06d}.mesh", checkpoint)
    parameters = [f"--iterations={args.updates+1}" if x.startswith("--iterations=") else x
                  for x in preset["optimization_arguments"]]
    command = [str(binary), *parameters, f"--mesh={checkpoint}"]
    env = dict(os.environ, OMP_NUM_THREADS="6", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               KELVIN_SYMMETRY="cubic", PYTHONDONTWRITEBYTECODE="1",
               LD_LIBRARY_PATH=build["runtime_library_path"])

    def save():
        temporary = args.output / "campaign.json.tmp"
        temporary.write_text(json.dumps(record, indent=2) + "\n")
        temporary.replace(args.output / "campaign.json")

    def terminate(child):
        os.killpg(child.pid, signal.SIGTERM)
        try:
            child.wait(timeout=60)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()

    checkpoints = sorted({step-args.restart_step for step in [10, 20, 40, 60, 80, 105]
                          if args.restart_step < step <= args.restart_step+args.updates} | {args.updates})
    stage = dict(name="optimization", command=command, started_at=time.time())
    record["stages"].append(stage)
    record.update(state="optimizing", local_bem_checkpoints=checkpoints)
    save()
    fem = bem = None
    code = 1
    try:
        with (optimization / "run.log").open("w") as log:
            fem = subprocess.Popen(command, cwd=optimization, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        stage["pid"] = fem.pid
        save()
        while True:
            history = optimization / "kelvin-ball.csv"
            rows = []
            if history.exists():
                with history.open() as stream:
                    rows = [row for row in csv.DictReader(stream)
                            if None not in row and row.get("omega_z") not in (None, "")]
            if rows:
                record.update(real_updates=int(rows[-1]["iteration"]), latest_fem=rows[-1])
            if fem.poll() is None and time.time()-stage["started_at"] > record["optimization_timeout_seconds"]:
                stage["timed_out"] = True
                terminate(fem)
            if "exit_code" not in stage and fem.poll() is not None:
                stage.update(exit_code=fem.returncode, finished_at=time.time())
                (optimization / "exit-code.txt").write_text(str(fem.returncode) + "\n")
                record["state"] = "independent_bem" if fem.returncode == 0 else "fem_failed_evaluating_saved_states"
            if bem is not None and (bem.poll() is not None or
                    time.time()-record["bem"][-1]["started_at"] > record["bem_checkpoint_timeout_seconds"]):
                item = record["bem"][-1]
                if bem.poll() is None:
                    item["timed_out"] = True
                    terminate(bem)
                item.update(exit_code=bem.returncode, finished_at=time.time())
                (Path(item["output"]).parent / "exit-code.txt").write_text(str(bem.returncode) + "\n")
                result_path = Path(item["output"]) / "campaign.json"
                if result_path.exists():
                    item["evaluation"] = json.loads(result_path.read_text())
                bem = None
            attempted = {item["local_step"] for item in record["bem"]}
            ready = [step for step in checkpoints if step <= record["real_updates"] and step not in attempted]
            # On native failure retain and assess the last fully evaluated geometry too.
            if fem.poll() is not None and record["real_updates"] not in attempted and record["real_updates"] > 0:
                ready = sorted(set(ready + [record["real_updates"]]))
            save()
            if bem is None and ready:
                step = ready[0]
                launch = args.output / "bem" / f"global-{step+args.restart_step:06d}"
                launch.mkdir(parents=True)
                output = launch / "evaluation"
                bem_command = [os.sys.executable, "-B", str(EXAMPLE / "scripts/evaluate-fem-bem.py"),
                    "--commit", args.commit, "--fem-campaign", str(args.output / "campaign.json"),
                    "--verifier", str(args.verifier), "--output", str(output), "--steps", str(step),
                    "--allow-running", "--backend", "auto", "--blas-threads", "2"]
                item = dict(local_step=step, global_step=step+args.restart_step, output=str(output),
                            command=bem_command, started_at=time.time())
                record["bem"].append(item)
                save()
                with (launch / "run.log").open("w") as log:
                    bem = subprocess.Popen(bem_command, stdout=log, stderr=subprocess.STDOUT,
                                           env=env, start_new_session=True)
                item["pid"] = bem.pid
            if fem.poll() is not None and bem is None and not ready:
                break
            save()
            time.sleep(30)
        good_bem = all(item.get("exit_code") == 0 and not item.get("timed_out") and
                       item.get("evaluation", {}).get("state") == "completed" for item in record["bem"])
        code = 0 if stage.get("exit_code") == 0 and not stage.get("timed_out") and record["real_updates"] == args.updates and good_bem else 1
        record["state"] = "completed" if code == 0 else "failed"
    except Exception as error:
        record.update(state="failed", error=str(error))
    finally:
        for child in [bem, fem]:
            if child is not None and child.poll() is None:
                terminate(child)
        record.update(exit_code=code, finished_at=time.time())
        save()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
