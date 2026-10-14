#!/usr/bin/env python3
"""Pinned historical FEM updates and independent CPU BEM measurements on 28h."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import time

EXAMPLE = Path(__file__).resolve().parents[1]
ROOT = EXAMPLE.parents[2]


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--build-record", type=Path, required=True)
    parser.add_argument("--verifier", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--historical-record", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--reuse-build-dir", action="store_true")
    args = parser.parse_args()
    if platform.node() != "28h":
        raise ValueError("This experiment must run entirely on 28h.")
    preset = json.loads((EXAMPLE / "reproduction/fem-domain-L15-legacy.json").read_text())
    for source, commit in [(ROOT, args.commit), (args.verifier, preset["historical_bem_commit"])]:
        if git(source, "rev-parse", "HEAD") != commit or git(source, "status", "--porcelain", "--untracked-files=all"):
            raise ValueError(f"Require exact clean pinned source: {source}")
    for path in [args.output, args.build_dir]:
        if (path.exists() and not (path == args.build_dir and args.reuse_build_dir)) or path.resolve().is_relative_to(ROOT):
            raise ValueError(f"Require a new directory outside source: {path}")
    args.output.mkdir(parents=True)
    build = json.loads(args.build_record.read_text())
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", KELVIN_SYMMETRY="cubic", PYTHONDONTWRITEBYTECODE="1",
               LD_LIBRARY_PATH=build["runtime_library_path"])
    record = dict(source=str(ROOT), commit=args.commit, host=platform.node(), preset=preset,
                  verifier=dict(path=str(args.verifier), commit=preset["historical_bem_commit"]),
                  build_record_sha256=digest(args.build_record), stages=[], state="building",
                  omp_threads=args.threads, started_at=time.time(), bem=[])

    def save():
        (args.output / "campaign.json").write_text(json.dumps(record, indent=2) + "\n")

    def run(name, command, directory, timeout):
        directory.mkdir(parents=True, exist_ok=False)
        stage = dict(name=name, command=list(map(str, command)), started_at=time.time())
        record["stages"].append(stage)
        save()
        with (directory / "run.log").open("w") as log:
            try:
                child = subprocess.run(stage["command"], cwd=directory, env=env,
                                       stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
                stage["exit_code"] = child.returncode
            except subprocess.TimeoutExpired:
                stage.update(exit_code=124, timed_out=True)
        stage["finished_at"] = time.time()
        (directory / "exit-code.txt").write_text(str(stage["exit_code"]) + "\n")
        save()
        if stage["exit_code"]:
            raise RuntimeError(f"{name} failed, actual exit {stage['exit_code']}: {directory / 'run.log'}")

    code = 1
    try:
        configure = list(build["commands"][0])
        configure[configure.index("-S") + 1] = str(ROOT)
        configure[configure.index("-B") + 1] = str(args.build_dir)
        run("configure", configure, args.output / "configure", 1800)
        run("build", ["cmake", "--build", args.build_dir, "--target", "KelvinBall",
                       "KelvinBallGradedChamberTest", "KelvinBallSeedTest", "-j", "2"], args.output / "build", 7200)
        binary = args.build_dir / "examples/ShapeOptimization/KelvinBall/KelvinBall"
        record["binary_sha256"] = digest(binary)
        run("graded_core_test", [binary.parent / "KelvinBallGradedChamberTest"], args.output / "graded-core-test", 60)
        run("seed_test", [binary.parent / "KelvinBallSeedTest"], args.output / "seed-test", 60)
        original = json.loads(args.historical_record.read_text())
        parameters = original["command"][1:]
        parameters = ["--iterations=1" if x.startswith("--iterations=") else x for x in parameters]
        run("historical_L2_initial", [binary, *parameters], args.output / "historical-L2-initial", 3600)
        first = next(csv.DictReader((args.output / "historical-L2-initial/kelvin-ball.csv").open()))
        control = original["summary"]["first"]
        comparison = {k:dict(actual=float(first[k]), historical=control[k])
                      for k in ["k", "c", "q", "rho", "volume"]}
        record["historical_initial_comparison"] = comparison
        # The restored example is run on the current core and compiler. Record
        # differences rather than treating legacy FEM coefficients as BEM truth.
        record["historical_initial_matches"] = all(math.isclose(v["actual"], v["historical"],
            rel_tol=1e-8, abs_tol=1e-10) for v in comparison.values())
        save()
        parameters = preset["optimization_arguments"]
        probe = ["--iterations=2" if x.startswith("--iterations=") else
                 "--step=0.005" if x.startswith("--step=") else x for x in parameters]
        record["state"] = "small_step_check"
        run("small_step_update", [binary, *probe], args.output / "small-step-update", 3600)
        rows = list(csv.DictReader((args.output / "small-step-update/kelvin-ball.csv").open()))
        if len(rows) != 2 or float(rows[0]["advection_increment"]) <= 0:
            raise ValueError("Small-step transport/reconstruction did not execute.")
        if any(abs(float(x["interface_shift_max"])) > 1e-12 for x in rows):
            raise ValueError("Distance projection moved Gamma.")
        record["small_step_check"] = dict(predicted=float(rows[1]["predicted_delta_rho"]),
                                          actual=float(rows[1]["actual_delta_rho"]))
        record["state"] = "optimizing"
        run("optimization", [binary, *parameters], args.output / "optimization", 14400)
        rows = list(csv.DictReader((args.output / "optimization/kelvin-ball.csv").open()))
        if len(rows) != preset["evaluated_states"]:
            raise ValueError("Wrong number of fully evaluated designs.")
        for row in rows[:-1]:
            if float(row["advection_increment"]) <= 0 or float(row["stage_9_seconds"]) <= 0:
                raise ValueError("A claimed real update did not transport and reconstruct.")
        record.update(real_updates=len(rows)-1, fem_history=rows, state="independent_bem")
        save()
        destination = args.output / "independent-bem"
        run("independent_bem_double", [args.python, "-B", EXAMPLE / "scripts/evaluate-fem-bem.py",
            "--commit", args.commit, "--fem-campaign", args.output / "campaign.json",
            "--verifier", args.verifier, "--output", destination], args.output / "bem-launch", 172800)
        measured = json.loads((destination / "campaign.json").read_text())
        record["bem"] = measured["evaluations"]
        record["bem_campaign_sha256"] = digest(destination / "campaign.json")
        save()
        record["state"] = "completed"
        code = 0
    except Exception as error:
        record.update(state="failed", error=str(error))
    finally:
        record.update(exit_code=code, finished_at=time.time())
        save()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
