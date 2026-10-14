#!/usr/bin/env python3
"""Frozen CPU-only 28h build, regression checks and 20-update FEM campaign."""

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


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--baseline-build-record", type=Path, required=True)
    parser.add_argument("--baseline-campaign", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if platform.node() != "28h":
        raise ValueError("This campaign is explicitly restricted to 28h.")
    if git("rev-parse", "HEAD") != args.commit or git("status", "--porcelain", "--untracked-files=all"):
        raise ValueError("Require the exact clean pushed source commit.")
    modules = git("submodule", "status", "third-party/mmg", "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in modules.splitlines()):
        raise ValueError("Require the exact initialized dependency submodules.")
    for path in [args.output, args.build_dir]:
        if path.resolve().is_relative_to(ROOT) or path.exists():
            raise ValueError(f"Require new directories outside source: {path}")
    args.output.mkdir(parents=True)
    preset = json.loads((EXAMPLE / "reproduction/fem-domain-L15.json").read_text())
    inherited = json.loads(args.baseline_build_record.read_text())
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", KELVIN_SYMMETRY="cubic")
    env["LD_LIBRARY_PATH"] = inherited["runtime_library_path"]
    record = dict(commit=args.commit, source=str(ROOT), hostname=platform.node(),
                  submodules=modules, preset=preset, threads=args.threads,
                  baseline_build_sha256=sha(args.baseline_build_record),
                  stages=[], state="building", started_at=time.time())
    def save():
        (args.output / "campaign.json").write_text(json.dumps(record, indent=2) + "\n")
    def run(name, command, cwd, timeout):
        cwd.mkdir(parents=True, exist_ok=True)
        stage = dict(name=name, command=[str(x) for x in command], started_at=time.time())
        record["stages"].append(stage)
        save()
        with (cwd / "run.log").open("w") as log:
            try:
                result = subprocess.run(stage["command"], cwd=cwd, env=env,
                                        stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
                stage["exit_code"] = result.returncode
            except subprocess.TimeoutExpired:
                stage.update(exit_code=124, timed_out=True)
        stage["finished_at"] = time.time()
        (cwd / "exit-code.txt").write_text(str(stage["exit_code"]) + "\n")
        save()
        if stage["exit_code"]:
            raise RuntimeError(f"{name} failed with actual exit {stage['exit_code']}; see {cwd / 'run.log'}")
    code = 1
    try:
        configure = list(inherited["commands"][0])
        configure[configure.index("-S") + 1] = str(ROOT)
        configure[configure.index("-B") + 1] = str(args.build_dir)
        run("configure", configure, args.output / "build-configure", 1800)
        run("build", ["cmake", "--build", args.build_dir, "--target", "KelvinBall",
                       "KelvinBallGradedChamberTest", "-j", "2"], args.output / "build", 7200)
        binary = args.build_dir / "examples/ShapeOptimization/KelvinBall/KelvinBall"
        test = binary.parent / "KelvinBallGradedChamberTest"
        record["binary_sha256"] = sha(binary)
        run("graded_grid_test", [test], args.output / "graded-grid-test", 60)
        run("original_regression", [binary, "--n=13", "--outer-radius=2", "--thickness-min=0",
            "--iterations=2", "--reconstruction=mmg", "--save-mesh"], args.output / "original-regression", 1800)
        history = list(csv.DictReader((args.output / "original-regression/kelvin-ball.csv").open()))
        baseline = json.loads(args.baseline_campaign.read_text())
        control = next(case for case in baseline["cases"] if case["name"] == "mmg")["summary"]
        if len(history) != 2:
            raise ValueError("Original regression must contain one real update.")
        for row, key in [(history[0], "first"), (history[1], "last")]:
            for metric in ["k", "c", "q", "rho", "volume"]:
                if not math.isclose(float(row[metric]), control[key][metric], rel_tol=1e-10, abs_tol=1e-12):
                    raise ValueError(f"Original baseline changed: {key} {metric}")
        record["original_regression"] = "matched initial/final k,c,q,rho,volume"
        run("large_box_geometry", [binary, *preset["optimization_arguments"], "--geometry-only"],
            args.output / "large-box-geometry", 1800)
        record["state"] = "optimizing"
        save()
        run("optimization", [binary, *preset["optimization_arguments"]], args.output / "optimization", 43200)
        rows = list(csv.DictReader((args.output / "optimization/kelvin-ball.csv").open()))
        if len(rows) != 21:
            raise ValueError("Twenty true updates require 21 evaluated states.")
        for row in rows[:-1]:
            if not (float(row["advection_increment"]) > 0 and float(row["stage_9_seconds"]) > 0):
                raise ValueError("A claimed update did not transport and reconstruct.")
        record.update(state="completed", evaluated_states=len(rows), real_updates=20,
                      initial={k: rows[0][k] for k in ["k", "c", "q", "rho", "volume"]},
                      final={k: rows[-1][k] for k in ["k", "c", "q", "rho", "volume"]})
        code = 0
    except Exception as error:
        record.update(state="failed", error=str(error))
    finally:
        record.update(finished_at=time.time(), exit_code=code)
        save()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
