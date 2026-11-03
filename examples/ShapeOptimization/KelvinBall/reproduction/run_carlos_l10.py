#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Run Carlos's November 3 L10/n60 experiment, then assess its final solved shape."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time
import traceback

UPSTREAM = "d948427a50c00264c81ca40b2f616220581c7a26"
UTILITIES = "b8f3c8c5929c21aa9bd72236c0c1e843e0a72a64"
VERIFIER = "12f38a9a33dbbb0a1b7c0d0affcfff91620b927a"
STATES = 200
PARAMETERS = ["--n=60", "--iterations=200", "--outer-radius=10",
              "--hmin-factor=0.1", "--hmax-factor=10", "--regularization=20",
              "--step=0.2", "--reconstruction=swift", "--mmg-adapt",
              "--thickness-min=0.05", "--swift-linear-solver=mumps",
              "--swift-linear-threads=4", "--swift-convergence-iterations-outer=30",
              "--swift-convergence-iterations-inner=15", "--save-mesh"]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def pinned(root, commit):
    if git(root, "rev-parse", "HEAD") != commit or git(root, "status", "--porcelain"):
        raise RuntimeError(f"Require clean frozen source: {root}")


def complete_history(path):
    """Ignore a partially flushed final line; reject corrupt completed states."""
    if not path.exists():
        return []
    raw = path.read_bytes()
    raw = raw[:raw.rfind(b"\n") + 1]
    rows = list(csv.DictReader(raw.decode().splitlines()))
    h = 10 / 59
    expected = {"outer_radius": 10, "initial_grid_spacing": h,
                "dt": .2*h, "regularization_length": 20*h, "thickness_min": .05}
    for iteration, row in enumerate(rows):
        if None in row or any(v is None or v == "" for v in row.values()):
            raise ValueError("Malformed complete CSV row")
        if int(row["iteration"]) != iteration or iteration >= STATES:
            raise ValueError("CSV iterations are not the requested contiguous states")
        for key in ("k", "c", "q", "rho", "volume", "target_volume"):
            if not math.isfinite(float(row[key])):
                raise ValueError(f"Nonfinite FEM diagnostic: {key}")
        if any(float(row[key]) <= 0 for key in ("k", "q", "volume", "target_volume")):
            raise ValueError("Nonpositive resistance or volume")
        for key, value in expected.items():
            if not math.isclose(float(row[key]), value, rel_tol=1e-12, abs_tol=1e-14):
                raise ValueError(f"Unexpected actual FEM parameter: {key}")
    return rows


def final_snapshot(directory, rows):
    if not rows:
        raise ValueError("No complete FEM state available for independent BEM")
    step = int(rows[-1]["iteration"])
    path = directory / f"KelvinBall.Chamber.mesh.{step:06d}.h5"
    if not path.is_file():
        raise ValueError("Final evaluated state lacks its original HDF5 snapshot")
    return step, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "build", "output", "prefix", "utilities", "verifier"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    if platform.node() != "28h":
        raise RuntimeError("All experiment stages are authorized only on 28h")
    for root, commit in ((args.source, args.commit), (args.utilities, UTILITIES),
                         (args.verifier, VERIFIER)):
        pinned(root, commit)
    if git(args.source, "merge-base", UPSTREAM, args.commit) != UPSTREAM:
        raise RuntimeError("Source does not descend from the approved upstream")
    # The only permitted C++ difference is a required dependent-type keyword.
    source_diff = git(args.source, "diff", "--name-only", UPSTREAM, args.commit,
                      "--", "src", "examples/ShapeOptimization/KelvinBall")
    numerical = "examples/ShapeOptimization/KelvinBall/KelvinBallOptimization.cpp"
    for name in source_diff.splitlines():
        if name != numerical and not name.startswith("examples/ShapeOptimization/KelvinBall/reproduction/"):
            raise RuntimeError(f"Unexpected upstream difference: {name}")
    original = git(args.source, "show", f"{UPSTREAM}:{numerical}")
    compiled = (args.source / numerical).read_text().strip()
    if compiled != original.replace("decltype(classifier)::Parameters parameters;",
                                    "typename decltype(classifier)::Parameters parameters;"):
        raise RuntimeError("Numerical source differs beyond the dependent-type compilation fix")
    modules = git(args.source, "submodule", "status", "third-party/mmg", "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in modules.splitlines()):
        raise RuntimeError("Dependencies are not at their pinned versions")
    for path in (args.build, args.output):
        if path.exists() or any(path.resolve().is_relative_to(root.resolve())
                               for root in (args.source, args.utilities, args.verifier)):
            raise RuntimeError(f"Require fresh external output/build: {path}")
    args.output.mkdir(parents=True)
    sys.path.insert(0, str(args.utilities / "scripts"))
    from reproduce_smallbox_original import export_snapshot
    from evaluate_saved_smallbox import validate_scan
    env = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", KELVIN_SYMMETRY="cubic", PYTHONDONTWRITEBYTECODE="1")
    lib = args.prefix / "lib/x86_64-linux-gnu"
    env["LD_LIBRARY_PATH"] = ":".join(map(str, (lib, lib / "openmpi/lib", lib / "hdf5/serial"))) + ":" + env.get("LD_LIBRARY_PATH", "")
    record = dict(status="preparing", started_at=time.time(), host=platform.node(), pid=os.getpid(),
                  source=str(args.source), commit=args.commit, upstream=UPSTREAM,
                  utilities=UTILITIES, verifier=VERIFIER, submodules=modules,
                  initial_shape="upstream unit sphere; no extra perturbation",
                  parameters=PARAMETERS, fem_threads=4, blas_threads=1,
                  requested_states=STATES, requested_real_updates=STATES-1,
                  smallbox_experiments=False, graded_initial_grid=False,
                  configuration_basis="Carlos Nov3 email outer-radius10/n60; other settings from upstream README",
                  stages=[], bem=[], limitations=[
                      "FEM rho is an internal diagnostic; independent BEM uses fixed epsilon factor .6",
                      "Original geometric seams are preserved by full-panel export",
                      "A completed iteration budget does not establish convergence"])

    def checkpoint():
        write(args.output / "campaign.json", record)

    def stop(signum, frame):
        raise KeyboardInterrupt(f"Received signal {signum}")

    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    def execute(name, command, budget, stage_env=None):
        pinned(args.source, args.commit)
        directory = args.output / name
        directory.mkdir(parents=True, exist_ok=False)
        item = dict(name=name, command=list(map(str, command)), started_at=time.time(),
                    timeout_seconds=budget, status="running")
        record["stages"].append(item)
        record["status"] = name
        checkpoint()
        interrupted = None
        with (directory / "run.log").open("w") as log:
            process = subprocess.Popen(item["command"], cwd=directory, env=stage_env or env,
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            item["pid"] = process.pid
            checkpoint()
            try:
                code = process.wait(timeout=budget)
            except BaseException as error:
                interrupted = error
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                code = process.returncode
                item["interrupted_by"] = type(error).__name__
        item.update(exit_code=code, finished_at=time.time(),
                    status="completed" if code == 0 and interrupted is None else "failed")
        (directory / "exit-code.txt").write_text(str(code) + "\n")
        checkpoint()
        if interrupted is not None and not isinstance(interrupted, subprocess.TimeoutExpired):
            raise interrupted
        return code

    code = 1
    try:
        checkpoint()
        mpi = lib / "openmpi"
        command = ["cmake", "-S", args.source / "examples/ShapeOptimization/KelvinBall/reproduction",
                   "-B", args.build, "-DCMAKE_BUILD_TYPE=Release",
                   f"-DCMAKE_PREFIX_PATH={args.prefix}", f"-DMUMPS_ROOT={args.prefix}",
                   f"-DMPI_C_COMPILER_INCLUDE_DIRS={mpi}/include;{mpi}/include/openmpi",
                   f"-DMPI_C_HEADER_DIR={mpi}/include", "-DMPI_C_LIB_NAMES=mpi",
                   f"-DMPI_mpi_LIBRARY={mpi}/lib/libmpi.so"]
        for name, command in (("configure", command), ("compile", ["cmake", "--build", args.build,
                                              "--target", "KelvinBall", "-j", "4"])):
            if execute(name, command, 4*3600):
                raise RuntimeError(f"{name} failed; see its actual exit and log")
        binary = args.build / "bin/KelvinBall"
        record["binary_sha256"] = sha(binary)
        checkpoint()
        native = execute("optimization", [binary, *PARAMETERS], 72*3600)
        record["native_exit_code"] = native
        simulation = args.output / "optimization"
        history_path = simulation / "kelvin-ball.csv"
        rows = complete_history(history_path)
        step, snapshot = final_snapshot(simulation, rows)
        record.update(fully_evaluated_states=len(rows), completed_real_updates=len(rows)-1,
                      budget_completed=(native == 0 and len(rows) == STATES),
                      final_fem={key: float(rows[-1][key]) for key in
                                 ("k", "c", "q", "rho", "volume", "target_volume")},
                      selected_bem_step=step, selected_bem_geometry_fem_evaluated=True,
                      history_sha256=sha(history_path))
        checkpoint()
        output = args.output / "bem" / f"step-{step:06d}"
        output.mkdir(parents=True)
        certificate = export_snapshot(snapshot, output / "original-double.mesh")
        write(output / "original-double.json", certificate)
        item = dict(step=step, input=certificate, status="prepared", solver_rtol=1e-6)
        record["bem"].append(item)
        checkpoint()
        pinned(args.utilities, UTILITIES)
        pinned(args.verifier, VERIFIER)
        command = [sys.executable, "-B", args.verifier / "scripts/bem_convergence_scan.py",
                   output / "original-double.mesh", "--panels", "full", "--subdivision-levels", "0",
                   "--epsilon-factors", "0.6", "--viscosity", "1", "--center", "0", "0", "0",
                   "--backend", "quotient-cg", "--solver-rtol", "1e-6", "--solver-maxiter", "6000",
                   "--preconditioner-panel-block-size", "256", "--matrix-free-row-chunk", "32",
                   "--solver-progress", "--matrix-free-device", "cpu",
                   "--json", output / "scan.json", "--csv", output / "scan.csv"]
        bem_env = dict(env, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        bem_exit = execute(f"bem/step-{step:06d}/solver", command, 72*3600, bem_env)
        item["solver_exit_code"] = bem_exit
        if bem_exit:
            item["status"] = "failed"
            raise RuntimeError(f"Final independent BEM returned {bem_exit}")
        row = validate_scan(output, args.verifier, certificate, solver_rtol=1e-6)
        if sha(snapshot) != certificate["source_sha256"] or sha(history_path) != record["history_sha256"] or sha(output / "original-double.mesh") != certificate["output_sha256"]:
            raise RuntimeError("Immutable FEM/BEM input changed")
        item.update(status="accepted", result={key: row[key] for key in
                     ("k", "c", "q", "rho", "linear_relative_residual", "rhs_relative_residuals")})
        record["status"] = "completed_with_valid_bem" if record["budget_completed"] else "fem_stopped_with_valid_final_bem"
        code = 0 if record["budget_completed"] else 1
    except BaseException as error:
        record.update(status="failed", error=str(error), error_type=type(error).__name__)
        traceback.print_exc()
    finally:
        record.update(exit_code=code, finished_at=time.time())
        checkpoint()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    sys.exit(main())
