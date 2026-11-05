#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Run the November 5 L10/n50 command on 28h without recurring monitoring."""
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

UPSTREAM = "aed431413caa8da37d2df1e813ca32befa931b00"
STATES = 200
THREADS = 12
PARAMETERS = ["--n=50", "--iterations=200", "--outer-radius=10",
              "--hmin-factor=0.05", "--hmax-factor=10", "--regularization=20",
              "--step=0.2", "--reconstruction=swift", "--mmg-adapt",
              "--thickness-min=0.05", "--swift-linear-solver=mumps",
              "--swift-linear-threads=12", "--swift-model-fit=10",
              "--swift-convergence-iterations-outer=30",
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


def complete_history(path, states=STATES):
    """Count only complete, contiguous solved states with the approved parameters."""
    if not path.exists():
        return []
    raw = path.read_bytes()
    raw = raw[:raw.rfind(b"\n") + 1]
    rows = list(csv.DictReader(raw.decode().splitlines()))
    h = 10 / 49
    expected = {"outer_radius": 10, "initial_grid_spacing": h,
                "hmin_factor": .05, "hmax_factor": 10,
                "dt": .2*h, "regularization_length": 20*h, "thickness_min": .05}
    for iteration, row in enumerate(rows):
        if None in row or any(v is None or v == "" for v in row.values()):
            raise ValueError("Malformed complete CSV row")
        if int(row["iteration"]) != iteration or iteration >= states:
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


def verify_blas(binary, env, library_directory):
    """Require the native executable to resolve LP64 BLAS/LAPACK to the selected OpenBLAS."""
    directory = library_directory.resolve()
    linked = subprocess.check_output(["ldd", str(binary)], env=env, text=True)
    libraries = {}
    for line in linked.splitlines():
        if "=>" not in line:
            continue
        name, target = line.split("=>", 1)
        name = name.strip()
        if name not in ("libblas.so.3", "liblapack.so.3", "libopenblas.so.0"):
            continue
        path = Path(target.strip().split()[0]).resolve()
        if not path.is_file() or not path.is_relative_to(directory):
            raise RuntimeError(f"Native executable would load an unintended math library: {line}")
        libraries[name] = dict(path=str(path), sha256=sha(path))
    if set(libraries) != {"libblas.so.3", "liblapack.so.3", "libopenblas.so.0"}:
        raise RuntimeError("Native executable does not resolve all three selected math libraries")
    probe = """import ctypes,json,sys
lib=ctypes.CDLL(sys.argv[1])
lib.openblas_get_config.restype=ctypes.c_char_p
lib.openblas_get_num_threads.restype=ctypes.c_int
print(json.dumps(dict(config=lib.openblas_get_config().decode(),threads=lib.openblas_get_num_threads())))
"""
    info = json.loads(subprocess.check_output(
        [sys.executable, "-B", "-c", probe, str(directory / "libopenblas.so.0")],
        env=env, text=True))
    if "USE64BITINT" in info["config"] or info["threads"] != THREADS:
        raise RuntimeError(f"Unexpected BLAS integer interface or thread count: {info}")
    return dict(libraries=libraries, runtime=info, ldd=linked)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "build", "output", "prefix"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--blas-library-dir", type=Path, required=True)
    parser.add_argument("--states", type=int, default=STATES)
    parser.add_argument("--reuse-campaign", type=Path,
                        help="Reuse a stopped campaign's verified binary without rebuilding it")
    args = parser.parse_args()
    if platform.node() != "28h":
        raise RuntimeError("All experiment stages are authorized only on 28h")
    if not 1 <= args.states <= STATES:
        raise ValueError("Require between 1 and 200 solved states")
    pinned(args.source, args.commit)
    if git(args.source, "merge-base", UPSTREAM, args.commit) != UPSTREAM:
        raise RuntimeError("Source does not descend from the approved upstream")
    numerical = "examples/ShapeOptimization/KelvinBall/KelvinBallOptimization.cpp"
    changed = git(args.source, "diff", "--name-only", UPSTREAM, args.commit,
                  "--", "src", "examples")
    for name in changed.splitlines():
        if name != numerical and not name.startswith("examples/ShapeOptimization/KelvinBall/reproduction/"):
            raise RuntimeError(f"Unexpected upstream difference: {name}")
    original = git(args.source, "show", f"{UPSTREAM}:{numerical}")
    compiled = (args.source / numerical).read_text().strip()
    if compiled != original.replace("decltype(classifier)::Parameters parameters;",
                                    "typename decltype(classifier)::Parameters parameters;"):
        raise RuntimeError("Numerical source differs beyond the compilation fix")
    modules = git(args.source, "submodule", "status", "third-party/mmg", "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in modules.splitlines()):
        raise RuntimeError("Dependencies are not at their pinned versions")
    for path in ((args.output,) if args.reuse_campaign else (args.build, args.output)):
        if path.exists() or path.resolve().is_relative_to(args.source.resolve()):
            raise RuntimeError(f"Require fresh external output/build: {path}")
    args.output.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(THREADS), OPENBLAS_NUM_THREADS=str(THREADS),
               MKL_NUM_THREADS=str(THREADS), KELVIN_SYMMETRY="cubic", PYTHONDONTWRITEBYTECODE="1")
    lib = args.prefix / "lib/x86_64-linux-gnu"
    env["LD_LIBRARY_PATH"] = ":".join(map(str, (args.blas_library_dir, lib, lib / "openmpi/lib", lib / "hdf5/serial"))) + ":" + env.get("LD_LIBRARY_PATH", "")
    parameters = [f"--iterations={args.states}" if p.startswith("--iterations=") else p
                  for p in PARAMETERS]
    record = dict(status="preparing", started_at=time.time(), host=platform.node(), pid=os.getpid(),
                  source=str(args.source), commit=args.commit, upstream=UPSTREAM,
                  submodules=modules, initial_shape="upstream unit sphere; no extra perturbation",
                  parameters=parameters, fem_threads=THREADS, blas_threads=THREADS,
                  blas_library_directory=str(args.blas_library_dir),
                  requested_states=args.states, requested_real_updates=args.states-1,
                  smallbox_experiments=False, graded_initial_grid=False,
                  configuration_basis="Carlos Nov5 11:51 Europe/Paris email; 12 threads on 28h",
                  stages=[], independent_bem=False, recurring_monitor=False,
                  native_wall_time_limit_seconds=None, automatic_restart=False)

    def checkpoint():
        write(args.output / "campaign.json", record)

    def stop(signum, frame):
        raise KeyboardInterrupt(f"Received signal {signum}")

    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    def execute(name, command, budget=None):
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
            process = subprocess.Popen(item["command"], cwd=directory, env=env,
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            item["pid"] = process.pid
            checkpoint()
            try:
                code = process.wait(timeout=budget)
            except BaseException as error:
                interrupted = error
                if process.poll() is None:
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
        if interrupted is not None:
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
        if args.reuse_campaign:
            previous = json.loads(args.reuse_campaign.read_text())
            previous_source = Path(previous["source"])
            pinned(previous_source, previous["commit"])
            if previous["status"] != "interrupted" or previous["upstream"] != UPSTREAM:
                raise RuntimeError("Require the preserved stopped campaign with the same upstream")
            stop_record = json.loads((args.reuse_campaign.parent / "user-stop.json").read_text())
            if stop_record["status"] != "stopped_by_user" or stop_record["remaining_processes"]:
                raise RuntimeError("Prior campaign has not been verified stopped")
            numerical_diff = git(args.source, "diff", "--name-only", previous["commit"],
                                 args.commit, "--", "src", "examples")
            if any(not p.startswith("examples/ShapeOptimization/KelvinBall/reproduction/")
                   for p in numerical_diff.splitlines()):
                raise RuntimeError("Reused binary has different numerical source")
            provenance = dict(source=str(previous_source), commit=previous["commit"],
                              campaign=str(args.reuse_campaign), stages={})
            for name in ("configure", "compile"):
                stage = next(s for s in previous["stages"] if s["name"] == name)
                actual = int((args.reuse_campaign.parent / name / "exit-code.txt").read_text())
                if actual != 0 or stage["exit_code"] != 0:
                    raise RuntimeError("Reused binary lacks a successful build record")
                provenance["stages"][name] = stage
            expected_binary = Path(provenance["stages"]["compile"]["command"][2]) / "bin/KelvinBall"
            if args.build / "bin/KelvinBall" != expected_binary or sha(expected_binary) != previous["binary_sha256"]:
                raise RuntimeError("Reused build path or binary hash differs")
            provenance["binary_sha256"] = previous["binary_sha256"]
            record["compiled_build_provenance"] = provenance
            checkpoint()
        else:
            for name, command in (("configure", command), ("compile", ["cmake", "--build", args.build,
                                              "--target", "KelvinBall", "-j", "4"])):
                if execute(name, command, 4*3600):
                    raise RuntimeError(f"{name} failed; see its actual exit and log")
        binary = args.build / "bin/KelvinBall"
        record["binary_sha256"] = sha(binary)
        record["math_libraries"] = verify_blas(binary, env, args.blas_library_dir)
        checkpoint()
        native = execute("optimization", [binary, *parameters])
        record["native_exit_code"] = native
        simulation = args.output / "optimization"
        history_path = simulation / "kelvin-ball.csv"
        rows = complete_history(history_path, args.states)
        record.update(fully_evaluated_states=len(rows), completed_real_updates=max(0, len(rows)-1),
                      budget_completed=(native == 0 and len(rows) == args.states))
        after = verify_blas(binary, env, args.blas_library_dir)
        if (after["libraries"] != record["math_libraries"]["libraries"] or
                after["runtime"] != record["math_libraries"]["runtime"]):
            raise RuntimeError("Math libraries changed during the run")
        if rows:
            step = int(rows[-1]["iteration"])
            snapshot = simulation / f"KelvinBall.Chamber.mesh.{step:06d}.h5"
            record.update(final_fem={key: float(rows[-1][key]) for key in
                                    ("k", "c", "q", "rho", "volume", "target_volume")},
                          last_solved_step=step, history_sha256=sha(history_path))
            if snapshot.is_file():
                record.update(last_solved_snapshot=str(snapshot), snapshot_sha256=sha(snapshot))
        record["status"] = "completed" if record["budget_completed"] else "native_ended_early"
        code = 0 if record["budget_completed"] else 1
    except BaseException as error:
        record.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                      error=str(error), error_type=type(error).__name__)
        traceback.print_exc()
    finally:
        record.update(exit_code=code, finished_at=time.time())
        checkpoint()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    sys.exit(main())
