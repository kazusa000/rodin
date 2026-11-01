#!/usr/bin/env python3
"""Frozen, 28h-only L15 SWIFT pilot, checkpoint verification and gated continuation."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
import traceback

import h5py
import numpy as np

UPSTREAM = "a9f1115ae71525978f086917a73958f95d77d1bf"
UTILITIES = "b8f3c8c5929c21aa9bd72236c0c1e843e0a72a64"
VERIFIER = "12f38a9a33dbbb0a1b7c0d0affcfff91620b927a"
PHYSICS = ["--graded-chamber", "--h=0.16666666666666666", "--outer-radius=15",
           "--near-field-radius=2", "--chamber-gradation=1.3", "--hmin-factor=0.1",
           "--hmax-factor=10", "--regularization=20", "--step=0.2", "--thickness-min=0.05",
           "--reconstruction=swift", "--mmg-adapt", "--swift-linear-solver=mumps",
           "--swift-linear-threads=4", "--swift-convergence-iterations-outer=30",
           "--swift-convergence-iterations-inner=15", "--save-mesh", "--save-checkpoints"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def pinned(root, commit):
    if git(root, "rev-parse", "HEAD") != commit or git(root, "status", "--porcelain"):
        raise RuntimeError(f"Source is not the clean pinned commit: {root}")


def datasets(path):
    data = {}
    with h5py.File(path, "r") as f:
        f.visititems(lambda k, obj: data.__setitem__(k, obj[()])
                     if isinstance(obj, h5py.Dataset) else None)
    if not data:
        raise ValueError(f"Empty HDF5: {path}")
    return data


def exact_hdf(left, right):
    a, b = datasets(left), datasets(right)
    if a.keys() != b.keys():
        raise ValueError(f"Checkpoint dataset keys differ: {left} / {right}")
    for key in a:
        if a[key].dtype != b[key].dtype or not np.array_equal(a[key], b[key]):
            raise ValueError(f"Checkpoint changed dataset {key}: {left}")
    return list(a)


def history(directory, count):
    path = directory / "kelvin-ball.csv"
    raw = path.read_text()
    if not raw.endswith("\n"):
        raise ValueError("Incomplete final CSV line")
    rows = list(csv.DictReader(raw.splitlines()))
    if len(rows) != count or [int(r["iteration"]) for r in rows] != list(range(count)):
        raise ValueError("Native run did not complete the requested evaluated states")
    for r in rows:
        if None in r or any(v is None or v == "" for v in r.values()):
            raise ValueError("Incomplete CSV row")
        for k in ("k", "c", "q", "rho", "volume", "target_volume", "mesh_quality_min"):
            if not math.isfinite(float(r[k])):
                raise ValueError(f"Nonfinite FEM diagnostic: {k}")
        if any(float(r[k]) <= 0 for k in ("k", "q", "volume", "mesh_quality_min")):
            raise ValueError("Invalid resistance, volume or mesh quality")
        for k, v in {"outer_radius": 15., "initial_grid_spacing": 1/6,
                     "dt": .2/6, "regularization_length": 20/6, "thickness_min": .05}.items():
            if not math.isclose(float(r[k]), v, rel_tol=1e-13, abs_tol=1e-15):
                raise ValueError(f"Numerical parameter changed: {k}")
    return rows


def compare_restart(pilot, restart):
    a, b = history(pilot, 6)[5], history(restart, 1)[0]
    for key in ("k", "c", "q", "rho", "volume", "target_volume"):
        if not math.isclose(float(a[key]), float(b[key]), rel_tol=1e-8, abs_tol=1e-10):
            raise ValueError(f"Restart failed numerical consistency: {key}: {a[key]} / {b[key]}")
    result = {"numerical_state": True, "checkpoint_datasets": {}, "fields": {}}
    for name in ("chamber.h5", "background.h5", "state.h5"):
        result["checkpoint_datasets"][name] = exact_hdf(
            pilot / "checkpoint-000005" / name, restart / "checkpoint-000005" / name)
    # XDMF attribute sidecars carry the complete P1 fields. Match by XDMF names,
    # not by a presumed attribute ordinal which can change with upstream output.
    import xml.etree.ElementTree as ET
    def field_path(directory, step, name):
        tree = ET.parse(directory / "KelvinBall.xdmf")
        for grid in tree.iter("Grid"):
            if grid.get("Name") != "Chamber":
                continue
            for attr in grid.findall("Attribute"):
                if attr.get("Name") == name:
                    text = attr.find("DataItem").text.strip()
                    filename, dataset = text.split(":", 1)
                    if f".{step:06d}." in filename:
                        return directory / filename, dataset
        raise ValueError(f"Cannot locate full {name} field for step {step}")
    for name, rtol, atol in (("Theta", 1e-7, 1e-9), ("Distance", 1e-12, 1e-12)):
        pa, ka = field_path(pilot, 5, name)
        pb, kb = field_path(restart, 0, name)
        with h5py.File(pa, "r") as fa, h5py.File(pb, "r") as fb:
            va, vb = fa[ka][()], fb[kb][()]
        if va.shape != vb.shape or not np.allclose(va, vb, rtol=rtol, atol=atol):
            raise ValueError(f"Restart changed full {name} field")
        result["fields"][name] = {"rtol": rtol, "atol": atol,
                                  "max_absolute_difference": float(np.max(np.abs(va-vb)))}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "build", "output", "prefix", "utilities", "verifier"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    if platform.node() != "28h":
        raise RuntimeError("This experiment is authorized only on 28h")
    pinned(args.source, args.commit)
    pinned(args.utilities, UTILITIES)
    pinned(args.verifier, VERIFIER)
    dependencies = git(args.source, "submodule", "status", "third-party/mmg", "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in dependencies.splitlines()):
        raise RuntimeError("Pinned submodules are missing or changed")
    if git(args.source, "merge-base", UPSTREAM, args.commit) != UPSTREAM:
        raise RuntimeError("Source does not descend from the approved upstream baseline")
    if git(args.source, "diff", "--name-only", UPSTREAM, args.commit, "--", "src"):
        raise RuntimeError("This experiment must not change upstream framework src/")
    for p in (args.build, args.output):
        if p.exists() or p.resolve().is_relative_to(args.source.resolve()):
            raise RuntimeError(f"Require fresh external path: {p}")
    args.output.mkdir(parents=True)
    sys.path.insert(0, str(args.utilities / "scripts"))
    from reproduce_smallbox_original import export_snapshot
    from evaluate_saved_smallbox import validate_scan
    env = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               KELVIN_SYMMETRY="cubic", PYTHONDONTWRITEBYTECODE="1")
    lib = args.prefix / "lib/x86_64-linux-gnu"
    env["LD_LIBRARY_PATH"] = str(lib) + ":" + str(lib / "openmpi/lib") + ":" + env.get("LD_LIBRARY_PATH", "")
    record = {"status": "preparing", "started_at": time.time(), "host": platform.node(),
              "pid": os.getpid(), "source": str(args.source), "commit": args.commit,
              "upstream": UPSTREAM, "utilities": UTILITIES, "verifier": VERIFIER,
              "submodules": dependencies,
              "parameters": PHYSICS, "fem_threads": 4, "blas_threads": 1,
              "initial_shape": "unit sphere", "pilot_real_updates": 5, "total_real_updates": 20,
              "smallbox_experiments": False, "stages": [], "bem": [], "validation": {},
              "scientific_success": False,
              "limitations": ["FEM rho is internal diagnostic", "BEM uses fixed regularization and full original panels; not a continuum or closure proof", "Pilot continuation gate is a diagnostic heuristic, not proof of reliable optimization"]}
    def checkpoint():
        write(args.output / "campaign.json", record)
    active = None
    def stop(signum, frame):
        raise KeyboardInterrupt(f"Received signal {signum}")
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    def execute(name, command, budget, stage_env=None, expected=0):
        nonlocal active
        pinned(args.source, args.commit)
        directory = args.output / name
        directory.mkdir(parents=True, exist_ok=False)
        item = {"name": name, "command": list(map(str, command)), "started_at": time.time(),
                "timeout_seconds": budget, "status": "running"}
        record["stages"].append(item)
        record["status"] = name
        checkpoint()
        with (directory / "run.log").open("w") as log:
            active = subprocess.Popen(item["command"], cwd=directory, env=stage_env or env,
                                      stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            item["pid"] = active.pid
            checkpoint()
            try:
                code = active.wait(timeout=budget)
            except BaseException as error:
                os.killpg(active.pid, signal.SIGTERM)
                try:
                    active.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(active.pid, signal.SIGKILL)
                    active.wait()
                item["interrupted_by"] = type(error).__name__
                item["exit_code"] = active.returncode
                item["finished_at"] = time.time()
                item["status"] = "interrupted"
                (directory / "exit-code.txt").write_text(str(active.returncode) + "\n")
                checkpoint()
                active = None
                raise
        active = None
        item.update(exit_code=code, finished_at=time.time(), status="completed" if code == expected else "failed")
        (directory / "exit-code.txt").write_text(str(code) + "\n")
        checkpoint()
        if (expected == 0 and code != 0) or (expected != 0 and code == 0):
            raise RuntimeError(f"Stage {name} returned {code}, expected {'nonzero' if expected else '0'}")
        return directory

    def bem(directory, local_step, global_step):
        snapshot = directory / f"KelvinBall.Chamber.mesh.{local_step:06d}.h5"
        output = args.output / "bem" / f"step-{global_step:06d}"
        output.mkdir(parents=True)
        certificate = export_snapshot(snapshot, output / "original-double.mesh")
        write(output / "original-double.json", certificate)
        item = {"global_step": global_step, "input": certificate, "status": "prepared",
                "solver_rtol": 1e-6}
        record["bem"].append(item)
        checkpoint()
        pinned(args.verifier, VERIFIER)
        if sha(snapshot) != certificate["source_sha256"]:
            raise ValueError("Input snapshot changed")
        command = [sys.executable, "-B", args.verifier / "scripts/bem_convergence_scan.py",
                   output / "original-double.mesh", "--panels", "full", "--subdivision-levels", "0",
                   "--epsilon-factors", "0.6", "--viscosity", "1", "--center", "0", "0", "0",
                   "--backend", "quotient-cg", "--solver-rtol", "1e-6", "--solver-maxiter", "6000",
                   "--preconditioner-panel-block-size", "256", "--matrix-free-row-chunk", "32",
                   "--solver-progress", "--matrix-free-device", "cpu",
                   "--json", output / "scan.json", "--csv", output / "scan.csv"]
        bem_env = dict(env, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        execute(f"bem/step-{global_step:06d}/solver", command, 72*3600, bem_env)
        row = validate_scan(output, args.verifier, certificate, solver_rtol=1e-6)
        if sha(snapshot) != certificate["source_sha256"] or sha(output / "original-double.mesh") != certificate["output_sha256"]:
            raise ValueError("BEM input changed during evaluation")
        # Sanitise unrelated NaN diagnostic columns, preserving the accepted scientific row on disk.
        item.update(status="accepted", result={k: row[k] for k in ("k", "c", "q", "rho")})
        checkpoint()

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
        execute("configure", command, 4*3600)
        execute("compile", ["cmake", "--build", args.build, "--target", "KelvinBall",
                            "KelvinL15CheckpointTest", "-j", "2"], 4*3600)
        binary = args.build / "bin/KelvinBall"
        test = args.build / "bin/KelvinL15CheckpointTest"
        record["binary_sha256"] = sha(binary)
        fixture = execute("validation/checkpoint", [test], 3600)
        record["validation"]["synthetic_roundtrip"] = exact_hdf(
            fixture / "fixture-original.h5", fixture / "fixture-restored.h5")
        corrupt = fixture / "fixture-corrupt.h5"
        shutil.copyfile(fixture / "fixture-original.h5", corrupt)
        with h5py.File(corrupt, "r+") as f:
            f["KelvinMMG/Tags0"][0] = np.uint64(2**63)
        execute("validation/reject-invalid-tag", [test, "--load", corrupt], 60, expected=1)
        record["validation"]["invalid_tag_rejected"] = True
        checkpoint()
        pilot = execute("pilot", [binary, *PHYSICS, "--iterations=6"], 8*3600)
        rows = history(pilot, 6)
        record["pilot_completed_real_updates"] = 5
        record["pilot_fem_rho"] = [float(row["rho"]) for row in rows]
        checkpoint()
        restart = execute("validation/restart", [binary, *PHYSICS, "--iterations=1",
                          f"--restart-checkpoint={pilot}/checkpoint-000005"], 2*3600)
        record["validation"]["restart"] = compare_restart(pilot, restart)
        checkpoint()
        bem(pilot, 5, 5)
        record["pilot_gate"] = {"restart_consistent": True, "bem_accepted": True,
                                "fem_net_improvement": float(rows[5]["rho"]) > float(rows[0]["rho"])}
        checkpoint()
        if not record["pilot_gate"]["fem_net_improvement"]:
            raise RuntimeError("Pilot requires review: no net FEM diagnostic improvement; no continuation launched")
        continuation = execute("continuation", [binary, *PHYSICS, "--iterations=16",
                               f"--restart-checkpoint={pilot}/checkpoint-000005"], 8*3600)
        final = history(continuation, 16)
        record["completed_real_updates"] = 20
        record["final_fem_rho"] = float(final[-1]["rho"])
        checkpoint()
        bem(continuation, 15, 20)
        record["status"] = "completed_budget"
        code = 0
    except BaseException as error:
        record["status"] = "stopped" if isinstance(error, KeyboardInterrupt) else "failed"
        record["error"] = str(error)
        traceback.print_exc()
    finally:
        record.update(exit_code=code, finished_at=time.time())
        checkpoint()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    sys.exit(main())
