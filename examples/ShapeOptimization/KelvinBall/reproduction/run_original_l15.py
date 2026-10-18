#!/usr/bin/env python3
"""Build original Rodin core, check the L2 control, then run L15 FEM from a sphere."""
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
import time

BASE = "9b3dfd7e1a49036f4f1074f3f7fc99f72bbb1f04"
H = 2 / 19


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "prefix", "build", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    if platform.node() != "28h":
        raise ValueError("Run only on the authorized 28h node")
    def git(*rest):
        return subprocess.check_output(["git", "-C", str(args.source), *rest], text=True).strip()
    if git("rev-parse", "HEAD") != args.commit or git("status", "--porcelain", "--untracked-files=all"):
        raise ValueError("Require clean pinned source")
    example = "examples/ShapeOptimization/KelvinBall/"
    unchanged = ["src", example + "KelvinBallOptimization.cpp", example + "Common.cpp",
                 example + "Metrics.cpp", example + "RotatedNitscheIntegrator.cpp", example + "SewedOutput.cpp"]
    if git("diff", BASE, "HEAD", "--", *unchanged):
        raise ValueError("Original core or update algorithm was changed")
    modules = git("submodule", "status", "third-party/mmg", "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in modules.splitlines()):
        raise ValueError("Require original pinned dependencies")
    for path in (args.build, args.output):
        if path.exists() or path.resolve().is_relative_to(args.source.resolve()):
            raise ValueError("Use fresh external output/build directories")
    args.output.mkdir(parents=True)
    record = dict(host=platform.node(), commit=args.commit, original_core=BASE, source=str(args.source),
                  source_status="", submodules=modules, unchanged_core_and_update_verified=True,
                  initial_shape="original unit sphere; no added perturbation", outer_halfwidth=15,
                  h=H, physical_step=.5 * H, regularization_length=100 * H, thickness_minimum=2 * H,
                  requested_states=20, requested_real_updates=19, bem_enabled=False,
                  state="building", started_at=time.time(), stages=[])
    def save():
        p=args.output / "campaign.json"
        pending=p.with_suffix(".json.pending")
        pending.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
        pending.replace(p)
    env = dict(os.environ, OMP_NUM_THREADS="6", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1")
    lib = args.prefix / "lib/x86_64-linux-gnu"
    env["LD_LIBRARY_PATH"] = ":".join(str(p) for p in (lib, lib / "openmpi/lib", lib / "hdf5/serial"))
    record["environment"] = {k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "LD_LIBRARY_PATH")}
    def run(name, command, hours):
        directory = args.output / name
        directory.mkdir(parents=True)
        item = dict(name=name, command=[str(x) for x in command], started_at=time.time())
        record["stages"].append(item)
        with (directory / "run.log").open("w") as log:
            child = subprocess.Popen(command, cwd=directory, env=env, stdout=log,
                                     stderr=subprocess.STDOUT, start_new_session=True)
            item["pid"] = child.pid
            save()
            try:
                code = child.wait(timeout=hours * 3600)
            except subprocess.TimeoutExpired:
                item["timed_out"] = True
                os.killpg(child.pid, signal.SIGTERM)
                try: code = child.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    code = child.wait()
        item.update(exit_code=code, finished_at=time.time())
        (directory / "exit-code.txt").write_text(str(code) + "\n")
        save()
        if code:
            raise RuntimeError(f"{name} exited {code}; preserve its log and raw output")
    code = 1
    try:
        run("configure", ["cmake", "-S", str(args.source / example / "reproduction"), "-B", str(args.build),
            "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_PREFIX_PATH={args.prefix}", f"-DMUMPS_ROOT={args.prefix}",
            f"-DMPI_C_HEADER_DIR={lib / 'openmpi/include'}",
            f"-DMPI_C_COMPILER_INCLUDE_DIRS={lib / 'openmpi/include'};{lib / 'openmpi/include/openmpi'}",
            "-DMPI_C_LIB_NAMES=mpi", f"-DMPI_mpi_LIBRARY={lib / 'openmpi/lib/libmpi.so'}",
            f"-DCARLOS_MUMPS_SUPPORT_LIBRARIES={lib / 'openmpi/lib/libmpi_cxx.so'}"], 1)
        run("compile", ["cmake", "--build", str(args.build), "--target", "KelvinBall", "KelvinBallGradedDomainTest", "-j", "2"], 4)
        binary = args.build / "bin/KelvinBall"
        record["binary_sha256"] = sha(binary)
        run("validation/configuration", [str(args.build / "bin/KelvinBallGradedDomainTest")], 1)
        record["state"] = "validating_original_smallbox"
        common = ["--regularization=100", "--step=0.5", "--advection-quadrature=3", "--reconstruction=mmg",
                  "--mmg-adapt", "--mmg-snap=0.005", "--thickness-min=2", "--thickness-weight=1",
                  "--motion-every=100", "--save-mesh"]
        run("validation/l2-uniform", [str(binary), "--n=20", "--outer-radius=2", *common, "--iterations=1"], 1)
        run("validation/l2-graded", [str(binary), "--graded-chamber", f"--h={H:.17g}",
                                    "--outer-radius=2", *common, "--iterations=1"], 1)
        expected = dict(k=90.26243247201937, c=-0.00016624440373000674,
                        q=27.630038343536466, volume=0.17390471514326827)
        controls = []
        for name in ("l2-uniform", "l2-graded"):
            with (args.output / f"validation/{name}/kelvin-ball.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            if len(rows) != 1 or any(not math.isclose(float(rows[0][k]), v, rel_tol=1e-8, abs_tol=1e-10)
                                     for k, v in expected.items()):
                raise ValueError(f"Original small-box numerical control failed: {name}")
            controls.append({k: rows[0][k] for k in expected})
        import h5py
        import numpy as np
        paths = [args.output / f"validation/{n}/KelvinBall.Chamber.mesh.000000.h5" for n in ("l2-uniform", "l2-graded")]
        with h5py.File(paths[0]) as a, h5py.File(paths[1]) as b:
            for key in ["Mesh/Geometry/Vertices", "Mesh/XDMF/Topology", *[f"Mesh/Attributes/{i}" for i in range(4)]]:
                if not np.array_equal(a[key][:], b[key][:]):
                    raise ValueError(f"The L2 graded limiting control changed the mesh: {key}")
        record["validation"] = dict(configuration_exit=0, smallbox_controls=controls,
                                     numerical_controls_passed=True, exact_l2_mesh_control_passed=True)
        record["state"] = "optimizing"
        run("optimization", [str(binary), "--graded-chamber", f"--h={H:.17g}", "--outer-radius=15",
            "--near-field-radius=2", "--chamber-gradation=1.3", *common, "--iterations=20"], 8)
        record["state"] = "completed_native_budget"
        code = 0
    except Exception as error:
        record.update(state="failed", error=str(error))
    finally:
        csv_path = args.output / "optimization/kelvin-ball.csv"
        if csv_path.exists():
            with csv_path.open() as stream:
                rows = [r for r in csv.DictReader(stream) if None not in r and r.get("omega_z") not in (None, "")]
            record["fully_evaluated_states"] = len(rows)
            record["completed_real_updates"] = max(0, len(rows) - 1)
            record["latest_fem_diagnostic"] = rows[-1] if rows else None
        record.update(exit_code=code, finished_at=time.time())
        save()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
