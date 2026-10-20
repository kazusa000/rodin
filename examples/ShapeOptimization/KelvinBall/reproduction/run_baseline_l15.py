#!/usr/bin/env python3
"""Validate the October 11 baseline before running its L15 graded initial background on 28h."""
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

BASE = "3ae2c0414ffebf11f5814a46510ddf4d8375ba0a"
H = 1 / 6


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "prefix", "build", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--reuse-build-record", type=Path)
    parser.add_argument("--restart-from", type=Path)
    parser.add_argument("--updates", type=int, default=20)
    args = parser.parse_args()
    if platform.node() != "28h":
        raise ValueError("Run only on the authorized 28h node")
    def git(*rest):
        return subprocess.check_output(["git", "-C", str(args.source), *rest], text=True).strip()
    if git("rev-parse", "HEAD") != args.commit or git("status", "--porcelain", "--untracked-files=all"):
        raise ValueError("Require clean pinned source")
    example = "examples/ShapeOptimization/KelvinBall/"
    unchanged = ["src", *[example + name for name in ("KelvinBallOptimization.h",
                 "Common.cpp", "Common.h", "Metrics.cpp", "Metrics.h", "RotatedNitscheIntegrator.cpp",
                 "RotatedNitscheIntegrator.h", "SewedOutput.cpp", "SewedOutput.h", "NodalTransport.h",
                 "Thickness.h", "RotatedCharacteristicContinuation.h", "HilbertIdentification.h")]]
    if git("diff", BASE, "HEAD", "--", *unchanged):
        raise ValueError("October 11 core or update algorithm was changed")
    original_optimizer = git("show", BASE + ":" + example + "KelvinBallOptimization.cpp")
    current_optimizer = (args.source / example / "KelvinBallOptimization.cpp").read_text().strip()
    update_start = "  Optional<MMG::Mesh> nextMesh;"
    if original_optimizer[original_optimizer.index(update_start):] != current_optimizer[current_optimizer.index(update_start):]:
        raise ValueError("The original numerical update loop changed")
    if args.updates <= 0:
        raise ValueError("Require a positive bounded update budget")
    modules = git("submodule", "status", "third-party/mmg", "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in modules.splitlines()):
        raise ValueError("Require original pinned dependencies")
    controller_source = Path(__file__).resolve().parents[4]
    controller_commit = subprocess.check_output(
        ["git", "-C", str(controller_source), "rev-parse", "HEAD"], text=True).strip()
    if subprocess.check_output(["git", "-C", str(controller_source), "status", "--porcelain", "--untracked-files=all"], text=True).strip():
        raise ValueError("Require clean pinned controller")
    if controller_commit != args.commit and subprocess.check_output(
        ["git", "-C", str(controller_source), "diff", args.commit, controller_commit, "--", "src", example], text=True).strip():
        changes = subprocess.check_output(["git", "-C", str(controller_source), "diff", "--name-only",
            args.commit, controller_commit, "--", "src", example], text=True).splitlines()
        if any(not name.startswith(example + "reproduction/") for name in changes):
            raise ValueError("Prebuilt numerical source differs from the controller's numerical source")
    if args.output.exists() or args.output.resolve().is_relative_to(args.source.resolve()):
        raise ValueError("Use a fresh external output directory")
    if args.build.resolve().is_relative_to(args.source.resolve()):
        raise ValueError("Keep the build outside the numerical source")
    prebuilt = None
    if args.reuse_build_record:
        prebuilt = json.loads(args.reuse_build_record.read_text())
        if (prebuilt.get("commit") != args.commit or
            Path(prebuilt["source"]).resolve() != args.source.resolve() or
            prebuilt.get("binary_sha256") != sha(args.build / "bin/KelvinBall") or
            not any(stage.get("name") == "compile" and stage.get("exit_code") == 0
                    for stage in prebuilt.get("stages", []))):
            raise ValueError("Require matching pinned source and a successfully compiled recorded binary")
        cached = (args.build / "CMakeCache.txt").read_text()
        if f"CMAKE_HOME_DIRECTORY:INTERNAL={args.source / example / 'reproduction'}" not in cached:
            raise ValueError("The recorded build uses another source directory")
    elif args.build.exists():
        raise ValueError("Use a fresh build or an explicitly validated existing build record")
    args.output.mkdir(parents=True)
    record = dict(host=platform.node(), commit=args.commit, controller_commit=controller_commit,
                  controller_source=str(controller_source), baseline=BASE, source=str(args.source),
                  source_status="", submodules=modules, unchanged_core_and_update_loop_verified=True,
                  initial_shape="original unit sphere; no added perturbation", outer_halfwidth=15,
                  h=H, physical_step=.2 * H, regularization_length=4 * H, thickness_minimum=0,
                  initial_grading="core2; exterior interval ratio1.3; maximum10h",
                  mmg_adapt=False, grading_policy="initial grid; baseline edge-length-driven MMG thereafter",
                  requested_states=20, requested_real_updates=19, bem_enabled=False,
                  state="building", started_at=time.time(), stages=[])
    def save():
        p=args.output / "campaign.json"
        pending=p.with_suffix(".json.pending")
        pending.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
        pending.replace(p)
    env = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
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
            except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                item["timed_out"] = isinstance(error, subprocess.TimeoutExpired)
                os.killpg(child.pid, signal.SIGTERM)
                try: code = child.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    code = child.wait()
        item.update(exit_code=code, finished_at=time.time())
        (directory / "exit-code.txt").write_text(str(code) + "\n")
        save()
        if code or "termination_signal" in record:
            raise RuntimeError(f"{name} exited {code}; preserve its log and raw output")
    def interrupted(signum, frame):
        record["termination_signal"] = signum
        raise KeyboardInterrupt(f"Controller received signal {signum}")
    signal.signal(signal.SIGTERM, interrupted)
    code = 1
    try:
        if prebuilt is None:
            run("configure", ["cmake", "-S", str(args.source / example / "reproduction"), "-B", str(args.build),
            "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_PREFIX_PATH={args.prefix}", f"-DMUMPS_ROOT={args.prefix}",
            f"-DMPI_C_HEADER_DIR={lib / 'openmpi/include'}",
            f"-DMPI_C_COMPILER_INCLUDE_DIRS={lib / 'openmpi/include'};{lib / 'openmpi/include/openmpi'}",
            "-DMPI_C_LIB_NAMES=mpi", f"-DMPI_mpi_LIBRARY={lib / 'openmpi/lib/libmpi.so'}",
            f"-DCARLOS_MUMPS_SUPPORT_LIBRARIES={lib / 'openmpi/lib/libmpi_cxx.so'}"], 1)
            run("compile", ["cmake", "--build", str(args.build), "--target", "KelvinBall", "KelvinBallGradedDomainTest", "-j", "2"], 4)
        else:
            record["compiled_build_provenance"] = dict(campaign=str(args.reuse_build_record),
                campaign_sha256=sha(args.reuse_build_record), numerical_source_commit=args.commit,
                stages=[stage for stage in prebuilt["stages"] if stage["name"] in ("configure", "compile")])
        binary = args.build / "bin/KelvinBall"
        record["binary_sha256"] = sha(binary)
        run("validation/configuration", [str(args.build / "bin/KelvinBallGradedDomainTest")], 1)
        record["state"] = "validating_original_smallbox"
        common = ["--regularization=4", "--thickness-min=0", "--reconstruction=mmg",
                  "--hmin-factor=0.1", "--hmax-factor=10", "--mmg-snap=0", "--save-mesh"]
        run("validation/l2-uniform", [str(binary), "--n=13", "--outer-radius=2", *common,
                                     "--step=0.1", "--iterations=2"], 1)
        run("validation/l2-graded", [str(binary), "--graded-chamber", f"--h={H:.17g}",
                                    "--outer-radius=2", *common, "--step=0.1", "--iterations=2"], 1)
        expected = [dict(k=90.8147229985936, c=0.002882121600478904,
                         q=28.078362671335412, volume=0.17302546463183208),
                    dict(k=94.58385022124627, c=0.002050327425223212,
                         q=29.44817736501169, volume=0.17372678650619647)]
        controls = []
        for name in ("l2-uniform", "l2-graded"):
            with (args.output / f"validation/{name}/kelvin-ball.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            if len(rows) != 2 or any(not math.isclose(float(row[k]), value, rel_tol=1e-8, abs_tol=1e-10)
                                    for row, reference in zip(rows, expected) for k, value in reference.items()):
                raise ValueError(f"October 11 small-box two-state regression failed: {name}")
            controls.append([{k: row[k] for k in expected[0]} for row in rows])
        import h5py
        import numpy as np
        for step in (0, 1):
            paths = [args.output / f"validation/{n}/KelvinBall.Chamber.mesh.{step:06d}.h5"
                     for n in ("l2-uniform", "l2-graded")]
            with h5py.File(paths[0]) as a, h5py.File(paths[1]) as b:
                for key in ["Mesh/Geometry/Vertices", "Mesh/XDMF/Topology", *[f"Mesh/Attributes/{i}" for i in range(4)]]:
                    if not np.array_equal(a[key][:], b[key][:]):
                        raise ValueError(f"L2 limiting mesh changed at step {step}: {key}")
        record["validation"] = dict(configuration_exit=0, smallbox_controls=controls,
                                     numerical_controls_passed=True, exact_l2_mesh_control_passed=True,
                                     reference_control_step_factor=0.1)
        command = [str(binary), "--graded-chamber", f"--h={H:.17g}", "--outer-radius=15",
                   "--near-field-radius=2", "--chamber-gradation=1.3", *common, "--step=0.2"]
        if args.restart_from:
            from restart_snapshot import export_restart
            parent = json.loads((args.restart_from / "campaign.json").read_text())
            if parent.get("exit_code") != 0 or parent.get("completed_real_updates") != 19:
                raise ValueError("Require the completed original 19-update L15 campaign")
            with (args.restart_from / "optimization/kelvin-ball.csv").open() as stream:
                previous = list(csv.DictReader(stream))
            if len(previous) != 20 or int(previous[-1]["iteration"]) != 19:
                raise ValueError("The parent CSV is incomplete")
            snapshot = args.restart_from / "optimization/KelvinBall.Chamber.mesh.000019.h5"
            checkpoint = args.output / "restart-original-double.mesh"
            record["restart_input"] = export_restart(snapshot, checkpoint, 15)
            volume_reference = previous[0]["target_volume"]
            restart = [f"--restart-mesh={checkpoint}", f"--original-target-volume={volume_reference}"]
            record.update(parent_campaign=str(args.restart_from), restart_global_step=19,
                          requested_states=args.updates + 1, requested_real_updates=args.updates,
                          requested_final_global_step=19 + args.updates,
                          original_target_volume=float(volume_reference),
                          initial_shape="original final global19 volume mesh; no geometric modification")
            record["state"] = "validating_restart_state"
            run("validation/restart", [*command, *restart, "--iterations=1"], 2)
            with (args.output / "validation/restart/kelvin-ball.csv").open() as stream:
                restarted = list(csv.DictReader(stream))
            if len(restarted) != 1 or any(not math.isclose(float(restarted[0][k]), float(previous[-1][k]),
                        rel_tol=1e-8, abs_tol=1e-10) for k in ("k", "c", "q", "volume")):
                raise ValueError("Restart state differs from the original final state")
            if float(restarted[0]["target_volume"]) != float(volume_reference):
                raise ValueError("Original volume reporting reference was lost")
            loaded_path = args.output / "validation/restart/KelvinBall.Chamber.mesh.000000.h5"
            with h5py.File(snapshot) as a, h5py.File(loaded_path) as b:
                for key in ("Mesh/Geometry/Vertices", "Mesh/XDMF/Topology", "Mesh/Attributes/3"):
                    if not np.array_equal(a[key][:], b[key][:]):
                        raise ValueError(f"The restart changed original geometry/connectivity/materials: {key}")
                for dim in (0, 1, 2):
                    if not np.array_equal(np.unique(a[f"Mesh/Attributes/{dim}"][:], return_counts=True),
                                          np.unique(b[f"Mesh/Attributes/{dim}"][:], return_counts=True)):
                        raise ValueError(f"Restart label histogram changed: dimension {dim}")
            record["validation"]["restart_state_passed"] = True
            record["validation"]["exact_restart_geometry_tets_materials_passed"] = True
            record["validation"]["restart_fem_diagnostics"] = restarted[0]
            command.extend([*restart, f"--iterations={args.updates + 1}"])
        else:
            command.append("--iterations=20")
        record["state"] = "optimizing"
        run("optimization", command, 8)
        record["state"] = "completed_native_budget"
        code = 0
    except (Exception, KeyboardInterrupt) as error:
        record.update(state="failed", error=str(error))
    finally:
        csv_path = args.output / "optimization/kelvin-ball.csv"
        if csv_path.exists():
            with csv_path.open() as stream:
                rows = [r for r in csv.DictReader(stream) if None not in r and all(value not in (None, "") for value in r.values())]
            record["fully_evaluated_states"] = len(rows)
            record["completed_real_updates"] = max(0, len(rows) - 1)
            if args.restart_from:
                record["latest_complete_global_step"] = 19 + max(0, len(rows) - 1)
            record["latest_fem_diagnostic"] = rows[-1] if rows else None
        record.update(exit_code=code, finished_at=time.time())
        save()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
