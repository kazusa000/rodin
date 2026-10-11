#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Build and exercise the unchanged Carlos KelvinBall baseline in isolation."""

import argparse
import json
import math
import os
from pathlib import Path
import platform
import sys

from . import carlos_reference as common

REFERENCE = common.EXAMPLE / "reproduction/clement-swift-bem.json"


def selected_commit(args, reference):
    return reference["compiler_fixed_commit"] if args.compiler_fixed else reference["baseline_commit"]


def audit_compile_fix(rodin, reference):
    diff = common.git(rodin, "diff", reference["baseline_commit"], reference["compiler_fixed_commit"])
    changed = [line for line in diff.splitlines()
               if line.startswith(("+", "-")) and not line.startswith(("+++", "---"))]
    expected = ["-        decltype(classifier)::Parameters parameters;",
                "+        typename decltype(classifier)::Parameters parameters;"]
    if changed != expected:
        raise ValueError("compiler-fixed baseline has changes beyond the audited typename fix")
    return dict(upstream_commit=reference["baseline_commit"], patch=diff,
                numerical_changes=False)


def configure_command(args, rodin, dest):
    return [args.cmake, "-S", str(rodin), "-B", str(dest),
            "-DCMAKE_BUILD_TYPE=Release", "-DRODIN_USE_MUMPS=ON",
            "-DRODIN_USE_PETSC=OFF", "-DRODIN_BUILD_EXAMPLES=ON",
            "-DRODIN_BUILD_UNIT_TESTS=OFF", "-DRODIN_BUILD_MANUFACTURED_TESTS=OFF",
            "-DRODIN_BUILD_CONVERGENCE_TESTS=OFF", "-DRODIN_BUILD_BENCHMARKS=OFF",
            "-DRODIN_BUILD_DOC=OFF", "-DRODIN_WITH_PY=OFF",
            "-DRODIN_INSTALL_RESOURCES=OFF", *args.cmake_arg]


def binary_path(build_dir, name):
    return Path(build_dir) / "examples/ShapeOptimization/KelvinBall" / name


def case_commands(build_dir, reference):
    sphere = [str(binary_path(build_dir, "KelvinBallSphere")), "--n=13",
              "--outer-radius=2", "--save-mesh"]
    return [("sphere", sphere), *[(method, [str(binary_path(build_dir, "KelvinBall")),
             *reference["baseline_arguments"], f"--reconstruction={method}"])
             for method in ("mmg", "swift")]]


def update_audit(summary):
    if summary["evaluated_designs"] != 2:
        raise ValueError("baseline needs the initial state and one evaluated real update")
    first, last = summary["first"], summary["last"]
    if last["actual_delta_rho"] is None or last["actual_delta_volume"] is None:
        raise ValueError("missing incoming update diagnostics")
    if not math.isclose(last["actual_delta_rho"], last["rho"]-first["rho"],
                        rel_tol=1e-9, abs_tol=1e-14):
        raise ValueError("inconsistent actual rho change")
    if first["advection_increment"] is None or first["advection_increment"] <= 0:
        raise ValueError("no real level-set transport took place")
    if first["stage_9_seconds"] is None or first["stage_9_seconds"] <= 0:
        raise ValueError("no reconstruction took place")
    return dict(evaluated_updates=1, actual_delta_rho=last["actual_delta_rho"],
                predicted_delta_rho=last["predicted_delta_rho"],
                volume_error=last["volume_error"],
                first_minimum_quality=first["mesh_quality_min"],
                last_minimum_quality=last["mesh_quality_min"],
                direct_solver=last["direct_solver"], bem_verified=False)


def build(args, reference):
    state = common.source_state(common.ROOT)
    selected = selected_commit(args, reference)
    rodin = common.source_state(args.rodin, selected)
    compiler_fix = audit_compile_fix(args.rodin, reference) if args.compiler_fixed else None
    modules = common.git(args.rodin, "submodule", "status", "third-party/mmg",
                         "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in modules.splitlines()):
        raise ValueError("initialize the exact MMG and termcolor submodules first")
    dest = common.external(args.build_dir, args.rodin)
    dest.mkdir(parents=True, exist_ok=False)
    commands = [configure_command(args, args.rodin, dest),
                [args.cmake, "--build", str(dest), "--target",
                 *reference["baseline_targets"], "-j", str(args.jobs)]]
    record = dict(stage="Carlos baseline build", source=state, rodin=rodin,
                  compiler_fix=compiler_fix,
                  submodules=modules, started_at=common.utcnow(), host=platform.node(),
                  commands=commands, success=False,
                  runtime_library_path=os.environ.get("LD_LIBRARY_PATH", ""))
    code = 1
    try:
        for index, command in enumerate(commands):
            common.save_json(dest / "build-record.json", record)
            code, timeout = common.execute(command, cwd=dest, env=os.environ.copy(),
                log=dest / f"build-{index}.log", hours=args.timeout_hours)
            record.update(exit_code=code, timed_out=timeout)
            if code or timeout:
                break
        else:
            record["executables"] = {name: common.digest(binary_path(dest, name))
                                      for name in reference["baseline_targets"]}
            common.source_state(args.rodin, selected)
            common.source_state(common.ROOT, state["commit"])
            record["success"] = True
    except BaseException as error:
        record["error"] = str(error)
        raise
    finally:
        record["finished_at"] = common.utcnow()
        common.save_json(dest / "build-record.json", record)
    return 0 if record["success"] else 1


def run(args, reference):
    built = json.loads((args.build_dir / "build-record.json").read_text())
    selected = selected_commit(args, reference)
    if not built.get("success") or built["rodin"]["commit"] != selected:
        raise ValueError("require a successful pinned Carlos baseline build")
    common.source_state(built["rodin"]["path"], selected)
    if args.compiler_fixed:
        audit_compile_fix(built["rodin"]["path"], reference)
    for name, digest in built["executables"].items():
        if common.digest(binary_path(args.build_dir, name)) != digest:
            raise ValueError(f"baseline binary changed: {name}")
    dest = common.external(args.output, built["rodin"]["path"])
    dest.mkdir(parents=True, exist_ok=False)
    record = dict(stage="Carlos baseline", source=common.source_state(common.ROOT),
                  compiler_fix=built.get("compiler_fix"),
                  rodin=built["rodin"], build_sha256=common.digest(args.build_dir / "build-record.json"),
                  started_at=common.utcnow(), host=platform.node(), cases=[], success=False)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", KELVIN_SYMMETRY="cubic")
    try:
        for name, command in case_commands(args.build_dir, reference):
            case_dir = dest / name
            case_dir.mkdir()
            case = dict(name=name, command=command, started_at=common.utcnow(), success=False)
            record["cases"].append(case)
            common.save_json(dest / "campaign.json", record)
            code, timeout = common.execute(command, cwd=case_dir, env=env,
                log=case_dir / "run.log", hours=args.timeout_hours)
            case.update(exit_code=code, timed_out=timeout, finished_at=common.utcnow())
            if code == 0 and not timeout:
                if name != "sphere":
                    case["summary"] = common.summary(case_dir / "kelvin-ball.csv")
                    case["update_audit"] = update_audit(case["summary"])
                case["success"] = True
            # Collect both backend outcomes even if one fails. Never restart it.
            common.save_json(case_dir / "run-record.json", case)
        common.source_state(built["rodin"]["path"], selected)
        record["success"] = all(case["success"] for case in record["cases"])
    except BaseException as error:
        record["error"] = str(error)
        raise
    finally:
        record["finished_at"] = common.utcnow()
        common.save_json(dest / "campaign.json", record)
        common.save_json(dest / "exit.json", dict(exit_code=0 if record["success"] else 1))
    return 0 if record["success"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)
    b = sub.add_parser("build")
    b.add_argument("--rodin", type=Path, required=True)
    b.add_argument("--cmake", default="cmake")
    b.add_argument("--cmake-arg", action="append", default=[])
    b.add_argument("--jobs", type=int, default=2)
    r = sub.add_parser("run")
    r.add_argument("--output", type=Path, required=True)
    r.add_argument("--threads", type=int, default=2)
    for p in (b, r):
        p.add_argument("--build-dir", type=Path, required=True)
        p.add_argument("--timeout-hours", type=float, default=4)
        p.add_argument("--compiler-fixed", action="store_true",
                       help="explicitly use the audited single-line GCC type-disambiguation fix")
    args = parser.parse_args()
    return globals()[args.stage](args, json.loads(REFERENCE.read_text()))


if __name__ == "__main__":
    raise SystemExit(main())
