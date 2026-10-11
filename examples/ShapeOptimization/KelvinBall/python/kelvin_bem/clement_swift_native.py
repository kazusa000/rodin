#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Pinned native CPU worker for the opt-in Carlos SWIFT/BEM interface."""

import argparse
import json
import os
from pathlib import Path
import platform

from . import carlos_reference as common
from . import clement_swift_baseline as baseline

REFERENCE = common.EXAMPLE / "reproduction/clement-swift-bem.json"


def command(args, binary, reference):
    fixed = ["--n=13", "--outer-radius=2", "--thickness-min=0"]
    if args.stage == "initial":
        return [str(binary), *fixed, "--geometry-only", "--save-mesh",
                "--seed-amplitude=0.02", "--reconstruction=mmg"]
    if args.stage == "step":
        return [str(binary), *fixed, "--iterations=1", "--boundary-step",
                f"--mesh={args.mesh}", f"--outer-values={args.outer_values}",
                f"--bem-values={args.bem_values}", f"--reconstruction={args.method}",
                f"--step={args.step}", "--regularization=4",
                f"--thickness-absolute={reference['thickness_absolute']}",
                f"--thickness-weight={reference['thickness_weight']}"]
    if args.stage == "reconstruct":
        return [str(binary), *fixed, "--reconstruct-only", f"--mesh={args.mesh}",
                f"--target={args.target}", f"--reconstruction={args.method}"]
    if args.stage == "thickness":
        return [str(binary), f"--mesh={args.mesh}",
                f"--minimum={reference['thickness_absolute']}", f"--fd={args.fd}"]
    raise ValueError(args.stage)


def build(args, reference):
    project = common.source_state(common.ROOT)
    rodin = common.source_state(args.rodin, reference["rodin_commit"])
    modules = common.git(args.rodin, "submodule", "status", "third-party/mmg",
                         "third-party/termcolor")
    if any(line.startswith(("-", "+", "U")) for line in modules.splitlines()):
        raise ValueError("require exact initialized submodules")
    output = common.external(args.build_dir, args.rodin)
    output.mkdir(parents=True, exist_ok=False)
    commands = [baseline.configure_command(args, args.rodin, output),
                [args.cmake, "--build", str(output), "--target",
                 *reference["native_targets"], "-j", str(args.jobs)]]
    record = dict(project=project, rodin=rodin, submodules=modules,
                  host=platform.node(), commands=commands, started_at=common.utcnow(),
                  runtime_library_path=os.environ.get("LD_LIBRARY_PATH", ""), success=False)
    try:
        for index, cmd in enumerate(commands):
            common.save_json(output / "build-record.json", record)
            code, timeout = common.execute(cmd, cwd=output, env=os.environ.copy(),
                log=output / f"build-{index}.log", hours=args.timeout_hours)
            record.update(exit_code=code, timed_out=timeout)
            if code or timeout:
                break
        else:
            record["executables"] = {name: common.digest(baseline.binary_path(output, name))
                                      for name in reference["native_targets"]}
            common.source_state(args.rodin, reference["rodin_commit"])
            common.source_state(common.ROOT, project["commit"])
            record["success"] = True
    except BaseException as error:
        record["error"] = str(error)
        raise
    finally:
        record["finished_at"] = common.utcnow()
        common.save_json(output / "build-record.json", record)
    return 0 if record["success"] else 1


def run(args, reference):
    project = common.source_state(common.ROOT)
    built = json.loads((args.build_dir / "build-record.json").read_text())
    if not built["success"] or built["rodin"]["commit"] != reference["rodin_commit"]:
        raise ValueError("require the matching successful native build")
    common.source_state(built["rodin"]["path"], reference["rodin_commit"])
    for name, digest in built["executables"].items():
        if common.digest(baseline.binary_path(args.build_dir, name)) != digest:
            raise ValueError("native executable changed")
    output = common.external(args.output, built["rodin"]["path"])
    output.mkdir(parents=True, exist_ok=False)
    inputs = {str(getattr(args, key)): common.digest(getattr(args, key))
              for key in ("mesh", "outer_values", "bem_values", "target", "direction")
              if getattr(args, key, None) is not None}
    record = dict(project=project, rodin=built["rodin"], host=platform.node(),
                  inputs=inputs, build_sha256=common.digest(args.build_dir/"build-record.json"),
                  started_at=common.utcnow(), stage=args.stage, success=False, cases=[])
    env = dict(os.environ, LD_LIBRARY_PATH=built["runtime_library_path"],
               OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS="1")
    if args.stage == "fd":
        import numpy as np
        from .clement_rho import fd_mesh
        fd_mesh(args.mesh, np.loadtxt(args.direction), args.amount, output/"diagnostic.mesh")
        commands = []
    elif args.stage == "geometry":
        from .clement_rho import geometry_diagnostics
        common.save_json(output/"geometry.json", geometry_diagnostics(args.mesh, "cubic"))
        commands = []
    elif args.stage == "checks":
        commands = [(name, [str(baseline.binary_path(args.build_dir, name))])
                    for name in reference["native_targets"] if name.endswith(
                        ("Test", "Diagnostic"))]
    else:
        name = "KelvinBallSoftThicknessDiagnostic" if args.stage == "thickness" else "KelvinBall"
        commands = [(args.stage, command(args, baseline.binary_path(args.build_dir, name), reference))]
    try:
        for name, cmd in commands:
            cwd = output/name if args.stage == "checks" else output
            if cwd != output:
                cwd.mkdir()
            case = dict(name=name, command=cmd)
            record["cases"].append(case)
            common.save_json(output/"run-record.json", record)
            code, timeout = common.execute(cmd, cwd=cwd, env=env,
                log=cwd/"run.log", hours=args.timeout_hours)
            case.update(exit_code=code, timed_out=timeout, success=code == 0 and not timeout)
        for path, digest in inputs.items():
            if common.digest(path) != digest:
                raise ValueError("native input changed")
        common.source_state(common.ROOT, project["commit"])
        record["success"] = all(case["success"] for case in record["cases"])
    except BaseException as error:
        record["error"] = str(error)
        raise
    finally:
        record["finished_at"] = common.utcnow()
        common.save_json(output/"run-record.json", record)
        common.save_json(output/"exit.json", dict(exit_code=0 if record["success"] else 1))
    return 0 if record["success"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)
    b = sub.add_parser("build")
    b.add_argument("--rodin", type=Path, default=common.ROOT)
    b.add_argument("--cmake", default="cmake")
    b.add_argument("--cmake-arg", action="append", default=[])
    b.add_argument("--jobs", type=int, default=2)
    for name in ("initial", "step", "reconstruct", "checks", "thickness", "geometry", "fd"):
        p = sub.add_parser(name)
        p.add_argument("--output", type=Path, required=True)
        p.add_argument("--threads", type=int, default=2)
        if name in ("step", "reconstruct", "thickness", "geometry", "fd"):
            p.add_argument("--mesh", type=Path, required=True)
        if name in ("step", "reconstruct"):
            p.add_argument("--method", choices=("mmg", "swift"), required=True)
        if name == "step":
            p.add_argument("--outer-values", type=Path, required=True)
            p.add_argument("--bem-values", type=Path, required=True)
            p.add_argument("--step", type=float, default=0.015625)
        if name == "reconstruct":
            p.add_argument("--target", type=Path, required=True)
        if name == "thickness":
            p.add_argument("--fd", type=float, default=0)
        if name == "fd":
            p.add_argument("--direction", type=Path, required=True)
            p.add_argument("--amount", type=float, required=True)
    for p in sub.choices.values():
        p.add_argument("--build-dir", type=Path, required=True)
        p.add_argument("--timeout-hours", type=float, default=1)
    args = parser.parse_args()
    reference = json.loads(REFERENCE.read_text())
    reference["rodin_commit"] = common.source_state(common.ROOT)["commit"]
    if args.stage == "build" and args.rodin.resolve() != common.ROOT:
        parser.error("native source must be this single Rodin checkout")
    return build(args, reference) if args.stage == "build" else run(args, reference)


if __name__ == "__main__":
    raise SystemExit(main())
