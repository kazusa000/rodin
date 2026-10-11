#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Mac controller: pinned CPU workers on 28h and all-panel GPU BEM on 2060.

Transfers only exact meshes and small boundary/certificate files. It never
transfers HDF5 fields and never writes remote source trees.
"""

import argparse
import csv
import json
import math
import os
import re
from pathlib import Path, PurePosixPath
import shlex
import subprocess

from . import carlos_reference as common
from .clement_rho import armijo, fragmentation_gate, thickness_gate, step_factors


def shell(words):
    return shlex.join(map(str, words))


def known_rejection(stage, code, timeout, text):
    """Only explicit candidate-local numerical failures are backtracked."""
    if not timeout and stage == "native" and code in (-6, 134) and (
            "Mesh reconstruction changed the fixed boundary geometry." in text):
        return "invalid fixed-boundary geometry"
    if not timeout and stage == "bem" and code == 1 and re.search(
            r"quotient-CG did not converge in [1-9][0-9]* iterations;", text):
        return "candidate BEM reached the configured iteration limit"
    return None


class CommandFailure(RuntimeError):
    def __init__(self, code, timeout, log):
        super().__init__(f"command exit={code}, timeout={timeout}, log={log}")
        self.code, self.timeout, self.log = code, timeout, Path(log)


class CandidateFailure(RuntimeError):
    def __init__(self, stage, code, timeout, log):
        super().__init__(f"{stage} candidate failed: exit={code}, log={log}")
        self.stage, self.code, self.timeout, self.log = stage, code, timeout, Path(log)

    def rejection(self):
        return known_rejection(self.stage,self.code,self.timeout,
                               self.log.read_text(errors="replace"))

    def evidence(self):
        return dict(stage=self.stage, exit_code=self.code, timed_out=self.timeout,
                    log=str(self.log))


def decide(initial_bem, candidate_bem, initial_thickness, candidate_thickness,
           initial_geometry, candidate_geometry, direction):
    thickness_gate(initial_thickness, candidate_thickness, 0.06, "cubic")
    if not fragmentation_gate(initial_geometry, candidate_geometry, "cubic"):
        return False, "new disconnected solid components"
    if not math.isclose(direction["penalty"], initial_thickness["penalty"],
                        rel_tol=1e-10, abs_tol=1e-14):
        raise ValueError("native/independent discrete penalty mismatch")
    if not math.isclose(direction["rho"], initial_bem["rho"], rel_tol=1e-12, abs_tol=1e-18):
        raise ValueError("native did not consume the current BEM coefficients")
    initial = initial_bem["rho"] - initial_thickness["penalty"]
    candidate = candidate_bem["rho"] - candidate_thickness["penalty"]
    slope = direction["merit_slope"]
    accepted = slope > 0 and armijo(initial, candidate, direction["dt"], slope)
    return accepted, "actual BEM rho-P Armijo" if accepted else "BEM merit insufficient"


class Controller:
    def __init__(self, args):
        self.args = args
        self.output = common.external(args.output)
        self.output.mkdir(parents=True, exist_ok=False)
        self.source = common.source_state(common.ROOT)
        self.cpu = PurePosixPath(args.cpu_output)
        self.gpu = PurePosixPath(args.gpu_output)
        self.logs = self.output/"commands"
        self.logs.mkdir()
        self.counter = 0

    def execute(self, words):
        self.counter += 1
        log = self.logs/f"{self.counter:04d}.log"
        code, timeout = common.execute(words, cwd=self.output, env=os.environ.copy(),
            log=log, hours=1.2)
        if code or timeout:
            raise CommandFailure(code, timeout, log)

    def ssh(self, host, words):
        self.execute(["ssh", "-o", "BatchMode=yes", host, shell(words)])

    def get(self, host, remote, local):
        Path(local).parent.mkdir(parents=True, exist_ok=True)
        self.execute(["scp", "-q", f"{host}:{remote}", str(local)])
        return Path(local)

    def put(self, host, local, remote):
        self.ssh(host, ["mkdir", "-p", str(PurePosixPath(remote).parent)])
        self.execute(["scp", "-q", str(local), f"{host}:{remote}"])

    def native(self, label, stage, *options):
        output = self.cpu/label
        launch_error = None
        try:
            self.ssh(self.args.cpu_host, ["env", "KELVIN_SYMMETRY=cubic",
            self.args.cpu_python, str(PurePosixPath(self.args.cpu_source)/
                "examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py"), "native", stage,
            "--build-dir", self.args.build, "--output", output, *options])
        except CommandFailure as error:
            launch_error = error
        try:
            local = self.get(self.args.cpu_host, output/"run-record.json",
                             self.output/label/"run-record.json")
        except CommandFailure:
            if launch_error:
                raise launch_error
            raise
        record = json.loads(local.read_text())
        if not record["success"]:
            log = self.get(self.args.cpu_host,output/"run.log",self.output/label/"run.log")
            if len(record.get("cases",[])) == 1:
                case = record["cases"][0]
                raise CandidateFailure("native",case["exit_code"],case["timed_out"],log)
            raise RuntimeError(f"native stage failed before a candidate solve: {label}")
        return output

    def diagnostic(self, label, mesh, stage):
        remote = self.native(label, stage, "--mesh", mesh)
        filename = "geometry.json" if stage == "geometry" else "thickness.json"
        return json.loads(self.get(self.args.cpu_host, remote/filename,
                          self.output/label/filename).read_text())

    def bem(self, label, host, mesh):
        local_mesh = self.get(host, mesh, self.output/label/"input.mesh")
        remote_mesh = self.gpu/"inputs"/(label+".mesh")
        self.put(self.args.gpu_host, local_mesh, remote_mesh)
        output = self.gpu/label
        try:
            self.ssh(self.args.gpu_host, ["env", "KELVIN_SYMMETRY=cubic",
            "OMP_NUM_THREADS=2", "OPENBLAS_NUM_THREADS=1", self.args.gpu_python,
            str(PurePosixPath(self.args.gpu_source)/"examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py"), "bem",
            "--mesh", remote_mesh, "--output", output,
            "--cuda-library", self.args.cuda_library])
        except CommandFailure as error:
            # A worker certificate distinguishes solver limits from import,
            # version, credential and transport failures. Never restart it.
            try:
                worker = json.loads(self.get(self.args.gpu_host,output/"worker.json",
                    self.output/label/"worker.json").read_text())
            except CommandFailure:
                raise error
            if not worker.get("success") and known_rejection("bem",error.code,error.timeout,
                                                             worker.get("error","")):
                raise CandidateFailure("bem",error.code,error.timeout,error.log) from error
            raise error
        for name in ("bem.json", "worker.json", "outer-values.txt"):
            self.get(self.args.gpu_host, output/name, self.output/label/name)
        result = json.loads((self.output/label/"bem.json").read_text())
        if result["mesh_sha256"] != common.digest(local_mesh):
            raise ValueError("BEM input identity mismatch")
        return result

    def prepare(self):
        # All nodes must be visible online before opening SSH; no credential changes.
        status = json.loads(subprocess.check_output(
            [self.args.tailscale, "status", "--json"], text=True))
        online = {p["TailscaleIPs"][0] for p in status["Peer"].values() if p["Online"]}
        if status["BackendState"] != "Running" or not status["Self"]["Online"] or not {
                "100.89.230.120", "100.86.213.89"}.issubset(online):
            raise RuntimeError("Mac, 28h and 2060 must be Tailscale-online")
        for host, source in ((self.args.cpu_host, self.args.cpu_source),
                             (self.args.gpu_host, self.args.gpu_source)):
            state = subprocess.check_output(["ssh", "-o", "BatchMode=yes", host,
                shell(["git", "-C", source, "rev-parse", "HEAD"])], text=True).strip()
            dirty = subprocess.check_output(["ssh", "-o", "BatchMode=yes", host,
                shell(["git", "-C", source, "status", "--porcelain"])], text=True).strip()
            if dirty or state != self.args.worker_commit:
                raise ValueError("remote worker source must be pinned and clean")

    def run(self):
        record = dict(source=self.source, started_at=common.utcnow(), success=False,
                      native_source=self.args.cpu_source, gpu_source=self.args.gpu_source,
                      worker_commit=self.args.worker_commit, trials=[])
        def save():
            common.save_json(self.output/"campaign.json", record)
        try:
            self.prepare()
            if self.args.mesh:
                source = self.cpu/"inputs"/"accepted.mesh"
                self.put(self.args.cpu_host, self.args.mesh, source)
            else:
                initial = self.native("initial", "initial")
                source = initial/"KelvinBallInitial.mesh"
            initial_bem = self.bem("initial-bem", self.args.cpu_host, source)
            initial_thickness = self.diagnostic("initial-thickness", source, "thickness")
            initial_geometry = self.diagnostic("initial-geometry", source, "geometry")
            record.update(initial=initial_bem, initial_thickness=initial_thickness,
                          initial_geometry=initial_geometry)
            coefficients = self.output/"coefficients.txt"
            coefficients.write_text(" ".join(format(initial_bem[k], ".17g")
                                             for k in ("k", "c", "q"))+"\n")
            outer = self.cpu/"inputs"/"outer-values.txt"
            values = self.cpu/"inputs"/"coefficients.txt"
            self.put(self.args.cpu_host, self.output/"initial-bem"/"outer-values.txt", outer)
            self.put(self.args.cpu_host, coefficients, values)
            save()
            for method in ("mmg", "swift"):
                for trial, step in enumerate(step_factors(self.args.trials)):
                    label = f"{method}-trial-{trial}"
                    try:
                        remote = self.native(label, "step", "--mesh", source,
                            "--outer-values", outer, "--bem-values", values,
                            "--method", method, "--step", step)
                    except CandidateFailure as error:
                        reason = error.rejection()
                        if not reason:
                            raise
                        record["trials"].append(dict(method=method,trial=trial,
                            step_factor=step,accepted=False,rejection=reason,
                            failure=error.evidence()))
                        save()
                        continue
                    direction = json.loads(self.get(self.args.cpu_host,
                        remote/"boundary-step.json", self.output/label/"boundary-step.json").read_text())
                    thickness = self.diagnostic(label+"-thickness", remote/"candidate.mesh", "thickness")
                    geometry = self.diagnostic(label+"-geometry", remote/"candidate.mesh", "geometry")
                    entry = dict(method=method, trial=trial, step_factor=step,
                                 direction=direction, thickness=thickness, geometry=geometry,
                                 accepted=False)
                    record["trials"].append(entry)
                    if not fragmentation_gate(initial_geometry, geometry, "cubic"):
                        entry["rejection"] = "new disconnected solid components before BEM"
                    else:
                        try:
                            measured = self.bem(label+"-bem", self.args.cpu_host, remote/"candidate.mesh")
                        except CandidateFailure as error:
                            reason = error.rejection()
                            if not reason:
                                raise
                            entry.update(rejection=reason,failure=error.evidence())
                            save()
                            continue
                        accepted, reason = decide(initial_bem, measured, initial_thickness,
                            thickness, initial_geometry, geometry, direction)
                        entry.update(accepted=accepted, acceptance=reason, bem=measured,
                            predicted_delta=direction["dt"]*direction["merit_slope"],
                            actual_delta=measured["rho"]-initial_bem["rho"])
                    # Pair reconstruction against a SINGLE unchanged source/target.
                    if method == "mmg" and trial == 0:
                        history = self.get(self.args.cpu_host, remote/"kelvin-ball.csv",
                                           self.output/label/"kelvin-ball.csv")
                        with history.open() as stream:
                            h = float(next(csv.DictReader(stream))["h"])
                        finite_differences = []
                        for factor in (0.01, 0.005):
                            signed = []
                            for sign in (1, -1):
                                name = f"fd-{factor}-{sign}"
                                perturbation = self.native(name, "fd", "--mesh", source,
                                    "--direction", remote/"direction.txt", "--amount", sign*factor*h)
                                result = self.bem(name+"-bem", self.args.cpu_host,
                                                  perturbation/"diagnostic.mesh")
                                signed.append(result["c"]/math.sqrt(result["k"]*result["q"]))
                            measured = (signed[0]-signed[1])/(2*factor*h)
                            predicted = math.copysign(1, initial_bem["c"])*direction["d_rho_theta"]
                            finite_differences.append(dict(epsilon=factor*h,
                                signed_bem_derivative=measured, predicted_signed_derivative=predicted,
                                same_sign=measured*predicted > 0,
                                relative_gap=abs(measured-predicted)/max(abs(measured),abs(predicted),1e-20),
                                diagnostic_only=True))
                        record["finite_differences"] = finite_differences
                        pair = []
                        for backend in ("mmg", "swift"):
                            try:
                                target = self.native("same-target-"+backend, "reconstruct",
                                    "--mesh", remote/"transport.mesh", "--target", remote/"target.sol",
                                    "--method", backend)
                            except CandidateFailure as error:
                                reason = error.rejection()
                                if not reason:
                                    raise
                                pair.append(dict(method=backend,success=False,rejection=reason,
                                    failure=error.evidence(),
                                    run=json.loads((self.output/("same-target-"+backend)/
                                        "run-record.json").read_text())))
                                continue
                            stats = self.diagnostic("same-target-"+backend+"-geometry",
                                                    target/"candidate.mesh", "geometry")
                            pair.append(dict(method=backend, success=True, geometry=stats,
                                run=json.loads((self.output/("same-target-"+backend)/
                                    "run-record.json").read_text())))
                        if pair[0]["run"]["inputs"] != pair[1]["run"]["inputs"]:
                            raise ValueError("paired reconstruction did not consume identical inputs")
                        record["same_target"] = pair
                    save()
                    if entry["accepted"]:
                        break
            record["success"] = all(any(e["accepted"] for e in record["trials"]
                if e["method"] == method) for method in ("mmg", "swift"))
            record["status"] = "both_backends_accepted" if record["success"] else "line_search_stalled"
        except BaseException as error:
            record.update(status="failed", error=str(error))
            raise
        finally:
            record["finished_at"] = common.utcnow()
            save()
        return 0 if record["success"] else 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nodes-config", type=Path,
                   help="JSON with deployment paths; numerical settings remain in the frozen preset")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--mesh", type=Path)
    p.add_argument("--trials", type=int, default=6)
    p.add_argument("--cpu-host", default="p90318@28h")
    p.add_argument("--gpu-host", default="2060")
    p.add_argument("--tailscale", default="tailscale")
    for key in ("cpu-source", "gpu-source", "worker-commit", "cpu-python",
                "gpu-python", "build", "cpu-output", "gpu-output", "cuda-library"):
        p.add_argument("--"+key)
    args = p.parse_args()
    deployment = json.loads(args.nodes_config.read_text()) if args.nodes_config else {}
    required = ("cpu_source", "gpu_source", "cpu_python", "gpu_python",
                "build", "cpu_output", "gpu_output", "cuda_library")
    allowed = set(required) | {"cpu_host", "gpu_host", "tailscale", "worker_commit"}
    if set(deployment) - allowed:
        p.error("nodes configuration contains unknown or numerical settings")
    for key, value in deployment.items():
        flag = "--" + key.replace("_", "-")
        if not any(arg == flag or arg.startswith(flag + "=") for arg in os.sys.argv[1:]):
            setattr(args, key, value)
    if args.worker_commit is None:
        args.worker_commit = common.source_state(common.ROOT)["commit"]
    if args.worker_commit != common.source_state(common.ROOT)["commit"]:
        p.error("controller and both workers must use the same Rodin commit")
    if any(getattr(args, key) is None for key in required):
        p.error("supply all deployment paths through --nodes-config or explicit flags")
    return Controller(args).run()


if __name__ == "__main__":
    raise SystemExit(main())
