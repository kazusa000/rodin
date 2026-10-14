#!/usr/bin/env python3
"""Assess preserved FEM states with original HDF5 coordinates and independent CPU BEM."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

import h5py
import numpy as np

EXAMPLE = Path(__file__).resolve().parents[1]
ROOT = EXAMPLE.parents[2]
BEM_COMMIT = "12f38a9a33dbbb0a1b7c0d0affcfff91620b927a"


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def prepare(snapshot, native, target):
    lines = native.read_text().splitlines()

    def section(name):
        start = lines.index(name) + 1
        while not lines[start].strip():
            start += 1
        count = int(lines[start])
        return start + 1, count

    start, count = section("Vertices")
    original = np.array([[float(x) for x in line.split()[:3]] for line in lines[start:start+count]])
    with h5py.File(snapshot, "r") as data:
        coordinates = np.asarray(data["/Mesh/Geometry/Vertices"], dtype=float)
        topology = np.asarray(data["/Mesh/XDMF/Topology"], dtype=np.int64).reshape(-1, 5)
        labels = np.asarray(data["/Mesh/Attributes/3"], dtype=np.int64).reshape(-1)
    if (coordinates.shape != original.shape or not np.isfinite(coordinates).all() or
            not np.allclose(coordinates, original, rtol=6e-6, atol=1e-9) or
            not np.all(topology[:, 0] == 6)):
        raise ValueError("Native mesh and HDF5 snapshot do not identify the same coordinates.")
    tet_start, tet_count = section("Tetrahedra")
    native_tets = np.array([[int(x) for x in line.split()] for line in lines[tet_start:tet_start+tet_count]])
    native_tets[:, :4] -= 1
    hdf_tets = np.column_stack([topology[:, 1:], labels])
    native_tets[:, :4].sort(axis=1)
    hdf_tets[:, :4].sort(axis=1)
    for tets in [native_tets, hdf_tets]:
        order = np.lexsort(tuple(tets[:, i] for i in reversed(range(5))))
        tets[:] = tets[order]
    if not np.array_equal(native_tets, hdf_tets):
        raise ValueError("Native mesh and HDF5 connectivity or material labels differ.")
    tri_start, tri_count = section("Triangles")
    gamma = np.array([[int(x)-1 for x in line.split()[:3]]
                      for line in lines[tri_start:tri_start+tri_count] if int(line.split()[-1]) == 13])
    if len(gamma) == 0:
        raise ValueError("No original Gamma panels.")
    for index, point in enumerate(coordinates):
        label = lines[start+index].split()[3]
        lines[start+index] = " ".join(format(float(x), ".17g") for x in point) + " " + label
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n")
    certificate = dict(snapshot=str(snapshot), snapshot_sha256=digest(snapshot),
        native_ascii=str(native), native_ascii_sha256=digest(native), output=str(target),
        output_sha256=digest(target), topology_unchanged=True, geometry_source="original HDF5 doubles",
        gamma_panels=len(gamma), geometry_modified=False,
        native_ascii_max_coordinate_error=float(np.max(np.abs(original-coordinates))),
        native_ascii_gamma_coordinate_error=float(np.max(np.abs(original[gamma]-coordinates[gamma]))))
    target.with_suffix(".json").write_text(json.dumps(certificate, indent=2) + "\n")
    return certificate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--fem-campaign", type=Path, required=True)
    parser.add_argument("--verifier", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--backend", choices=("quotient", "quotient-cg"), default="quotient")
    parser.add_argument("--blas-threads", type=int, default=6)
    args = parser.parse_args()
    if platform.node() != "28h":
        raise ValueError("Independent evaluation is explicitly CPU-only on 28h.")
    for source, commit in [(ROOT, args.commit), (args.verifier, BEM_COMMIT)]:
        if git(source, "rev-parse", "HEAD") != commit or git(source, "status", "--porcelain", "--untracked-files=all"):
            raise ValueError("Require exact clean pinned controller and verifier.")
    fem_record = args.fem_campaign.read_bytes()
    fem = json.loads(fem_record)
    stage = next(s for s in fem["stages"] if s["name"] == "optimization")
    if stage.get("exit_code") != 0 or fem.get("real_updates") != 5:
        raise ValueError("Require five completed native FEM updates before evaluation.")
    simulation = args.fem_campaign.parent / "optimization"
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "native-fem-campaign.json").write_bytes(fem_record)
    if args.blas_threads < 1:
        raise ValueError("BLAS threads must be positive.")
    threads = str(args.blas_threads if args.backend == "quotient" else 1)
    env = dict(os.environ, KELVIN_SYMMETRY="cubic", OPENBLAS_NUM_THREADS=threads,
               MKL_NUM_THREADS=threads, OMP_NUM_THREADS="6", PYTHONDONTWRITEBYTECODE="1")
    record = dict(controller_commit=args.commit, native_fem_commit=fem["commit"],
                  fem_campaign_sha256=hashlib.sha256(fem_record).hexdigest(), source=str(ROOT),
                  verifier_commit=BEM_COMMIT, device="cpu", backend=args.backend,
                  blas_threads=int(threads), host=platform.node(),
                  started_at=time.time(), evaluations=[], state="independent_bem")

    def save():
        (args.output / "campaign.json").write_text(json.dumps(record, indent=2) + "\n")

    code = 1
    try:
        for step in [0, 1, 3, 5]:
            directory = args.output / f"step-{step:06d}"
            directory.mkdir()
            mesh = directory / "original-double.mesh"
            certificate = prepare(simulation / f"KelvinBall.Chamber.mesh.{step:06d}.h5",
                                  simulation / f"KelvinBallState-{step:06d}.mesh", mesh)
            command = [os.sys.executable, "-B", str(args.verifier / "scripts/bem_convergence_scan.py"), str(mesh),
                "--panels", "full", "--epsilon-factors", "0.6", "--subdivision-levels", "0",
                "--viscosity", "1", "--center", "0", "0", "0", "--backend", args.backend,
                "--max-matrix-gib", "12",
                "--solver-rtol", "1e-10", "--solver-maxiter", "6000", "--preconditioner-panel-block-size", "256",
                "--matrix-free-row-chunk", "32", "--solver-progress", "--matrix-free-device", "cpu",
                "--json", str(directory / "scan.json"), "--csv", str(directory / "scan.csv")]
            item = dict(step=step, input=certificate, command=command, started_at=time.time())
            record["evaluations"].append(item)
            save()
            with (directory / "run.log").open("w") as log:
                child = subprocess.run(command, cwd=directory, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, timeout=43200)
            item.update(exit_code=child.returncode, finished_at=time.time())
            (directory / "exit-code.txt").write_text(str(child.returncode) + "\n")
            if child.returncode:
                raise RuntimeError(f"BEM step {step} actual exit {child.returncode}")
            result = json.loads((directory / "scan.json").read_text())[0]
            if (result.get("status") != "ok" or not result.get("solver_converged") or
                    result.get("linear_relative_residual", float("inf")) > 1e-10):
                raise RuntimeError(f"BEM step {step} did not converge: {result}")
            item["result"] = result
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
