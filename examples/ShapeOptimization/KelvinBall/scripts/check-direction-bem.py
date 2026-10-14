#!/usr/bin/env python3
"""Fixed-panel signed BEM directional differences; no optimization update."""
import argparse
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import subprocess

import numpy as np

EXAMPLE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("independent_bem", EXAMPLE / "scripts/evaluate-fem-bem.py")
bem = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bem)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--fem-campaign", type=Path, required=True)
    parser.add_argument("--verifier", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if platform.node() != "28h":
        raise ValueError("This diagnostic must run on 28h.")
    for source, commit in [(bem.ROOT, args.commit), (args.verifier, bem.BEM_COMMIT)]:
        if bem.git(source, "rev-parse", "HEAD") != commit or bem.git(source, "status", "--porcelain", "--untracked-files=all"):
            raise ValueError("Require exact clean diagnostic and BEM source.")
    native = json.loads(args.fem_campaign.read_text())
    simulation = args.fem_campaign.parent / "optimization"
    args.output.mkdir(parents=True, exist_ok=False)
    base = args.output / "base-original-double.mesh"
    certificate = bem.prepare(simulation / "KelvinBall.Chamber.mesh.000000.h5",
                              simulation / "KelvinBallState-000000.mesh", base)
    lines = base.read_text().splitlines()
    start = lines.index("Vertices") + 1
    while not lines[start].strip():
        start += 1
    count = int(lines[start])
    start += 1
    points = np.array([[float(v) for v in line.split()[:3]] for line in lines[start:start+count]])
    direction_file = simulation / "direction-000000.txt"
    direction = np.loadtxt(direction_file)
    if direction.shape != points.shape:
        raise ValueError("Direction and mesh numbering differ.")
    row = next(csv.DictReader((simulation / "kelvin-ball.csv").open()))
    predicted = math.copysign(1, float(row["c"])) * float(row["d_rho_theta"])
    record = dict(controller_commit=args.commit, native_fem_commit=native["commit"],
        verifier_commit=bem.BEM_COMMIT, direction_sha256=bem.digest(direction_file),
        base=certificate, diagnostic="fixed-connectivity signed eta central differences",
        includes_transport_or_mmg=False, fem_signed_prediction=predicted, differences=[], state="running")
    env = dict(os.environ, KELVIN_SYMMETRY="cubic", OPENBLAS_NUM_THREADS="2",
               MKL_NUM_THREADS="2", OMP_NUM_THREADS="2", PYTHONDONTWRITEBYTECODE="1")

    def save():
        (args.output / "directional-fd.json").write_text(json.dumps(record, indent=2) + "\n")

    code = 1
    try:
        for factor in [0.01, 0.005]:
            delta = factor * float(row["h"])
            values = []
            for sign in [1, -1]:
                directory = args.output / f"factor-{factor}-sign-{sign}"
                directory.mkdir()
                mesh = directory / "perturbed.mesh"
                altered = list(lines)
                for index, point in enumerate(points + sign * delta * direction):
                    label = lines[start+index].split()[3]
                    altered[start+index] = " ".join(format(float(v), ".17g") for v in point) + " " + label
                mesh.write_text("\n".join(altered) + "\n")
                command = [os.sys.executable, "-B", str(args.verifier / "scripts/bem_convergence_scan.py"), str(mesh),
                    "--panels", "full", "--epsilon-factors", "0.6", "--subdivision-levels", "0", "--backend", "quotient",
                    "--max-matrix-gib", "12", "--json", str(directory / "scan.json"), "--csv", str(directory / "scan.csv")]
                (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
                with (directory / "run.log").open("w") as log:
                    child = subprocess.run(command, cwd=directory, env=env, stdout=log,
                                           stderr=subprocess.STDOUT, timeout=600)
                (directory / "exit-code.txt").write_text(str(child.returncode) + "\n")
                if child.returncode:
                    raise RuntimeError("Directional BEM solve failed.")
                result = json.loads((directory / "scan.json").read_text())[0]
                if result.get("status") != "ok" or result.get("linear_relative_residual", float("inf")) > 1e-10:
                    raise RuntimeError("Directional BEM solve did not meet residual tolerance.")
                values.append(result["c"] / math.sqrt(result["k"] * result["q"]))
            slope = (values[0]-values[1]) / (2*delta)
            record["differences"].append(dict(delta=delta, plus_signed_eta=values[0], minus_signed_eta=values[1],
                bem_signed_slope=slope, absolute_prediction_difference=abs(slope-predicted)))
            save()
        record["state"] = "completed"
        code = 0
    except Exception as error:
        record.update(state="failed", error=str(error))
    finally:
        record["exit_code"] = code
        save()
        (args.output / "exit-code.txt").write_text(str(code) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
