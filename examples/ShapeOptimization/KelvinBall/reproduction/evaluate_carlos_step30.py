#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Independently assess the displayed L10 step30; never launch or alter FEM."""
import argparse
import json
import os
from pathlib import Path
import platform
import signal
import sys
import time
import traceback

from run_carlos_l10 import complete_history, git, pinned, sha, write, UTILITIES, VERIFIER

FEM = "2511f3bf0b1f655927c65e577ac4b40eceb4ce34"
SNAPSHOT_SHA = "84283b5c488a22a2c3cb85f7132705b798a563e704173958bdbeda125d9ff7c1"
BINARY_SHA = "f5c81a35f062f666ff578eca7891090b57be4b317c84f414ad2b4b19de54f4a5"


def validate_input(campaign):
    parent = json.loads(campaign.read_text())
    if parent.get("commit") != FEM or parent.get("binary_sha256") != BINARY_SHA:
        raise ValueError("Require the frozen Carlos L10 FEM experiment")
    simulation = campaign.parent / "optimization"
    rows = complete_history(simulation / "kelvin-ball.csv")
    row = next((row for row in rows if int(row["iteration"]) == 30), None)
    if row is None:
        raise ValueError("Displayed step30 has no complete FEM row")
    snapshot = simulation / "KelvinBall.Chamber.mesh.000030.h5"
    if not snapshot.is_file() or sha(snapshot) != SNAPSHOT_SHA:
        raise ValueError("Step30 HDF5 differs from the displayed shape")
    return row, snapshot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("fem-campaign", "output", "utilities", "verifier"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if platform.node() != "28h":
        raise RuntimeError("Assessment is authorized only on 28h")
    source = Path(__file__).resolve().parents[4]
    for path, commit in ((source, args.commit), (args.utilities, UTILITIES), (args.verifier, VERIFIER)):
        pinned(path, commit)
    row, snapshot = validate_input(args.fem_campaign)
    if args.output.exists() or any(args.output.resolve().is_relative_to(path.resolve())
            for path in (source, args.utilities, args.verifier, args.fem_campaign.parent)):
        raise ValueError("Require fresh output separate from sources and original FEM")
    args.output.mkdir(parents=True)
    sys.path.insert(0, str(args.utilities / "scripts"))
    from reproduce_smallbox_original import export_snapshot, execute
    from evaluate_saved_smallbox import validate_scan
    record = dict(status="preparing", started_at=time.time(), host=platform.node(), pid=os.getpid(),
                  source=str(source), commit=args.commit, fem_commit=FEM,
                  fem_campaign=str(args.fem_campaign), step=30, fem_row=row,
                  source_snapshot=str(snapshot), source_snapshot_sha256=SNAPSHOT_SHA,
                  utilities=UTILITIES, verifier=VERIFIER, solver_rtol=1e-6,
                  full_original_panels=True, optimization_interrupted=False,
                  optimization_restarted=False, prepare_only=args.prepare_only,
                  limitations=["Original geometric seams and full panel topology preserved",
                               "Fixed BEM regularization epsilon factor .6"])

    def checkpoint():
        write(args.output / "campaign.json", record)

    def stop(signum, frame):
        raise KeyboardInterrupt(f"Received signal {signum}")

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    code = 1
    try:
        checkpoint()
        certificate = export_snapshot(snapshot, args.output / "original-double.mesh")
        write(args.output / "original-double.json", certificate)
        record["input"] = certificate
        if certificate["source_sha256"] != SNAPSHOT_SHA:
            raise ValueError("Export source does not match the displayed shape")
        if args.prepare_only:
            record["status"] = "prepared_only"
            code = 0
        else:
            record["status"] = "independent_bem_running"
            checkpoint()
            env = dict(os.environ, KELVIN_SYMMETRY="cubic", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
            command = [sys.executable, "-B", str(args.verifier / "scripts/bem_convergence_scan.py"),
                       str(args.output / "original-double.mesh"), "--panels", "full",
                       "--subdivision-levels", "0", "--epsilon-factors", "0.6", "--viscosity", "1",
                       "--center", "0", "0", "0", "--backend", "quotient-cg", "--solver-rtol", "1e-6",
                       "--solver-maxiter", "6000", "--preconditioner-panel-block-size", "256",
                       "--matrix-free-row-chunk", "32", "--solver-progress", "--matrix-free-device", "cpu",
                       "--json", str(args.output / "scan.json"), "--csv", str(args.output / "scan.csv")]
            stage = dict(status="running")
            record["solver"] = stage
            actual = execute(command, args.output / "solver", env, 72, stage, checkpoint)
            stage["status"] = "completed" if actual == 0 else "failed"
            if actual:
                raise RuntimeError(f"Independent step30 BEM exited {actual}")
            result = validate_scan(args.output, args.verifier, certificate, solver_rtol=1e-6)
            fresh, _ = validate_input(args.fem_campaign)
            if fresh != row or sha(args.output / "original-double.mesh") != certificate["output_sha256"]:
                raise ValueError("Original FEM row or BEM geometry changed")
            record["result"] = {key: result[key] for key in
                    ("k", "c", "q", "rho", "linear_relative_residual", "rhs_relative_residuals")}
            record["status"] = "completed_with_valid_bem"
            code = 0
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
