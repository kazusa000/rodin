#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Panel/regularization convergence scan for the selected-symmetry BEM scorer."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Callable, Iterable

SCRIPT_DIR = Path(__file__).resolve().parent

from . import bem_kcq_score as bem


def parse_int_list(raw: str) -> list[int]:
    values = [0 if part.lower() == "full" else int(part) for part in raw.replace(",", " ").split()]
    if not values or any(value < 0 for value in values):
        raise ValueError("integer list must contain nonnegative values (0/full means all panels)")
    return values


def parse_float_list(raw: str) -> list[float]:
    values = [float(part) for part in raw.replace(",", " ").split()]
    if not values or any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("float list must contain finite positive values")
    return values


def run_scan(
    mesh: Path,
    *,
    panels: Iterable[int],
    epsilon_factors: Iterable[float],
    viscosity: float,
    center: bem.Point,
    backend: str = "full",
    subdivision_levels: Iterable[int] = (0,),
    max_matrix_gib: float = 2.0,
    checkpoint: Callable[[list[dict[str, object]]], None] | None = None,
    solver_options: dict[str, object] | None = None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    panel_counts, epsilons = list(panels), list(epsilon_factors)
    if backend not in {"full", "quotient", "quotient-cg"}:
        raise ValueError("backend must be full, quotient, or quotient-cg")
    if not panel_counts or not epsilons:
        raise ValueError("panel and epsilon lists must be nonempty")
    if any(not math.isfinite(value) or value <= 0.0 for value in epsilons):
        raise ValueError("epsilon factors must be finite and positive")
    if not math.isfinite(max_matrix_gib) or max_matrix_gib <= 0:
        raise ValueError("max_matrix_gib must be finite and positive")
    original_count = len(bem.read_gamma_triangles(mesh))
    for level in subdivision_levels:
        if level < 0:
            raise ValueError("subdivision levels must be nonnegative")
        raw_count = original_count * 4**level
        for panel_count in panel_counts:
            if panel_count < 0:
                raise ValueError("panel count must be nonnegative")
            used_count = min(panel_count, raw_count) if panel_count else raw_count
            dimension = used_count * (9 if backend != "full" else 3 * bem.symmetry.ORDER)
            matrix_bytes = 0 if backend == "quotient-cg" else dimension**2 * 8
            for epsilon_factor in epsilons:
                row: dict[str, object] = {
                    "mesh": str(mesh), "symmetry": bem.symmetry.NAME,
                    "subdivision_levels": level, "requested_panels": panel_count,
                    "panels": used_count, "raw_panels": raw_count,
                    "original_panels": original_count,
                    "full_panels": used_count * bem.symmetry.ORDER,
                    "epsilon_factor": epsilon_factor, "backend": backend,
                    "system_dimension": dimension, "matrix_bytes": matrix_bytes,
                    "equivalent_dense_matrix_bytes": dimension**2 * 8,
                    "status": "skipped_matrix_memory_limit",
                }
                if matrix_bytes <= max_matrix_gib * 1024**3:
                    row.update(_score_row(mesh, panel_count, epsilon_factor, viscosity,
                                          center, backend, level, solver_options or {}))
                rows.append(row)
                if checkpoint is not None:
                    checkpoint(rows)
    return rows


def _score_row(mesh, panel_count, epsilon_factor, viscosity, center, backend, level, solver_options):
    result = bem.score_kcq(
        mesh,
        max_chamber_panels=panel_count,
        epsilon_factor=epsilon_factor,
        viscosity=viscosity,
        center=center,
        backend=backend,
        subdivision_levels=level,
        **solver_options,
    )
    return {
        "status": "ok",
        "epsilon": float(result["epsilon"]),
        "k": float(result["k"]),
        "c": float(result["c"]),
        "q": float(result["q"]),
        "rho": float(result["rho"]),
        "eK": float(result["eK"]),
        "eC": float(result["eC"]),
        "eQ": float(result["eQ"]),
        "reciprocity_error": float(result["reciprocity_error"]),
        "solver": str(result["solver"]),
        "linear_relative_residual": float(result["linear_relative_residual"]),
        "elapsed_seconds": float(result["elapsed_seconds"]),
        "solver_converged": result.get("solver_converged", result["solver"] == "solve"),
        "solver_rtol": result.get("solver_rtol"),
        "rhs_relative_residuals": result.get("rhs_relative_residuals"),
        "iterations_by_rhs": result.get("iterations_by_rhs"),
        "operator_applications": result.get("operator_applications"),
        "preconditioner_panel_block_size": result.get("preconditioner_panel_block_size"),
        "preconditioner_bytes": result.get("preconditioner_bytes"),
        "matrix_free_device": result.get("matrix_free_device"),
        "gpu_rhs_relative_residuals": result.get("gpu_rhs_relative_residuals"),
        "cpu_rhs_relative_residuals": result.get("cpu_rhs_relative_residuals"),
        "cpu_residual_verified": result.get("cpu_residual_verified"),
        "cuda": result.get("cuda"),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "mesh",
        "symmetry",
        "subdivision_levels",
        "requested_panels",
        "panels",
        "raw_panels",
        "original_panels",
        "full_panels",
        "epsilon_factor",
        "epsilon",
        "k",
        "c",
        "q",
        "rho",
        "eK",
        "eC",
        "eQ",
        "reciprocity_error",
        "solver",
        "backend",
        "system_dimension",
        "matrix_bytes",
        "equivalent_dense_matrix_bytes",
        "linear_relative_residual",
        "elapsed_seconds",
        "status",
        "solver_converged",
        "solver_rtol",
        "rhs_relative_residuals",
        "iterations_by_rhs",
        "operator_applications",
        "preconditioner_panel_block_size",
        "preconditioner_bytes",
        "matrix_free_device",
        "gpu_rhs_relative_residuals",
        "cpu_rhs_relative_residuals",
        "cpu_residual_verified",
        "cuda",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_name(path.name + ".pending")
    with pending.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    pending.replace(path)


def write_json(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n", encoding="ascii")
    pending.replace(path)


def print_table(rows: list[dict[str, object]], *, header: bool = True) -> None:
    if header:
        print("subdivision panels epsilon_factor k c q rho eK eC eQ", flush=True)
    for row in rows:
        if row["status"] != "ok":
            print(f"{row['subdivision_levels']} {row['panels']} {row['epsilon_factor']} "
                  f"{row['status']} matrix_bytes={row['matrix_bytes']}", flush=True)
            continue
        print(
            f"{row['subdivision_levels']} "
            f"{row['panels']} "
            f"{row['epsilon_factor']:.8g} "
            f"{row['k']:.12g} "
            f"{row['c']:.12g} "
            f"{row['q']:.12g} "
            f"{row['rho']:.12g} "
            f"{row['eK']:.4g} "
            f"{row['eC']:.4g} "
            f"{row['eQ']:.4g}", flush=True)


def provenance(mesh: Path) -> dict[str, object]:
    root = SCRIPT_DIR.parent
    def git(*args):
        try:
            return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unavailable"
    return {
        "mesh_sha256": hashlib.sha256(mesh.read_bytes()).hexdigest(),
        "source_commit": git("rev-parse", "HEAD"),
        "source_status": git("status", "--porcelain"),
        "scorer_sha256": hashlib.sha256((SCRIPT_DIR / "bem_kcq_score.py").read_bytes()).hexdigest(),
        "python": sys.version, "numpy": bem.np.__version__, "platform": platform.platform(),
        "thread_environment": {name: os.environ.get(name) for name in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }


def convergence_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    """Descriptive spread only: not a continuum-accuracy certificate."""
    valid = [row for row in rows if row["status"] == "ok"]
    all_panel_rows = [row for row in valid if row["panels"] == row["raw_panels"]]
    result: dict[str, object] = {
        "completed_evaluations": len(valid), "skipped_evaluations": len(rows) - len(valid),
        "all_panels_evaluations": len(all_panel_rows),
        "least_squares_fallbacks": sum(row["solver"] == "lstsq" for row in valid),
        "matrix_free_solves": sum(row["solver"] == "pcg" for row in valid),
        "note": "Panel and regularization sensitivity on the fixed piecewise-planar mesh; not a continuum error bound.",
    }
    if all_panel_rows:
        finest = max(row["subdivision_levels"] for row in all_panel_rows)
        fine_rows = [row for row in all_panel_rows if row["subdivision_levels"] == finest]
        rho_values = [float(row["rho"]) for row in fine_rows]
        result.update({"finest_all_panels_subdivision": finest,
                       "finest_all_panels_rho_min": min(rho_values),
                       "finest_all_panels_rho_max": max(rho_values),
                       "finest_all_panels_rho_relative_spread":
                           (max(rho_values) - min(rho_values)) / max(max(rho_values), 1e-300)})
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--panels", default="256,512,full", help="0 or full uses every (possibly subdivided) Gamma panel.")
    parser.add_argument("--epsilon-factors", default="0.4,0.6,0.8")
    parser.add_argument("--viscosity", type=float, default=1.0)
    parser.add_argument("--center", type=float, nargs=3, default=(0.0, 0.0, 0.0))
    parser.add_argument("--backend", choices=("full", "quotient", "quotient-cg"), default="full")
    parser.add_argument("--subdivision-levels", default="0", help="Comma-separated nonnegative levels; each level quadruples source panels.")
    parser.add_argument("--max-matrix-gib", type=float, default=2.0,
                        help="Skip larger matrices explicitly. Factorization peak memory is greater than this allocation.")
    parser.add_argument("--solver-rtol", type=float, default=1e-10)
    parser.add_argument("--solver-maxiter", type=int, default=500)
    parser.add_argument("--matrix-free-row-chunk", type=int, default=32)
    parser.add_argument("--solver-progress", action="store_true")
    parser.add_argument("--preconditioner-panel-block-size", type=int, default=16)
    parser.add_argument("--matrix-free-device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--cuda-library")
    parser.add_argument("--cuda-device-index", type=int, default=0)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--json", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = provenance(args.mesh)
    metadata.update({"arguments": vars(args) | {"mesh": str(args.mesh), "csv": str(args.csv), "json": str(args.json)},
                     "completed": False})
    def checkpoint(rows):
        if args.csv is not None:
            write_csv(args.csv, rows)
        if args.json is not None:
            write_json(args.json, rows)
            metadata["convergence"] = convergence_summary(rows)
            args.json.with_suffix(".metadata.json").write_text(
                json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="ascii",
            )
    def progress(rows):
        checkpoint(rows)
        print_table(rows[-1:], header=len(rows) == 1)
    checkpoint([])
    rows = run_scan(
        args.mesh,
        panels=parse_int_list(args.panels),
        epsilon_factors=parse_float_list(args.epsilon_factors),
        viscosity=args.viscosity,
        center=tuple(args.center),
        backend=args.backend,
        subdivision_levels=parse_int_list(args.subdivision_levels),
        max_matrix_gib=args.max_matrix_gib,
        checkpoint=progress,
        solver_options={
            "solver_rtol": args.solver_rtol,
            "solver_maxiter": args.solver_maxiter,
            "matrix_free_row_chunk": args.matrix_free_row_chunk,
            "solver_progress": args.solver_progress,
            "preconditioner_panel_block_size": args.preconditioner_panel_block_size,
            "matrix_free_device": args.matrix_free_device,
            "cuda_library": args.cuda_library,
            "cuda_device_index": args.cuda_device_index,
        },
    )
    metadata["completed"] = True
    checkpoint(rows)


if __name__ == "__main__":
    main()
