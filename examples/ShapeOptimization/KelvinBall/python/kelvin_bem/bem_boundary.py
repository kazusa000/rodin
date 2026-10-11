#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""All-panel BEM measurement and six outer-boundary velocity traces."""
import copy
import math
from pathlib import Path
import time

import numpy as np
from . import run_support as common

REFERENCE = common.EXAMPLE / "reproduction/clement-swift-bem.json"


def read_mesh(path):
    # Existing reader, with one-based Medit indices converted exactly once.
    from .mesh_io import read_medit_mesh
    x, faces, cells = read_medit_mesh(Path(path))
    faces, cells = np.asarray(faces, int), np.asarray(cells, int)
    return np.asarray(x), faces[:, :3]-1, faces[:, 3], cells[:, :4]-1, cells[:, 4]


def exterior_velocity(points, centers, areas, density, epsilon, chunk=16):
    """Same regularized single-layer kernel as the pinned BEM, all panels."""
    out = np.empty((len(points), 3, 6))
    weighted = density * areas[:, None, None]
    for start in range(0, len(points), chunk):
        d = points[start:start+chunk, None, :] - centers[None, :, :]
        r2 = np.einsum("mni,mni->mn", d, d)
        factor = 1 / (8*np.pi*(r2+epsilon**2)**1.5)
        iso = (r2+2*epsilon**2)*factor
        dot = np.einsum("mni,nil->mnl", d, weighted)
        out[start:start+chunk] = (np.einsum("mn,nil->mil", iso, weighted)
            + np.einsum("mn,mni,mnl->mil", factor, d, dot))
    return out


def resistance(matrix):
    k, q = float(np.trace(matrix[:3, :3])/3), float(np.trace(matrix[3:, 3:])/3)
    c = float(np.trace(matrix[3:, :3])/3)
    ct = float(np.trace(matrix[:3, 3:])/3)
    if not np.isfinite(matrix).all() or min(k, q) <= 0:
        raise ValueError("non-finite or non-positive body resistance")
    return dict(k=k, c=c, c_reciprocal=ct, q=q, rho=abs(c)/math.sqrt(k*q),
                reciprocal_gap=abs(c-ct)/math.sqrt(k*q), matrix=matrix.tolist())


def bem_trace(args, reference, mesh, dest):
    started = time.perf_counter()
    rtol, maxiter = bem_solver_settings(args, reference)
    source = common.source_state(common.ROOT)
    # Import the actual pinned verifier; never fall back to the older local scorer.
    from . import bem_solver as bem
    if Path(bem.__file__).resolve() != Path(__file__).with_name("bem_solver.py"):
        raise RuntimeError("wrong BEM verifier module")
    panels = bem.read_gamma_panels(mesh)
    rotations = bem.symmetry.rotations()
    if len(rotations) != (12 if args.symmetry == "terra" else 24):
        raise RuntimeError("BEM group mismatch")
    epsilon = reference["epsilon_factor"] * math.sqrt(np.mean([p.area for p in panels]))
    cpu_operator = bem.QuotientStokesletOperator(panels, rotations, epsilon=epsilon,
                                            viscosity=1, row_chunk=32)
    operator = cpu_operator
    cuda_library = getattr(args, "cuda_library", None)
    if cuda_library:
        from .cuda_operator import make_cuda_operator
        operator = make_cuda_operator(cpu_operator, library=cuda_library)
    rhs = bem._quotient_rhs(panels)
    try:
        solution, diagnostics = bem.solve_quotient_cg(operator, rhs,
            rtol=rtol, maxiter=maxiter,
            panel_block_size=reference["panel_block_size"], progress=True)
        if cuda_library:
            relative = np.linalg.norm(rhs-cpu_operator.matmat(solution), axis=0)/np.linalg.norm(rhs, axis=0)
            if not np.isfinite(relative).all() or np.any(relative > rtol):
                raise RuntimeError("CUDA solution failed independent CPU residual check")
            diagnostics.update(cpu_rhs_relative_residuals=relative.tolist(),
                               cuda=operator.cuda.metadata, cpu_residual_verified=True)
    finally:
        if cuda_library:
            operator.cuda.close()
    full = bem._expand_panels(panels, rotations)
    density = bem._expand_quotient_tractions(panels, rotations, solution, (0, 0, 0))
    centers = np.array([p.centroid for p in full])
    areas = np.array([p.area for p in full])
    matrix = np.empty((6, 6))
    for load in range(6):
        force, torque = bem.integrated_force_torque(full, density[:, load], center=(0, 0, 0))
        matrix[:, load] = np.r_[force, torque]
    density = density.reshape(len(full), 3, 6)
    x, faces, labels, _, _ = read_mesh(mesh)
    outer_ids = np.unique(faces[labels == 7])
    if len(outer_ids) == 0:
        raise ValueError("mesh has no Outer=7 faces")
    points = x[outer_ids]
    values = exterior_velocity(points, centers, areas, density, epsilon)
    table = np.c_[points, values.transpose(0, 2, 1).reshape(-1, 18)]
    np.savetxt(dest / "outer-values.txt", table, fmt="%.17g", header=str(len(table)), comments="")
    np.savez_compressed(dest / "density.npz", centers=centers, areas=areas,
                        density=density, epsilon=epsilon, quotient_solution=solution)
    result = dict(measurement="all original panels, free-space regularized BEM diagnostic",
        source=source, mesh_sha256=common.digest(mesh), group_order=len(rotations),
        chamber_panels=len(panels), full_panels=len(full), subsampled=False,
        epsilon=epsilon, epsilon_factor=reference["epsilon_factor"],
        device="CUDA FP64 with CPU residual check" if cuda_library else "CPU FP64",
        diagnostics=diagnostics, outer_nodes=len(points),
        solver_rtol=rtol, solver_maxiter=maxiter,
        elapsed_seconds=time.perf_counter()-started,
        table_sha256=common.digest(dest / "outer-values.txt"), **resistance(matrix))
    common.source_state(common.ROOT, source["commit"])
    common.save_json(dest / "bem.json", result)
    weak_rtol = getattr(args, "bem_weak_rtol", None)
    if weak_rtol is not None and rtol > weak_rtol and result["rho"] < 1e-4:
        # The optimization-grade residual can overwhelm tiny sphere coupling.
        # One bounded refinement; never accept the coarse observable silently.
        common.save_json(dest / "bem-fast.json", result)
        fine = copy.copy(args)
        fine.bem_rtol = weak_rtol
        refined = bem_trace(fine, reference, mesh, dest)
        refined["weak_coupling_refinement"] = dict(threshold=1e-4,
            requested_rtol=rtol, effective_rtol=weak_rtol,
            coarse_rho=result["rho"], rho_difference=abs(refined["rho"]-result["rho"]))
        refined["elapsed_seconds"] = time.perf_counter()-started
        common.save_json(dest / "bem.json", refined)
        return refined
    return result


def bem_solver_settings(args, reference):
    """Explicit optimization overrides do not modify the frozen reference."""
    rtol = getattr(args, "bem_rtol", None)
    maxiter = getattr(args, "bem_maxiter", None)
    rtol = reference["solver_rtol"] if rtol is None else rtol
    maxiter = reference["solver_maxiter"] if maxiter is None else maxiter
    if not math.isfinite(rtol) or not 0 < rtol <= 1e-3:
        raise ValueError("BEM rtol must be finite and in (0, 1e-3]")
    if not isinstance(maxiter, int) or maxiter <= 0:
        raise ValueError("BEM maxiter must be a positive integer")
    return rtol, maxiter
