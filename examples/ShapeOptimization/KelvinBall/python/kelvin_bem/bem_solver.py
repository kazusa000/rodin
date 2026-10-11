#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Infinite-domain regularized-Stokeslet KCQ scorer for chamber Gamma meshes."""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from itertools import permutations, product
from pathlib import Path

import sys
from . import symmetry as symmetry

import numpy as np


GAMMA = 13
Point = tuple[float, float, float]
Matrix = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
Triangle = tuple[Point, Point, Point]


@dataclass(frozen=True)
class Panel:
    centroid: Point
    area: float
    normal: Point


def det3(matrix: Matrix) -> int:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def proper_octahedral_rotations() -> list[Matrix]:
    rotations: list[Matrix] = []
    for perm in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            matrix = [[0, 0, 0] for _ in range(3)]
            for row, col in enumerate(perm):
                matrix[row][col] = signs[row]
            rotation = tuple(tuple(row) for row in matrix)  # type: ignore[assignment]
            if det3(rotation) == 1:
                rotations.append(rotation)
    return sorted(set(rotations))


def mat_vec(matrix: Matrix, point: Point) -> Point:
    return tuple(
        float(sum(matrix[row][col] * point[col] for col in range(3)))
        for row in range(3)
    )  # type: ignore[return-value]


def read_gamma_triangles(mesh: Path) -> list[Triangle]:
    lines = [line.strip() for line in mesh.read_text(encoding="ascii", errors="replace").splitlines()]
    vertices: list[Point] = []
    triangles: list[tuple[int, int, int]] = []
    i = 0
    while i < len(lines):
        token = lines[i]
        if token == "Vertices":
            count = int(lines[i + 1])
            for line in lines[i + 2 : i + 2 + count]:
                x, y, z, *_ = line.split()
                vertices.append((float(x), float(y), float(z)))
            i += 2 + count
            continue
        if token == "Triangles":
            count = int(lines[i + 1])
            for line in lines[i + 2 : i + 2 + count]:
                a, b, c, attr = line.split()[:4]
                if int(attr) == GAMMA:
                    triangles.append((int(a) - 1, int(b) - 1, int(c) - 1))
            i += 2 + count
            continue
        i += 1
    if not vertices:
        raise ValueError(f"mesh has no vertices: {mesh}")
    if not triangles:
        raise ValueError(f"mesh has no Gamma triangles: {mesh}")
    return [(vertices[a], vertices[b], vertices[c]) for a, b, c in triangles]


def subdivide_triangles(triangles: list[Triangle], levels: int) -> list[Triangle]:
    """Refine quadrature geometry without changing the piecewise-planar surface."""
    if levels < 0:
        raise ValueError("subdivision levels must be nonnegative")
    for _ in range(levels):
        refined = []
        for a, b, c in triangles:
            ab = tuple((np.asarray(a) + np.asarray(b)) / 2.0)
            bc = tuple((np.asarray(b) + np.asarray(c)) / 2.0)
            ca = tuple((np.asarray(c) + np.asarray(a)) / 2.0)
            refined.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
        triangles = refined
    return triangles


def read_gamma_panels(mesh: Path) -> list[Panel]:
    return [triangle_panel(*triangle) for triangle in read_gamma_triangles(mesh)]


def triangle_panel(a: Point, b: Point, c: Point) -> Panel:
    ab = np.array(b) - np.array(a)
    ac = np.array(c) - np.array(a)
    cross = np.cross(ab, ac)
    norm = float(np.linalg.norm(cross))
    area = 0.5 * norm
    normal = tuple((cross / norm).tolist()) if norm > 0.0 else (0.0, 0.0, 0.0)
    centroid = tuple(((np.array(a) + np.array(b) + np.array(c)) / 3.0).tolist())
    return Panel(centroid=centroid, area=area, normal=normal)  # type: ignore[arg-type]


def downsample_panel_indices(panels: list[Panel], target: int) -> list[int]:
    if target <= 0 or len(panels) <= target:
        return list(range(len(panels)))
    points = np.array([panel.centroid for panel in panels], dtype=float)
    areas = np.array([panel.area for panel in panels], dtype=float)
    centroid = np.average(points, axis=0, weights=np.maximum(areas, 1e-300))
    selected: list[int] = [int(np.argmax(np.sum((points - centroid) ** 2, axis=1)))]
    min_dist2 = np.sum((points - points[selected[0]]) ** 2, axis=1)
    while len(selected) < target:
        score = min_dist2 * np.sqrt(np.maximum(areas, 0.0))
        score[selected] = -1.0
        index = int(np.argmax(score))
        selected.append(index)
        dist2 = np.sum((points - points[index]) ** 2, axis=1)
        min_dist2 = np.minimum(min_dist2, dist2)
    return sorted(selected)


def downsample_panels(panels: list[Panel], target: int, *, rescale_area: bool = False) -> list[Panel]:
    if target <= 0 or len(panels) <= target:
        return panels
    indices = downsample_panel_indices(panels, target)
    sampled = [panels[int(index)] for index in indices]
    if not rescale_area:
        return sampled
    original_area = sum(panel.area for panel in panels)
    sampled_area = sum(panel.area for panel in sampled)
    if sampled_area <= 0.0:
        return sampled
    scale = original_area / sampled_area
    return [
        Panel(centroid=panel.centroid, area=panel.area * scale, normal=panel.normal)
        for panel in sampled
    ]


def _expand_panels(chamber_panels: list[Panel], rotations: list[Matrix]) -> list[Panel]:
    full: list[Panel] = []
    for rotation in rotations:
        for panel in chamber_panels:
            full.append(
                Panel(
                    centroid=mat_vec(rotation, panel.centroid),
                    area=panel.area,
                    normal=mat_vec(rotation, panel.normal),
                )
            )
    return full


def regularized_stokeslet_block(
    x: Point,
    y: Point,
    *,
    epsilon: float,
    viscosity: float,
) -> list[list[float]]:
    r = np.array(x, dtype=float) - np.array(y, dtype=float)
    r2 = float(np.dot(r, r))
    eps2 = epsilon * epsilon
    denom = (r2 + eps2) ** 1.5
    block = ((r2 + 2.0 * eps2) * np.eye(3) + np.outer(r, r)) / denom
    block *= 1.0 / (8.0 * math.pi * viscosity)
    return block.tolist()


def assemble_regularized_stokeslet_matrix(
    panels: list[Panel],
    *,
    epsilon: float,
    viscosity: float,
) -> np.ndarray:
    points = np.array([panel.centroid for panel in panels], dtype=float)
    areas = np.array([panel.area for panel in panels], dtype=float)
    diff = points[:, None, :] - points[None, :, :]
    r2 = np.sum(diff * diff, axis=2)
    eps2 = epsilon * epsilon
    denom = np.power(r2 + eps2, 1.5)
    coeff = 1.0 / (8.0 * math.pi * viscosity)
    n = len(panels)
    matrix = np.empty((3 * n, 3 * n), dtype=float)
    for row in range(3):
        for col in range(3):
            block = diff[:, :, row] * diff[:, :, col] / denom
            if row == col:
                block = block + (r2 + 2.0 * eps2) / denom
            matrix[row::3, col::3] = coeff * block * areas[None, :]
    return matrix


def rigid_rhs(panels: list[Panel], load: int, *, center: Point) -> np.ndarray:
    points = np.array([panel.centroid for panel in panels], dtype=float)
    rhs = np.zeros((len(panels), 3), dtype=float)
    if load < 3:
        rhs[:, load] = 1.0
    else:
        omega = np.zeros(3, dtype=float)
        omega[load - 3] = 1.0
        rhs = np.cross(np.broadcast_to(omega, points.shape), points - np.array(center, dtype=float))
    return rhs.reshape(3 * len(panels))


def integrated_force_torque(
    panels: list[Panel],
    traction: np.ndarray,
    *,
    center: Point,
) -> tuple[np.ndarray, np.ndarray]:
    forces = traction.reshape(len(panels), 3) * np.array([panel.area for panel in panels])[:, None]
    points = np.array([panel.centroid for panel in panels], dtype=float)
    force = np.sum(forces, axis=0)
    torque = np.sum(np.cross(points - np.array(center, dtype=float), forces), axis=0)
    return force, torque


def cross_matrix(point: Point) -> np.ndarray:
    """Matrix representing point cross vector."""
    x, y, z = point
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def assemble_quotient_stokeslet_matrix(
    panels: list[Panel],
    rotations: list[Matrix],
    *,
    epsilon: float,
    viscosity: float,
    row_chunk: int = 32,
) -> np.ndarray:
    """Exact restriction of the full discrete operator to rigid-load responses.

    For either three translation or three rotation loads, the response tensor
    obeys F(g x) = g F(x) g.T.  Each chamber panel consequently has nine
    unknowns, not 3 * len(rotations).  Unknowns and equations are weighted by
    sqrt(panel area), making this restriction symmetric.  No surface averaging,
    additional quadrature approximation, or regularization change is introduced.
    """
    if row_chunk <= 0:
        raise ValueError("row_chunk must be positive")
    n = len(panels)
    points = np.asarray([panel.centroid for panel in panels], dtype=float)
    sqrt_areas = np.sqrt([panel.area for panel in panels])
    matrix = np.zeros((9 * n, 9 * n), dtype=float)
    blocks = matrix.reshape(n, 3, 3, n, 3, 3)
    identity = np.eye(3)
    coeff = 1.0 / (8.0 * math.pi * viscosity)
    for rotation in rotations:
        g = np.asarray(rotation, dtype=float)
        rotated_points = points @ g.T
        for start in range(0, n, row_chunk):
            stop = min(n, start + row_chunk)
            diff = points[start:stop, None, :] - rotated_points[None, :, :]
            r2 = np.sum(diff * diff, axis=-1)
            kernel = diff[..., :, None] * diff[..., None, :]
            kernel += (r2 + 2.0 * epsilon**2)[..., None, None] * identity
            kernel *= (coeff / (r2 + epsilon**2)**1.5)[..., None, None]
            kernel = kernel @ g
            kernel *= (sqrt_areas[start:stop, None] * sqrt_areas[None, :])[..., None, None]
            # Both supported groups consist of signed permutation matrices.
            # Keep the load-column indices explicit to avoid a 9x temporary.
            for load in range(3):
                column = int(np.argmax(np.abs(g[load])))
                blocks[start:stop, :, load, :, :, column] += (
                    g[load, column] * kernel.transpose(0, 2, 1, 3)
                )
    return matrix


class QuotientStokesletOperator:
    """Matrix-free application of exactly the same weighted quotient matrix.

    Each call shares kernel evaluations across both three-load response
    tensors.  Only row_chunk x n pairwise distances are stored at one time.
    The block-Jacobi preconditioner includes all self-orbit image interactions,
    not merely the central diagonal regularized Stokeslet.
    """

    def __init__(self, panels, rotations, *, epsilon, viscosity, row_chunk=32):
        if row_chunk <= 0:
            raise ValueError("row_chunk must be positive")
        self.points = np.asarray([panel.centroid for panel in panels], dtype=float)
        self.panels = panels
        self.sqrt_areas = np.sqrt([panel.area for panel in panels])
        self.rotations = [np.asarray(g, dtype=float) for g in rotations]
        self.rotated_points = [self.points @ g.T for g in self.rotations]
        self.epsilon = epsilon
        self.viscosity = viscosity
        self.coeff = 1.0 / (8.0 * math.pi * viscosity)
        self.row_chunk = row_chunk
        self.n = len(panels)
        self.applications = 0

    def matmat(self, values: np.ndarray) -> np.ndarray:
        if values.ndim != 2 or values.shape[0] != 9 * self.n:
            raise ValueError("quotient operator expects a (9*n, nrhs) array")
        self.applications += 1
        nrhs = values.shape[1]
        response = values.reshape(self.n, 3, 3, nrhs)
        output = np.zeros((self.n, 3, 3 * nrhs), dtype=float)
        eps2 = self.epsilon**2
        for g, rotated_points in zip(self.rotations, self.rotated_points):
            rotated = np.einsum("ab,nbck,dc->nadk", g, response, g, optimize=True)
            rotated = (rotated * self.sqrt_areas[:, None, None, None]).reshape(self.n, 3, 3 * nrhs)
            flat = np.ascontiguousarray(rotated.reshape(self.n, -1))
            components = [np.ascontiguousarray(rotated[:, axis, :]) for axis in range(3)]
            for start in range(0, self.n, self.row_chunk):
                stop = min(self.n, start + self.row_chunk)
                diff = self.points[start:stop, None, :] - rotated_points[None, :, :]
                r2 = np.sum(diff * diff, axis=-1)
                beta = self.coeff / (r2 + eps2)**1.5
                isotropic = (r2 + 2.0 * eps2) * beta
                chunk = (isotropic @ flat).reshape(stop - start, 3, 3 * nrhs)
                # The dyadic kernel is symmetric; reuse its six unique entries.
                for row in range(3):
                    for col in range(row, 3):
                        factor = beta * diff[:, :, row] * diff[:, :, col]
                        chunk[:, row, :] += factor @ components[col]
                        if row != col:
                            chunk[:, col, :] += factor @ components[row]
                output[start:stop] += chunk
        output *= self.sqrt_areas[:, None, None]
        return output.reshape(9 * self.n, nrhs)

    def inverse_diagonal_blocks(self) -> np.ndarray:
        blocks = np.zeros((self.n, 3, 3, 3, 3), dtype=float)
        eps2 = self.epsilon**2
        for g, rotated_points in zip(self.rotations, self.rotated_points):
            diff = self.points - rotated_points
            r2 = np.sum(diff * diff, axis=-1)
            kernel = diff[:, :, None] * diff[:, None, :]
            kernel += (r2 + 2.0 * eps2)[:, None, None] * np.eye(3)
            kernel *= (self.coeff / (r2 + eps2)**1.5)[:, None, None]
            kernel = kernel @ g
            kernel *= (self.sqrt_areas**2)[:, None, None]
            blocks += np.einsum("nab,lm->nalbm", kernel, g)
        blocks = blocks.reshape(self.n, 9, 9)
        # Symmetrization affects only the preconditioner, not the operator.
        blocks = 0.5 * (blocks + blocks.transpose(0, 2, 1))
        try:
            np.linalg.cholesky(blocks)
            return np.linalg.inv(blocks)
        except np.linalg.LinAlgError as error:
            raise RuntimeError("quotient-CG self-orbit block preconditioner is not positive definite") from error

    def make_preconditioner(self, panel_block_size=16):
        if panel_block_size <= 0:
            raise ValueError("preconditioner panel block size must be positive")
        if panel_block_size == 1:
            inverse_blocks = self.inverse_diagonal_blocks()
            def apply(values):
                return (inverse_blocks @ values.reshape(self.n, 9, -1)).reshape(values.shape)
            return apply, inverse_blocks.nbytes
        # Nearest neighbors in the quotient metric can straddle identified
        # chamber faces. These blocks remain principal submatrices of A.
        unassigned = np.ones(self.n, dtype=bool)
        blocks = []
        while np.any(unassigned):
            seed = int(np.flatnonzero(unassigned)[0])
            distance = np.full(self.n, np.inf)
            for rotated_points in self.rotated_points:
                difference = rotated_points - self.points[seed]
                distance = np.minimum(distance, np.sum(difference * difference, axis=1))
            available = np.flatnonzero(unassigned)
            count = min(panel_block_size, len(available))
            selected = available[np.argsort(distance[available], kind="stable")[:count]]
            unassigned[selected] = False
            block = assemble_quotient_stokeslet_matrix(
                [self.panels[index] for index in selected], self.rotations,
                epsilon=self.epsilon, viscosity=self.viscosity,
            )
            block = 0.5 * (block + block.T)
            try:
                np.linalg.cholesky(block)
                inverse = np.linalg.inv(block)
            except np.linalg.LinAlgError as error:
                raise RuntimeError("quotient-CG geometric block preconditioner is not positive definite") from error
            indices = (9 * selected[:, None] + np.arange(9)).ravel()
            blocks.append((indices, inverse))
        def apply(values):
            output = np.empty_like(values)
            for indices, inverse in blocks:
                output[indices] = inverse @ values[indices]
            return output
        return apply, sum(inverse.nbytes + indices.nbytes for indices, inverse in blocks)


def solve_quotient_cg(
    operator: QuotientStokesletOperator,
    rhs: np.ndarray,
    *,
    rtol: float = 1e-10,
    maxiter: int = 500,
    progress: bool = False,
    panel_block_size: int = 16,
) -> tuple[np.ndarray, dict[str, object]]:
    """Independent PCG recurrences sharing the two-RHS operator application.

    Check the true residual every 50 iterations and whenever a recursive
    residual appears converged. Preserve conjugacy unless the residual gap is
    material (10% of the true residual) or recursive convergence was spurious;
    restarting at every diagnostic check can prevent difficult SPD systems
    from converging. Only a fresh true residual certifies convergence. No
    least-squares or loose-tolerance fallback is permitted.
    """
    if not math.isfinite(rtol) or rtol <= 0.0 or maxiter <= 0:
        raise ValueError("CG rtol and maxiter must be positive")
    if rhs.ndim != 2 or rhs.shape[0] != 9 * operator.n:
        raise ValueError("CG right-hand sides must have shape (9*n, nrhs)")
    precondition, preconditioner_bytes = operator.make_preconditioner(panel_block_size)
    rhs_norm = np.linalg.norm(rhs, axis=0)
    denominator = np.where(rhs_norm > 0, rhs_norm, 1.0)
    solution = np.zeros_like(rhs)
    residual = rhs.copy()
    z = precondition(residual)
    direction = z.copy()
    rz = np.sum(residual * z, axis=0)
    active = rhs_norm > 0
    iterations_by_rhs = np.zeros(rhs.shape[1], dtype=int)
    history = []
    for iteration in range(1, maxiter + 1):
        if not np.any(active):
            break
        if np.any(~np.isfinite(rz[active])) or np.any(rz[active] <= 0.0):
            raise RuntimeError(f"quotient-CG invalid preconditioned residual at iteration {iteration}")
        direction[:, ~active] = 0.0
        action = operator.matmat(direction)
        curvature = np.sum(direction * action, axis=0)
        if np.any(~np.isfinite(curvature[active])) or np.any(curvature[active] <= 0.0):
            raise RuntimeError(f"quotient-CG nonpositive/nonfinite curvature at iteration {iteration}")
        alpha = np.zeros_like(rz)
        alpha[active] = rz[active] / curvature[active]
        solution += direction * alpha
        residual -= action * alpha
        recursive_relative = np.linalg.norm(residual, axis=0) / denominator
        restart = np.zeros(rhs.shape[1], dtype=bool)
        refresh = (iteration % 50 == 0 or iteration == maxiter
                   or np.any(active & (recursive_relative <= rtol)))
        if refresh:
            true_residual = rhs - operator.matmat(solution)
            relative = np.linalg.norm(true_residual, axis=0) / denominator
            if np.any(~np.isfinite(relative)):
                raise RuntimeError(f"quotient-CG nonfinite true residual at iteration {iteration}")
            gap = np.linalg.norm(true_residual - residual, axis=0) / denominator
            previously_active = active.copy()
            newly_converged = active & (relative <= rtol)
            iterations_by_rhs[newly_converged] = iteration
            active = relative > rtol
            restart = active & (~previously_active | (gap > 0.1 * relative)
                                | (recursive_relative <= rtol))
            residual[:, restart] = true_residual[:, restart]
            history.append({"iteration": iteration, "true_relative_residuals": relative.tolist(),
                            "relative_residual_gaps": gap.tolist(),
                            "restarted_rhs": np.flatnonzero(restart).tolist()})
            if progress:
                print(f"quotient-cg iteration={iteration} true_relative_residuals={relative.tolist()} "
                      f"restarted_rhs={np.flatnonzero(restart).tolist()}",
                      file=sys.stderr, flush=True)
            if not np.any(active):
                break
        z = precondition(residual)
        next_rz = np.sum(residual * z, axis=0)
        beta = np.zeros_like(rz)
        continuing = active & ~restart
        beta[continuing] = next_rz[continuing] / rz[continuing]
        direction = z + direction * beta
        rz = next_rz
    true_residual = rhs - operator.matmat(solution)
    relative = np.linalg.norm(true_residual, axis=0) / denominator
    if np.any(~np.isfinite(relative)) or np.any(relative > rtol):
        raise RuntimeError(f"quotient-CG did not converge in {maxiter} iterations; "
                           f"true relative residuals={relative.tolist()}, required={rtol}")
    return solution, {
        "solver_converged": True,
        "solver_rtol": rtol,
        "solver_maxiter": maxiter,
        "iterations_by_rhs": iterations_by_rhs.tolist(),
        "rhs_relative_residuals": relative.tolist(),
        "linear_relative_residual": float(np.linalg.norm(true_residual) / max(np.linalg.norm(rhs), 1e-300)),
        "operator_applications": operator.applications,
        "residual_history": history,
        "residual_restart_policy": "gap_above_10_percent_or_false_recursive_convergence",
        "preconditioner": "quotient-geometric-block-jacobi",
        "preconditioner_panel_block_size": panel_block_size,
        "preconditioner_bytes": preconditioner_bytes,
        "matrix_free_row_chunk": operator.row_chunk,
    }


def _quotient_rhs(panels: list[Panel]) -> np.ndarray:
    sqrt_areas = np.sqrt([panel.area for panel in panels])
    rhs = np.empty((len(panels), 3, 3, 2), dtype=float)
    rhs[:, :, :, 0] = np.eye(3)
    rhs[:, :, :, 1] = [-cross_matrix(panel.centroid) for panel in panels]
    rhs *= sqrt_areas[:, None, None, None]
    return rhs.reshape(9 * len(panels), 2)


def _expand_quotient_tractions(panels, rotations, solution, center):
    sqrt_areas = np.sqrt([panel.area for panel in panels])
    response = solution.reshape(len(panels), 3, 3, 2) / sqrt_areas[:, None, None, None]
    # omega x (x-center) = omega x x + center x omega.
    shift = cross_matrix(center)
    full_response = []
    for rotation in rotations:
        g = np.asarray(rotation, dtype=float)
        translation = g @ response[:, :, :, 0] @ g.T
        rotation_response = g @ response[:, :, :, 1] @ g.T + translation @ shift
        full_response.append(np.concatenate((translation, rotation_response), axis=-1))
    return np.concatenate(full_response, axis=0).reshape(-1, 6)


def solve_quotient_tractions(
    panels: list[Panel],
    rotations: list[Matrix],
    *,
    epsilon: float,
    viscosity: float,
    center: Point,
) -> tuple[np.ndarray, str, float]:
    """Return full-orbit six-load tractions, in the original panel ordering."""
    matrix = assemble_quotient_stokeslet_matrix(
        panels, rotations, epsilon=epsilon, viscosity=viscosity,
    )
    rhs = _quotient_rhs(panels)
    try:
        solution = np.linalg.solve(matrix, rhs)
        solver = "solve"
    except np.linalg.LinAlgError:
        solution = np.linalg.lstsq(matrix, rhs, rcond=1e-10)[0]
        solver = "lstsq"
    relative_residual = float(np.linalg.norm(matrix @ solution - rhs) / np.linalg.norm(rhs))
    return _expand_quotient_tractions(panels, rotations, solution, center), solver, relative_residual


def anisotropy_error(matrix: np.ndarray) -> float:
    trace_average = float(np.trace(matrix) / 3.0)
    scale = max(1.0, abs(trace_average))
    return float(np.linalg.norm(matrix - trace_average * np.eye(3), ord="fro") / scale)


def score_kcq(
    mesh: Path,
    *,
    max_chamber_panels: int,
    epsilon_factor: float,
    viscosity: float,
    center: Point,
    backend: str = "full",
    subdivision_levels: int = 0,
    solver_rtol: float = 1e-10,
    solver_maxiter: int = 500,
    matrix_free_row_chunk: int = 32,
    solver_progress: bool = False,
    preconditioner_panel_block_size: int = 16,
    matrix_free_device: str = "cpu",
    cuda_library: str | None = None,
    cuda_device_index: int = 0,
) -> dict[str, object]:
    triangles = read_gamma_triangles(mesh)
    original_count = len(triangles)
    chamber_panels = [triangle_panel(*triangle) for triangle in subdivide_triangles(triangles, subdivision_levels)]
    result = score_panels(
        chamber_panels,
        mesh_label=str(mesh),
        max_chamber_panels=max_chamber_panels,
        epsilon_factor=epsilon_factor,
        viscosity=viscosity,
        center=center,
        backend=backend,
        solver_rtol=solver_rtol,
        solver_maxiter=solver_maxiter,
        matrix_free_row_chunk=matrix_free_row_chunk,
        solver_progress=solver_progress,
        preconditioner_panel_block_size=preconditioner_panel_block_size,
        matrix_free_device=matrix_free_device,
        cuda_library=cuda_library,
        cuda_device_index=cuda_device_index,
    )
    result["chamber_panels_original"] = original_count
    result["surface_subdivision_levels"] = subdivision_levels
    return result


def score_panels(
    chamber_panels: list[Panel],
    *,
    mesh_label: str,
    max_chamber_panels: int,
    epsilon_factor: float,
    viscosity: float,
    center: Point,
    backend: str = "full",
    solver_rtol: float = 1e-10,
    solver_maxiter: int = 500,
    matrix_free_row_chunk: int = 32,
    solver_progress: bool = False,
    preconditioner_panel_block_size: int = 16,
    matrix_free_device: str = "cpu",
    cuda_library: str | None = None,
    cuda_device_index: int = 0,
) -> dict[str, object]:
    started = time.perf_counter()
    if backend not in {"full", "quotient", "quotient-cg"}:
        raise ValueError("backend must be full, quotient, or quotient-cg")
    if matrix_free_device not in {"cpu", "cuda"}:
        raise ValueError("matrix-free device must be cpu or cuda")
    if matrix_free_device == "cuda" and backend != "quotient-cg":
        raise ValueError("CUDA matmat is only available with quotient-cg")
    if matrix_free_device == "cpu" and (cuda_library is not None or cuda_device_index != 0):
        raise ValueError("CUDA settings require explicit --matrix-free-device cuda")
    if not chamber_panels or any(panel.area <= 0.0 for panel in chamber_panels):
        raise ValueError("BEM requires nonempty, positive-area panels")
    if not math.isfinite(epsilon_factor) or epsilon_factor <= 0.0:
        raise ValueError("epsilon_factor must be finite and positive")
    if not math.isfinite(viscosity) or viscosity <= 0.0:
        raise ValueError("viscosity must be finite and positive")
    sampled = downsample_panels(chamber_panels, max_chamber_panels, rescale_area=True)
    rotations = (proper_tetrahedral_rotations() if symmetry.IS_TERRA
                 else proper_octahedral_rotations())
    panels = _expand_panels(sampled, rotations)
    mean_area = sum(panel.area for panel in panels) / len(panels)
    epsilon = epsilon_factor * math.sqrt(mean_area)
    diagnostics = {}
    if backend == "quotient-cg":
        cpu_operator = QuotientStokesletOperator(
            sampled, rotations, epsilon=epsilon, viscosity=viscosity,
            row_chunk=matrix_free_row_chunk,
        )
        operator = cpu_operator
        if matrix_free_device == "cuda":
            from .cuda_operator import make_cuda_operator
            operator = make_cuda_operator(cpu_operator, library=cuda_library, device_index=cuda_device_index)
        try:
            rhs = _quotient_rhs(sampled)
            solution, diagnostics = solve_quotient_cg(
                operator, rhs, rtol=solver_rtol,
                maxiter=solver_maxiter, progress=solver_progress,
                panel_block_size=preconditioner_panel_block_size,
            )
            diagnostics["matrix_free_device"] = matrix_free_device
            if matrix_free_device == "cuda":
                # Independent reference implementation, not another GPU residual.
                cpu_residual = rhs - cpu_operator.matmat(solution)
                cpu_relative = np.linalg.norm(cpu_residual, axis=0) / np.linalg.norm(rhs, axis=0)
                diagnostics["gpu_rhs_relative_residuals"] = diagnostics["rhs_relative_residuals"]
                diagnostics["cpu_rhs_relative_residuals"] = cpu_relative.tolist()
                if np.any(~np.isfinite(cpu_relative)) or np.any(cpu_relative > solver_rtol):
                    raise RuntimeError(f"CUDA solution failed independent CPU true residual: {cpu_relative.tolist()}")
                diagnostics["rhs_relative_residuals"] = cpu_relative.tolist()
                diagnostics["linear_relative_residual"] = float(np.linalg.norm(cpu_residual) / np.linalg.norm(rhs))
                diagnostics["cpu_residual_verified"] = True
                diagnostics["cpu_verification_operator_applications"] = cpu_operator.applications
                diagnostics["cuda"] = operator.cuda.metadata
        finally:
            if matrix_free_device == "cuda":
                operator.cuda.close()
        tractions = _expand_quotient_tractions(sampled, rotations, solution, center)
        solver = "pcg"
        residual = diagnostics["linear_relative_residual"]
        system_dimension = 9 * len(sampled)
    elif backend == "quotient":
        tractions, solver, residual = solve_quotient_tractions(
            sampled, rotations, epsilon=epsilon, viscosity=viscosity, center=center,
        )
        system_dimension = 9 * len(sampled)
    else:
        matrix = assemble_regularized_stokeslet_matrix(panels, epsilon=epsilon, viscosity=viscosity)
        rhs = np.column_stack([rigid_rhs(panels, load, center=center) for load in range(6)])
        try:
            tractions = np.linalg.solve(matrix, rhs)
            solver = "solve"
        except np.linalg.LinAlgError:
            tractions = np.linalg.lstsq(matrix, rhs, rcond=1e-10)[0]
            solver = "lstsq"
        residual = float(np.linalg.norm(matrix @ tractions - rhs) / np.linalg.norm(rhs))
        system_dimension = len(matrix)

    K = np.zeros((3, 3), dtype=float)
    C = np.zeros((3, 3), dtype=float)
    Q = np.zeros((3, 3), dtype=float)
    force_due_rotation = np.zeros((3, 3), dtype=float)
    for load in range(6):
        force, torque = integrated_force_torque(panels, tractions[:, load], center=center)
        if load < 3:
            K[:, load] = force
            C[:, load] = torque
        else:
            force_due_rotation[:, load - 3] = force
            Q[:, load - 3] = torque

    k = float(np.trace(K) / 3.0)
    c = float(np.trace(C) / 3.0)
    q = float(np.trace(Q) / 3.0)
    rho = abs(c) / math.sqrt(k * q) if k > 0.0 and q > 0.0 else math.nan
    return {
        "mesh": mesh_label,
        "symmetry": symmetry.NAME,
        "symmetry_order": symmetry.ORDER,
        "solver": solver,
        "backend": backend,
        "system_dimension": system_dimension,
        "matrix_bytes": 0 if backend == "quotient-cg" else system_dimension**2 * np.dtype(float).itemsize,
        "equivalent_dense_matrix_bytes": system_dimension**2 * np.dtype(float).itemsize,
        "linear_relative_residual": residual,
        "elapsed_seconds": time.perf_counter() - started,
        "viscosity": viscosity,
        "epsilon": epsilon,
        "epsilon_factor": epsilon_factor,
        "chamber_panels_raw": len(chamber_panels),
        "chamber_panels_used": len(sampled),
        "full_panels": len(panels),
        "surface_area": sum(panel.area for panel in panels),
        "K": K.tolist(),
        "C": C.tolist(),
        "Q": Q.tolist(),
        "force_due_rotation": force_due_rotation.tolist(),
        "k": k,
        "c": c,
        "q": q,
        "rho": rho,
        "eK": anisotropy_error(K),
        "eC": anisotropy_error(C),
        "eQ": anisotropy_error(Q),
        "reciprocity_error": float(np.linalg.norm(force_due_rotation - C.T, ord="fro") / max(1.0, np.linalg.norm(C, ord="fro"))),
        **diagnostics,
    }


def write_summary(path: Path, result: dict[str, object]) -> None:
    lines: list[str] = []
    for key in ("k", "c", "q", "rho", "epsilon", "epsilon_factor", "surface_area"):
        if key in result:
            lines.append(f"{key}={result[key]}")
    for name in ("K", "C", "Q"):
        matrix = result.get(name)
        if not isinstance(matrix, list):
            continue
        for i, row in enumerate(matrix):
            if not isinstance(row, list):
                continue
            for j, value in enumerate(row):
                lines.append(f"{name}{i}{j}={value}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--max-chamber-panels", type=int, default=32)
    parser.add_argument("--backend", choices=("full", "quotient", "quotient-cg"), default="full")
    parser.add_argument("--subdivision-levels", type=int, default=0,
                        help="Split every source Gamma triangle into 4^levels coplanar triangles before sampling.")
    parser.add_argument("--solver-rtol", type=float, default=1e-10, help="Each quotient-CG RHS must pass an explicit relative residual check.")
    parser.add_argument("--solver-maxiter", type=int, default=500)
    parser.add_argument("--matrix-free-row-chunk", type=int, default=32)
    parser.add_argument("--solver-progress", action="store_true", help="Report true quotient-CG residuals to stderr.")
    parser.add_argument("--preconditioner-panel-block-size", type=int, default=16)
    parser.add_argument("--matrix-free-device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--cuda-library", help="Explicit fixed-source FP64 shared library; no fallback.")
    parser.add_argument("--cuda-device-index", type=int, default=0)
    parser.add_argument("--epsilon-factor", type=float, default=0.6)
    parser.add_argument("--viscosity", type=float, default=1.0)
    parser.add_argument("--center", type=float, nargs=3, default=(0.0, 0.0, 0.0))
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON only.")
    parser.add_argument("--summary-path", type=Path, help="Write key=value metrics for C++ optimizer integration.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = score_kcq(
        args.mesh,
        max_chamber_panels=args.max_chamber_panels,
        epsilon_factor=args.epsilon_factor,
        viscosity=args.viscosity,
        center=tuple(args.center),
        backend=args.backend,
        subdivision_levels=args.subdivision_levels,
        solver_rtol=args.solver_rtol,
        solver_maxiter=args.solver_maxiter,
        matrix_free_row_chunk=args.matrix_free_row_chunk,
        solver_progress=args.solver_progress,
        preconditioner_panel_block_size=args.preconditioner_panel_block_size,
        matrix_free_device=args.matrix_free_device,
        cuda_library=args.cuda_library,
        cuda_device_index=args.cuda_device_index,
    )
    if args.summary_path is not None:
        write_summary(args.summary_path, result)
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return
    print(f"BEM KCQ mesh={result['mesh']}")
    print(
        "panels "
        f"chamber={result['chamber_panels_used']}/{result['chamber_panels_raw']} "
        f"full={result['full_panels']} epsilon={result['epsilon']:.8g} solver={result['solver']}"
    )
    print(f"k={result['k']:.12g} c={result['c']:.12g} q={result['q']:.12g} rho={result['rho']:.12g}")
    print(f"eK={result['eK']:.4g} eC={result['eC']:.4g} eQ={result['eQ']:.4g} reciprocity={result['reciprocity_error']:.4g}")
    print("K=" + json.dumps(result["K"]))
    print("C=" + json.dumps(result["C"]))
    print("Q=" + json.dumps(result["Q"]))


def proper_tetrahedral_rotations() -> list[Matrix]:
    rotations: list[Matrix] = []
    for perm in permutations(range(3)):
        inversions = sum(
            perm[i] > perm[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        if inversions % 2:
            continue
        for signs in product((-1, 1), repeat=3):
            if signs[0] * signs[1] * signs[2] != 1:
                continue
            matrix = [[0, 0, 0] for _ in range(3)]
            for row, col in enumerate(perm):
                matrix[row][col] = signs[row]
            rotation = tuple(tuple(row) for row in matrix)  # type: ignore[assignment]
            if det3(rotation) == 1:
                rotations.append(rotation)
    unique = sorted(set(rotations))
    if len(unique) != 12:
        raise RuntimeError(f"expected 12 proper tetrahedral rotations, got {len(unique)}")
    return unique


if __name__ == "__main__":
    main()
