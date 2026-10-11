# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Unchanged geometric, merit and fixed-connectivity diagnostic contracts."""
import math
from pathlib import Path
import numpy as np
from . import clement_boundary as boundary
from . import carlos_reference as common


def geometry_diagnostics(mesh, symmetry="terra"):
    from dataclasses import asdict
    from . import export_carlos_interface as exporter
    from .tetra_solid_connectivity import solid_connectivity
    x, faces, labels, cells, cell_labels = boundary.read_mesh(mesh)
    gamma, volume = exporter.material_interface(x, cells, cell_labels)
    saved = {tuple(sorted(f)) for f in faces[labels == 13]}
    if saved != {tuple(sorted(f)) for f in gamma}:
        raise ValueError("Gamma labels differ from the solid/fluid interface")
    rotations = (exporter.geometry.proper_tetrahedral_rotations() if symmetry == "terra"
                 else exporter.geometry.proper_octahedral_rotations())
    stats = exporter.geometry.full24_watertightness(
        chamber_vertices=x.tolist(), chamber_triangles=gamma.tolist(),
        rotations=rotations, tolerance=1e-8)
    return dict(volume_full=len(rotations)*volume, chamber_panels=len(gamma),
                **{f"full{len(rotations)}": asdict(stats)}, mesh_sha256=common.digest(mesh),
                solid_connectivity=solid_connectivity(x, cells, cell_labels, rotations,
                    symmetry=symmetry))


def armijo(initial, candidate, dt, slope):
    if not all(math.isfinite(v) for v in (initial, candidate, dt, slope)) or min(dt, slope) <= 0:
        raise ValueError("Armijo requires finite measurements and a positive ascent slope")
    return candidate >= initial + max(1e-12, 1e-4*dt*slope)


def thickness_gate(initial, candidate, minimum, symmetry="terra"):
    """Validate thickness diagnostics; rho-P Armijo handles the soft cost.

    In particular P=0 is not a hard prohibition on entering the thin region.
    Neither sampled thickness nor this check guarantees pointwise thickness.
    """
    for result in (initial, candidate):
        if result.get("group_order") != (12 if symmetry == "terra" else 24) or not result.get("passed"):
            raise ValueError("thickness requires a valid symmetry-matched diagnostic")
        if not all(math.isfinite(result[k]) and result[k] >= 0 for k in ("penalty", "deepest")):
            raise ValueError("invalid thickness diagnostic")
        if not math.isclose(result["minimum"], minimum, rel_tol=0, abs_tol=1e-14):
            raise ValueError("thickness minimum changed")
    return True


def fragmentation_gate(initial, candidate, symmetry="terra"):
    """Reject additional full-group solid components, not open Gamma seams."""
    for geometry in (initial, candidate):
        diagnostic = geometry["solid_connectivity"]
        if diagnostic["group_order"] != (12 if symmetry == "terra" else 24) or diagnostic["full_components"] < 1:
            raise ValueError("invalid symmetry-matched solid connectivity diagnostic")
    return candidate["solid_connectivity"]["full_components"] <= initial["solid_connectivity"]["full_components"]


def step_factors(count):
    if not 1 <= count <= 12:
        raise ValueError("backtracking count must be between 1 and 12")
    return tuple(.5 / 2**i for i in range(count))


def fd_mesh(source, direction, amount, target):
    """Fixed-connectivity derivative diagnostic ONLY, never an optimization update."""
    lines = Path(source).read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == "Vertices") + 1
    while not lines[start].strip():
        start += 1
    count = int(lines[start])
    start += 1
    nodes = np.array([[float(v) for v in lines[start+i].split()[:3]] for i in range(count)])
    if direction.shape != (count, 6) or not np.allclose(nodes, direction[:, :3], rtol=0, atol=1e-13):
        raise ValueError("direction/mesh node ordering mismatch")
    for i, point in enumerate(nodes + amount*direction[:, 3:]):
        label = lines[start+i].split()[3]
        lines[start+i] = " ".join(format(v, ".17g") for v in point) + " " + label
    Path(target).write_text("\n".join(lines) + "\n")
