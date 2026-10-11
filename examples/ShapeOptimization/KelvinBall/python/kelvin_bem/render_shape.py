#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Render full 24-copy views from a 1/24 chamber Gamma surface."""

from __future__ import annotations

import argparse
from pathlib import Path


from . import mesh_io as periodic_gmsh
from . import periodic_geometry as geometry
from . import render_support


def write_full24_mesh(
    source: Path,
    target: Path,
    *,
    tolerance: float | None = None,
) -> geometry.Watertightness:
    vertices, triangles, _tetrahedra = periodic_gmsh.read_medit_mesh(source)
    gamma = [
        (a - 1, b - 1, c - 1)
        for a, b, c, attr in triangles
        if attr == periodic_gmsh.GAMMA
    ]
    extent = max(abs(coordinate) for point in vertices for coordinate in point)
    weld_tolerance = (
        tolerance
        if tolerance is not None
        else max(1e-12, extent * 1e-8)
    )
    welded = geometry.weld_full24_surface(
        chamber_vertices=vertices,
        chamber_triangles=gamma,
        rotations=geometry.proper_octahedral_rotations(),
        tolerance=weld_tolerance,
    )
    full_triangles = [
        (a + 1, b + 1, c + 1, periodic_gmsh.GAMMA)
        for a, b, c in welded.triangles
    ]

    target.parent.mkdir(parents=True, exist_ok=True)
    periodic_gmsh.write_medit_mesh(target, welded.vertices, full_triangles, [])
    return welded.watertightness


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tolerance", type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    full_mesh = args.out_dir / "full24_gamma.mesh"
    stats = write_full24_mesh(
        args.mesh,
        full_mesh,
        tolerance=args.tolerance,
    )
    print(
        "full24 welding"
        f" merged_vertices={stats.merged_vertices}"
        f" triangles={stats.triangles}"
        f" boundary_edges={stats.boundary_edges}"
        f" nonmanifold_edges={stats.nonmanifold_edges}"
        f" degenerate_triangles={stats.degenerate_triangles}"
    )
    for path in render_support.render_views(full_mesh, args.out_dir):
        print(path)


if __name__ == "__main__":
    main()
