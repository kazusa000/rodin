#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Checks for basic GammaShape mesh view rendering."""

from __future__ import annotations

import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


MESH = """MeshVersionFormatted 2
Dimension
3
Vertices
5
0 0 0 0
1 0 0 0
0 1 0 0
0 0 1 0
10 10 10 0
Triangles
2
1 2 3 13
1 2 5 7
Tetrahedra
1
1 2 3 4 2
End
"""


CHAMBER_PATCH = """MeshVersionFormatted 2
Dimension
3
Vertices
4
1 0 0 0
1 1 1 0
1 1 -1 0
0 0 0 0
Triangles
1
1 2 3 13
Tetrahedra
1
1 2 3 4 2
End
"""


def test_parse_medit_surface_keeps_gamma_triangles_only() -> None:
    from kelvin_bem import render_support as render

    with tempfile.TemporaryDirectory() as tmp:
        mesh = Path(tmp) / "tiny.mesh"
        mesh.write_text(MESH, encoding="ascii")
        vertices, triangles = render.parse_medit_surface(mesh)

    if len(vertices) != 5:
        raise AssertionError("all vertices should be parsed")
    if triangles != [(0, 1, 2)]:
        raise AssertionError("only GammaShape=13 triangles should be rendered")


def test_surface_axis_limits_ignore_unused_outer_vertices() -> None:
    from kelvin_bem import render_support as render

    with tempfile.TemporaryDirectory() as tmp:
        mesh = Path(tmp) / "tiny.mesh"
        mesh.write_text(MESH, encoding="ascii")
        vertices, triangles = render.parse_medit_surface(mesh)
        limits = render.surface_axis_limits(vertices, triangles)

    for lo, hi in limits:
        if hi - lo > 2.0:
            raise AssertionError("view limits should be based on GammaShape vertices only")


def test_render_views_writes_four_pngs() -> None:
    from kelvin_bem import render_support as render

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        mesh = root / "tiny.mesh"
        out = root / "views"
        mesh.write_text(MESH, encoding="ascii")
        written = render.render_views(mesh, out)

        if [path.name for path in written] != ["front.png", "side.png", "top.png", "iso.png"]:
            raise AssertionError("renderer should emit the four standard views")
        for path in written:
            if not path.exists() or path.stat().st_size <= 0:
                raise AssertionError(f"rendered PNG should be non-empty: {path}")


def test_full24_writer_welds_rotated_vertex_orbits() -> None:
    from kelvin_bem import render_support as render
    from kelvin_bem import render_shape as full24

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / "chamber.mesh"
        target = root / "full24.mesh"
        source.write_text(CHAMBER_PATCH, encoding="ascii")
        stats = full24.write_full24_mesh(source, target, tolerance=1e-10)
        vertices, triangles = render.parse_medit_surface(target)

    if len(vertices) != stats.merged_vertices:
        raise AssertionError(f"writer did not emit only welded vertices: {stats}")
    if len(triangles) != 24:
        raise AssertionError(f"writer should emit all 24 rotated triangles: {stats}")
    if stats.boundary_edges or stats.nonmanifold_edges or stats.degenerate_triangles:
        raise AssertionError(f"welded cube patch should be a valid closed surface: {stats}")


if __name__ == "__main__":
    test_parse_medit_surface_keeps_gamma_triangles_only()
    test_surface_axis_limits_ignore_unused_outer_vertices()
    test_render_views_writes_four_pngs()
    test_full24_writer_welds_rotated_vertex_orbits()
    print("ok")
