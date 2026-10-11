#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Report periodic cut-face and selected-symmetry Gamma diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

from .carlos_reference import ROOT, EXAMPLE

from . import kelvin_symmetry as symmetry
from . import mesh_io as periodic_gmsh
from .mesh_io import GAMMA, SIGMA_MINUS, SIGMA_PLUS, SIGMA_XY, SIGMA_TETRA_DEPENDENT
from . import periodic_interface_geometry as geometry


ROT_XY_Z_NEG: geometry.Matrix = (
    (0.0, 1.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 0.0, -1.0),
)
ROT_X_PLUS_90: geometry.Matrix = (
    (1.0, 0.0, 0.0),
    (0.0, 0.0, -1.0),
    (0.0, 1.0, 0.0),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--phi", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--tolerance", type=float)
    return parser.parse_args()


def main() -> None:
    if symmetry.IS_TERRA:
        return _terra_main()
    args = parse_args()
    vertices, triangles, tetrahedra = periodic_gmsh.read_medit_mesh(args.mesh)
    surfaces = {
        attr: [
            (a - 1, b - 1, c - 1)
            for a, b, c, triangle_attr in triangles
            if triangle_attr == attr
        ]
        for attr in (GAMMA, SIGMA_PLUS, SIGMA_MINUS, SIGMA_XY)
    }
    extent = max(abs(coordinate) for point in vertices for coordinate in point)
    tolerance = args.tolerance or max(1e-12, extent * 1e-8)
    if tolerance <= 0.0:
        raise ValueError("diagnostic tolerance must be positive")

    report: dict[str, object] = {
        "mesh": str(args.mesh.resolve()),
        "tolerance": tolerance,
        "full24_watertightness": asdict(
            geometry.full24_watertightness(
                chamber_vertices=vertices,
                chamber_triangles=surfaces[GAMMA],
                rotations=geometry.proper_octahedral_rotations(),
                tolerance=tolerance,
            )
        ),
        "paired_gamma_seams": {},
        "paired_partitions": {},
        "paired_traces": {},
    }

    partition_by_face = _boundary_partition_by_face(tetrahedra)
    pair_specs = (
        (
            "sigma_minus_to_plus",
            surfaces[SIGMA_PLUS],
            surfaces[SIGMA_MINUS],
            ROT_X_PLUS_90,
        ),
        (
            "sigma_xy_self",
            surfaces[SIGMA_XY],
            surfaces[SIGMA_XY],
            ROT_XY_Z_NEG,
        ),
    )
    gamma_edges = _surface_edges(surfaces[GAMMA])
    seam_segments = {
        attr: _edge_segments(
            vertices,
            gamma_edges.intersection(_surface_edges(surfaces[attr])),
        )
        for attr in (SIGMA_PLUS, SIGMA_MINUS, SIGMA_XY)
    }
    report["paired_gamma_seams"] = {
        "sigma_minus_to_plus": _seam_report(
            seam_segments[SIGMA_PLUS],
            seam_segments[SIGMA_MINUS],
            ROT_X_PLUS_90,
        ),
        "sigma_xy_self": _seam_report(
            seam_segments[SIGMA_XY],
            seam_segments[SIGMA_XY],
            ROT_XY_Z_NEG,
        ),
    }
    partition_reports: dict[str, object] = {}
    for name, master_triangles, dependent_triangles, dependent_to_master in pair_specs:
        if not master_triangles or not dependent_triangles:
            continue
        master_regions = [
            (triangle, _face_region(triangle, partition_by_face))
            for triangle in master_triangles
        ]
        dependent_regions = [
            (triangle, _face_region(triangle, partition_by_face))
            for triangle in dependent_triangles
        ]
        partition_reports[name] = asdict(
            geometry.compare_paired_partition(
                master_vertices=vertices,
                master_regions=master_regions,
                dependent_vertices=vertices,
                dependent_regions=dependent_regions,
                dependent_to_master=dependent_to_master,
            )
        )
    report["paired_partitions"] = partition_reports

    if args.phi is not None:
        values = _read_values(args.phi)
        if len(values) != len(vertices):
            raise ValueError(
                "phi value count does not match mesh vertices: "
                f"values={len(values)} vertices={len(vertices)}"
            )
        trace_reports: dict[str, object] = {}
        for name, master_triangles, dependent_triangles, dependent_to_master in pair_specs:
            if not master_triangles or not dependent_triangles:
                continue
            trace_reports[name] = {
                "comparison": asdict(
                    geometry.compare_paired_trace(
                        master_vertices=vertices,
                        master_triangles=master_triangles,
                        master_values=values,
                        dependent_vertices=vertices,
                        dependent_triangles=dependent_triangles,
                        dependent_values=values,
                        dependent_to_master=dependent_to_master,
                    )
                ),
                "zero_contour_hausdorff": geometry.mapped_zero_contour_hausdorff(
                    master_vertices=vertices,
                    master_triangles=master_triangles,
                    master_values=values,
                    dependent_vertices=vertices,
                    dependent_triangles=dependent_triangles,
                    dependent_values=values,
                    dependent_to_master=dependent_to_master,
                ),
            }
        report["paired_traces"] = trace_reports

    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    sys.stdout.write(output)
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(output, encoding="ascii")


def _terra_main() -> None:
    args = parse_args()
    vertices, triangles, tetrahedra = periodic_gmsh.read_medit_mesh(args.mesh)
    surfaces = {
        attr: [
            (a - 1, b - 1, c - 1)
            for a, b, c, triangle_attr in triangles
            if triangle_attr == attr
        ]
        for attr in (
            GAMMA, SIGMA_PLUS, SIGMA_MINUS, SIGMA_XY,
            SIGMA_TETRA_DEPENDENT,
        )
    }
    extent = max(abs(coordinate) for point in vertices for coordinate in point)
    tolerance = args.tolerance or max(1e-12, extent * 1e-8)
    if tolerance <= 0.0:
        raise ValueError("diagnostic tolerance must be positive")
    report: dict[str, object] = {
        "mesh": str(args.mesh.resolve()),
        "tolerance": tolerance,
        "full12_watertightness": asdict(
            geometry.full12_watertightness(
                chamber_vertices=vertices,
                chamber_triangles=surfaces[GAMMA],
                rotations=geometry.proper_tetrahedral_rotations(),
                tolerance=tolerance,
            )
        ),
        "paired_gamma_seams": {},
        "paired_partitions": {},
        "paired_traces": {},
    }

    partition_by_face = _boundary_partition_by_face(tetrahedra)
    pair_specs = (
        (
            "tetra_plus",
            surfaces[SIGMA_PLUS],
            surfaces[SIGMA_MINUS],
            geometry.transpose_matrix(geometry.ROT_TETRA_PLUS),
        ),
        (
            "tetra_minus",
            surfaces[SIGMA_XY],
            surfaces[SIGMA_TETRA_DEPENDENT],
            geometry.transpose_matrix(geometry.ROT_TETRA_MINUS),
        ),
    )
    gamma_edges = _surface_edges(surfaces[GAMMA])
    seam_segments = {
        "plus_master": _edge_segments(
            vertices, gamma_edges.intersection(_surface_edges(surfaces[SIGMA_PLUS]))),
        "plus_dependent": _edge_segments(
            vertices, gamma_edges.intersection(_surface_edges(surfaces[SIGMA_MINUS]))),
        "minus_master": _edge_segments(
            vertices, gamma_edges.intersection(_surface_edges(surfaces[SIGMA_XY]))),
        "minus_dependent": _edge_segments(
            vertices,
            gamma_edges.intersection(_surface_edges(surfaces[SIGMA_TETRA_DEPENDENT]))),
    }
    report["paired_gamma_seams"] = {
        "tetra_plus": _seam_report(
            seam_segments["plus_master"],
            seam_segments["plus_dependent"],
            geometry.transpose_matrix(geometry.ROT_TETRA_PLUS),
        ),
        "tetra_minus": _seam_report(
            seam_segments["minus_master"],
            seam_segments["minus_dependent"],
            geometry.transpose_matrix(geometry.ROT_TETRA_MINUS),
        ),
    }
    partition_reports: dict[str, object] = {}
    for name, master_triangles, dependent_triangles, dependent_to_master in pair_specs:
        if not master_triangles or not dependent_triangles:
            continue
        master_regions = [
            (triangle, _face_region(triangle, partition_by_face))
            for triangle in master_triangles
        ]
        dependent_regions = [
            (triangle, _face_region(triangle, partition_by_face))
            for triangle in dependent_triangles
        ]
        partition_reports[name] = asdict(
            geometry.compare_paired_partition(
                master_vertices=vertices,
                master_regions=master_regions,
                dependent_vertices=vertices,
                dependent_regions=dependent_regions,
                dependent_to_master=dependent_to_master,
            )
        )
    report["paired_partitions"] = partition_reports

    if args.phi is not None:
        values = _read_values(args.phi)
        if len(values) != len(vertices):
            raise ValueError(
                "phi value count does not match mesh vertices: "
                f"values={len(values)} vertices={len(vertices)}"
            )
        trace_reports: dict[str, object] = {}
        for name, master_triangles, dependent_triangles, dependent_to_master in pair_specs:
            if not master_triangles or not dependent_triangles:
                continue
            trace_reports[name] = {
                "comparison": asdict(
                    geometry.compare_paired_trace(
                        master_vertices=vertices,
                        master_triangles=master_triangles,
                        master_values=values,
                        dependent_vertices=vertices,
                        dependent_triangles=dependent_triangles,
                        dependent_values=values,
                        dependent_to_master=dependent_to_master,
                    )
                ),
                "zero_contour_hausdorff": geometry.mapped_zero_contour_hausdorff(
                    master_vertices=vertices,
                    master_triangles=master_triangles,
                    master_values=values,
                    dependent_vertices=vertices,
                    dependent_triangles=dependent_triangles,
                    dependent_values=values,
                    dependent_to_master=dependent_to_master,
                ),
            }
        report["paired_traces"] = trace_reports

    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    sys.stdout.write(output)
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(output, encoding="ascii")


def _read_values(path: Path) -> list[float]:
    return [
        float(line.strip())
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip()
    ]


def _boundary_partition_by_face(
    tetrahedra: list[periodic_gmsh.Tet],
) -> dict[tuple[int, int, int], list[int]]:
    partition: dict[tuple[int, int, int], list[int]] = {}
    for a, b, c, d, attr in tetrahedra:
        vertices = (a - 1, b - 1, c - 1, d - 1)
        for face in (
            (vertices[0], vertices[1], vertices[2]),
            (vertices[0], vertices[1], vertices[3]),
            (vertices[0], vertices[2], vertices[3]),
            (vertices[1], vertices[2], vertices[3]),
        ):
            partition.setdefault(tuple(sorted(face)), []).append(attr)
    return partition


def _surface_edges(
    triangles: list[geometry.Triangle],
) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for triangle in triangles:
        for first, second in ((0, 1), (1, 2), (2, 0)):
            edges.add(tuple(sorted((triangle[first], triangle[second]))))
    return edges


def _edge_segments(
    vertices: list[geometry.Point],
    edges: set[tuple[int, int]],
) -> list[geometry.Segment]:
    return [(vertices[first], vertices[second]) for first, second in sorted(edges)]


def _seam_report(
    master_segments: list[geometry.Segment],
    dependent_segments: list[geometry.Segment],
    dependent_to_master: geometry.Matrix,
) -> dict[str, float | int]:
    return {
        "master_segments": len(master_segments),
        "dependent_segments": len(dependent_segments),
        "hausdorff": geometry.mapped_segment_hausdorff(
            master_segments=master_segments,
            dependent_segments=dependent_segments,
            dependent_to_master=dependent_to_master,
        ),
    }


def _face_region(
    triangle: geometry.Triangle,
    partition_by_face: dict[tuple[int, int, int], list[int]],
) -> int:
    adjacent = partition_by_face.get(tuple(sorted(triangle)), [])
    if len(adjacent) != 1:
        raise ValueError(
            "chamber cut triangle must have exactly one adjacent tetrahedron: "
            f"triangle={triangle} adjacent={adjacent}"
        )
    return adjacent[0]


if __name__ == "__main__":
    main()
