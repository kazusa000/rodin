#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Tests for rotation-aware diagnostics on non-matching chamber faces."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

from kelvin_bem import periodic_interface_geometry as geometry


ROT_X_MINUS_90: geometry.Matrix = (
    (1.0, 0.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, -1.0, 0.0),
)
ROT_X_PLUS_90: geometry.Matrix = (
    (1.0, 0.0, 0.0),
    (0.0, 0.0, -1.0),
    (0.0, 1.0, 0.0),
)


def _apply(matrix: geometry.Matrix, point: geometry.Point) -> geometry.Point:
    return tuple(
        sum(matrix[row][col] * point[col] for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _nonmatching_face_data() -> tuple[
    list[geometry.Point],
    list[geometry.Triangle],
    list[geometry.Point],
    list[geometry.Triangle],
]:
    master_vertices: list[geometry.Point] = [
        (1.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (2.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
    ]
    master_triangles: list[geometry.Triangle] = [(0, 1, 2), (0, 2, 3)]
    dependent_vertices = [
        _apply(ROT_X_MINUS_90, point)
        for point in master_vertices
    ]
    dependent_vertices.append(_apply(ROT_X_MINUS_90, (1.5, 0.5, 0.5)))
    dependent_triangles: list[geometry.Triangle] = [
        (0, 1, 4),
        (1, 2, 4),
        (2, 3, 4),
        (3, 0, 4),
    ]
    return master_vertices, master_triangles, dependent_vertices, dependent_triangles


def _master_affine(point: geometry.Point) -> float:
    return point[0] + 2.0 * point[1] - 0.3


def test_nonmatching_affine_traces_agree_after_rotation() -> None:
    master_vertices, master_triangles, dependent_vertices, dependent_triangles = (
        _nonmatching_face_data()
    )
    master_values = [_master_affine(point) for point in master_vertices]
    dependent_values = [
        _master_affine(_apply(ROT_X_PLUS_90, point))
        for point in dependent_vertices
    ]

    stats = geometry.compare_paired_trace(
        master_vertices=master_vertices,
        master_triangles=master_triangles,
        master_values=master_values,
        dependent_vertices=dependent_vertices,
        dependent_triangles=dependent_triangles,
        dependent_values=dependent_values,
        dependent_to_master=ROT_X_PLUS_90,
    )

    if stats.sample_count <= 0:
        raise AssertionError("trace comparison should visit quadrature samples")
    if stats.uncovered_samples != 0:
        raise AssertionError(f"mapped face should be fully covered: {stats}")
    if stats.max_abs_difference >= 1e-12:
        raise AssertionError(f"mapped affine trace should agree pointwise: {stats}")
    if stats.rms_difference >= 1e-12:
        raise AssertionError(f"mapped affine trace should agree in RMS: {stats}")
    if stats.sign_mismatch_area >= 1e-12:
        raise AssertionError(f"mapped affine trace signs should agree: {stats}")


def test_nonmatching_trace_detects_a_perturbed_dependent_value() -> None:
    master_vertices, master_triangles, dependent_vertices, dependent_triangles = (
        _nonmatching_face_data()
    )
    master_values = [_master_affine(point) for point in master_vertices]
    dependent_values = [
        _master_affine(_apply(ROT_X_PLUS_90, point))
        for point in dependent_vertices
    ]
    dependent_values[-1] += 0.25

    stats = geometry.compare_paired_trace(
        master_vertices=master_vertices,
        master_triangles=master_triangles,
        master_values=master_values,
        dependent_vertices=dependent_vertices,
        dependent_triangles=dependent_triangles,
        dependent_values=dependent_values,
        dependent_to_master=ROT_X_PLUS_90,
    )

    if stats.max_abs_difference <= 1e-3:
        raise AssertionError(f"trace perturbation should be detected: {stats}")
    if stats.rms_difference <= 1e-4:
        raise AssertionError(f"trace perturbation should affect RMS: {stats}")


def test_skinny_exact_partition_is_not_reported_as_uncovered() -> None:
    master_vertices: list[geometry.Point] = [
        (0.6977364772734991, 0.6977364772734991, 0.16015203501588082),
        (0.6977663658274733, 0.6977663658274733, 0.15964442703752624),
        (0.701708702504665, 0.701708702504665, 0.08703882797796147),
    ]
    dependent_vertices = [
        _apply(geometry.ROT_XY_Z_NEG, point)
        for point in master_vertices
    ]

    stats = geometry.compare_paired_partition(
        master_vertices=master_vertices,
        master_regions=[((0, 1, 2), 1)],
        dependent_vertices=dependent_vertices,
        dependent_regions=[((2, 1, 0), 1)],
        dependent_to_master=geometry.ROT_XY_Z_NEG,
    )

    if stats.uncovered_samples != 0:
        raise AssertionError(f"exact skinny partition must remain covered: {stats}")
    if stats.sample_count != 6:
        raise AssertionError(f"all bidirectional samples must be located: {stats}")


def test_mapped_zero_contours_agree_on_nonmatching_faces() -> None:
    master_vertices, master_triangles, dependent_vertices, dependent_triangles = (
        _nonmatching_face_data()
    )
    master_values = [point[0] - 1.4 for point in master_vertices]
    dependent_values = [
        _apply(ROT_X_PLUS_90, point)[0] - 1.4
        for point in dependent_vertices
    ]

    distance = geometry.mapped_zero_contour_hausdorff(
        master_vertices=master_vertices,
        master_triangles=master_triangles,
        master_values=master_values,
        dependent_vertices=dependent_vertices,
        dependent_triangles=dependent_triangles,
        dependent_values=dependent_values,
        dependent_to_master=ROT_X_PLUS_90,
    )

    if distance >= 1e-12:
        raise AssertionError(f"mapped affine zero contours should agree: {distance}")


def test_mapped_zero_contour_detects_a_shift() -> None:
    master_vertices, master_triangles, dependent_vertices, dependent_triangles = (
        _nonmatching_face_data()
    )
    master_values = [point[0] - 1.4 for point in master_vertices]
    dependent_values = [
        _apply(ROT_X_PLUS_90, point)[0] - 1.4 + 0.1
        for point in dependent_vertices
    ]

    distance = geometry.mapped_zero_contour_hausdorff(
        master_vertices=master_vertices,
        master_triangles=master_triangles,
        master_values=master_values,
        dependent_vertices=dependent_vertices,
        dependent_triangles=dependent_triangles,
        dependent_values=dependent_values,
        dependent_to_master=ROT_X_PLUS_90,
    )

    if distance <= 0.09:
        raise AssertionError(f"shifted zero contour should be detected: {distance}")


def test_full24_cube_patch_is_watertight() -> None:
    vertices: list[geometry.Point] = [
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
    ]
    triangles: list[geometry.Triangle] = [(0, 1, 2)]
    rotations = geometry.proper_octahedral_rotations()

    stats = geometry.full24_watertightness(
        chamber_vertices=vertices,
        chamber_triangles=triangles,
        rotations=rotations,
        tolerance=1e-10,
    )

    if len(rotations) != 24:
        raise AssertionError(f"expected 24 proper rotations, got {len(rotations)}")
    if stats.boundary_edges != 0 or stats.nonmanifold_edges != 0:
        raise AssertionError(f"full cube surface should be watertight: {stats}")
    if stats.degenerate_triangles != 0:
        raise AssertionError(f"full cube surface should be non-degenerate: {stats}")


def test_full24_cube_patch_returns_welded_connectivity() -> None:
    vertices: list[geometry.Point] = [
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
    ]
    welded = geometry.weld_full24_surface(
        chamber_vertices=vertices,
        chamber_triangles=[(0, 1, 2)],
        rotations=geometry.proper_octahedral_rotations(),
        tolerance=1e-10,
    )

    if len(welded.vertices) != welded.watertightness.merged_vertices:
        raise AssertionError(f"welded vertex count is inconsistent: {welded}")
    if len(welded.triangles) != 24:
        raise AssertionError(f"every rotation should contribute one triangle: {welded}")
    if max(index for triangle in welded.triangles for index in triangle) >= len(welded.vertices):
        raise AssertionError(f"welded triangle index is out of range: {welded}")


def test_aabb_point_welding_crosses_round_key_boundaries() -> None:
    tolerance = 1e-8
    points: list[geometry.Point] = [
        (0.49 * tolerance, 0.0, 0.0),
        (0.51 * tolerance, 0.0, 0.0),
    ]

    representatives, labels = geometry._weld_points_by_distance(points, tolerance)

    if len(representatives) != 1 or labels != [0, 0]:
        raise AssertionError(
            "sub-tolerance points across a round-key boundary must be welded: "
            f"representatives={representatives} labels={labels}"
        )


def test_aabb_point_welding_rejects_transitive_overmerge() -> None:
    tolerance = 1.0
    points: list[geometry.Point] = [
        (0.0, 0.0, 0.0),
        (0.75, 0.0, 0.0),
        (1.5, 0.0, 0.0),
    ]

    try:
        geometry._weld_points_by_distance(points, tolerance)
    except ValueError as error:
        if "transitive point-welding cluster" not in str(error):
            raise AssertionError(f"unexpected transitive-cluster error: {error}") from error
    else:
        raise AssertionError("transitive welding beyond the tolerance must be rejected")


def test_full24_patch_reports_a_missing_rotation() -> None:
    vertices: list[geometry.Point] = [
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
    ]
    triangles: list[geometry.Triangle] = [(0, 1, 2)]
    rotations = geometry.proper_octahedral_rotations()[:-1]

    stats = geometry.full24_watertightness(
        chamber_vertices=vertices,
        chamber_triangles=triangles,
        rotations=rotations,
        tolerance=1e-10,
    )

    if stats.boundary_edges <= 0:
        raise AssertionError(f"missing rotated patch should create open edges: {stats}")


def test_geometry_diagnostic_cli_emits_json() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        mesh = Path(tmp) / "chamber.mesh"
        mesh.write_text(
            """MeshVersionFormatted 2
Dimension 3

Vertices
4
0 0 0 0
1 0 0 0
1 1 1 0
1 1 -1 0

Triangles
4
2 3 4 13
1 2 3 31
1 4 2 32
1 3 4 33

Tetrahedra
1
1 2 3 4 2

End
""",
            encoding="ascii",
        )
        command = [
            sys.executable,
            str(ROOT / "scripts" / "kelvin-bem.py"), "geometry-check",
            str(mesh),
        ]
        result = subprocess.run(command, check=False, text=True, capture_output=True)
        if result.returncode != 0:
            raise AssertionError(f"geometry CLI failed: {result.stdout}\n{result.stderr}")
        report = json.loads(result.stdout)
        required = {
            "mesh",
            "full24_watertightness",
            "paired_gamma_seams",
            "paired_partitions",
            "paired_traces",
        }
        if not required.issubset(report):
            raise AssertionError(f"geometry CLI report is incomplete: {report}")
        if report["full24_watertightness"]["boundary_edges"] != 0:
            raise AssertionError(f"cube chamber fixture should close: {report}")
        for name, seam in report["paired_gamma_seams"].items():
            if seam["master_segments"] <= 0 or seam["dependent_segments"] <= 0:
                raise AssertionError(f"fixture seam should contain segments for {name}: {report}")
            if seam["hausdorff"] >= 1e-12:
                raise AssertionError(f"fixture seam should map exactly for {name}: {report}")


if __name__ == "__main__":
    test_nonmatching_affine_traces_agree_after_rotation()
    test_nonmatching_trace_detects_a_perturbed_dependent_value()
    test_skinny_exact_partition_is_not_reported_as_uncovered()
    test_mapped_zero_contours_agree_on_nonmatching_faces()
    test_mapped_zero_contour_detects_a_shift()
    test_full24_cube_patch_is_watertight()
    test_full24_cube_patch_returns_welded_connectivity()
    test_aabb_point_welding_crosses_round_key_boundaries()
    test_aabb_point_welding_rejects_transitive_overmerge()
    test_full24_patch_reports_a_missing_rotation()
    test_geometry_diagnostic_cli_emits_json()
    print("ok")
