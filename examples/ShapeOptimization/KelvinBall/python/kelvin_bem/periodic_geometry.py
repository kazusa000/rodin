#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Geometry diagnostics for non-matching periodic chamber interfaces."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import permutations, product
import math


Point = tuple[float, float, float]
Triangle = tuple[int, int, int]
Segment = tuple[Point, Point]
Matrix = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]


@dataclass(frozen=True)
class TraceComparison:
    sample_count: int
    uncovered_samples: int
    max_abs_difference: float
    rms_difference: float
    sign_mismatch_area: float
    sampled_area: float


@dataclass(frozen=True)
class Watertightness:
    merged_vertices: int
    triangles: int
    boundary_edges: int
    nonmanifold_edges: int
    degenerate_triangles: int


@dataclass(frozen=True)
class WeldedSurface:
    vertices: list[Point]
    triangles: list[Triangle]
    watertightness: Watertightness


@dataclass
class _PointAABBNode:
    minimum: Point
    maximum: Point
    indices: tuple[int, ...] = ()
    left: "_PointAABBNode | None" = None
    right: "_PointAABBNode | None" = None


class _DisjointSet:
    def __init__(self, count: int) -> None:
        self.parent = list(range(count))
        self.rank = [0] * count

    def find(self, index: int) -> int:
        if self.parent[index] != index:
            self.parent[index] = self.find(self.parent[index])
        return self.parent[index]

    def join(self, first: int, second: int) -> None:
        first_root = self.find(first)
        second_root = self.find(second)
        if first_root == second_root:
            return
        if self.rank[first_root] < self.rank[second_root]:
            first_root, second_root = second_root, first_root
        self.parent[second_root] = first_root
        if self.rank[first_root] == self.rank[second_root]:
            self.rank[first_root] += 1


@dataclass(frozen=True)
class PartitionComparison:
    sample_count: int
    uncovered_samples: int
    mismatch_area: float
    sampled_area: float


ROT_XY_Z_NEG: Matrix = (
    (0.0, 1.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 0.0, -1.0),
)

ROT_X_MINUS_90: Matrix = (
    (1.0, 0.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, -1.0, 0.0),
)


ROT_TETRA_PLUS: Matrix = (
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
)

ROT_TETRA_MINUS: Matrix = (
    (0.0, 1.0, 0.0),
    (0.0, 0.0, -1.0),
    (-1.0, 0.0, 0.0),
)


_TRIANGLE_SAMPLES: tuple[tuple[float, float, float], ...] = (
    (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
    (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
    (1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0),
)


def apply_matrix(matrix: Matrix, point: Point) -> Point:
    return tuple(
        sum(matrix[row][col] * point[col] for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def transpose_matrix(matrix: Matrix) -> Matrix:
    return tuple(
        tuple(matrix[col][row] for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def proper_octahedral_rotations() -> list[Matrix]:
    rotations: set[Matrix] = set()
    for permutation in permutations(range(3)):
        for signs in product((-1.0, 1.0), repeat=3):
            rows = [[0.0, 0.0, 0.0] for _ in range(3)]
            for row, column in enumerate(permutation):
                rows[row][column] = signs[row]
            matrix: Matrix = tuple(tuple(row) for row in rows)  # type: ignore[assignment]
            if round(_determinant(matrix)) == 1:
                rotations.add(matrix)
    return sorted(rotations)


def full24_watertightness(
    *,
    chamber_vertices: list[Point],
    chamber_triangles: list[Triangle],
    rotations: list[Matrix],
    tolerance: float,
) -> Watertightness:
    return weld_full24_surface(
        chamber_vertices=chamber_vertices,
        chamber_triangles=chamber_triangles,
        rotations=rotations,
        tolerance=tolerance,
    ).watertightness


def weld_full24_surface(
    *,
    chamber_vertices: list[Point],
    chamber_triangles: list[Triangle],
    rotations: list[Matrix],
    tolerance: float,
) -> WeldedSurface:
    """Expand and weld a chamber surface using AABB candidates and exact distance."""
    if tolerance <= 0.0:
        raise ValueError("welding tolerance must be positive")
    if not chamber_triangles:
        raise ValueError("full-surface welding requires at least one Gamma triangle")
    used_vertices = sorted({index for triangle in chamber_triangles for index in triangle})
    if used_vertices[0] < 0 or used_vertices[-1] >= len(chamber_vertices):
        raise ValueError("Gamma triangle index is out of range")

    compact_index = {vertex: index for index, vertex in enumerate(used_vertices)}
    compact_triangles = [
        tuple(compact_index[vertex] for vertex in triangle)
        for triangle in chamber_triangles
    ]
    expanded_vertices: list[Point] = []
    expanded_triangles: list[Triangle] = []
    for rotation in rotations:
        offset = len(expanded_vertices)
        expanded_vertices.extend(
            apply_matrix(rotation, chamber_vertices[vertex])
            for vertex in used_vertices
        )
        expanded_triangles.extend(
            tuple(offset + vertex for vertex in triangle)
            for triangle in compact_triangles
        )

    welded_vertices, labels = _weld_points_by_distance(
        expanded_vertices,
        tolerance,
    )
    welded_triangles = [
        tuple(labels[vertex] for vertex in triangle)
        for triangle in expanded_triangles
    ]
    edge_incidence: Counter[tuple[int, int]] = Counter()
    degenerate_triangles = 0
    for triangle in welded_triangles:
        if len(set(triangle)) != 3:
            degenerate_triangles += 1
            continue
        for first, second in ((0, 1), (1, 2), (2, 0)):
            edge = tuple(sorted((triangle[first], triangle[second])))
            edge_incidence[edge] += 1

    watertightness = Watertightness(
        merged_vertices=len(welded_vertices),
        triangles=len(welded_triangles),
        boundary_edges=sum(count == 1 for count in edge_incidence.values()),
        nonmanifold_edges=sum(count > 2 for count in edge_incidence.values()),
        degenerate_triangles=degenerate_triangles,
    )
    return WeldedSurface(
        vertices=welded_vertices,
        triangles=welded_triangles,
        watertightness=watertightness,
    )


def _weld_points_by_distance(
    points: list[Point],
    tolerance: float,
) -> tuple[list[Point], list[int]]:
    if tolerance <= 0.0:
        raise ValueError("point-welding tolerance must be positive")
    if not points:
        return [], []

    tree = _build_point_aabb(points, tuple(range(len(points))))
    disjoint = _DisjointSet(len(points))
    tolerance_squared = math.nextafter(tolerance * tolerance, math.inf)
    for index, point in enumerate(points):
        query_minimum = tuple(coordinate - tolerance for coordinate in point)
        query_maximum = tuple(coordinate + tolerance for coordinate in point)
        candidates: list[int] = []
        _query_point_aabb(tree, query_minimum, query_maximum, candidates)
        for candidate in candidates:
            if candidate <= index:
                continue
            if _distance_squared(point, points[candidate]) <= tolerance_squared:
                disjoint.join(index, candidate)

    groups: dict[int, list[int]] = {}
    for index in range(len(points)):
        groups.setdefault(disjoint.find(index), []).append(index)
    ordered_groups = sorted(groups.values(), key=lambda members: members[0])
    labels = [0] * len(points)
    representatives: list[Point] = []
    for label, members in enumerate(ordered_groups):
        maximum_diameter_squared = 0.0
        for position, first in enumerate(members):
            for second in members[position + 1 :]:
                maximum_diameter_squared = max(
                    maximum_diameter_squared,
                    _distance_squared(points[first], points[second]),
                )
        if maximum_diameter_squared > tolerance_squared:
            raise ValueError(
                "transitive point-welding cluster exceeds tolerance: "
                f"diameter={math.sqrt(maximum_diameter_squared)} "
                f"tolerance={tolerance} members={members[:8]}"
            )
        representatives.append(points[members[0]])
        for member in members:
            labels[member] = label
    return representatives, labels


def _build_point_aabb(
    points: list[Point],
    indices: tuple[int, ...],
    leaf_size: int = 16,
) -> _PointAABBNode:
    minimum = tuple(min(points[index][axis] for index in indices) for axis in range(3))
    maximum = tuple(max(points[index][axis] for index in indices) for axis in range(3))
    if len(indices) <= leaf_size:
        return _PointAABBNode(minimum=minimum, maximum=maximum, indices=indices)
    axis = max(range(3), key=lambda item: maximum[item] - minimum[item])
    ordered = tuple(sorted(indices, key=lambda index: (points[index][axis], index)))
    middle = len(ordered) // 2
    return _PointAABBNode(
        minimum=minimum,
        maximum=maximum,
        left=_build_point_aabb(points, ordered[:middle], leaf_size),
        right=_build_point_aabb(points, ordered[middle:], leaf_size),
    )


def _query_point_aabb(
    node: _PointAABBNode,
    query_minimum: Point,
    query_maximum: Point,
    candidates: list[int],
) -> None:
    if any(
        node.maximum[axis] < query_minimum[axis]
        or query_maximum[axis] < node.minimum[axis]
        for axis in range(3)
    ):
        return
    if node.indices:
        candidates.extend(node.indices)
        return
    if node.left is not None:
        _query_point_aabb(node.left, query_minimum, query_maximum, candidates)
    if node.right is not None:
        _query_point_aabb(node.right, query_minimum, query_maximum, candidates)


def compare_paired_trace(
    *,
    master_vertices: list[Point],
    master_triangles: list[Triangle],
    master_values: list[float],
    dependent_vertices: list[Point],
    dependent_triangles: list[Triangle],
    dependent_values: list[float],
    dependent_to_master: Matrix,
) -> TraceComparison:
    _validate_trace_input(master_vertices, master_triangles, master_values, "master")
    _validate_trace_input(
        dependent_vertices,
        dependent_triangles,
        dependent_values,
        "dependent",
    )

    extent = max(
        1.0,
        _coordinate_extent(master_vertices),
        _coordinate_extent(dependent_vertices),
    )
    distance_tolerance = extent * 1e-9
    accum = _TraceAccumulator()
    _accumulate_trace_direction(
        source_vertices=dependent_vertices,
        source_triangles=dependent_triangles,
        source_values=dependent_values,
        target_vertices=master_vertices,
        target_triangles=master_triangles,
        target_values=master_values,
        source_to_target=dependent_to_master,
        distance_tolerance=distance_tolerance,
        accum=accum,
    )
    _accumulate_trace_direction(
        source_vertices=master_vertices,
        source_triangles=master_triangles,
        source_values=master_values,
        target_vertices=dependent_vertices,
        target_triangles=dependent_triangles,
        target_values=dependent_values,
        source_to_target=transpose_matrix(dependent_to_master),
        distance_tolerance=distance_tolerance,
        accum=accum,
    )
    rms = (
        math.sqrt(accum.weighted_squared_error / accum.sampled_area)
        if accum.sampled_area > 0.0
        else 0.0
    )
    return TraceComparison(
        sample_count=accum.sample_count,
        uncovered_samples=accum.uncovered_samples,
        max_abs_difference=accum.max_abs_difference,
        rms_difference=rms,
        sign_mismatch_area=accum.sign_mismatch_area,
        sampled_area=accum.sampled_area,
    )


def compare_paired_partition(
    *,
    master_vertices: list[Point],
    master_regions: list[tuple[Triangle, int]],
    dependent_vertices: list[Point],
    dependent_regions: list[tuple[Triangle, int]],
    dependent_to_master: Matrix,
) -> PartitionComparison:
    master_triangles = [triangle for triangle, _region in master_regions]
    dependent_triangles = [triangle for triangle, _region in dependent_regions]
    _validate_triangles(master_vertices, master_triangles, "master")
    _validate_triangles(dependent_vertices, dependent_triangles, "dependent")
    extent = max(
        1.0,
        _coordinate_extent(master_vertices),
        _coordinate_extent(dependent_vertices),
    )
    distance_tolerance = extent * 1e-9
    accum = _PartitionAccumulator()
    _accumulate_partition_direction(
        source_vertices=dependent_vertices,
        source_regions=dependent_regions,
        target_vertices=master_vertices,
        target_regions=master_regions,
        source_to_target=dependent_to_master,
        distance_tolerance=distance_tolerance,
        accum=accum,
    )
    _accumulate_partition_direction(
        source_vertices=master_vertices,
        source_regions=master_regions,
        target_vertices=dependent_vertices,
        target_regions=dependent_regions,
        source_to_target=transpose_matrix(dependent_to_master),
        distance_tolerance=distance_tolerance,
        accum=accum,
    )
    return PartitionComparison(
        sample_count=accum.sample_count,
        uncovered_samples=accum.uncovered_samples,
        mismatch_area=accum.mismatch_area,
        sampled_area=accum.sampled_area,
    )


def mapped_zero_contour_hausdorff(
    *,
    master_vertices: list[Point],
    master_triangles: list[Triangle],
    master_values: list[float],
    dependent_vertices: list[Point],
    dependent_triangles: list[Triangle],
    dependent_values: list[float],
    dependent_to_master: Matrix,
) -> float:
    _validate_trace_input(master_vertices, master_triangles, master_values, "master")
    _validate_trace_input(
        dependent_vertices,
        dependent_triangles,
        dependent_values,
        "dependent",
    )
    master_segments = _extract_zero_segments(
        master_vertices,
        master_triangles,
        master_values,
    )
    dependent_segments = [
        (
            apply_matrix(dependent_to_master, segment[0]),
            apply_matrix(dependent_to_master, segment[1]),
        )
        for segment in _extract_zero_segments(
            dependent_vertices,
            dependent_triangles,
            dependent_values,
        )
    ]
    if not master_segments and not dependent_segments:
        return 0.0
    if not master_segments or not dependent_segments:
        return math.inf
    return max(
        _directed_endpoint_hausdorff(master_segments, dependent_segments),
        _directed_endpoint_hausdorff(dependent_segments, master_segments),
    )


def mapped_segment_hausdorff(
    *,
    master_segments: list[Segment],
    dependent_segments: list[Segment],
    dependent_to_master: Matrix,
) -> float:
    mapped_dependent = [
        (
            apply_matrix(dependent_to_master, segment[0]),
            apply_matrix(dependent_to_master, segment[1]),
        )
        for segment in dependent_segments
    ]
    if not master_segments and not mapped_dependent:
        return 0.0
    if not master_segments or not mapped_dependent:
        return math.inf
    return max(
        _directed_endpoint_hausdorff(master_segments, mapped_dependent),
        _directed_endpoint_hausdorff(mapped_dependent, master_segments),
    )


@dataclass
class _TraceAccumulator:
    sample_count: int = 0
    uncovered_samples: int = 0
    max_abs_difference: float = 0.0
    weighted_squared_error: float = 0.0
    sign_mismatch_area: float = 0.0
    sampled_area: float = 0.0


@dataclass
class _PartitionAccumulator:
    sample_count: int = 0
    uncovered_samples: int = 0
    mismatch_area: float = 0.0
    sampled_area: float = 0.0


def _accumulate_trace_direction(
    *,
    source_vertices: list[Point],
    source_triangles: list[Triangle],
    source_values: list[float],
    target_vertices: list[Point],
    target_triangles: list[Triangle],
    target_values: list[float],
    source_to_target: Matrix,
    distance_tolerance: float,
    accum: _TraceAccumulator,
) -> None:
    target_locator = _TraceLocator(
        target_vertices,
        target_triangles,
        distance_tolerance,
    )
    for triangle in source_triangles:
        source_points = tuple(source_vertices[index] for index in triangle)
        source_area = _triangle_area(*source_points)
        if source_area <= 0.0:
            continue
        weight = source_area / len(_TRIANGLE_SAMPLES)
        source_nodal_values = tuple(source_values[index] for index in triangle)
        for barycentric in _TRIANGLE_SAMPLES:
            source_point = _interpolate_point(source_points, barycentric)
            target_point = apply_matrix(source_to_target, source_point)
            located = target_locator.locate(target_point)
            if located is None:
                accum.uncovered_samples += 1
                continue
            target_triangle, target_barycentric = located
            source_value = _interpolate_scalar(source_nodal_values, barycentric)
            target_value = _interpolate_scalar(
                tuple(target_values[index] for index in target_triangle),
                target_barycentric,
            )
            difference = source_value - target_value
            accum.sample_count += 1
            accum.sampled_area += weight
            accum.max_abs_difference = max(
                accum.max_abs_difference,
                abs(difference),
            )
            accum.weighted_squared_error += weight * difference * difference
            if _scalar_sign(source_value) != _scalar_sign(target_value):
                accum.sign_mismatch_area += weight


def _accumulate_partition_direction(
    *,
    source_vertices: list[Point],
    source_regions: list[tuple[Triangle, int]],
    target_vertices: list[Point],
    target_regions: list[tuple[Triangle, int]],
    source_to_target: Matrix,
    distance_tolerance: float,
    accum: _PartitionAccumulator,
) -> None:
    target_triangles = [triangle for triangle, _region in target_regions]
    target_locator = _TraceLocator(
        target_vertices,
        target_triangles,
        distance_tolerance,
    )
    target_region_by_triangle = {
        triangle: region for triangle, region in target_regions
    }
    for triangle, source_region in source_regions:
        source_points = tuple(source_vertices[index] for index in triangle)
        source_area = _triangle_area(*source_points)
        if source_area <= 0.0:
            continue
        weight = source_area / len(_TRIANGLE_SAMPLES)
        for barycentric in _TRIANGLE_SAMPLES:
            target_point = apply_matrix(
                source_to_target,
                _interpolate_point(source_points, barycentric),
            )
            located = target_locator.locate(target_point)
            if located is None:
                accum.uncovered_samples += 1
                continue
            target_triangle, _target_barycentric = located
            accum.sample_count += 1
            accum.sampled_area += weight
            if source_region != target_region_by_triangle[target_triangle]:
                accum.mismatch_area += weight


def _validate_trace_input(
    vertices: list[Point],
    triangles: list[Triangle],
    values: list[float],
    label: str,
) -> None:
    if len(values) != len(vertices):
        raise ValueError(
            f"{label} trace value count does not match vertices: "
            f"values={len(values)} vertices={len(vertices)}"
        )
    _validate_triangles(vertices, triangles, label)


def _validate_triangles(
    vertices: list[Point],
    triangles: list[Triangle],
    label: str,
) -> None:
    for triangle in triangles:
        if len(set(triangle)) != 3:
            raise ValueError(f"{label} trace contains a degenerate triangle index")
        for index in triangle:
            if index < 0 or index >= len(vertices):
                raise ValueError(f"{label} trace triangle index is out of range")


def _extract_zero_segments(
    vertices: list[Point],
    triangles: list[Triangle],
    values: list[float],
) -> list[Segment]:
    coordinate_tolerance = max(1.0, _coordinate_extent(vertices)) * 1e-12
    value_scale = max([1.0, *(abs(value) for value in values)])
    value_tolerance = value_scale * 1e-14
    segments: list[Segment] = []
    for triangle in triangles:
        points = tuple(vertices[index] for index in triangle)
        nodal_values = tuple(values[index] for index in triangle)
        intersections: list[Point] = []
        for first, second in ((0, 1), (1, 2), (2, 0)):
            a = points[first]
            b = points[second]
            va = nodal_values[first]
            vb = nodal_values[second]
            a_zero = abs(va) <= value_tolerance
            b_zero = abs(vb) <= value_tolerance
            if a_zero:
                _append_unique_point(intersections, a, coordinate_tolerance)
            if b_zero:
                _append_unique_point(intersections, b, coordinate_tolerance)
            if not a_zero and not b_zero and va * vb < 0.0:
                fraction = va / (va - vb)
                point = _add(a, _scale(fraction, _sub(b, a)))
                _append_unique_point(intersections, point, coordinate_tolerance)
        if len(intersections) < 2:
            continue
        first, second = _farthest_pair(intersections)
        if _distance_squared(first, second) > coordinate_tolerance * coordinate_tolerance:
            segments.append((first, second))
    return segments


def _append_unique_point(points: list[Point], point: Point, tolerance: float) -> None:
    tolerance_squared = tolerance * tolerance
    if all(_distance_squared(point, existing) > tolerance_squared for existing in points):
        points.append(point)


def _farthest_pair(points: list[Point]) -> tuple[Point, Point]:
    best = points[0], points[1]
    best_distance = _distance_squared(*best)
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            distance = _distance_squared(points[i], points[j])
            if distance > best_distance:
                best = points[i], points[j]
                best_distance = distance
    return best


def _directed_endpoint_hausdorff(
    source: list[Segment],
    target: list[Segment],
) -> float:
    maximum = 0.0
    for segment in source:
        for point in segment:
            distance = min(_point_segment_distance(point, other) for other in target)
            maximum = max(maximum, distance)
    return maximum


def _point_segment_distance(point: Point, segment: Segment) -> float:
    a, b = segment
    direction = _sub(b, a)
    denominator = _dot(direction, direction)
    if denominator <= 0.0:
        return math.sqrt(_distance_squared(point, a))
    fraction = max(0.0, min(1.0, _dot(_sub(point, a), direction) / denominator))
    projection = _add(a, _scale(fraction, direction))
    return math.sqrt(_distance_squared(point, projection))


def _distance_squared(a: Point, b: Point) -> float:
    difference = _sub(a, b)
    return _dot(difference, difference)


def _coordinate_extent(vertices: list[Point]) -> float:
    if not vertices:
        return 0.0
    return max(abs(coordinate) for point in vertices for coordinate in point)


class _TraceLocator:
    def __init__(
        self,
        vertices: list[Point],
        triangles: list[Triangle],
        distance_tolerance: float,
    ) -> None:
        self.vertices = vertices
        self.triangles = triangles
        self.distance_tolerance_squared = distance_tolerance * distance_tolerance
        referenced_vertices = sorted({vertex for triangle in triangles for vertex in triangle})
        surface_points = [vertices[vertex] for vertex in referenced_vertices]
        spans = [
            max((point[axis] for point in surface_points), default=0.0)
            - min((point[axis] for point in surface_points), default=0.0)
            for axis in range(3)
        ]
        self.axes = tuple(sorted(range(3), key=lambda axis: spans[axis], reverse=True)[:2])
        self.minimum = tuple(
            min((point[axis] for point in surface_points), default=0.0)
            for axis in self.axes
        )
        self.maximum = tuple(
            max((point[axis] for point in surface_points), default=0.0)
            for axis in self.axes
        )
        self.cells_per_axis = max(1, int(math.sqrt(max(1, len(triangles)))))
        self.grid: dict[tuple[int, int], list[int]] = {}
        for index, triangle in enumerate(triangles):
            points = [vertices[vertex] for vertex in triangle]
            lower = tuple(min(point[axis] for point in points) for axis in self.axes)
            upper = tuple(max(point[axis] for point in points) for axis in self.axes)
            first = tuple(self._cell(coordinate, dimension) for dimension, coordinate in enumerate(lower))
            last = tuple(self._cell(coordinate, dimension) for dimension, coordinate in enumerate(upper))
            for i in range(first[0], last[0] + 1):
                for j in range(first[1], last[1] + 1):
                    self.grid.setdefault((i, j), []).append(index)

    def locate(
        self,
        point: Point,
    ) -> tuple[Triangle, tuple[float, float, float]] | None:
        center = tuple(
            self._cell(point[axis], dimension)
            for dimension, axis in enumerate(self.axes)
        )
        candidates: set[int] = set()
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                key = (center[0] + di, center[1] + dj)
                candidates.update(self.grid.get(key, ()))
        barycentric_tolerance = 1e-9
        best: tuple[Triangle, tuple[float, float, float]] | None = None
        best_distance = math.inf
        for index in candidates:
            triangle = self.triangles[index]
            a, b, c = (self.vertices[vertex] for vertex in triangle)
            result = _triangle_barycentric(point, a, b, c)
            if result is None:
                continue
            barycentric, distance_squared = result
            if distance_squared > self.distance_tolerance_squared:
                continue
            if any(
                coordinate < -barycentric_tolerance
                or coordinate > 1.0 + barycentric_tolerance
                for coordinate in barycentric
            ):
                continue
            if distance_squared < best_distance:
                best = triangle, barycentric
                best_distance = distance_squared
        return best

    def _cell(self, coordinate: float, dimension: int) -> int:
        lower = self.minimum[dimension]
        upper = self.maximum[dimension]
        if upper <= lower:
            return 0
        scaled = (coordinate - lower) / (upper - lower)
        return max(
            0,
            min(self.cells_per_axis - 1, int(math.floor(scaled * self.cells_per_axis))),
        )


def _triangle_barycentric(
    point: Point,
    a: Point,
    b: Point,
    c: Point,
) -> tuple[tuple[float, float, float], float] | None:
    v0 = _sub(b, a)
    v1 = _sub(c, a)
    v2 = _sub(point, a)
    normal = _cross(v0, v1)
    dropped_axis = max(range(3), key=lambda axis: abs(normal[axis]))
    first_axis, second_axis = (
        axis for axis in range(3) if axis != dropped_axis
    )
    denominator = (
        v0[first_axis] * v1[second_axis]
        - v0[second_axis] * v1[first_axis]
    )
    if abs(denominator) <= 1e-30:
        return None
    second = (
        v2[first_axis] * v1[second_axis]
        - v2[second_axis] * v1[first_axis]
    ) / denominator
    third = (
        v0[first_axis] * v2[second_axis]
        - v0[second_axis] * v2[first_axis]
    ) / denominator
    first = 1.0 - second - third
    projection = _add(a, _add(_scale(second, v0), _scale(third, v1)))
    distance_squared = _dot(_sub(point, projection), _sub(point, projection))
    return (first, second, third), distance_squared


def _interpolate_point(
    points: tuple[Point, Point, Point],
    barycentric: tuple[float, float, float],
) -> Point:
    return tuple(
        sum(barycentric[i] * points[i][component] for i in range(3))
        for component in range(3)
    )  # type: ignore[return-value]


def _interpolate_scalar(
    values: tuple[float, float, float],
    barycentric: tuple[float, float, float],
) -> float:
    return sum(barycentric[i] * values[i] for i in range(3))


def _triangle_area(a: Point, b: Point, c: Point) -> float:
    cross = _cross(_sub(b, a), _sub(c, a))
    return 0.5 * math.sqrt(_dot(cross, cross))


def _scalar_sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def _add(a: Point, b: Point) -> Point:
    return a[0] + b[0], a[1] + b[1], a[2] + b[2]


def _sub(a: Point, b: Point) -> Point:
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def _scale(factor: float, point: Point) -> Point:
    return factor * point[0], factor * point[1], factor * point[2]


def _dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Point, b: Point) -> Point:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _determinant(matrix: Matrix) -> float:
    return (
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def proper_tetrahedral_rotations() -> list[Matrix]:
    """Return the 12 rotations preserving the tetrahedron xyz-sign product +1."""
    rotations: set[Matrix] = set()
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        if inversions % 2:
            continue
        for signs in product((-1.0, 1.0), repeat=3):
            if signs[0] * signs[1] * signs[2] != 1.0:
                continue
            rows = [[0.0, 0.0, 0.0] for _ in range(3)]
            for row, column in enumerate(permutation):
                rows[row][column] = signs[row]
            matrix: Matrix = tuple(tuple(row) for row in rows)  # type: ignore[assignment]
            if round(_determinant(matrix)) == 1:
                rotations.add(matrix)
    ordered = sorted(rotations)
    if len(ordered) != 12:
        raise RuntimeError(f"expected 12 proper tetrahedral rotations, got {len(ordered)}")
    return ordered

full12_watertightness = full24_watertightness
weld_full12_surface = weld_full24_surface
