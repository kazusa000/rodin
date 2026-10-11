# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Original Medit reader/writer, isolated from historical meshing drivers."""
from __future__ import annotations
from pathlib import Path
Point = tuple[float, float, float]
Triangle = tuple[int, int, int, int]
Tet = tuple[int, int, int, int, int]
GAMMA = 13
SIGMA_PLUS = 31
SIGMA_MINUS = 32
SIGMA_XY = 33
SIGMA_TETRA_DEPENDENT = 34


def read_medit_mesh(path: Path) -> tuple[list[Point], list[Triangle], list[Tet]]:
    vertices: list[Point] = []
    triangles: list[Triangle] = []
    tetrahedra: list[Tet] = []

    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        lines = iter(handle)
        for raw in lines:
            text = raw.strip()
            if not text or text.startswith("#"):
                continue
            key = text.lower()
            if key == "vertices":
                count = int(next(lines).strip())
                for _ in range(count):
                    values = next(lines).split()
                    vertices.append(tuple(float(value) for value in values[:3]))  # type: ignore[arg-type]
            elif key == "triangles":
                count = int(next(lines).strip())
                for _ in range(count):
                    a, b, c, attr = (int(value) for value in next(lines).split()[:4])
                    triangles.append((a, b, c, attr))
            elif key == "tetrahedra":
                count = int(next(lines).strip())
                for _ in range(count):
                    a, b, c, d, attr = (int(value) for value in next(lines).split()[:5])
                    tetrahedra.append((a, b, c, d, attr))
            elif key == "end":
                break

    if not vertices:
        raise RuntimeError(f"MEDIT mesh has no vertices: {path}")
    if not triangles:
        raise RuntimeError(f"MEDIT mesh has no triangles: {path}")
    if not tetrahedra:
        raise RuntimeError(f"MEDIT mesh has no tetrahedra: {path}")
    return vertices, triangles, tetrahedra


def write_medit_mesh(
    path: Path,
    vertices: list[Point],
    triangles: list[Triangle],
    tetrahedra: list[Tet],
) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("MeshVersionFormatted 2\n")
        handle.write("Dimension 3\n\n")
        handle.write("Vertices\n")
        handle.write(f"{len(vertices)}\n")
        for p in vertices:
            handle.write(f"{p[0]:.17g} {p[1]:.17g} {p[2]:.17g} 0\n")
        handle.write("\nTriangles\n")
        handle.write(f"{len(triangles)}\n")
        for a, b, c, attr in triangles:
            handle.write(f"{a} {b} {c} {attr}\n")
        handle.write("\nTetrahedra\n")
        handle.write(f"{len(tetrahedra)}\n")
        for a, b, c, d, attr in tetrahedra:
            handle.write(f"{a} {b} {c} {d} {attr}\n")
        handle.write("\nEnd\n")
