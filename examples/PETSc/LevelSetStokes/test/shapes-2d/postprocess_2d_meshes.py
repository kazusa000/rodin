from __future__ import annotations

from pathlib import Path
from typing import Callable


BASE = Path(__file__).resolve().parents[2] / "2dshapes"


def inside_annulus(x: float, y: float) -> bool:
    dx = x - 0.5
    dy = y - 0.5
    r2 = dx * dx + dy * dy
    return (0.08 ** 2) <= r2 <= (0.18 ** 2)


def inside_eccentric_ring(x: float, y: float) -> bool:
    outer = ((x - 0.5) / 0.19) ** 2 + ((y - 0.5) / 0.16) ** 2 <= 1.0
    inner = ((x - 0.56) / 0.10) ** 2 + ((y - 0.48) / 0.07) ** 2 <= 1.0
    return outer and not inner


def inside_double_hole(x: float, y: float) -> bool:
    outer = ((x - 0.5) / 0.22) ** 2 + ((y - 0.5) / 0.18) ** 2 <= 1.0
    hole1 = ((x - 0.43) / 0.07) ** 2 + ((y - 0.5) / 0.05) ** 2 <= 1.0
    hole2 = ((x - 0.57) / 0.07) ** 2 + ((y - 0.5) / 0.05) ** 2 <= 1.0
    return outer and not hole1 and not hole2


def is_outer_boundary(x: float, y: float, tol: float = 5e-3) -> bool:
    return (
        abs(x) < tol
        or abs(y) < tol
        or abs(x - 1.0) < tol
        or abs(y - 1.0) < tol
    )


def rewrite_mesh(path: Path, inside: Callable[[float, float], bool]) -> None:
    lines = path.read_text().splitlines()
    out: list[str] = []
    vertices: dict[int, tuple[float, float]] = {}

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        out.append(lines[i])
        i += 1

        if line == "Vertices":
            count = int(lines[i].strip())
            out.append(lines[i])
            i += 1
            for vid in range(1, count + 1):
                raw = lines[i]
                parts = raw.split()
                x, y = float(parts[0]), float(parts[1])
                vertices[vid] = (x, y)
                out.append(raw)
                i += 1
            continue

        if line == "Edges":
            count = int(lines[i].strip())
            out.append(lines[i])
            i += 1
            for _ in range(count):
                parts = lines[i].split()
                a, b = int(parts[0]), int(parts[1])
                xa, ya = vertices[a]
                xb, yb = vertices[b]
                mx = 0.5 * (xa + xb)
                my = 0.5 * (ya + yb)
                ref = 7 if is_outer_boundary(mx, my) else 13
                out.append(f"{a} {b} {ref}")
                i += 1
            continue

        if line == "Triangles":
            count = int(lines[i].strip())
            out.append(lines[i])
            i += 1
            for _ in range(count):
                parts = lines[i].split()
                a, b, c = int(parts[0]), int(parts[1]), int(parts[2])
                xa, ya = vertices[a]
                xb, yb = vertices[b]
                xc, yc = vertices[c]
                mx = (xa + xb + xc) / 3.0
                my = (ya + yb + yc) / 3.0
                ref = 2 if inside(mx, my) else 3
                out.append(f"{a} {b} {c} {ref}")
                i += 1
            continue

    path.write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    rewrite_mesh(BASE / "annulus_init_2d.mesh", inside_annulus)
    rewrite_mesh(BASE / "eccentric_ring_init_2d.mesh", inside_eccentric_ring)
    rewrite_mesh(BASE / "double_hole_init_2d.mesh", inside_double_hole)
