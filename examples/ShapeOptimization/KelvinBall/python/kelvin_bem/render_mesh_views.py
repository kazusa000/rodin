#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Render basic PNG views of a MEDIT GammaShape surface."""

from __future__ import annotations

import argparse
from pathlib import Path


Point = tuple[float, float, float]
Triangle = tuple[int, int, int]
GAMMA_SHAPE = 13
VIEW_SPECS = (
    ("front", 0.0, -90.0),
    ("side", 0.0, 0.0),
    ("top", 90.0, -90.0),
    ("iso", 25.0, -45.0),
)


def parse_medit_surface(path: Path) -> tuple[list[Point], list[Triangle]]:
    lines = [line.strip() for line in path.read_text(encoding="ascii", errors="replace").splitlines()]
    vertices: list[Point] = []
    triangles: list[Triangle] = []
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
                if int(attr) == GAMMA_SHAPE:
                    triangles.append((int(a) - 1, int(b) - 1, int(c) - 1))
            i += 2 + count
            continue
        i += 1
    if not vertices:
        raise ValueError(f"No vertices found in {path}")
    if not triangles:
        raise ValueError(f"No GammaShape={GAMMA_SHAPE} triangles found in {path}")
    return vertices, triangles


def axis_limits(vertices: list[Point]) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    xs = [p[0] for p in vertices]
    ys = [p[1] for p in vertices]
    zs = [p[2] for p in vertices]
    center = ((min(xs) + max(xs)) / 2.0, (min(ys) + max(ys)) / 2.0, (min(zs) + max(zs)) / 2.0)
    radius = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)) / 2.0
    radius = max(radius, 1e-12) * 1.08
    return tuple((c - radius, c + radius) for c in center)  # type: ignore[return-value]


def surface_axis_limits(
    vertices: list[Point],
    triangles: list[Triangle],
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    used = sorted({idx for tri in triangles for idx in tri})
    return axis_limits([vertices[idx] for idx in used])


def render_views(mesh: Path, out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    vertices, triangles = parse_medit_surface(mesh)
    polys = [[vertices[a], vertices[b], vertices[c]] for a, b, c in triangles]
    limits = surface_axis_limits(vertices, triangles)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for name, elev, azim in VIEW_SPECS:
        fig = plt.figure(figsize=(6, 6), dpi=180)
        ax = fig.add_subplot(111, projection="3d")
        collection = Poly3DCollection(polys, linewidths=0.08, alpha=1.0)
        collection.set_facecolor((0.76, 0.82, 0.90, 1.0))
        collection.set_edgecolor((0.18, 0.20, 0.24, 0.20))
        ax.add_collection3d(collection)
        ax.set_xlim(*limits[0])
        ax.set_ylim(*limits[1])
        ax.set_zlim(*limits[2])
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        try:
            ax.set_box_aspect((1, 1, 1))
        except AttributeError:
            pass
        fig.tight_layout(pad=0)
        target = out_dir / f"{name}.png"
        fig.savefig(target, transparent=False, facecolor="white")
        plt.close(fig)
        written.append(target)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for path in render_views(args.mesh, args.out_dir):
        print(path)


if __name__ == "__main__":
    main()
