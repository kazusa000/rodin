#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

BASE = Path("/home/wjj/Project/stage/rodin/examples/PETSc/LevelSetStokes")
BASE_3D = BASE / "3dshapes"

EPS = 1e-8


def parse_mesh(path: Path):
    lines = path.read_text().splitlines()
    i = 0
    data = {}
    while i < len(lines):
        key = lines[i].strip()
        if key in {"Vertices", "Triangles", "Tetrahedra"}:
            n = int(lines[i + 1].strip())
            data[key] = lines[i + 2 : i + 2 + n]
            i += 2 + n
        else:
            i += 1
    return lines, data


def write_mesh(path: Path, lines, vertices, triangles, tetrahedra):
    out = []
    i = 0
    while i < len(lines):
        key = lines[i].strip()
        if key == "Vertices":
            out.extend([lines[i], lines[i + 1], *vertices])
            n = int(lines[i + 1].strip())
            i += 2 + n
        elif key == "Triangles":
            out.extend([lines[i], lines[i + 1], *triangles])
            n = int(lines[i + 1].strip())
            i += 2 + n
        elif key == "Tetrahedra":
            out.extend([lines[i], lines[i + 1], *tetrahedra])
            n = int(lines[i + 1].strip())
            i += 2 + n
        else:
            out.append(lines[i])
            i += 1
    path.write_text("\n".join(out) + "\n")


def on_outer_boundary(pt):
    return any(abs(c) < EPS or abs(c - 1.0) < EPS for c in pt)


def inside_ellipsoid(pt):
    x, y, z = pt
    return ((x - 0.5) / 0.275) ** 2 + ((y - 0.5) / 0.187) ** 2 + ((z - 0.5) / 0.187) ** 2 <= 1.0


def inside_oblate(pt):
    x, y, z = pt
    return ((x - 0.5) / 0.242) ** 2 + ((y - 0.5) / 0.242) ** 2 + ((z - 0.5) / 0.121) ** 2 <= 1.0


def inside_prolate(pt):
    x, y, z = pt
    return ((x - 0.5) / 0.187) ** 2 + ((y - 0.5) / 0.187) ** 2 + ((z - 0.5) / 0.297) ** 2 <= 1.0


def inside_capsule(pt):
    x, y, z = pt
    ax, ay, az = 0.34, 0.5, 0.5
    bx, by, bz = 0.66, 0.5, 0.5
    px, py, pz = x - ax, y - ay, z - az
    vx, vy, vz = bx - ax, by - ay, bz - az
    vv = vx * vx + vy * vy + vz * vz
    t = max(0.0, min(1.0, (px * vx + py * vy + pz * vz) / vv))
    qx = ax + t * vx
    qy = ay + t * vy
    qz = az + t * vz
    dx, dy, dz = x - qx, y - qy, z - qz
    return dx * dx + dy * dy + dz * dz <= 0.12 * 0.12


def inside_double_sphere(pt):
    x, y, z = pt
    s1 = (x - 0.42) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.18 * 0.18
    s2 = (x - 0.58) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.18 * 0.18
    return s1 or s2


def inside_triaxial(pt):
    x, y, z = pt
    return ((x - 0.5) / 0.26) ** 2 + ((y - 0.5) / 0.20) ** 2 + ((z - 0.5) / 0.16) ** 2 <= 1.0


def inside_peanut(pt):
    x, y, z = pt
    s1 = (x - 0.41) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.17 * 0.17
    s2 = (x - 0.59) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.17 * 0.17
    return s1 or s2


def inside_three_lobe(pt):
    x, y, z = pt
    c0 = (x - 0.52) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.11 * 0.11
    c1 = (x - 0.62) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.14 * 0.14
    c2 = (x - 0.45) ** 2 + (y - 0.60) ** 2 + (z - 0.5) ** 2 <= 0.14 * 0.14
    c3 = (x - 0.45) ** 2 + (y - 0.40) ** 2 + (z - 0.5) ** 2 <= 0.14 * 0.14
    return c0 or c1 or c2 or c3


def inside_dumbbell(pt):
    x, y, z = pt
    s1 = (x - 0.38) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.15 * 0.15
    s2 = (x - 0.62) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.15 * 0.15
    ax, ay, az = 0.38, 0.5, 0.5
    bx, by, bz = 0.62, 0.5, 0.5
    px, py, pz = x - ax, y - ay, z - az
    vx, vy, vz = bx - ax, by - ay, bz - az
    vv = vx * vx + vy * vy + vz * vz
    t = max(0.0, min(1.0, (px * vx + py * vy + pz * vz) / vv))
    qx = ax + t * vx
    qy = ay + t * vy
    qz = az + t * vz
    dx, dy, dz = x - qx, y - qy, z - qz
    neck = dx * dx + dy * dy + dz * dz <= 0.07 * 0.07
    return s1 or s2 or neck


def inside_sphere_bump(pt):
    x, y, z = pt
    sphere = (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 <= 0.19 * 0.19
    bump = (x - 0.62) ** 2 + (y - 0.54) ** 2 + (z - 0.5) ** 2 <= 0.05 * 0.05
    return sphere or bump


def inside_torus(pt):
    x, y, z = pt
    rho = math.hypot(x - 0.5, y - 0.5)
    return (rho - 0.14) ** 2 + (z - 0.5) ** 2 <= 0.08 * 0.08


def inside_rounded_three_lobe(pt):
    x, y, z = pt
    c0 = (x - 0.50) ** 2 + (y - 0.50) ** 2 + (z - 0.50) ** 2 <= 0.12 * 0.12
    c1 = (x - 0.62) ** 2 + (y - 0.50) ** 2 + (z - 0.50) ** 2 <= 0.13 * 0.13
    c2 = (x - 0.44) ** 2 + (y - 0.60) ** 2 + (z - 0.50) ** 2 <= 0.13 * 0.13
    c3 = (x - 0.44) ** 2 + (y - 0.40) ** 2 + (z - 0.50) ** 2 <= 0.13 * 0.13
    return c0 or c1 or c2 or c3


def inside_rounded_four_lobe(pt):
    x, y, z = pt
    c0 = (x - 0.50) ** 2 + (y - 0.50) ** 2 + (z - 0.50) ** 2 <= 0.11 * 0.11
    c1 = (x - 0.50) ** 2 + (y - 0.66) ** 2 + (z - 0.50) ** 2 <= 0.12 * 0.12
    c2 = (x - 0.50) ** 2 + (y - 0.34) ** 2 + (z - 0.50) ** 2 <= 0.12 * 0.12
    c3 = (x - 0.66) ** 2 + (y - 0.50) ** 2 + (z - 0.50) ** 2 <= 0.12 * 0.12
    c4 = (x - 0.34) ** 2 + (y - 0.50) ** 2 + (z - 0.50) ** 2 <= 0.12 * 0.12
    return c0 or c1 or c2 or c3 or c4


def inside_seven_lobe_star(pt):
    x, y, z = pt
    bulbs = [
        (0.50, 0.50, 0.50, 0.095),
        (0.50, 0.69, 0.50, 0.095),
        (0.65, 0.62, 0.50, 0.095),
        (0.69, 0.47, 0.50, 0.095),
        (0.60, 0.33, 0.50, 0.095),
        (0.43, 0.30, 0.50, 0.095),
        (0.31, 0.41, 0.50, 0.095),
        (0.34, 0.59, 0.50, 0.095),
    ]
    return any((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r * r for cx, cy, cz, r in bulbs)


SHAPES = {
    "ellipsoid_init.mesh": inside_ellipsoid,
    "oblate_init.mesh": inside_oblate,
    "prolate_init.mesh": inside_prolate,
    "capsule_init.mesh": inside_capsule,
    "double_sphere_init.mesh": inside_double_sphere,
    "triaxial_init.mesh": inside_triaxial,
    "peanut_init.mesh": inside_peanut,
    "three_lobe_init.mesh": inside_three_lobe,
    "dumbbell_init.mesh": inside_dumbbell,
    "sphere_bump_init.mesh": inside_sphere_bump,
    "torus_init.mesh": inside_torus,
    "rounded_three_lobe_init.mesh": inside_rounded_three_lobe,
    "rounded_four_lobe_init.mesh": inside_rounded_four_lobe,
    "seven_lobe_star_init.mesh": inside_seven_lobe_star,
}


for name, inside_fn in SHAPES.items():
    path = BASE_3D / name
    lines, data = parse_mesh(path)

    vertices = []
    coords = [None]
    for line in data["Vertices"]:
        x, y, z, _ref = line.split()
        coords.append((float(x), float(y), float(z)))
        vertices.append(f"{x} {y} {z} 0")

    triangles = []
    for line in data["Triangles"]:
        a, b, c, _ref = line.split()
        ia, ib, ic = int(a), int(b), int(c)
        cx = (coords[ia][0] + coords[ib][0] + coords[ic][0]) / 3.0
        cy = (coords[ia][1] + coords[ib][1] + coords[ic][1]) / 3.0
        cz = (coords[ia][2] + coords[ib][2] + coords[ic][2]) / 3.0
        ref = 7 if on_outer_boundary((cx, cy, cz)) else 13
        triangles.append(f"{a} {b} {c} {ref}")

    tetrahedra = []
    for line in data["Tetrahedra"]:
        a, b, c, d, _ref = line.split()
        ids = [int(a), int(b), int(c), int(d)]
        cx = sum(coords[i][0] for i in ids) / 4.0
        cy = sum(coords[i][1] for i in ids) / 4.0
        cz = sum(coords[i][2] for i in ids) / 4.0
        ref = 2 if inside_fn((cx, cy, cz)) else 3
        tetrahedra.append(f"{a} {b} {c} {d} {ref}")

    write_mesh(path, lines, vertices, triangles, tetrahedra)
