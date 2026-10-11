# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Solid connectivity through positive-area faces, including nonmatching T12 cuts.

This checks the tetrahedral solid, not the open exported Gamma triangle soup.
Point/edge contact does not connect volumes. No mesh is altered or welded.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def triangle_overlap_area(a, b):
    """Convex clipping in a common plane, independent of triangle orientation."""
    def cross(u, v):
        return float(u[0]*v[1] - u[1]*v[0])
    orientation = np.sign(cross(b[1]-b[0], b[2]-b[0]))
    polygon = list(a)
    for p, q in zip(b, np.roll(b, -1, axis=0)):
        if not polygon:
            return 0.
        clipped = []
        for s, t in zip(polygon, polygon[1:]+polygon[:1]):
            ds = orientation*cross(q-p, s-p)
            dt = orientation*cross(q-p, t-p)
            if ds >= 0:
                clipped.append(s)
            if (ds < 0) != (dt < 0):
                clipped.append(s+(t-s)*ds/(ds-dt))
        polygon = clipped
    if len(polygon) < 3:
        return 0.
    p = np.asarray(polygon)
    return abs(sum(cross(s, t) for s, t in zip(p, np.roll(p, -1, axis=0))))/2


def solid_connectivity(x, cells, labels, rotations, *, symmetry="terra"):
    """Count solid components after gluing positive-area rotated cut overlaps."""
    t = np.asarray(cells)[np.asarray(labels) == 2]
    if not len(t):
        raise ValueError("empty solid")
    f = np.sort(t[:, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]].reshape(-1, 3), axis=1)
    unique, inverse, counts = np.unique(f, axis=0, return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError("nonmanifold solid tetrahedral face")
    owners = np.repeat(np.arange(len(t)), 4)
    order = np.argsort(inverse, kind="stable")
    offsets = np.r_[0, np.cumsum(counts)]
    shared = np.flatnonzero(counts == 2)
    a, b = owners[order[offsets[shared]]], owners[order[offsets[shared]+1]]
    n, components = connected_components(
        coo_matrix((np.ones(len(a)), (a, b)), shape=(len(t), len(t))), directed=False)
    parent = list(range(n*len(rotations)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    exterior = np.flatnonzero(counts == 1)
    triangles = x[unique[exterior]]
    face_components = components[owners[order[offsets[exterior]]]]
    scale = max(float(np.max(np.abs(x))), 1.)
    tol = 1e-8*scale
    if symmetry not in ("terra", "cubic") or len(rotations) != (12 if symmetry == "terra" else 24):
        raise ValueError("solid connectivity symmetry/group mismatch")
    # Cubic chamber x>=y>=|z| has three planes; x=y is self-paired
    # across its positive/negative-z halves by a proper rotation.
    normals = np.array(([[-1, 0, 1], [-1, 0, -1], [0, -1, 1], [0, -1, -1]]
                        if symmetry == "terra" else [[-1, 1, 0], [0, -1, 1], [0, -1, -1]]), float)/np.sqrt(2)
    planes = {}
    for normal in normals:
        mask = np.max(np.abs(triangles @ normal), axis=1) <= tol
        for r, rotation in enumerate(rotations):
            rotation = np.asarray(rotation)
            outward = rotation @ normal
            first = np.flatnonzero(np.abs(outward) > 1e-8)[0]
            side = 1 if outward[first] > 0 else -1
            canonical = outward*side
            key = tuple(np.round(canonical, 8))
            projection = np.delete(np.eye(3), np.argmax(np.abs(canonical)), axis=0)
            points = (triangles[mask] @ rotation.T) @ projection.T
            ids = face_components[mask]+r*n
            planes.setdefault(key, {1: [], -1: []})[side].extend(zip(points, ids))
    joins = 0
    for sides in planes.values():
        if not sides[1] or not sides[-1]:
            continue
        right = np.array([item[0] for item in sides[-1]])
        ids = [int(item[1]) for item in sides[-1]]
        lower, upper = right.min(axis=1), right.max(axis=1)
        for triangle, owner in sides[1]:
            candidates = np.flatnonzero(np.all(upper > triangle.min(axis=0), axis=1)
                & np.all(lower < triangle.max(axis=0), axis=1))
            for j in candidates:
                ra, rb = root(int(owner)), root(ids[j])
                if ra != rb and triangle_overlap_area(triangle, right[j]) > 1e-12*scale**2:
                    parent[rb] = ra
                    joins += 1
    xyz = x[t]
    volumes = np.abs(np.linalg.det(xyz[:, 1:]-xyz[:, :1]))/6
    chamber_volumes = np.bincount(components, weights=volumes, minlength=n)
    full = {}
    for r in range(len(rotations)):
        for i, volume in enumerate(chamber_volumes):
            full[root(r*n+i)] = full.get(root(r*n+i), 0.)+float(volume)
    return dict(chamber_components=int(n), full_components=len(full),
        component_volumes_full=sorted(full.values(), reverse=True),
        positive_area_cut_joins=joins, group_order=len(rotations),
        method="solid tetrahedra shared faces and positive-area rotated cut overlaps",
        length_tolerance=tol, projected_area_tolerance=1e-12*scale**2)
