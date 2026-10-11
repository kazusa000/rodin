#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Independent linear-tetrahedron quality and unchanged-target fitting audit."""
import argparse
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from . import carlos_reference as common
from .clement_boundary import read_mesh
from .clement_rho import geometry_diagnostics


def quality(vertices, cells):
    tetra = vertices[cells]
    volume = np.abs(np.linalg.det(np.transpose(tetra[:,1:]-tetra[:,:1],(0,2,1))))/6
    edges = sum(np.sum((tetra[:,i]-tetra[:,j])**2, axis=1)
                for i in range(4) for j in range(i+1,4))
    return 12*np.power(3*volume,2/3)/edges


def evaluate(vertices, cells, values, points):
    tetra = vertices[cells]
    origin = tetra[:,0]
    inverse = np.linalg.inv(np.transpose(tetra[:,1:]-tetra[:,:1], (0,2,1)))
    low, high = tetra.min(axis=1), tetra.max(axis=1)
    tree = cKDTree(tetra.mean(axis=1))
    result = np.full(len(points), np.nan)
    for i, point in enumerate(points):
        _, near = tree.query(point, k=min(32,len(cells)))
        for candidates in (np.atleast_1d(near),
                           np.flatnonzero(np.all((point >= low-1e-10)&(point <= high+1e-10),axis=1))):
            if len(candidates) == 0:
                continue
            local = np.einsum("nij,nj->ni", inverse[candidates], point-origin[candidates])
            bary = np.column_stack((1-local.sum(axis=1),local))
            inside = np.flatnonzero((bary.min(axis=1) >= -1e-9)&(bary.max(axis=1) <= 1+1e-9))
            if len(inside):
                k = inside[0]
                result[i] = bary[k]@values[cells[candidates[k]]]
                break
    # A nearest cell accelerates the search, never supplies extrapolated values.
    return result


def audit(source, target, candidate):
    x, _, _, cells, _ = read_mesh(source)
    y, faces, labels, new_cells, _ = read_mesh(candidate)
    tokens = Path(target).read_text().split()
    start = tokens.index("SolAtVertices")
    count, fields, kind = map(int, tokens[start+1:start+4])
    values = np.array(list(map(float,tokens[start+4:start+4+count])))
    if count != len(x) or fields != 1 or kind != 1 or not np.isfinite(values).all():
        raise ValueError("target must be scalar P1 on the unchanged source mesh")
    gamma = faces[labels == 13]
    nodes = np.unique(gamma)
    samples = np.vstack((y[nodes], y[gamma].mean(axis=1)))
    fitted = evaluate(x, cells, values, samples)
    finite = fitted[np.isfinite(fitted)]
    q = quality(y, new_cells)
    fixed = {}
    for label in (7,31,32,33,35):
        ids = np.unique(faces[labels == label])
        if not len(ids):
            raise ValueError("a fixed boundary label is missing")
        points = y[ids]
        residual = (points[:,0]-2 if label == 7 else points[:,1]-points[:,2] if label == 31
                    else points[:,1]+points[:,2] if label == 32 else points[:,0]-points[:,1])
        fixed[str(label)] = float(np.max(np.abs(residual)))
    return dict(source_sha256=common.digest(source), target_sha256=common.digest(target),
        candidate_sha256=common.digest(candidate), target_recomputed=False,
        quality_definition="12*(3*abs(tetra_volume))^(2/3)/sum_squared_edges",
        quality_min=float(q.min()), quality_mean=float(q.mean()),
        fixed_plane_residuals=fixed, fit_samples=len(samples),
        fit_unlocated=int(np.count_nonzero(~np.isfinite(fitted))),
        level_set_residual_max=float(np.abs(finite).max()) if len(finite) else None,
        level_set_residual_rms=float(np.sqrt(np.mean(finite**2))) if len(finite) else None,
        geometry=geometry_diagnostics(candidate,"cubic"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("source","target","candidate","output"):
        p.add_argument("--"+name,type=Path,required=True)
    args = p.parse_args()
    output = common.external(args.output)
    if output.exists():
        raise FileExistsError(output)
    common.source_state(common.ROOT)
    common.save_json(output,audit(args.source,args.target,args.candidate))


if __name__ == "__main__":
    main()
