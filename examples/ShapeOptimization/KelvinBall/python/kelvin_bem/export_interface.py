#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Export the unchanged material interface from a Carlos chamber HDF5 snapshot."""

import argparse
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np

from .run_support import ROOT
from . import mesh_io as mesh_io
from . import periodic_geometry as geometry


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load_snapshot(path):
    with h5py.File(path, "r") as data:
        vertices = np.asarray(data["/Mesh/Geometry/Vertices"], dtype=float)
        stream = np.asarray(data["/Mesh/XDMF/Topology"], dtype=np.int64)
        labels = np.asarray(data["/Mesh/Attributes/3"], dtype=np.int64).reshape(-1)
    if stream.ndim != 1 or len(stream) != 5 * len(labels):
        raise ValueError("expected upstream linear tetrahedral mixed topology")
    cells = stream.reshape(-1, 5)
    if not np.all(cells[:, 0] == 6):
        raise ValueError("snapshot contains non-tetrahedral elements")
    tets = cells[:, 1:]
    if vertices.ndim != 2 or vertices.shape[1] != 3 or not np.isfinite(vertices).all():
        raise ValueError("invalid coordinates")
    if len(tets) == 0 or tets.min() < 0 or tets.max() >= len(vertices):
        raise ValueError("invalid tetrahedral connectivity")
    if set(labels.tolist()) != {2, 3}:
        raise ValueError("expected a chamber containing both solid=2 and fluid=3")
    return vertices, tets, labels


def material_interface(vertices, tets, labels):
    """Find solid/fluid common faces, excluding artificial chamber boundaries."""
    templates = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
    faces = np.sort(tets[:, templates].reshape(-1, 3), axis=1)
    owners = np.repeat(np.arange(len(tets)), 4)
    _, inverse, counts = np.unique(faces, axis=0, return_inverse=True, return_counts=True)
    if counts.max() > 2:
        raise ValueError("nonmanifold tetrahedral face")
    solid_count = np.bincount(inverse, weights=(labels[owners] == 2).astype(int))
    interface = (counts == 2) & (solid_count == 1)
    selection = interface[inverse] & (labels[owners] == 2)
    triangles = faces[selection].copy()
    if not len(triangles):
        raise ValueError("no shared solid/fluid interface")
    xyz = vertices[triangles]
    normals = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    if np.any(np.linalg.norm(normals, axis=1) == 0):
        raise ValueError("degenerate interface triangle")
    inside = vertices[tets[owners[selection]]].mean(axis=1) - xyz.mean(axis=1)
    flip = np.einsum("ij,ij->i", normals, inside) > 0
    triangles[flip] = triangles[flip][:, [0, 2, 1]]
    tetra_xyz = vertices[tets]
    det = np.einsum("ij,ij->i", tetra_xyz[:, 1] - tetra_xyz[:, 0],
                    np.cross(tetra_xyz[:, 2] - tetra_xyz[:, 0], tetra_xyz[:, 3] - tetra_xyz[:, 0]))
    if np.any(det == 0):
        raise ValueError("zero-volume tetrahedron")
    return triangles, float(np.abs(det[labels == 2]).sum() / 6)


def seam_report(vertices, gamma, tolerance=1e-8, symmetry="cubic"):
    """Distinguish unmatched triangulation from displaced paired seam curves."""
    edges = np.sort(np.concatenate((gamma[:, [0, 1]], gamma[:, [1, 2]], gamma[:, [2, 0]])), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    cuts = {name: [] for name in ("plus", "minus", "xy_plus", "xy_minus")}
    unclassified = 0
    for edge in unique[counts == 1]:
        pair = vertices[edge]
        if symmetry == "terra":
            name = None
            for candidate, axis, sign in (("plus", 1, -1), ("minus", 0, -1),
                                           ("xy_plus", 0, 1), ("xy_minus", 1, 1)):
                if np.max(np.abs(pair[:, axis] + sign * pair[:, 2])) <= tolerance:
                    name = candidate
                    break
            if name is None:
                unclassified += 1
                continue
        elif np.max(np.abs(pair[:, 0] - pair[:, 1])) <= tolerance:
            name = "xy_minus" if pair[:, 2].mean() < 0 else "xy_plus"
        elif np.max(np.abs(pair[:, 1] - pair[:, 2])) <= tolerance:
            name = "plus"
        elif np.max(np.abs(pair[:, 1] + pair[:, 2])) <= tolerance:
            name = "minus"
        else:
            unclassified += 1
            continue
        cuts[name].append(tuple(map(tuple, pair)))
    result = {"unclassified_boundary_edges": unclassified, "plane_tolerance": tolerance}
    pairs = (
        ("y_minus_to_y_plus", "plus", "minus", ((1, 0, 0), (0, 0, -1), (0, 1, 0))),
        ("xy_minus_to_xy_plus", "xy_plus", "xy_minus", ((0, 1, 0), (1, 0, 0), (0, 0, -1))),
    ) if symmetry == "cubic" else (
        ("y_z_to_x_z", "minus", "plus", ((0, 0, 1), (1, 0, 0), (0, 1, 0))),
        ("x_neg_z_to_y_neg_z", "xy_minus", "xy_plus", ((0, 1, 0), (0, 0, -1), (-1, 0, 0))),
    )
    for name, master, slave, rotation in pairs:
        distance = geometry.mapped_segment_hausdorff(master_segments=cuts[master],
            dependent_segments=cuts[slave], dependent_to_master=rotation)
        result[name] = {"master_edges": len(cuts[master]), "slave_edges": len(cuts[slave]),
                        "bidirectional_endpoint_to_polyline_distance": distance if np.isfinite(distance) else None}
    return result


def export(source, output, symmetry="cubic"):
    if symmetry not in ("cubic", "terra"):
        raise ValueError("unknown symmetry")
    order = 12 if symmetry == "terra" else 24
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.is_relative_to(ROOT):
        raise ValueError("write generated meshes outside the source repository")
    report_path = output.with_suffix(".json")
    if output.exists() or report_path.exists():
        raise FileExistsError(output)
    vertices, tets, labels = load_snapshot(source)
    gamma, volume = material_interface(vertices, tets, labels)
    stats = geometry.full24_watertightness(
        chamber_vertices=vertices.tolist(), chamber_triangles=gamma.tolist(),
        rotations=(geometry.proper_tetrahedral_rotations() if symmetry == "terra"
                   else geometry.proper_octahedral_rotations()), tolerance=1e-8)
    from dataclasses import asdict
    report = {"source": str(source), "source_sha256": sha256(source),
              "symmetry": symmetry, "group_order": order, "vertices": len(vertices),
              "tetrahedra": len(tets), "interface_triangles": len(gamma),
              "solid_volume_chamber": volume, "solid_volume_full": order * volume,
              "geometry_modified": False, f"full{order}_watertightness": asdict(stats),
              "paired_seams": seam_report(vertices, gamma, symmetry=symmetry),
              "interface_method": "shared solid=2/fluid=3 tetrahedral faces; outward orientation"}
    output.parent.mkdir(parents=True, exist_ok=True)
    mesh_io.write_medit_mesh(output, vertices.tolist(),
        [(*map(lambda x: int(x) + 1, face), 13) for face in gamma],
        [(*map(lambda x: int(x) + 1, cell), int(label)) for cell, label in zip(tets, labels)])
    report["output_sha256"] = sha256(output)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--symmetry", choices=("cubic", "terra"), default="cubic")
    args = parser.parse_args()
    if os.environ.get("KELVIN_SYMMETRY", args.symmetry) != args.symmetry:
        parser.error("snapshot symmetry and KELVIN_SYMMETRY disagree")
    print(json.dumps(export(args.snapshot, args.output, args.symmetry), indent=2))


if __name__ == "__main__":
    main()
