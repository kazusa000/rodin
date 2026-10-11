#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Combine exact XDMF coordinates with the same saved chamber's Medit labels."""

import argparse
from pathlib import Path

import numpy as np

from . import carlos_reference as common
from . import export_carlos_interface as interface


def match(snapshot, skeleton):
    vertices, tets, labels = snapshot
    rounded, faces, cells = skeleton
    rounded, faces, cells = np.asarray(rounded), np.asarray(faces, int), np.asarray(cells, int)
    if rounded.shape != vertices.shape or not np.array_equal(cells[:, :4]-1, tets):
        raise ValueError("skeleton has different nodes or ordered tetrahedral connectivity")
    if not np.array_equal(cells[:, 4], labels):
        raise ValueError("skeleton material labels differ from the exact snapshot")
    if not np.isfinite(rounded).all() or not np.allclose(rounded, vertices, rtol=5e-6, atol=5e-12):
        raise ValueError("skeleton is not the rounded serialization of this snapshot")
    gamma, _ = interface.material_interface(vertices, tets, labels)
    expected = {tuple(sorted(face)) for face in gamma}
    saved = {tuple(sorted(face[:3]-1)) for face in faces if face[3] == 13}
    if saved != expected:
        raise ValueError("skeleton Gamma differs from the exact material interface")
    return vertices, faces, cells


def export(snapshot, skeleton, output):
    snapshot, skeleton, output = map(lambda p: Path(p).resolve(), (snapshot, skeleton, output))
    output = common.external(output)
    certificate = output.with_suffix(".json")
    if output.exists() or certificate.exists():
        raise FileExistsError(output)
    sources = {str(p): common.digest(p) for p in (snapshot, skeleton)}
    vertices, faces, cells = match(interface.load_snapshot(snapshot),
                                  interface.mesh_io.read_medit_mesh(skeleton))
    output.parent.mkdir(parents=True, exist_ok=True)
    interface.mesh_io.write_medit_mesh(output, vertices.tolist(), faces.tolist(), cells.tolist())
    restored, restored_faces, restored_cells = interface.mesh_io.read_medit_mesh(output)
    if not (np.array_equal(restored, vertices) and np.array_equal(restored_faces, faces)
            and np.array_equal(restored_cells, cells)):
        raise ValueError("exact chamber did not survive the serialization round trip")
    for path, digest in sources.items():
        if common.digest(path) != digest:
            raise ValueError("source changed during exact chamber export")
    record = dict(sources=sources, output_sha256=common.digest(output),
                  vertices=len(vertices), tetrahedra=len(cells),
                  geometry_modified=False, connectivity_modified=False,
                  original_double_coordinates=True, face_labels_preserved=True,
                  method="exact HDF5 coordinates and identical ordered Medit skeleton")
    common.save_json(certificate, record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "skeleton", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    export(args.snapshot, args.skeleton, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
