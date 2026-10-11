#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Export an existing linear XDMF nodal target without recomputing its values."""
import argparse
from pathlib import Path
import h5py
import numpy as np
from . import run_support as common
from . import export_chamber as chamber
from . import export_interface as interface


def export(snapshot, skeleton, field, output):
    snapshot, skeleton, field, output = map(lambda p: Path(p).resolve(),
                                           (snapshot, skeleton, field, output))
    output = common.external(output)
    if output.exists():
        raise FileExistsError(output)
    sources = {str(p): common.digest(p) for p in (snapshot, skeleton, field)}
    vertices, faces, cells = chamber.match(interface.load_snapshot(snapshot),
                                          interface.mesh_io.read_medit_mesh(skeleton))
    with h5py.File(field, "r") as data:
        values = np.asarray(data["/GridFunction/Values/Data"], dtype=float)
        dimension = np.asarray(data["/GridFunction/Meta/Dimension"]).item()
    if dimension != 1 or values.shape != (len(vertices),) or not np.isfinite(values).all():
        raise ValueError("expected one scalar XDMF value per linear visualization vertex")
    output.mkdir(parents=True)
    interface.mesh_io.write_medit_mesh(output/"source.mesh", vertices.tolist(),
                                       faces.tolist(), cells.tolist())
    text = ("MeshVersionFormatted 2\nDimension 3\nSolAtVertices\n"+str(len(values))+
            "\n1 1\n"+"\n".join(format(v, ".17g") for v in values)+"\nEnd\n")
    (output/"target.sol").write_text(text)
    restored = np.array([float(v) for v in text.splitlines()[5:-1]])
    if not np.array_equal(restored, values):
        raise ValueError("nodal target changed in serialization")
    for path, digest in sources.items():
        if common.digest(path) != digest:
            raise ValueError("original stopping-case evidence changed")
    record = dict(sources=sources, mesh_sha256=common.digest(output/"source.mesh"),
        target_sha256=common.digest(output/"target.sol"), vertices=len(values),
        geometry_modified=False, target_recomputed=False,
        nodal_order="linear XDMF writer evaluates each vertex by its mesh index")
    common.save_json(output/"certificate.json", record)
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "skeleton", "field", "output"):
        p.add_argument("--"+name, type=Path, required=True)
    args = p.parse_args()
    export(args.snapshot, args.skeleton, args.field, args.output)


if __name__ == "__main__":
    main()
