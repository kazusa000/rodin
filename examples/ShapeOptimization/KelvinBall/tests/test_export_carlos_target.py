from pathlib import Path
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
import sys
import h5py
import numpy as np
import pytest
from kelvin_bem import export_carlos_target as target
from test_export_carlos_exact_chamber import fixture


def inputs(tmp_path):
    (x, cells, labels), (rounded, faces, _) = fixture()
    snapshot, field, skeleton = [tmp_path/name for name in ("snapshot.h5","field.h5","input.mesh")]
    with h5py.File(snapshot, "w") as f:
        f["/Mesh/Geometry/Vertices"] = x
        f["/Mesh/XDMF/Topology"] = np.column_stack((np.full(len(cells),6), cells)).ravel()
        f["/Mesh/Attributes/3"] = labels
    values = np.array([.12345678901234567, -1e-12, 0, 3.2, -5.1])
    with h5py.File(field, "w") as f:
        f["/GridFunction/Values/Data"] = values
        f["/GridFunction/Meta/Dimension"] = [1]
    target.interface.mesh_io.write_medit_mesh(skeleton, rounded.tolist(), faces.tolist(),
                                             np.column_stack((cells+1, labels)).tolist())
    return snapshot, skeleton, field, values


def test_preserves_saved_target_and_exact_geometry(tmp_path):
    snapshot, skeleton, field, values = inputs(tmp_path)
    result = target.export(snapshot, skeleton, field, tmp_path/"export")
    text = (tmp_path/"export/target.sol").read_text().splitlines()
    np.testing.assert_array_equal([float(v) for v in text[5:-1]], values)
    assert not result["target_recomputed"] and not result["geometry_modified"]


def test_wrong_nodal_dimension_rejected(tmp_path):
    snapshot, skeleton, field, _ = inputs(tmp_path)
    with h5py.File(field, "a") as f:
        f["/GridFunction/Meta/Dimension"][0] = 3
    with pytest.raises(ValueError, match="scalar"):
        target.export(snapshot, skeleton, field, tmp_path/"export")
