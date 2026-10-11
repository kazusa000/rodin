from pathlib import Path
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
import sys

import numpy as np
import pytest

from kelvin_bem import export_carlos_exact_chamber as exact


def fixture():
    x = np.array([[0, 0, 0], [1.123456789012345, 0, 0], [0, 1, 0],
                  [0, 0, -1], [0, 0, 1]], float)
    tets = np.array([[0, 1, 2, 3], [0, 1, 2, 4]])
    labels = np.array([2, 3])
    rounded = np.array([[float(format(v, ".6g")) for v in point] for point in x])
    return (x, tets, labels), (rounded, np.array([[1, 2, 3, 13]]),
                              np.column_stack((tets+1, labels)))


def test_restores_exact_coordinates_not_rounded_values():
    snapshot, skeleton = fixture()
    vertices, faces, cells = exact.match(snapshot, skeleton)
    np.testing.assert_array_equal(vertices, snapshot[0])
    assert vertices[1, 0] != skeleton[0][1, 0]
    np.testing.assert_array_equal(faces, skeleton[1])
    np.testing.assert_array_equal(cells, skeleton[2])


def test_wrong_skeleton_is_never_used_for_another_shape():
    snapshot, skeleton = fixture()
    skeleton[0][1, 0] += .01
    with pytest.raises(ValueError, match="rounded serialization"):
        exact.match(snapshot, skeleton)


def test_reordered_cells_are_rejected_before_label_transfer():
    snapshot, skeleton = fixture()
    with pytest.raises(ValueError, match="ordered tetrahedral"):
        exact.match(snapshot, (skeleton[0], skeleton[1], skeleton[2][::-1]))


def test_false_gamma_is_rejected():
    snapshot, skeleton = fixture()
    skeleton[1][0, 2] = 4
    with pytest.raises(ValueError, match="material interface"):
        exact.match(snapshot, skeleton)
