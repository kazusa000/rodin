#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Unit checks for the report2 BEM acceptance metric."""

from __future__ import annotations

import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "bem_solver.py"


def load_module():
    from kelvin_bem import bem_solver as module
    return module


def test_proper_octahedral_rotations_are_24() -> None:
    module = load_module()
    rotations = module.proper_octahedral_rotations()
    assert len(rotations) == 24
    assert len({tuple(tuple(row) for row in rotation) for rotation in rotations}) == 24
    assert all(round(module.det3(rotation)) == 1 for rotation in rotations)


def test_regularized_stokeslet_block_is_symmetric_positive() -> None:
    module = load_module()
    block = module.regularized_stokeslet_block(
        (0.0, 0.0, 0.0),
        (0.3, -0.2, 0.5),
        epsilon=0.1,
        viscosity=1.0,
    )
    for i in range(3):
        assert block[i][i] > 0.0
        for j in range(3):
            assert abs(block[i][j] - block[j][i]) <= 1e-14


def test_panel_downsampling_is_geometry_order_independent() -> None:
    module = load_module()
    panels = [
        module.Panel(
            centroid=(float(i % 5), float(i // 5), 0.1 * float(i)),
            area=1.0 + 0.03 * i,
            normal=(1.0, 0.0, 0.0),
        )
        for i in range(15)
    ]
    permutation = [7, 2, 14, 0, 11, 5, 9, 3, 13, 1, 10, 4, 12, 6, 8]
    forward = module.downsample_panels(panels, target=6)
    reordered = module.downsample_panels([panels[i] for i in permutation], target=6)
    assert sorted(panel.centroid for panel in forward) == sorted(
        panel.centroid for panel in reordered
    )


def test_summary_writer_emits_kcq_and_rho() -> None:
    module = load_module()
    result = {
        "k": 1.0,
        "c": -0.25,
        "q": 2.0,
        "rho": 0.125,
        "K": [[1.0, 0.0, 0.0], [0.0, 1.1, 0.0], [0.0, 0.0, 0.9]],
        "C": [[-0.2, 0.0, 0.0], [0.0, -0.3, 0.0], [0.0, 0.0, -0.25]],
        "Q": [[2.0, 0.0, 0.0], [0.0, 2.1, 0.0], [0.0, 0.0, 1.9]],
    }
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "bem.summary"
        module.write_summary(path, result)
        text = path.read_text(encoding="ascii")
    for token in ("k=1", "c=-0.25", "q=2", "rho=0.125"):
        assert token in text


if __name__ == "__main__":
    test_proper_octahedral_rotations_are_24()
    test_regularized_stokeslet_block_is_symmetric_positive()
    test_panel_downsampling_is_geometry_order_independent()
    test_summary_writer_emits_kcq_and_rho()
    print("test_bem_solver: ok")
