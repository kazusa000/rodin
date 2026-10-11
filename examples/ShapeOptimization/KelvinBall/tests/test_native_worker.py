import argparse
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
from pathlib import Path
from kelvin_bem import native_worker as native


def test_same_mesh_and_coefficients_for_both_backends():
    args = argparse.Namespace(stage="step", mesh=Path("/input.mesh"),
        outer_values=Path("/outer.txt"), bem_values=Path("/bem.txt"),
        step=0.015625, method="mmg")
    reference = dict(thickness_absolute=0.06, thickness_weight=1)
    mmg = native.command(args, Path("/KelvinBall"), reference)
    args.method = "swift"
    swift = native.command(args, Path("/KelvinBall"), reference)
    assert [a for a in mmg if not a.startswith("--reconstruction=")] == [
        a for a in swift if not a.startswith("--reconstruction=")]
    assert "--boundary-step" in swift and "--iterations=1" in swift
    assert "--thickness-min=0" in swift


def test_common_fixed_seed_not_selected():
    args = argparse.Namespace(stage="initial")
    command = native.command(args, Path("/KelvinBall"), {})
    assert "--seed-amplitude=0.02" in command
    assert "--reconstruction=mmg" in command
    assert "--geometry-only" in command


def test_reconstruction_does_not_recompute_target():
    args = argparse.Namespace(stage="reconstruct", mesh=Path("/same.mesh"),
        target=Path("/same.sol"), method="swift")
    command = native.command(args, Path("/KelvinBall"), {})
    assert "--target=/same.sol" in command and "--mesh=/same.mesh" in command
    assert "--boundary-step" not in command
