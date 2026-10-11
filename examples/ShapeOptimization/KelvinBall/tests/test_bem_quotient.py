"""Exact discrete symmetry reduction, with physical and scaling checks."""
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.

from pathlib import Path
import sys

import numpy as np
import pytest

from kelvin_bem import bem_kcq_score as bem


def panels():
    # Unequal areas and a chiral point arrangement exercise the coupling blocks.
    return [
        bem.Panel((1.0, 0.11, 0.24), 0.10, (1.0, 0.0, 0.0)),
        bem.Panel((0.8, 0.40, 0.19), 0.08, (1.0, 0.0, 0.0)),
        bem.Panel((0.9, 0.12, -0.22), 0.12, (1.0, 0.0, 0.0)),
        bem.Panel((1.1, -0.18, 0.33), 0.06, (1.0, 0.0, 0.0)),
    ]


def rotations():
    return (bem.proper_tetrahedral_rotations() if bem.symmetry.IS_TERRA
            else bem.proper_octahedral_rotations())


def score(chamber, backend="quotient", center=(0.0, 0.0, 0.0), viscosity=1.0):
    return bem.score_panels(
        chamber, mesh_label="test", max_chamber_panels=0,
        epsilon_factor=0.6, viscosity=viscosity, center=center, backend=backend,
    )


@pytest.mark.parametrize("center", [(0.0, 0.0, 0.0), (0.1, 0.2, -0.1)])
def test_full_and_quotient_six_loads_agree(center):
    full = score(panels(), "full", center)
    quotient = score(panels(), "quotient", center)
    for key in ("K", "C", "Q", "force_due_rotation", "k", "c", "q", "rho"):
        np.testing.assert_allclose(quotient[key], full[key], rtol=2e-11, atol=2e-13)
    assert quotient["epsilon"] == full["epsilon"]
    assert quotient["surface_area"] == full["surface_area"]
    assert quotient["full_panels"] == full["full_panels"]
    assert quotient["linear_relative_residual"] < 1e-12
    assert quotient["system_dimension"] == 9 * len(panels())


def test_quotient_is_exact_operator_restriction_not_only_a_scalar_match():
    chamber, group = panels(), rotations()
    full_panels = bem._expand_panels(chamber, group)
    epsilon = 0.17
    full_matrix = bem.assemble_regularized_stokeslet_matrix(
        full_panels, epsilon=epsilon, viscosity=1.3,
    )
    quotient = bem.assemble_quotient_stokeslet_matrix(
        chamber, group, epsilon=epsilon, viscosity=1.3, row_chunk=2,
    )
    np.testing.assert_allclose(quotient, quotient.T, rtol=2e-14, atol=2e-15)
    response = np.random.default_rng(429).normal(size=(len(chamber), 3, 3))
    expanded = np.concatenate([np.asarray(g) @ response @ np.asarray(g).T for g in group])
    action = (full_matrix @ expanded.reshape(-1, 3)).reshape(-1, 3, 3)
    identity_index = group.index(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    base_action = action[identity_index * len(chamber):(identity_index + 1) * len(chamber)]
    sqrt_area = np.sqrt([panel.area for panel in chamber])[:, None, None]
    reduced_action = quotient @ (sqrt_area * response).ravel()
    np.testing.assert_allclose(reduced_action, (sqrt_area * base_action).ravel(),
                               rtol=2e-13, atol=2e-14)


def test_scale_and_viscosity_laws_and_handedness():
    original = score(panels())
    length_scale, viscosity_scale = 2.3, 1.7
    scaled = [bem.Panel(tuple(length_scale * np.asarray(p.centroid)),
                        p.area * length_scale**2, p.normal) for p in panels()]
    result = score(scaled, viscosity=viscosity_scale)
    for key, power in (("k", 1), ("c", 2), ("q", 3)):
        np.testing.assert_allclose(result[key], original[key] * viscosity_scale * length_scale**power,
                                   rtol=2e-11, atol=1e-12)
    np.testing.assert_allclose(result["rho"], original["rho"], rtol=2e-11, atol=1e-14)
    inverted = [bem.Panel(tuple(-np.asarray(p.centroid)), p.area,
                          tuple(-np.asarray(p.normal))) for p in panels()]
    mirror = score(inverted)
    np.testing.assert_allclose(mirror["c"], -original["c"], rtol=2e-11, atol=1e-13)
    np.testing.assert_allclose(mirror["rho"], original["rho"], rtol=2e-11, atol=1e-14)


def test_sphere_approaches_analytic_resistance():
    sphere = bem.read_gamma_panels(
        Path(__file__).resolve().parents[1] / "input/chamber_1_24_sphere_a8.mesh",
    )
    errors = []
    for count in (24, 128):
        chamber = bem.downsample_panels(sphere, count, rescale_area=True)
        if bem.symmetry.IS_TERRA:
            # T is an index-two subgroup of O; two O chambers form a T chamber.
            outside = next(g for g in bem.proper_octahedral_rotations()
                           if g not in bem.proper_tetrahedral_rotations())
            chamber = chamber + bem._expand_panels(chamber, [outside])
        result = score(chamber)
        error = (abs(result["k"] / (6 * np.pi) - 1),
                 abs(result["q"] / (8 * np.pi) - 1))
        errors.append(error)
        assert result["rho"] < 1e-12
    assert errors[1][0] < errors[0][0]
    assert errors[1][1] < errors[0][1]
    assert errors[1][0] < 0.01
    assert errors[1][1] < 0.025


def test_default_backend_remains_full():
    result = bem.score_panels(panels(), mesh_label="test", max_chamber_panels=0,
                             epsilon_factor=0.6, viscosity=1.0, center=(0.0, 0.0, 0.0))
    assert result["backend"] == "full"


@pytest.mark.parametrize("epsilon,viscosity", [(0.0, 1.0), (float("nan"), 1.0), (0.6, -1.0)])
def test_invalid_physical_parameters_are_rejected(epsilon, viscosity):
    with pytest.raises(ValueError):
        bem.score_panels(panels(), mesh_label="test", max_chamber_panels=0,
                         epsilon_factor=epsilon, viscosity=viscosity, center=(0.0, 0.0, 0.0))
