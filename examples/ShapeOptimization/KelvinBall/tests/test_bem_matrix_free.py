"""Matrix-free action, PCG accuracy, resource behavior, and hard failures."""
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.

from pathlib import Path
import sys

import numpy as np
import pytest

from kelvin_bem import bem_kcq_score as bem
from kelvin_bem import bem_convergence_scan as scan


def panels():
    return [
        bem.Panel((1.0, 0.11, 0.24), 0.10, (1.0, 0.0, 0.0)),
        bem.Panel((0.8, 0.40, 0.19), 0.08, (1.0, 0.0, 0.0)),
        bem.Panel((0.9, 0.12, -0.22), 0.12, (1.0, 0.0, 0.0)),
        bem.Panel((1.1, -0.18, 0.33), 0.06, (1.0, 0.0, 0.0)),
    ]


def group():
    return (bem.proper_tetrahedral_rotations() if bem.symmetry.IS_TERRA
            else bem.proper_octahedral_rotations())


def operator():
    return bem.QuotientStokesletOperator(panels(), group(), epsilon=0.17,
                                        viscosity=1.3, row_chunk=2)


def test_matrix_free_action_equals_assembled_operator_for_arbitrary_fields():
    op = operator()
    matrix = bem.assemble_quotient_stokeslet_matrix(panels(), group(), epsilon=0.17, viscosity=1.3)
    values = np.random.default_rng(431).normal(size=(9 * len(panels()), 5))
    np.testing.assert_allclose(op.matmat(values), matrix @ values, rtol=2e-13, atol=2e-15)
    inverses = op.inverse_diagonal_blocks()
    for index, inverse in enumerate(inverses):
        diagonal = matrix[9 * index:9 * index + 9, 9 * index:9 * index + 9]
        np.testing.assert_allclose(inverse @ diagonal, np.eye(9), rtol=1e-13, atol=2e-15)


@pytest.mark.parametrize("block_size", [1, 2, 4])
def test_pcg_solution_passes_independent_dense_residual_check(block_size):
    op = operator()
    rhs = bem._quotient_rhs(panels())
    result, metrics = bem.solve_quotient_cg(op, rhs, rtol=1e-12, maxiter=200,
                                          panel_block_size=block_size)
    matrix = bem.assemble_quotient_stokeslet_matrix(panels(), group(), epsilon=0.17, viscosity=1.3)
    relative = np.linalg.norm(matrix @ result - rhs, axis=0) / np.linalg.norm(rhs, axis=0)
    assert np.all(relative <= 1.01e-12)
    np.testing.assert_allclose(result, np.linalg.solve(matrix, rhs), rtol=2e-9, atol=2e-10)
    assert metrics["solver_converged"] is True
    assert len(metrics["rhs_relative_residuals"]) == 2
    assert all(iterations > 0 for iterations in metrics["iterations_by_rhs"])


@pytest.mark.parametrize("center", [(0.0, 0.0, 0.0), (0.1, 0.2, -0.1)])
def test_matrix_free_matches_direct_six_load_resistance(center):
    kwargs = dict(mesh_label="test", max_chamber_panels=0, epsilon_factor=0.6,
                  viscosity=1.0, center=center)
    direct = bem.score_panels(panels(), backend="quotient", **kwargs)
    iterative = bem.score_panels(panels(), backend="quotient-cg", solver_rtol=1e-12,
                                 preconditioner_panel_block_size=1, **kwargs)
    for key in ("K", "C", "Q", "force_due_rotation", "rho"):
        np.testing.assert_allclose(iterative[key], direct[key], rtol=2e-8, atol=2e-10)
    assert iterative["matrix_bytes"] == 0
    assert iterative["equivalent_dense_matrix_bytes"] == direct["matrix_bytes"]
    assert iterative["solver"] == "pcg"
    assert iterative["solver_converged"] is True


def test_nonconvergence_is_a_hard_failure_not_a_score():
    with pytest.raises(RuntimeError, match="did not converge"):
        bem.solve_quotient_cg(operator(), bem._quotient_rhs(panels()),
                              rtol=1e-14, maxiter=1, panel_block_size=1)


def test_no_least_squares_fallback_for_nonpositive_curvature(monkeypatch):
    op = operator()
    monkeypatch.setattr(op, "matmat", lambda values: -values)
    with pytest.raises(RuntimeError, match="nonpositive/nonfinite curvature"):
        bem.solve_quotient_cg(op, bem._quotient_rhs(panels()), panel_block_size=1)


def test_zero_rhs_requires_no_iteration():
    solution, metrics = bem.solve_quotient_cg(operator(), np.zeros((36, 2)))
    assert np.count_nonzero(solution) == 0
    assert metrics["iterations_by_rhs"] == [0, 0]
    assert metrics["rhs_relative_residuals"] == [0.0, 0.0]


def test_scan_accepts_matrix_free_without_dense_allocation_and_counts_pcg(tmp_path):
    mesh = Path(__file__).resolve().parents[1] / "input/chamber_1_24_sphere_a8.mesh"
    rows = scan.run_scan(mesh, panels=[8], epsilon_factors=[0.6], viscosity=1.0,
                         center=(0.0, 0.0, 0.0), backend="quotient-cg", max_matrix_gib=1e-9,
                         solver_options={"solver_maxiter": 100, "preconditioner_panel_block_size": 8})
    assert rows[0]["status"] == "ok"
    assert rows[0]["matrix_bytes"] == 0
    assert len(rows[0]["rhs_relative_residuals"]) == 2
    summary = scan.convergence_summary(rows)
    assert summary["least_squares_fallbacks"] == 0
    assert summary["matrix_free_solves"] == 1
    scan.write_csv(tmp_path / "cg.csv", rows)
