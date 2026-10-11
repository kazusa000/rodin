"""True residual checks must not discard healthy CG conjugacy."""
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
from pathlib import Path

import numpy as np

from kelvin_bem import bem_solver as bem


class DiagonalOperator:
    def __init__(self, n=40, corrupt_first=False):
        self.n = n
        self.diagonal = np.geomspace(1., 500., 9*n)
        self.applications = 0
        self.row_chunk = 32
        self.corrupt_first = corrupt_first

    def matmat(self, value):
        self.applications += 1
        result = self.diagonal[:, None] * value
        if self.corrupt_first and self.applications == 1:
            result = result + 0.01 * value
        return result

    def make_preconditioner(self, _):
        return lambda value: value.copy(), 0


def test_true_residual_monitor_keeps_conjugacy_past_fifty_iterations():
    op = DiagonalOperator()
    rhs = np.random.default_rng(218).normal(size=(9*op.n, 2))
    solution, result = bem.solve_quotient_cg(op, rhs, rtol=1e-11, maxiter=500)
    np.testing.assert_allclose(solution, rhs/op.diagonal[:, None], rtol=1e-8, atol=1e-10)
    checks = [row for row in result['residual_history'] if row['iteration'] % 50 == 0]
    assert len(checks) >= 2
    assert all(row['restarted_rhs'] == [] for row in checks)
    assert max(result['rhs_relative_residuals']) <= 1e-11


def test_material_recurrence_gap_restarts_without_false_success():
    op = DiagonalOperator(n=4, corrupt_first=True)
    rhs = np.random.default_rng(109).normal(size=(9*op.n, 2))
    solution, result = bem.solve_quotient_cg(op, rhs, rtol=1e-11, maxiter=500)
    assert any(row['restarted_rhs'] for row in result['residual_history'])
    residual = op.diagonal[:, None] * solution - rhs
    assert np.all(np.linalg.norm(residual, axis=0)/np.linalg.norm(rhs, axis=0) <= 1e-11)


def test_zero_rhs_alongside_active_rhs_keeps_exact_zero_solution():
    op = DiagonalOperator(n=8)
    rhs = np.column_stack((np.zeros(9*op.n), np.ones(9*op.n)))
    solution, result = bem.solve_quotient_cg(op, rhs, rtol=1e-11, maxiter=500)
    assert np.count_nonzero(solution[:, 0]) == 0
    assert result['iterations_by_rhs'][0] == 0
    assert result['rhs_relative_residuals'][0] == 0


def test_previously_converged_rhs_can_reactivate_without_false_success():
    class ReactivationProbe(DiagonalOperator):
        def matmat(self, values):
            result = super().matmat(values)
            # A deliberate one-call fault at the iteration-50 true check.
            # This probes robustness, not a changing physical SPD operator.
            if self.applications == 52:
                assert values[0, 0] == 1.0
                assert np.count_nonzero(values[:, 0]) == 1
                result[0, 0] += 1e-7
            return result

    op = ReactivationProbe()
    first = np.zeros(9*op.n)
    first[0] = 1.
    rhs = np.column_stack((first, np.random.default_rng(217).normal(size=9*op.n)))
    solution, result = bem.solve_quotient_cg(op, rhs, rtol=1e-11, maxiter=600)
    checks = {row['iteration']: row for row in result['residual_history']}
    assert checks[1]['true_relative_residuals'][0] == 0.
    assert checks[50]['true_relative_residuals'][0] > 1e-11
    assert 0 in checks[50]['restarted_rhs']
    assert result['iterations_by_rhs'][0] > 50
    residual = op.diagonal[:, None]*solution-rhs
    assert np.all(np.linalg.norm(residual, axis=0)/np.linalg.norm(rhs, axis=0) <= 1e-11)
