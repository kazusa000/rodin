"""CPU-side failure gates always run; GPU arithmetic is explicitly opt-in."""
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
from kelvin_bem import bem_kcq_score as bem
from kelvin_bem import bem_cuda_operator as cuda


def panels():
    return [bem.Panel((1., .11, .24), .10, (1., 0., 0.)),
            bem.Panel((.8, .40, .19), .08, (1., 0., 0.)),
            bem.Panel((.9, .12, -.22), .12, (1., 0., 0.)),
            bem.Panel((1.1, -.18, .33), .06, (1., 0., 0.))]


def group():
    return bem.proper_tetrahedral_rotations() if bem.symmetry.IS_TERRA else bem.proper_octahedral_rotations()


def score_options():
    return dict(mesh_label="cuda-regression", max_chamber_panels=0, epsilon_factor=.6,
                viscosity=1., center=(0., 0., 0.))


def test_cuda_requires_explicit_backend_and_library():
    with pytest.raises(ValueError, match="only available with quotient-cg"):
        bem.score_panels(panels(), matrix_free_device="cuda", **score_options())
    with pytest.raises(ValueError, match="explicit --matrix-free-device"):
        bem.score_panels(panels(), cuda_library="not-used.so", **score_options())
    with pytest.raises(ValueError, match="explicit --cuda-library"):
        bem.score_panels(panels(), backend="quotient-cg", matrix_free_device="cuda", **score_options())


def test_missing_cuda_library_is_not_a_cpu_fallback(tmp_path):
    with pytest.raises(FileNotFoundError):
        bem.score_panels(panels(), backend="quotient-cg", matrix_free_device="cuda",
                         cuda_library=str(tmp_path / "absent.so"), **score_options())


def test_independent_cpu_residual_rejects_a_self_consistent_wrong_gpu_operator(monkeypatch):
    closed = []
    def fake(cpu, **kwargs):
        class WrongOperator:
            n = cpu.n
            row_chunk = cpu.row_chunk
            applications = 0
            cuda = SimpleNamespace(metadata={"fake": True}, close=lambda: closed.append(True))
            def make_preconditioner(self, panel_block_size):
                return cpu.make_preconditioner(panel_block_size)
            def matmat(self, values):
                self.applications += 1
                return 2. * cpu.matmat(values)
        return WrongOperator()
    monkeypatch.setattr(cuda, "make_cuda_operator", fake)
    with pytest.raises(RuntimeError, match="independent CPU true residual"):
        bem.score_panels(panels(), backend="quotient-cg", matrix_free_device="cuda",
                         cuda_library="fake", preconditioner_panel_block_size=4, **score_options())
    assert closed == [True]


def test_cpu_default_does_not_construct_cuda_operator(monkeypatch):
    monkeypatch.setattr(cuda, "make_cuda_operator", lambda *a, **k: pytest.fail("unexpected CUDA selection"))
    result = bem.score_panels(panels(), backend="quotient-cg", preconditioner_panel_block_size=4, **score_options())
    assert result["matrix_free_device"] == "cpu"


native = pytest.mark.skipif(os.environ.get("KELVIN_RUN_CUDA") != "1", reason="explicit CUDA hardware test")


@native
@pytest.mark.parametrize("rhs_count", [1, 2])
@pytest.mark.parametrize("panel_count", [4, 33, 129])
def test_cuda_matmat_matches_dense_and_cpu_for_random_signed_fields(rhs_count, panel_count):
    rng = np.random.default_rng(919)
    shape = [bem.Panel(tuple(point), float(area), (1., 0., 0.)) for point, area in
             zip(rng.uniform((.6, -.3, -.2), (1.2, .4, .3), (panel_count, 3)),
                 rng.uniform(.01, .1, panel_count))]
    op = bem.QuotientStokesletOperator(shape, group(), epsilon=.17, viscosity=1.3)
    gpu = cuda.make_cuda_operator(op, library=os.environ["BEM_CUDA_LIBRARY"])
    values = np.random.default_rng(920).normal(size=(9*panel_count, rhs_count)).astype(np.float64)
    matrix = bem.assemble_quotient_stokeslet_matrix(shape, group(), epsilon=.17, viscosity=1.3)
    try:
        measured = gpu.matmat(values)
        np.testing.assert_allclose(measured, matrix @ values, rtol=3e-12, atol=3e-14)
        np.testing.assert_allclose(measured, op.matmat(values), rtol=3e-12, atol=3e-14)
        assert gpu.cuda.metadata["dtype"] == "float64"
        assert gpu.cuda.metadata["fast_math"] is False
        with pytest.raises(ValueError, match="float64"):
            gpu.matmat(values.astype(np.float32))
    finally:
        gpu.cuda.close()


@native
@pytest.mark.parametrize("center", [(0., 0., 0.), (.1, .2, -.1)])
def test_cuda_pcg_resistance_matches_direct_with_both_residual_certificates(center):
    options = dict(score_options(), center=center)
    direct = bem.score_panels(panels(), backend="quotient", **options)
    result = bem.score_panels(panels(), backend="quotient-cg", matrix_free_device="cuda",
        cuda_library=os.environ["BEM_CUDA_LIBRARY"], preconditioner_panel_block_size=2,
        solver_rtol=1e-12, **options)
    for key in ("K", "C", "Q", "force_due_rotation", "rho"):
        np.testing.assert_allclose(result[key], direct[key], rtol=2e-8, atol=2e-10)
    assert result["cpu_residual_verified"]
    assert max(result["cpu_rhs_relative_residuals"]) <= 1e-12
    assert max(result["gpu_rhs_relative_residuals"]) <= 1e-12
    assert result["cuda"]["source_sha256"]


@native
def test_default_tolerance_matches_cpu_pcg_without_changing_acceptance():
    options = dict(score_options(), backend="quotient-cg", preconditioner_panel_block_size=2)
    cpu = bem.score_panels(panels(), **options)
    gpu = bem.score_panels(panels(), matrix_free_device="cuda",
                          cuda_library=os.environ["BEM_CUDA_LIBRARY"], **options)
    assert gpu["solver_rtol"] == cpu["solver_rtol"] == 1e-10
    for key in ("K", "C", "Q", "rho"):
        np.testing.assert_allclose(gpu[key], cpu[key], rtol=2e-10, atol=2e-12)
    assert max(gpu["cpu_rhs_relative_residuals"]) <= 1e-10
    assert max(gpu["gpu_rhs_relative_residuals"]) <= 1e-10


@native
def test_cuda_scale_mirror_and_nonconvergence():
    options = dict(score_options(), backend="quotient-cg", matrix_free_device="cuda",
        cuda_library=os.environ["BEM_CUDA_LIBRARY"], preconditioner_panel_block_size=4)
    base = bem.score_panels(panels(), **options)
    scaled = [bem.Panel(tuple(1.7*np.array(p.centroid)), 1.7**2*p.area, p.normal) for p in panels()]
    mirror = [bem.Panel(tuple(-np.array(p.centroid)), p.area, tuple(-np.array(p.normal))) for p in panels()]
    for shape, factor, sign in ((scaled, 1.7, 1), (mirror, 1., -1)):
        result = bem.score_panels(shape, **options)
        np.testing.assert_allclose(result["rho"], base["rho"], rtol=2e-8, atol=2e-10)
        np.testing.assert_allclose(result["c"], sign*factor**2*base["c"], rtol=2e-8, atol=2e-10)
    options["preconditioner_panel_block_size"] = 1
    with pytest.raises(RuntimeError, match="did not converge"):
        bem.score_panels(panels(), solver_maxiter=1, solver_rtol=1e-14, **options)
