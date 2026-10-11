import json
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from kelvin_bem import clement_swift_baseline as baseline


def test_baseline_is_unmodified_carlos_with_two_backends():
    reference = json.loads(baseline.REFERENCE.read_text())
    assert reference["baseline_commit"] == "79718dc7d018cf41c727dc2df2634628b42612c4"
    cases = baseline.case_commands(Path("/build"), reference)
    assert [name for name, _ in cases] == ["sphere", "mmg", "swift"]
    for _, command in cases[1:]:
        assert "--iterations=2" in command
        assert "--thickness-min=0" in command
        assert not any("outer-values" in word or "perturb" in word for word in command)


def test_build_directly_targets_carlos_source_without_transplant():
    args = SimpleNamespace(cmake="cmake", cmake_arg=["-DMUMPS_ROOT=/deps"])
    command = baseline.configure_command(args, Path("/rodin"), Path("/build"))
    assert command[:5] == ["cmake", "-S", "/rodin", "-B", "/build"]
    assert "-DRODIN_USE_MUMPS=ON" in command
    assert "-DRODIN_INSTALL_RESOURCES=OFF" in command
    assert "-DMUMPS_ROOT=/deps" in command


def audit_summary():
    first = dict(rho=1e-6, advection_increment=.01, stage_9_seconds=2,
                 mesh_quality_min=.2)
    last = dict(rho=3e-6, actual_delta_rho=2e-6, actual_delta_volume=.001,
                predicted_delta_rho=1e-6, volume_error=.001,
                mesh_quality_min=.1, direct_solver="MUMPS")
    return dict(evaluated_designs=2, first=first, last=last)


def test_successful_update_needs_transport_and_reconstruction():
    summary = audit_summary()
    result = baseline.update_audit(summary)
    assert result["evaluated_updates"] == 1
    assert result["bem_verified"] is False
    summary["first"]["advection_increment"] = 0
    with pytest.raises(ValueError, match="real level-set"):
        baseline.update_audit(summary)


def test_initial_state_does_not_count_as_real_update():
    summary = audit_summary()
    summary["evaluated_designs"] = 1
    with pytest.raises(ValueError, match="one evaluated"):
        baseline.update_audit(summary)


def test_inconsistent_measured_change_is_rejected():
    summary = audit_summary()
    summary["last"]["actual_delta_rho"] = .01
    with pytest.raises(ValueError, match="inconsistent"):
        baseline.update_audit(summary)


def test_compile_fix_only_allows_type_disambiguation(monkeypatch):
    reference = json.loads(baseline.REFERENCE.read_text())
    exact = ("-        decltype(classifier)::Parameters parameters;\n"
             "+        typename decltype(classifier)::Parameters parameters;")
    monkeypatch.setattr(baseline.common, "git", lambda *args: exact)
    assert baseline.audit_compile_fix(Path("/rodin"), reference)["numerical_changes"] is False
    monkeypatch.setattr(baseline.common, "git", lambda *args: exact + "\n+  step = .5;")
    with pytest.raises(ValueError, match="beyond"):
        baseline.audit_compile_fix(Path("/rodin"), reference)


def test_original_baseline_stays_the_default():
    reference = json.loads(baseline.REFERENCE.read_text())
    assert baseline.selected_commit(SimpleNamespace(compiler_fixed=False), reference) == reference["baseline_commit"]
    assert baseline.selected_commit(SimpleNamespace(compiler_fixed=True), reference) == reference["compiler_fixed_commit"]
