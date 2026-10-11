# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Migration contracts: one checkout and unchanged numerical definitions."""
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from kelvin_bem import run_support as common
from kelvin_bem import cuda_operator as cuda

CLEANUP = json.loads((common.EXAMPLE / "reproduction/cleanup-provenance.json").read_text())


class StripImports(ast.NodeTransformer):
    def visit_Import(self, node):
        return None

    def visit_ImportFrom(self, node):
        return None

    def visit_Expr(self, node):
        if isinstance(node.value, ast.Call) and ast.unparse(node.value.func).startswith("sys.path."):
            return None
        return self.generic_visit(node)

    def visit_Constant(self, node):
        # Only the bundled-module filename guard changed; no numerical literal
        # or arbitrary string is normalized away.
        for old, new in CLEANUP["module_renames"].items():
            if node.value == new + ".py":
                node.value = old + ".py"
        return node


def test_retained_original_numerical_definitions_are_unchanged():
    data = json.loads((common.EXAMPLE / "reproduction/migration-provenance.json").read_text())
    checked = 0
    retired = 0
    for row in data["definitions"]:
        original = Path(row["destination"])
        old_name = original.stem
        renamed = original.with_name(CLEANUP["module_renames"].get(old_name, old_name) + ".py")
        tree = ast.parse((common.EXAMPLE / renamed).read_text())
        nodes = {n.name: n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
        for name, expected in row["unchanged_definitions"].items():
            if name in CLEANUP["removed_definitions"].get(old_name, []):
                assert name not in nodes
                retired += 1
                continue
            normalized = ast.dump(StripImports().visit(nodes[name]), include_attributes=False)
            assert hashlib.sha256(normalized.encode()).hexdigest() == expected
            checked += 1
    assert checked == CLEANUP["retained_definitions"] == 99
    assert checked + retired == CLEANUP["original_definitions"] == 109


def test_native_cuda_preset_and_frozen_migration_evidence_are_unchanged():
    example_relative = common.EXAMPLE.relative_to(common.ROOT)
    paths = list(common.EXAMPLE.glob("*.cpp")) + list(common.EXAMPLE.glob("*.h"))
    paths += [common.EXAMPLE / name for name in (
        "CMakeLists.txt", "cuda/bem_quotient_cuda.cu",
        "input/chamber_1_24_sphere_a8.mesh", "reproduction/clement-swift-bem.json",
        "reproduction/migration-provenance.json", "reproduction/migration-validation.json")]
    for path in paths:
        original = subprocess.check_output(["git", "-C", str(common.ROOT), "show",
            CLEANUP["baseline_commit"] + ":" + str(example_relative / path.relative_to(common.EXAMPLE))])
        assert path.read_bytes() == original, path


def test_source_roots_and_cuda_are_bundled():
    assert common.ROOT / "examples/ShapeOptimization/KelvinBall" == common.EXAMPLE
    reference = json.loads((common.EXAMPLE / "reproduction/clement-swift-bem.json").read_text())
    assert cuda.SOURCE == common.EXAMPLE / "cuda/bem_quotient_cuda.cu"
    assert common.digest(cuda.SOURCE) == reference["cuda_source_sha256"]
    for path in (common.ROOT / "out/results", common.EXAMPLE / "results"):
        with pytest.raises(ValueError, match="inside source"):
            common.external(path)


def test_geometry_contract_runs_on_the_versioned_sphere_mesh():
    from kelvin_bem.acceptance import geometry_diagnostics
    mesh = common.EXAMPLE / "input/chamber_1_24_sphere_a8.mesh"
    result = geometry_diagnostics(mesh, "cubic")
    assert result["mesh_sha256"] == common.digest(mesh)
    assert result["volume_full"] > 0 and result["chamber_panels"] > 0
    assert result["full24"]["triangles"] == 24 * result["chamber_panels"]
    assert result["solid_connectivity"]["group_order"] == 24


@pytest.mark.parametrize("command", [
    ("native",), ("bem",), ("compare",), ("build-cuda",), ("baseline",), ("viewer",),
    ("diagnostics", "audit"), ("diagnostics", "geometry"),
    ("diagnostics", "export-chamber"), ("diagnostics", "export-target"),
    ("diagnostics", "export-interface"), ("diagnostics", "score"),
    ("diagnostics", "scan"), ("diagnostics", "plot"), ("diagnostics", "render"),
    ("paired-check",), ("geometry-check",)])
def test_entry_from_an_unrelated_directory_needs_no_03_checkout(tmp_path, command):
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    result = subprocess.run([sys.executable, str(common.EXAMPLE / "scripts/kelvin-bem.py"),
                             *command, "--help"], cwd=tmp_path, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "--verifier" not in result.stdout


def test_python_package_has_no_absolute_personal_paths_or_old_tree_imports():
    package = common.EXAMPLE / "python/kelvin_bem"
    for path in package.glob("*.py"):
        text = path.read_text()
        assert "/Users/wjj/" not in text and "/home/wjj/" not in text
        assert ' / "src"' not in text and ' / "scripts"))' not in text


def test_nodes_config_cannot_override_numerical_settings(tmp_path, monkeypatch):
    from kelvin_bem import compare_backends as pair
    config = tmp_path / "nodes.json"
    config.write_text(json.dumps({"solver_rtol": 0.1}))
    monkeypatch.setattr(sys, "argv", ["pair", "--output", str(tmp_path / "out"),
                                     "--nodes-config", str(config)])
    with pytest.raises(SystemExit) as caught:
        pair.main()
    assert caught.value.code == 2


def test_explicit_deployment_flag_wins_over_nodes_config(tmp_path, monkeypatch):
    from kelvin_bem import compare_backends as pair
    config = tmp_path / "nodes.json"
    values = {key: "/configured" for key in ("cpu_source", "gpu_source", "cpu_python",
              "gpu_python", "build", "cpu_output", "gpu_output", "cuda_library")}
    values.update(cpu_tailscale_ip="100.1.2.3", gpu_tailscale_ip="100.4.5.6")
    config.write_text(json.dumps(values))
    monkeypatch.setattr(sys, "argv", ["pair", "--output", str(tmp_path / "out"),
               "--nodes-config", str(config), "--build=/explicit"])
    monkeypatch.setattr(pair.common, "source_state", lambda *a: {"commit": "same-sha"})
    def controller(args):
        assert args.build == "/explicit"
        assert args.cpu_tailscale_ip == "100.1.2.3" and args.gpu_tailscale_ip == "100.4.5.6"
        assert args.worker_commit == "same-sha"
        from types import SimpleNamespace
        return SimpleNamespace(run=lambda: 0)
    monkeypatch.setattr(pair, "Controller", controller)
    assert pair.main() == 0
