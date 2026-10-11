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
from kelvin_bem import carlos_reference as common
from kelvin_bem import bem_cuda_operator as cuda


class StripImports(ast.NodeTransformer):
    def visit_Import(self, node):
        return None

    def visit_ImportFrom(self, node):
        return None

    def visit_Expr(self, node):
        if isinstance(node.value, ast.Call) and ast.unparse(node.value.func).startswith("sys.path."):
            return None
        return self.generic_visit(node)


def test_109_original_numerical_definitions_are_unchanged():
    data = json.loads((common.EXAMPLE / "reproduction/migration-provenance.json").read_text())
    checked = 0
    for row in data["definitions"]:
        tree = ast.parse((common.EXAMPLE / row["destination"]).read_text())
        nodes = {n.name: n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
        for name, expected in row["unchanged_definitions"].items():
            normalized = ast.dump(StripImports().visit(nodes[name]), include_attributes=False)
            assert hashlib.sha256(normalized.encode()).hexdigest() == expected
            checked += 1
    assert checked == 109


def test_source_roots_and_cuda_are_bundled():
    assert common.ROOT / "examples/ShapeOptimization/KelvinBall" == common.EXAMPLE
    reference = json.loads((common.EXAMPLE / "reproduction/clement-swift-bem.json").read_text())
    assert cuda.SOURCE == common.EXAMPLE / "cuda/bem_quotient_cuda.cu"
    assert common.digest(cuda.SOURCE) == reference["cuda_source_sha256"]
    for path in (common.ROOT / "out/results", common.EXAMPLE / "results"):
        with pytest.raises(ValueError, match="inside source"):
            common.external(path)


@pytest.mark.parametrize("command", ["native", "bem", "paired-check", "build-cuda", "viewer"])
def test_entry_from_an_unrelated_directory_needs_no_03_checkout(tmp_path, command):
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    result = subprocess.run([sys.executable, str(common.EXAMPLE / "scripts/kelvin-bem.py"),
                             command, "--help"], cwd=tmp_path, env=env,
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
    from kelvin_bem import clement_swift_pair as pair
    config = tmp_path / "nodes.json"
    config.write_text(json.dumps({"solver_rtol": 0.1}))
    monkeypatch.setattr(sys, "argv", ["pair", "--output", str(tmp_path / "out"),
                                     "--nodes-config", str(config)])
    with pytest.raises(SystemExit) as caught:
        pair.main()
    assert caught.value.code == 2


def test_explicit_deployment_flag_wins_over_nodes_config(tmp_path, monkeypatch):
    from kelvin_bem import clement_swift_pair as pair
    config = tmp_path / "nodes.json"
    values = {key: "/configured" for key in ("cpu_source", "gpu_source", "cpu_python",
              "gpu_python", "build", "cpu_output", "gpu_output", "cuda_library")}
    config.write_text(json.dumps(values))
    monkeypatch.setattr(sys, "argv", ["pair", "--output", str(tmp_path / "out"),
               "--nodes-config", str(config), "--build=/explicit"])
    monkeypatch.setattr(pair.common, "source_state", lambda *a: {"commit": "same-sha"})
    def controller(args):
        assert args.build == "/explicit"
        assert args.worker_commit == "same-sha"
        from types import SimpleNamespace
        return SimpleNamespace(run=lambda: 0)
    monkeypatch.setattr(pair, "Controller", controller)
    assert pair.main() == 0
