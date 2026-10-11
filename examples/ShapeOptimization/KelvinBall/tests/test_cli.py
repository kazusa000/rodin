# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Command routing keeps workflow, setup and diagnostics separate."""
import sys
from types import SimpleNamespace

import pytest
from kelvin_bem import __main__ as cli


@pytest.mark.parametrize("arguments,module", [
    (["native"], "native_worker"), (["bem"], "bem_worker"),
    (["compare"], "compare_backends"), (["viewer"], "history_viewer"),
    (["build-cuda"], "build_cuda"), (["baseline"], "native_baseline"),
    *[(["diagnostics", name], module) for name, module in cli.DIAGNOSTIC_COMMANDS.items()],
    (["paired-check"], "compare_backends"), (["geometry-check"], "geometry_check"),
    (["audit"], "reconstruction_audit"),
])
def test_dispatch_preserves_worker_arguments_and_exit_status(monkeypatch, arguments, module):
    monkeypatch.setattr(sys, "argv", ["kelvin-bem.py", *arguments, "--output", "/outside/source"])
    def load(name, package):
        assert name == "." + module and package == "kelvin_bem"
        assert sys.argv == ["kelvin-bem.py", "--output", "/outside/source"]
        return SimpleNamespace(main=lambda: 7)
    monkeypatch.setattr(cli.importlib, "import_module", load)
    assert cli.main() == 7


@pytest.mark.parametrize("arguments", [[], ["--help"], ["diagnostics"], ["diagnostics", "--help"]])
def test_help_never_starts_a_worker(monkeypatch, capsys, arguments):
    monkeypatch.setattr(sys, "argv", ["kelvin-bem.py", *arguments])
    monkeypatch.setattr(cli.importlib, "import_module",
        lambda *a, **k: pytest.fail("help must not start a numerical stage"))
    assert cli.main() == 0
    output = capsys.readouterr().out
    assert "Usage:" in output
    if "diagnostics" not in arguments:
        assert "not a long optimization loop" in output
        assert "paired-check" not in output


@pytest.mark.parametrize("arguments", [["missing"], ["diagnostics", "missing"]])
def test_unknown_command_fails_without_importing_a_worker(monkeypatch, capsys, arguments):
    monkeypatch.setattr(sys, "argv", ["kelvin-bem.py", *arguments])
    monkeypatch.setattr(cli.importlib, "import_module",
        lambda *a, **k: pytest.fail("invalid command must not start a stage"))
    assert cli.main() == 2
    assert "Unknown" in capsys.readouterr().err
