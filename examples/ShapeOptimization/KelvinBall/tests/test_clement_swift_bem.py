import json
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from kelvin_bem import clement_swift_bem as worker


def test_worker_uses_retained_all_panel_speed_protocol(monkeypatch):
    monkeypatch.setenv("KELVIN_SYMMETRY", "cubic")
    reference = json.loads(worker.REFERENCE.read_text())
    args = worker.settings(SimpleNamespace(), reference)
    assert args.symmetry == "cubic"
    assert args.bem_rtol == 1e-6
    assert args.bem_weak_rtol == 1e-9
    assert args.bem_maxiter == 2000
    assert reference["panel_block_size"] == 256


def test_worker_never_switches_an_imported_tetra_mode(monkeypatch):
    monkeypatch.setenv("KELVIN_SYMMETRY", "terra")
    with pytest.raises(ValueError, match="fresh process"):
        worker.settings(SimpleNamespace(), json.loads(worker.REFERENCE.read_text()))
