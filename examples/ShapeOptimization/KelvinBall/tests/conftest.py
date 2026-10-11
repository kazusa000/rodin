# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
import os
from pathlib import Path
import sys

os.environ.setdefault("KELVIN_SYMMETRY", "cubic")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
