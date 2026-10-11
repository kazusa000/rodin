#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Run the bundled package from any working directory; no 03 checkout needed."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from kelvin_bem.__main__ import main

if __name__ == "__main__":
    raise SystemExit(main())
