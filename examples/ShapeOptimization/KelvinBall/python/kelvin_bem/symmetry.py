# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Process-wide symmetry selection shared by all Kelvin Python stages.

Set KELVIN_SYMMETRY before starting Python; never switch an imported numerical
pipeline in place. Independent cubic and terra processes may run concurrently.
"""

import os

NAME = os.environ.get("KELVIN_SYMMETRY", "cubic")
if NAME not in {"cubic", "terra"}:
    raise ValueError("KELVIN_SYMMETRY must be 'cubic' or 'terra'")
IS_TERRA = NAME == "terra"
ORDER = 12 if IS_TERRA else 24
FULL_PREFIX = f"full{ORDER}"


def rotations():
    # Lazy import avoids a cycle with the geometry implementation.
    from . import periodic_geometry as geometry

    return (geometry.proper_tetrahedral_rotations() if IS_TERRA
            else geometry.proper_octahedral_rotations())
