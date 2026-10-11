#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Measure one cubic mesh and export BEM coefficients and outer-boundary data."""

import argparse
import json
import os
from pathlib import Path
import platform

from . import run_support as common

REFERENCE = common.EXAMPLE / "reproduction/clement-swift-bem.json"


def settings(args, reference):
    if os.environ.get("KELVIN_SYMMETRY") != "cubic":
        raise ValueError("this worker requires KELVIN_SYMMETRY=cubic in a fresh process")
    args.symmetry = "cubic"
    args.bem_rtol = reference["solver_rtol"]
    args.bem_weak_rtol = reference["weak_coupling_rtol"]
    args.bem_maxiter = reference["solver_maxiter"]
    return args


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("mesh", "output", "cuda-library"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = settings(parser.parse_args(), json.loads(REFERENCE.read_text()))
    reference = json.loads(REFERENCE.read_text())
    source = common.source_state(common.ROOT)
    cuda_hash = common.digest(args.cuda_library)
    from .cuda_operator import SOURCE
    if common.digest(SOURCE) != reference["cuda_source_sha256"]:
        raise ValueError("CUDA source differs from the frozen all-panel BEM protocol")
    mesh_hash = common.digest(args.mesh)
    output = common.external(args.output)
    output.mkdir(parents=True, exist_ok=False)
    record = dict(source=source, host=platform.node(), started_at=common.utcnow(),
                  mesh=str(args.mesh.resolve()), mesh_sha256=mesh_hash,
                  cuda_sha256=cuda_hash, success=False, shape_updates=0)
    try:
        from . import bem_boundary as boundary
        result = boundary.bem_trace(args, reference, args.mesh, output)
        if (result["group_order"] != 24 or result["subsampled"]
                or result["full_panels"] != 24 * result["chamber_panels"]
                or result["mesh_sha256"] != mesh_hash):
            raise ValueError("incorrect panel or symmetry provenance")
        common.source_state(common.ROOT, source["commit"])
        if common.digest(args.mesh) != mesh_hash:
            raise ValueError("input mesh changed during BEM measurement")
        record.update(success=True, measurement=result["measurement"],
                      effective_rtol=result["solver_rtol"],
                      rho=result["rho"], bem_sha256=common.digest(output / "bem.json"))
    except BaseException as error:
        record["error"] = str(error)
        raise
    finally:
        record["finished_at"] = common.utcnow()
        common.save_json(output / "worker.json", record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
