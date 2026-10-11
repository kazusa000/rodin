#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Build the optional FP64-only operator outside the source checkout."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from .run_support import ROOT, EXAMPLE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nvcc", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error("build output must be outside the source checkout")
    if output.exists():
        parser.error("output already exists; use a new fixed-version build directory")
    source = EXAMPLE / "cuda/bem_quotient_cuda.cu"
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    pending = output.with_name(output.name + ".pending")
    command = [str(args.nvcc.resolve()), "-shared", "-Xcompiler=-fPIC", "-O3", "-std=c++17",
               "--gpu-architecture=compute_75", "--gpu-code=sm_75,compute_75", "--fmad=false",
               f'-DBEM_CUDA_SOURCE_SHA256="{digest}"', str(source), "-o", str(pending)]
    subprocess.run(command, check=True)
    pending.replace(output)
    metadata = {"command": command, "source_sha256": digest,
                "library_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
                "nvcc": subprocess.check_output([str(args.nvcc), "--version"], text=True)}
    output.with_suffix(output.suffix + ".build.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
