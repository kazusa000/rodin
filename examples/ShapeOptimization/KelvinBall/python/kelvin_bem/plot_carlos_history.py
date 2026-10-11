#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Plot measured Carlos reproduction history without mixing FEM and BEM values."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .carlos_reference import summary


def plot(history, output, bem_entries=()):
    summary(history)
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    with Path(history).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    field = lambda key: np.array([float(row[key]) for row in rows])
    iterations = field("iteration")
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    axes[0].semilogy(iterations, field("rho"), "o-", color="#185d88", ms=4,
                     label="Bounded FEM (Carlos reference)")
    axes[0].axhline(0.01, color="#777777", ls="--", lw=1, label="0.01 target")
    for index, (step, path) in enumerate(bem_entries):
        data = json.loads(Path(path).read_text())
        if isinstance(data, dict):
            data = [data]
        selected = [row for row in data if row.get("status", "ok") == "ok"
                    and row.get("requested_panels", row.get("panels")) == 0]
        if not selected:
            raise ValueError("plot BEM only when all source panels were used")
        y = [row["rho"] for row in selected]
        axes[0].semilogy([step] * len(y), y, "s", color="#c65f26", ms=5,
                         label="Free-space BEM diagnostic*" if index == 0 else None)
    axes[0].set(ylabel="Normalized coupling", title="Coupling (distinct measurement models)")
    axes[0].legend(fontsize=7, loc="best")
    volume_error = 100 * field("volume_error") / field("target_volume")
    axes[1].plot(iterations, volume_error, "o-", color="#185d88", ms=4)
    axes[1].axhline(0, color="#777777", lw=1)
    axes[1].set(ylabel="Volume drift (%)", title="Volume relative to initial sphere")
    axes[2].semilogy(iterations, field("mesh_quality_min"), "o-", color="#185d88", ms=4)
    axes[2].set(ylabel="Minimum tetrahedron quality", title="Mesh quality")
    for axis in axes:
        axis.set_xlabel("Design index (0 = initial sphere)")
        axis.grid(alpha=0.2)
    if bem_entries:
        fig.supxlabel("*Unmodified exported interface; consult paired-seam and BEM convergence diagnostics.", fontsize=8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bem", action="append", default=[], metavar="STEP=JSON")
    args = parser.parse_args()
    entries = [(int(item.split("=", 1)[0]), Path(item.split("=", 1)[1])) for item in args.bem]
    plot(args.history, args.output, entries)


if __name__ == "__main__":
    main()
