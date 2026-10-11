#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Build and run Carlos's pinned cubic implementation without changing report2."""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[5]
EXAMPLE = Path(__file__).resolve().parents[2]


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def source_state(root, expected=None):
    root = Path(root).resolve()
    sha = git(root, "rev-parse", "HEAD")
    status = git(root, "status", "--porcelain", "--untracked-files=all")
    if status or (expected is not None and sha != expected):
        raise ValueError(f"source must be clean and pinned: {root}, SHA={sha}, status={status!r}")
    return {"path": str(root), "commit": sha, "status": status}


def external(path, *roots):
    path = Path(path).resolve()
    for root in (ROOT, *roots):
        if path.is_relative_to(Path(root).resolve()):
            raise ValueError(f"generated output cannot be inside source: {path}")
    return path


def save_json(path, value):
    path = Path(path)
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def execute(command, *, cwd, env, log, hours):
    """Bound this process group only; always return its observed exit status."""
    with Path(log).open("wb") as stream:
        proc = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, start_new_session=True)
        timed_out = False
        try:
            code = proc.wait(timeout=hours * 3600)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                code = proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                code = proc.wait()
        except BaseException:
            os.killpg(proc.pid, signal.SIGTERM)
            proc.wait()
            raise
    return code, timed_out


def summary(path, symmetry="cubic"):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("history has no completed evaluated designs")
    values = []
    for row in rows:
        converted = {}
        for key, value in row.items():
            if key in ("assembly_backend", "direct_solver"):
                converted[key] = value
            else:
                number = float(value)
                converted[key] = number if math.isfinite(number) else None
        values.append(converted)
    for row in values:
        if not all(row[key] is not None and math.isfinite(row[key])
                   for key in ("iteration", "k", "c", "q", "rho")):
            raise ValueError("history has non-finite resistance values")
        if row["k"] <= 0 or row["q"] <= 0:
            raise ValueError("history has non-positive resistance")
        expected = abs(row["c"]) / math.sqrt(row["k"] * row["q"])
        if not math.isclose(row["rho"], expected, rel_tol=1e-9, abs_tol=1e-14):
            raise ValueError("history rho is inconsistent with k,c,q")
    iterations = [int(row["iteration"]) for row in values]
    if iterations != list(range(len(rows))):
        raise ValueError("history must contain consecutive designs from zero")
    if symmetry not in ("cubic", "terra"):
        raise ValueError("unknown symmetry")
    return {"measurement": "bounded-domain FEM", "symmetry": symmetry, "group_order": 12 if symmetry == "terra" else 24,
            "evaluated_designs": len(rows), "first": values[0], "last": values[-1],
            "max_rho": max(row["rho"] for row in values), "bem_verified": False}
