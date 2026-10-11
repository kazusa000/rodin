#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Build a full-resolution T12 or O24 viewer from retained Clément BEM history."""

import argparse
import base64
from collections import deque
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import struct
import sys

from .carlos_reference import ROOT, EXAMPLE


def collect(pilot_path, continuation_path, resume_path, reader_source):
    """Read original meshes on their owning host; do not write remote files."""
    from . import clement_boundary as boundary
    from . import export_carlos_interface as exporter

    campaigns = [json.loads((Path(p) / "campaign.json").read_text())
                 for p in (pilot_path, continuation_path, resume_path)]
    pilot, continuation, resume = campaigns
    if any(c.get("symmetry") != "terra" for c in campaigns):
        raise ValueError("all campaigns must use Terra T12")
    if (continuation["history"][-1]["bem"]["mesh_sha256"]
            != resume["history"][0]["bem"]["mesh_sha256"]):
        raise ValueError("resume seed differs from the retained checkpoint")
    if (pilot["final"]["mesh_sha256"]
            != continuation["history"][0]["bem"]["mesh_sha256"]):
        raise ValueError("continuation seed differs from the accepted pilot")
    for c in (continuation, resume):
        if (c["history"][-1]["update"] != c["accepted_updates"]
                or c["history"][-1]["bem"] != c["final"]
                or c["history"][-1]["mesh"] != c["final_mesh"]):
            raise ValueError("last accepted checkpoint is inconsistent")

    entries = [("初始扰动球", -1, Path(pilot_path) / "geometry/chamber.mesh",
                pilot["initial"], pilot["initial_geometry"])]
    for h in continuation["history"]:
        label = "单步试跑后 · 续跑基线" if h["update"] == 0 else f"续跑第 {h['update']} 步"
        entries.append((label, h["update"], Path(h["mesh"]), h["bem"], h["geometry"]))
    offset = continuation["accepted_updates"]
    for h in resume["history"][1:]:
        step = offset + h["update"]
        entries.append((f"续跑第 {step} 步", step, Path(h["mesh"]), h["bem"], h["geometry"]))
    return collect_frames(entries, campaigns, resume_path, boundary, exporter)


def chain_entries(campaigns, symmetry="terra"):
    """Join accepted histories only, removing the repeated continuation seeds."""
    entries, offset, previous = [], 0, None
    for campaign in campaigns:
        history = campaign["history"]
        if symmetry not in ("terra", "cubic") or campaign.get("symmetry") != symmetry:
            raise ValueError("all campaigns must use the selected symmetry")
        if (not history or [h["update"] for h in history] != list(range(len(history)))
                or history[-1]["update"] != campaign["accepted_updates"]
                or history[-1]["bem"] != campaign["final"]
                or history[-1]["mesh"] != campaign["final_mesh"]):
            raise ValueError("last accepted checkpoint is inconsistent")
        if previous is not None and history[0]["bem"]["mesh_sha256"] != previous:
            raise ValueError("continuation seed differs from the accepted checkpoint")
        for h in history[(1 if previous is not None else 0):]:
            step = offset + h["update"]
            label = "初始扰动球" if step == 0 else f"第 {step} 步"
            entries.append((label, step, Path(h["mesh"]), h["bem"], h["geometry"]))
        offset += campaign["accepted_updates"]
        previous = history[-1]["bem"]["mesh_sha256"]
    return entries


def collect_chain(paths, reader_source, symmetry="terra"):
    from . import clement_boundary as boundary
    from . import export_carlos_interface as exporter

    campaigns = [json.loads((Path(p) / "campaign.json").read_text()) for p in paths]
    return collect_frames(chain_entries(campaigns,symmetry), campaigns, paths[-1], boundary, exporter,symmetry)


def collect_frames(entries, campaigns, last_path, boundary, exporter, symmetry="terra"):
    if symmetry not in ("terra", "cubic"):
        raise ValueError("unknown viewer symmetry")
    group = 12 if symmetry == "terra" else 24
    frames = []
    for label, step, mesh, score, geometry in entries:
        digest = hashlib.sha256(mesh.read_bytes()).hexdigest()
        validate_score(score, digest,symmetry)
        if geometry["mesh_sha256"] != digest:
            raise ValueError("geometry certificate input hash differs")
        points, saved_faces, attributes, cells, materials = boundary.read_mesh(mesh)
        triangles, volume = exporter.material_interface(points, cells, materials)
        if ({tuple(sorted(f)) for f in triangles}
                != {tuple(sorted(f)) for f in saved_faces[attributes == 13]}):
            raise ValueError("saved Gamma differs from material interface")
        if (len(triangles) != score["chamber_panels"]
                or not math.isclose(group * volume, geometry["volume_full"], rel_tol=1e-12)):
            raise ValueError("original mesh count or volume differs from certificate")
        used = sorted({int(i) for face in triangles for i in face})
        remap = {i: j for j, i in enumerate(used)}
        frames.append(dict(step=step, label=label, mesh_sha256=digest,
            points=points[used].reshape(-1).tolist(),
            triangles=[remap[int(i)] for face in triangles for i in face],
            volume=group * volume, rho=score["rho"], k=score["k"], c=score["c"], q=score["q"],
            fullTriangles=score["full_panels"],
            solverRtol=score.get("solver_rtol", 1e-10),
            boundaryEdges=geometry[f"full{group}"]["boundary_edges"]))
    frames[-1]["label"] += " · 最后接受"
    return dict(runId=Path(last_path).name, symmetryLabel=("Clément · tetra T(12)" if symmetry=="terra" else "Clément · cubic O(24)"),
                scoreType="bem", rotations=(exporter.geometry.proper_tetrahedral_rotations() if symmetry=="terra"
                                            else exporter.geometry.proper_octahedral_rotations()),
                frames=frames, status=campaigns[-1]["status"], sourceCommits=[c["source"]["commit"] for c in campaigns],
                geometry_modified=False)


def validate_score(score, digest,symmetry="terra"):
    if symmetry not in ("terra", "cubic"):
        raise ValueError("unknown viewer symmetry")
    group = 12 if symmetry=="terra" else 24
    if score["mesh_sha256"] != digest:
        raise ValueError("BEM input hash differs from original mesh")
    if (score["group_order"] != group or score["subsampled"]
            or score["full_panels"] != group * score["chamber_panels"]):
        raise ValueError(f"viewer requires all original {'T12' if symmetry=='terra' else 'O24'} panels")
    diagnostics = score["diagnostics"]
    rtol = score.get("solver_rtol", 1e-10)
    if (rtol not in (1e-10, 1e-9, 1e-6)
            or diagnostics.get("solver_rtol", rtol) != rtol):
        raise ValueError("unknown or inconsistent BEM tolerance certificate")
    residuals = diagnostics["rhs_relative_residuals"] + diagnostics["cpu_rhs_relative_residuals"]
    if (not diagnostics["solver_converged"] or not diagnostics["cpu_residual_verified"]
            or len(diagnostics["rhs_relative_residuals"]) != 2
            or len(diagnostics["cpu_rhs_relative_residuals"]) != 2
            or not all(math.isfinite(v) and 0 <= v <= rtol for v in residuals)):
        raise ValueError("BEM checkpoint did not pass CPU/GPU residual acceptance")
    if (not all(math.isfinite(score[k]) for k in ("k", "c", "q", "rho"))
            or min(score["k"], score["q"]) <= 0
            or not math.isclose(score["rho"], abs(score["c"]) / math.sqrt(score["k"] * score["q"]), rel_tol=1e-13)):
        raise ValueError("BEM rho does not match k/c/q")


def render(data, output):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    if output.is_relative_to(ROOT):
        raise ValueError("write generated viewers outside the source repository")
    encoded = base64.b64encode(gzip.compress(pack_data(data), mtime=0)).decode()
    template = (EXAMPLE / "scripts/templates/carlos-shape.html").read_text()
    if template.count("__CARLOS_PAYLOAD__") != 1:
        raise ValueError("viewer template payload marker must occur once")
    html = template.replace("__CARLOS_PAYLOAD__", encoded)
    if len(html.encode()) >= 1_000_000:
        raise ValueError(f"all-panel viewer exceeds the inline size budget: {len(html.encode())} bytes")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)
    return dict(path=str(output), bytes=len(html.encode()), frames=len(data["frames"]),
                steps=[f["step"] for f in data["frames"]], geometry_modified=False,
                payload="lossless gzip of original float64 interface coordinates",
                initial_rho=data["frames"][0]["rho"], final_rho=data["frames"][-1]["rho"])


def pack_data(data):
    """Lossless float64/uint32 storage, without dropping any original panel."""
    metadata = dict(data, frames=[])
    buffers = []
    for frame in data["frames"]:
        points, triangles = frame["points"], frame["triangles"]
        if (len(points) % 3 or len(triangles) % 3
                or not all(math.isfinite(v) for v in points)
                or not all(isinstance(i, int) and 0 <= i < len(points)//3 for i in triangles)):
            raise ValueError("invalid original interface arrays")
        points, triangles = locality_order(points, triangles)
        # Oriented faces are cyclically reordered (never reflected), then sorted.
        # Store their exact connectivity as small nonnegative index differences.
        deltas, previous = [], 0
        for i in range(0, len(triangles), 3):
            a, b, c = triangles[i:i+3]
            deltas.extend((a-previous, b-a, c-a))
            previous = a
        index_width = 2 if len(points)//3 <= 65535 else 4
        info = {k: v for k, v in frame.items() if k not in ("points", "triangles")}
        info.update(pointValues=len(points), triangleValues=len(triangles), indexWidth=index_width)
        metadata["frames"].append(info)
        # Group byte planes to expose identical exponents/index high bytes to gzip.
        # This permutation is exactly reversible, not numeric quantization.
        for values, fmt, width in ((points, "d", 8),
                                   (deltas, "H" if index_width == 2 else "I", index_width)):
            raw = struct.pack(f"<{len(values)}{fmt}", *values)
            buffers.append(b"".join(raw[i::width] for i in range(width)))
    header = json.dumps(metadata, separators=(",", ":"), allow_nan=False).encode()
    return b"MSH3" + struct.pack("<I", len(header)) + header + b"".join(buffers)


def locality_order(points, triangles):
    """Reindex nodes for compression; retain every coordinate and oriented face."""
    adjacency = [set() for _ in range(len(points)//3)]
    faces = [triangles[i:i+3] for i in range(0, len(triangles), 3)]
    for a, b, c in faces:
        adjacency[a].update((b, c))
        adjacency[b].update((a, c))
        adjacency[c].update((a, b))
    order, seen = [], set()
    for seed in range(len(adjacency)):
        if seed in seen:
            continue
        queue = deque([seed])
        seen.add(seed)
        while queue:
            vertex = queue.popleft()
            order.append(vertex)
            for neighbor in sorted(adjacency[vertex]):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
    remap = {vertex: i for i, vertex in enumerate(order)}
    ordered_faces = []
    for face in faces:
        mapped = [remap[v] for v in face]
        index = mapped.index(min(mapped))
        ordered_faces.append(mapped[index:] + mapped[:index])
    ordered_faces.sort()
    return ([x for v in order for x in points[3*v:3*v+3]],
            [i for face in ordered_faces for i in face])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--remote")
    parser.add_argument("--symmetry", choices=("terra","cubic"),default="terra")
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--continuation", type=Path)
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--campaigns", type=Path, nargs="+")
    parser.add_argument("--reader-source", type=Path, default=None, help="remote Rodin checkout for --remote")
    parser.add_argument("--python", default="python3")
    parser.add_argument("--emit-data", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--record", type=Path)
    parser.add_argument("--data-cache", type=Path, help="retain compact Gamma data, not volume fields")
    args = parser.parse_args()
    old_paths = (args.pilot, args.continuation, args.resume)
    if args.campaigns:
        if any(old_paths):
            parser.error("--campaigns cannot be combined with the legacy history options")
        history_args = ["--campaigns", *map(str, args.campaigns)]
    else:
        if args.symmetry!="terra":
            parser.error("cubic viewers require --campaigns; legacy history is T12 only")
        if not all(old_paths):
            parser.error("supply --campaigns or all three legacy history options")
        history_args = ["--pilot", str(args.pilot), "--continuation", str(args.continuation),
                        "--resume", str(args.resume)]
    def gather():
        return (collect_chain(args.campaigns, args.reader_source,args.symmetry) if args.campaigns
                else collect(*old_paths, args.reader_source))
    if args.emit_data:
        print(json.dumps(gather(), separators=(",", ":")))
        return
    if args.output is None:
        parser.error("--output is required unless --emit-data is selected")
    if args.remote:
        # A read-only export through stdin; no remote checkout or result file is changed.
        import shlex
        if args.reader_source is None:
            parser.error("--remote requires --reader-source pointing to the remote Rodin checkout")
        entry = args.reader_source / "examples/ShapeOptimization/KelvinBall/scripts/kelvin-bem.py"
        command = ["env", f"KELVIN_SYMMETRY={args.symmetry}", "PYTHONDONTWRITEBYTECODE=1",
                   "OPENBLAS_NUM_THREADS=1", args.python, str(entry), "viewer", "--emit-data",
                   *history_args, "--symmetry", args.symmetry]
        result = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                                 args.remote, shlex.join(command)],
                                text=True, capture_output=True, check=True)
        data = json.loads(result.stdout)
    else:
        data = gather()
    if args.data_cache:
        if args.data_cache.exists():
            raise FileExistsError(args.data_cache)
        if args.data_cache.resolve().is_relative_to(ROOT):
            raise ValueError("write generated data outside the source repository")
        args.data_cache.parent.mkdir(parents=True, exist_ok=True)
        args.data_cache.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False))
    record = render(data, args.output)
    record.update(run_id=data["runId"], source_commits=data["sourceCommits"],
                  frame_mesh_sha256=[f["mesh_sha256"] for f in data["frames"]])
    if args.record:
        if args.record.exists():
            raise FileExistsError(args.record)
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
