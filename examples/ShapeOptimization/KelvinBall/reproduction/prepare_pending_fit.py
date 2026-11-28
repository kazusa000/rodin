"""Join archived KelvinBall scoring segments and reference the last unfinished fit.

Inputs are read-only. Snapshot numbering in each segment starts at zero.
"""
import argparse
import csv
import json
from pathlib import Path
import re


def prepare(destination, segments):
    records = []
    snapshots = []
    header = None
    for source in segments:
        source = source.resolve()
        data = (source / "kelvin-ball.csv").read_bytes()
        if not data.endswith(b"\n"):
            raise ValueError(f"Incomplete final CSV line in {source}")
        lines = data.decode().splitlines()
        table = csv.DictReader(lines)
        if header is None:
            header = lines[0]
        elif lines[0] != header:
            raise ValueError("Inconsistent CSV headers")
        local = 0
        for row in table:
            if int(row["iteration"]) != len(records) or None in row or None in row.values():
                raise ValueError("Source CSVs must cover consecutive scored states from zero")
            snapshots.append((source, local))
            records.append(lines[1 + local])
            local += 1
        if local == 0:
            raise ValueError(f"Empty score history in {source}")

    step = len(records)
    log = re.sub(r"\x1b\[[0-9;]*m", "", (source / "run.log").read_text(errors="replace"))
    iterations = list(re.finditer(r"Iteration\s+(\d+) of \d+", log))
    if not iterations or int(iterations[-1][1]) != step + 1:
        raise ValueError("Log does not end in the expected unfinished iteration")
    phase = log[iterations[-1].end():]
    if "Stage 9: Fitting the advected interface with SWIFT." not in phase:
        raise ValueError("Saved transport did not reach the SWIFT fit")
    match = re.search(r"MinSTCut classification.*?Solid cells:\s*(\d+).*?Fluid cells:\s*(\d+).*?Interface triangles:\s*(\d+)", phase, re.S)
    if not match:
        raise ValueError("Missing post-adaptation material counts")
    solid, fluid, interface = map(int, match.groups())
    snapshots.append((source, local))
    links = {}
    for global_step, (directory, local_index) in enumerate(snapshots):
        for field in ("mesh", "Advected"):
            stem = f"KelvinBall.Chamber.{field}."
            path = directory / f"{stem}{local_index:06}.h5"
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing {path}")
            links[f"{stem}{global_step:06}.h5"] = path
    destination.mkdir()
    for name, target in links.items():
        (destination / name).symlink_to(target)
    (destination / "kelvin-ball.csv").write_text("\n".join([header, *records]) + "\n")
    (destination / "pending-fit.txt").write_text(f"{step} {solid + fluid}\n")
    (destination / "sources.json").write_text(json.dumps({
        "segments": [str(x.resolve()) for x in segments], "last_complete_csv": step - 1,
        "pending_transport": step, "pending_input": str(snapshots[-1][0]),
        "pending_file_index": local, "after_adaptation": {
            "solid": solid, "fluid": fluid, "interface_triangles": interface}
    }, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("segments", nargs="+", type=Path)
    args = parser.parse_args()
    prepare(args.destination, args.segments)
