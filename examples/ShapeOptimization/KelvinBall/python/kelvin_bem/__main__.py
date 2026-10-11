# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Dispatch existing workflows without introducing a second solver path."""
import importlib
import sys

COMMANDS = {
    "baseline": "clement_swift_baseline",
    "native": "clement_swift_native",
    "bem": "clement_swift_bem",
    "paired-check": "clement_swift_pair",
    "build-cuda": "build_bem_cuda",
    "export-chamber": "export_carlos_exact_chamber",
    "export-target": "export_carlos_target",
    "export-interface": "export_carlos_interface",
    "audit": "audit_clement_reconstruction",
    "geometry-check": "check_periodic_interface_geometry",
    "viewer": "render_clement_history",
    "render": "render_periodic_chamber_full24_views",
    "plot": "plot_carlos_history",
    "score": "bem_kcq_score",
    "scan": "bem_convergence_scan",
}


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("KelvinBall single-checkout commands: " + ", ".join(COMMANDS))
        print("Use <command> --help for the retained workflow arguments.")
        return 0
    command = sys.argv.pop(1)
    if command not in COMMANDS:
        print(f"Unknown KelvinBall command: {command}", file=sys.stderr)
        return 2
    module = importlib.import_module("." + COMMANDS[command], package=__package__)
    return module.main()


if __name__ == "__main__":
    raise SystemExit(main())
