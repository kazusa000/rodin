# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""One entry point for KelvinBall workers, backend comparison and diagnostics."""
import importlib
import sys

WORKFLOW_COMMANDS = {
    "native": ("native_worker", "Build/check the solver, generate an initial mesh, or propose one step"),
    "bem": ("bem_worker", "Measure one mesh and export six BEM outer-boundary velocity traces"),
    "compare": ("compare_backends", "Compare one accepted MMG/SWIFT update from the same mesh"),
    "viewer": ("history_viewer", "Build a rotating, iteration-slider viewer from saved histories"),
}
SETUP_COMMANDS = {
    "build-cuda": ("build_cuda", "Build the optional FP64 CUDA operator"),
    "baseline": ("native_baseline", "Check the original Carlos baseline without BEM"),
}
DIAGNOSTIC_COMMANDS = {
    "audit": "reconstruction_audit",
    "geometry": "geometry_check",
    "export-chamber": "export_chamber",
    "export-target": "export_target",
    "export-interface": "export_interface",
    "score": "bem_solver",
    "scan": "convergence_scan",
    "plot": "plot_history",
    "render": "render_shape",
}
# Preserve published command lines without maintaining duplicate modules.
LEGACY_COMMANDS = {
    "paired-check": ("compare",),
    "geometry-check": ("diagnostics", "geometry"),
    **{name: ("diagnostics", name) for name in DIAGNOSTIC_COMMANDS if name != "geometry"},
}


def show_help():
    print("Usage: kelvin-bem.py <command> [arguments]")
    for title, commands in (("Main workflow", WORKFLOW_COMMANDS), ("Setup and baseline", SETUP_COMMANDS)):
        print("\n" + title + ":")
        for name, (_, description) in commands.items():
            print(f"  {name:12} {description}")
    print("\nAuxiliary tools:\n  diagnostics  Geometry, reconstruction, exports, score/scan, plots and renders")
    print("\nStart with native build/checks; use compare for the bounded one-update backend check.")
    print("compare is not a long optimization loop. Use <command> --help for arguments.")


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        show_help()
        return 0
    command = sys.argv.pop(1)
    if command in LEGACY_COMMANDS:
        replacement = LEGACY_COMMANDS[command]
        print(f"Legacy command '{command}': use '{' '.join(replacement)}'.", file=sys.stderr)
        command = replacement[0]
        sys.argv[1:1] = replacement[1:]
    if command == "diagnostics":
        if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
            print("Usage: kelvin-bem.py diagnostics <tool> [arguments]")
            print("Tools: " + ", ".join(DIAGNOSTIC_COMMANDS))
            return 0
        tool = sys.argv.pop(1)
        if tool not in DIAGNOSTIC_COMMANDS:
            print(f"Unknown diagnostic tool: {tool}", file=sys.stderr)
            return 2
        module_name = DIAGNOSTIC_COMMANDS[tool]
    else:
        commands = {**WORKFLOW_COMMANDS, **SETUP_COMMANDS}
        if command not in commands:
            print(f"Unknown KelvinBall command: {command}", file=sys.stderr)
            return 2
        module_name = commands[command][0]
    module = importlib.import_module("." + module_name, package=__package__)
    return module.main()


if __name__ == "__main__":
    raise SystemExit(main())
