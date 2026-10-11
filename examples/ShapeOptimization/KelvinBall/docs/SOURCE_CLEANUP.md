# Source cleanup after the single-checkout migration

Baseline: `7deabf12c31d56b8f553b259b880e427a51d68ca` on
`codex/clement-swift-bem`. The cleanup changes names, command organization,
unused imports and deployment configuration. It does not change native C++,
CUDA arithmetic, the numerical preset, acceptance thresholds or old experiments.

## Retired code

A whole-Rodin-tree reference search and a package/test call-graph review found
no consumer of the following operations outside their own implementation:

- Three projection result types and five associated functions in the old
  `periodic_interface_geometry.py`: paired scalar trace projection, H1
  correction extension, cubic/T12 periodic level-set projection and skeleton
  orbit correction. The eight definitions occupy 566 lines.
- `expand_full24_panels` and `expand_full12_panels`: uncalled wrappers;
  the actual shared `_expand_panels` operation remains unchanged.
- Unused imports carried over from the extracted historical drivers.

These operations were not part of the current optimization path and are
removed, not archived as an alternate solver. Their original implementation
is recoverable from the baseline Git commit. Rotation, seam, watertightness,
connectivity and sampled-thickness diagnostics remain.

Carlos's `Thickness.h` and the current `SoftThickness.h` serve different
native modes; both remain. Tests, the versioned sphere mesh and frozen
migration records also remain.

## Names and entry points

Modules now describe their responsibilities rather than experiment authors:
`native_worker`, `bem_worker`, `compare_backends`, `bem_solver`,
`bem_boundary`, `acceptance`, `run_support`, `periodic_geometry`,
`solid_connectivity`, `history_viewer` and named export/audit/render tools.
The exhaustive mapping is in
[cleanup-provenance.json](../reproduction/cleanup-provenance.json).
Author/source attribution stays in the provenance records and method documents.

The public command groups are:

- Main workflow: `native`, `bem`, `compare`, `viewer`.
- Setup/original baseline: `build-cuda`, `baseline`.
- Auxiliary tools: `diagnostics <tool>`.

Published `paired-check` and flat diagnostic commands forward to these
implementations with a notice. Old T12/cubic histories and generated HTML
viewers are not rewritten. Removed internal Python module names are not kept
as redundant wrapper files; external imports should use the new names.

The current controller still validates one accepted update per backend; this
cleanup does not add a long optimization loop. Tailscale IPs are configurable
deployment fields; their defaults and all numerical settings are unchanged.

The long Carlos native-method documentation is moved to
[CARLOS_METHOD.md](CARLOS_METHOD.md). The top-level README is the short shared
entry guide, with operational details in [BEM_WORKFLOW.md](../BEM_WORKFLOW.md).

## Verification

The existing migration AST check uses the rename mapping and verifies all
99 retained definitions against the original frozen hashes. Only import
declarations and the exact renamed module-filename guard are normalized;
no numerical constants or formulas are ignored. Ten explicitly retired
definitions must be absent.

Regression checks cover grouped/legacy command dispatch, worker arguments
and exit status, help from an unrelated directory, existing BEM solves,
geometry, thickness/merit decisions, saved-field IO and T12/cubic history
compatibility. A separate byte comparison verifies that native source,
CMake, CUDA, the sphere input, numerical preset and frozen migration evidence
are identical to the baseline.

Local verification on 11 October 2026: **141 checks passed, 10 CUDA checks
explicitly skipped**, independently with `KELVIN_SYMMETRY=cubic` and
`KELVIN_SYMMETRY=terra`. The latter checks the symmetry-aware BEM/geometry
and retained-history tools; it does not claim a native T12 optimization run.
All 26 Python modules are reachable from the grouped command entry point,
and local documentation links resolve. No new optimization, native rebuild
or GPU campaign was run for this source-only cleanup.
