# Migration validation, 11 October 2026

Executing source: `d625db418e4cf7359a044cc1f60913a2e54fdab4`, clean on Mac, 28h and 2060.
Later evidence/documentation commits do not change the runtime package or native
solver. Build and result records retain the executing commit, rather than
claiming that a documentation commit ran these numerical cases.

The shared source is now [kazusa000/rodin, codex/clement-swift-bem](https://github.com/kazusa000/rodin/tree/codex/clement-swift-bem).
Use [BEM_WORKFLOW.md](BEM_WORKFLOW.md) for installation and commands. No separate
Project03 or fixed external BEM checkout is needed.

## Validation

- Mac: 102 Python checks passed; 10 CUDA checks explicitly skipped.
- 28h: 102 Python checks passed; 10 CUDA checks explicitly skipped, and all four
  native boundary, pressure-gauge, seed and thickness checks passed.
- 2060: all 112 Python checks passed, including actual CUDA checks.
- The native generator reproduced the unchanged fixed perturbed sphere with
  mesh SHA256 `7f95825378b02bdc6086398a7ebf5d611a009bfb3cd924c923a3396b2cba0fc6`.
- Seven same-input BEM measurements (initial, one MMG candidate, one SWIFT
  candidate and four FD diagnostics) reproduced the old `k,c,q,rho` exactly;
  all seven outer-wall tables also have identical byte hashes.
- Candidate mesh hashes, direction/penalty/geometry diagnostics and actual
  BEM `rho-P` acceptance decisions match the previous paired campaign.
  CPU and GPU true residual certificates meet the recorded effective tolerance.
  Original Stokes residuals are `1.0461e-14` and `1.0889e-14`.
- Both backends reproduced the old saved difficult target: candidate mesh
  hashes and complete quality/fitting/geometry audits are identical.
- The bundled viewer read the retained cubic history without an old 03 reader:
  18 original full-resolution frames, steps 0–17, final rho
  `0.0015002521841784946`, with no geometry modification.

| Shared input / backend | Full panels | BEM rho | Soft penalty P | Accepted |
|---|---:|---:|---:|---|
| Initial perturbed sphere | 4776 | 4.1506460999709566e-8 | 0 | baseline |
| MMG, one update | 15384 | 5.803767952054189e-8 | 5.694304508306868e-9 | yes |
| SWIFT, one update | 5256 | 4.6554504998347627e-7 | 0 | yes |

These are bounded migration comparisons, not a new full optimization or evidence
that SWIFT always gives higher-quality cells. In the saved difficult case,
minimum tetrahedron quality is 0.0007620624054099174 for MMG and
0.005414016514944722 for SWIFT. The existing FEM/BEM derivative magnitude
mismatch and open seams remain unchanged.

## Evidence and preservation

The compact, public checksums and measurements are in
[reproduction/migration-validation.json](reproduction/migration-validation.json).
Detailed retained records use the external workspace-relative directory
`out/project-03/clement-swift/migration-20261011-d625db41` on the owning nodes.
No HDF5, density arrays or matrices
were downloaded, committed or copied into Rodin. Only compact JSON/logs,
the bounded controller's exact input meshes and lossless interface viewer
data were returned.

The independent auditor is versioned as `tests/audit_migration_outputs.py`.
It compares retained outputs; it does not execute either solver or need
Project03 source. For example, provide the old paired output, new paired
output, old/new saved-target audits, new native-check records and the executing
pin through its explicit `--old`, `--new`, `--old-replay`, `--replay`,
`--checks`, `--pin`, and external `--output` flags.

Project03 stays at `067f8e9b3a42cfe6155c87406733a21122f615e8` in its original
integration worktree. Its main checkout, Rodin's original main checkout and
the existing root-workspace report/README changes are untouched. The retained
cubic and T12 campaign/final-mesh hashes were checked again after all runs;
they match the pre-migration values recorded in the evidence JSON.

The first imported snapshot `6e64eb22` exposed a missing `common` import in
the extracted geometry helper. It was repaired in `d625db41`, with a real-mesh
regression; numerical definitions stayed unchanged. That failed migration
snapshot remains separate and is not marked as a successful campaign.

No algorithm, numerical tolerance, thickness setting, reconstruction backend
implementation or old experiment was changed by this migration. The 109
numerical-definition checks and unchanged CUDA arithmetic body are recorded in
[reproduction/migration-provenance.json](reproduction/migration-provenance.json).
